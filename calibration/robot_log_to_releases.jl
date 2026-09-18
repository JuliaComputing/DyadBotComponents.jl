# Cut release tests out of a robot log for calibrate_plant.jl.
#
# The runtime writes one row per 200 Hz tick (BalansRuntime DyadTickRecord) with
# the IMU signals already in the model frame, the wheel position in metres, the
# position reference, and the torque command. This script
#
#   1. rebuilds the tilt with the same complementary filter as the deployed
#      controller (TiltEstimator: alpha, angle_zero), because the log has the
#      raw IMU signals but not the estimate,
#   2. rescales the position by R / wheel_radius_in_runtime, because the runtime
#      converts counts to metres with the radius baked into the deploy, not the
#      measured one,
#   3. plots the whole log so you can pick the hands-off windows, and
#   4. writes one data/release_robot_<n>.csv per window.
#
# A release window starts when you let go and ends when the robot is still. Do
# not include stretches where a hand is on the robot: they look like good
# balancing and are not.
#
# Run:
#   julia +1.12 --project=. robot_log_to_releases.jl run.arrow             # overview plot only
#   julia +1.12 --project=. robot_log_to_releases.jl run.arrow 12.0:15.0 31.2:34.2
#
# Windows are start:end in seconds from the first log row. run.csv works too.

include(joinpath(@__DIR__, "common.jl"))
using CSV, DataFrames, Arrow, Plots
mkpath(PLOTS); mkpath(DATA)

const LOGCFG = CFG["robot_log"]

function read_log(path)
    df = endswith(path, ".arrow") ? DataFrame(Arrow.Table(path)) : CSV.read(path, DataFrame)
    for c in ("now", "accel_x", "accel_z", "gyro_y", "position", "torque")
        c in names(df) || error("log has no column '$c'; columns are $(names(df))")
    end
    df.t = (Float64.(df.now) .- Float64(df.now[1])) ./ 1e9
    return df
end

"""
The deployed TiltEstimator, replayed over the log. It must match
`BalansBotDeploy/dyad/estimator.dyad` exactly, or every fit is against a tilt
the robot never saw.

`T_bias` is the gyro-bias state added on 2026-09-18. A plain complementary
filter does NOT reject a gyro offset: its high pass acts on the integrated
angle, so the gyro path is `tau/(1 + tau*s)` on the RATE with DC gain
`tau = alpha*Ts/(1 - alpha)`, and an offset settles at `tau` times itself. At
`alpha = 0.995` that is 0.995 s, so 1 deg/s of bias was 1 degree of tilt error
and 31 cm of position. Logs written BEFORE that change should be replayed with
`T_bias = Inf`, which disables the state.
"""
function estimate_tilt(df; alpha = LOGCFG["alpha"], Ts = LOGCFG["Ts"],
                       angle_zero = deg2rad(LOGCFG["angle_zero_deg"]),
                       T_bias = get(LOGCFG, "T_bias", 20.0))
    acc = atan.(df.accel_x, df.accel_z)
    raw = similar(acc)
    bias = zero(eltype(acc))
    raw[1] = acc[1]
    for k in 2:length(acc)
        bias += (Ts / T_bias) * (df.gyro_y[k] - bias)
        raw[k] = alpha * (raw[k - 1] + Ts * (df.gyro_y[k] - bias)) + (1 - alpha) * acc[k]
    end
    return raw .- angle_zero
end

function convert_log(df)
    scale = MEAS.R / LOGCFG["wheel_radius_in_runtime"]
    out = DataFrame("timestamp" => df.t, "plant.theta" => estimate_tilt(df),
                    "plant.x" => df.position .* scale, "torque" => df.torque,
                    "pos_ref" => "pos_ref" in names(df) ? df.pos_ref .* scale : zeros(nrow(df)))
    maxref = maximum(abs, out.pos_ref)
    maxref > 1e-6 && println("the log has a position reference of up to ", round(maxref, digits = 3), " m: use calibrate_tracking.jl, not the release fit")
    return out
end

function overview(out, path)
    p1 = plot(out.timestamp, rad2deg.(out[!, "plant.theta"]); label = "tilt estimate", color = "#2a78d6",
              ylabel = "tilt [deg]", title = basename(path))
    p2 = plot(out.timestamp, 100 .* out[!, "plant.x"]; label = "position (rescaled)", color = "#1baf7a", ylabel = "x [cm]")
    p3 = plot(out.timestamp, out.torque; label = "torque command", color = "#eb6834", ylabel = "N m", xlabel = "time [s]")
    p = plot(p1, p2, p3; layout = (3, 1), size = (1100, 800), lw = 1.5, framestyle = :box, gridalpha = 0.25,
             foreground_color_legend = nothing, background_color_legend = nothing, dpi = 130, left_margin = 6Plots.mm)
    file = joinpath(PLOTS, "robot_log_overview.png")
    savefig(p, file)
    println("overview written to ", file)
end

function cut_windows(out, windows)
    rm.(filter(p -> startswith(basename(p), "release_robot_"), readdir(DATA; join = true)))
    for (i, (t0, t1)) in enumerate(windows)
        rows = (out.timestamp .>= t0) .& (out.timestamp .<= t1)
        any(rows) || error("window $t0:$t1 has no rows; the log runs 0 to $(out.timestamp[end]) s")
        w = out[rows, ["timestamp", "plant.theta", "plant.x"]]
        w.timestamp .-= w.timestamp[1]
        CSV.write(joinpath(DATA, "release_robot_$(i).csv"), w)
        println("release $i: ", t0, " to ", t1, " s, ", nrow(w), " rows, release tilt ",
                round(rad2deg(w[1, "plant.theta"]), digits = 2), " deg at x = ", round(100 * w[1, "plant.x"], digits = 2), " cm")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: robot_log_to_releases.jl <run.arrow|run.csv> [start:end ...]")
    df = read_log(ARGS[1])
    println(nrow(df), " rows, ", round(df.t[end], digits = 1), " s, mean tick ",
            round(1000 * df.t[end] / (nrow(df) - 1), digits = 2), " ms")
    out = convert_log(df)
    overview(out, ARGS[1])
    windows = [(parse(Float64, split(a, ":")[1]), parse(Float64, split(a, ":")[2])) for a in ARGS[2:end]]
    isempty(windows) ? println("no windows given; look at the overview and rerun with start:end pairs") : cut_windows(out, windows)
end
