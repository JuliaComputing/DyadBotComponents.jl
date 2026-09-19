# Plot one robot run: what it was asked to do, what it did, and what it cost.
#
# Run:  julia +1.12 --project=. plot_run.jl data/run_amp10.arrow [t0 t1]
#
# 4 panels, all on one time axis:
#   position   reference against the wheel position, the thing being controlled
#   tilt       the estimator output and the accelerometer angle it is built from
#   duty       the command, with the limit and the break-away band marked
#   heading    from the gyro and from the wheel difference, which disagree when
#              the robot slips

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))
using Arrow, DataFrames, Statistics, Printf, Plots
mkpath(PLOTS)

path = length(ARGS) >= 1 ? ARGS[1] : joinpath(DATA, "run_amp10.arrow")
raw = read_log(path); conv = convert_log(raw)
t = conv.timestamp
t0 = length(ARGS) >= 3 ? parse(Float64, ARGS[2]) : 0.0
t1 = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : t[end]
sel = (t .>= t0) .& (t .<= t1)

ts = t[sel]
x = 100 .* conv[!, "plant.x"][sel]
ref = 100 .* conv.pos_ref[sel]
tilt = rad2deg.(conv[!, "plant.theta"][sel])
acc = rad2deg.(atan.(raw.accel_x, raw.accel_z) .- deg2rad(CFG["robot_log"]["angle_zero_deg"]))[sel]
duty = Float64.(raw.duty)[sel]
gy = Float64.(raw.raw_gyro_y)
head = (cumsum(gy .- median(gy[1:min(400, end)])) .* 0.005)[sel]
d = Float64.(raw.wheel_a) .- Float64.(raw.wheel_b); d .-= d[1]
head_enc = rad2deg.(d ./ 1320 .* 2pi .* MEAS.R ./ 0.13)[sel]

const C_REF, C_ROBOT, C_ACC, C_DUTY, C_ENC = "#1baf7a", "#2a78d6", "#8a8985", "#eb6834", "#9b59b6"
base = (framestyle = :box, gridalpha = 0.25, lw = 1.3,
        foreground_color_legend = nothing, background_color_legend = nothing)

p1 = plot(ts, ref; label = "reference", color = C_REF, ls = :dash, base...)
plot!(p1, ts, x; label = "wheel position", color = C_ROBOT, base...)
plot!(p1; ylabel = "position [cm]", title = basename(path) * "  (" * string(round(t[end], digits = 0)) * " s)")

p2 = plot(ts, acc; label = "accelerometer angle", color = C_ACC, lw = 0.8,
          framestyle = :box, gridalpha = 0.25,
          foreground_color_legend = nothing, background_color_legend = nothing)
plot!(p2, ts, tilt; label = "estimator tilt", color = C_ROBOT, lw = 1.3)
plot!(p2; ylabel = "tilt [deg]", ylims = (-8, 8))

p3 = plot(ts, duty; label = "duty", color = C_DUTY, base...)
hline!(p3, [0.077, -0.077]; label = "break-away", color = "#c0392b", ls = :dot, lw = 1)
hline!(p3, [0.122, -0.122]; label = "", color = "#c0392b", ls = :dot, lw = 1)
hline!(p3, [1.0, -1.0]; label = "limit", color = "#000000", ls = :dash, lw = 1)
plot!(p3; ylabel = "duty", ylims = (-1.05, 1.05))

p4 = plot(ts, head; label = "heading, gyro", color = C_ROBOT, base...)
plot!(p4, ts, head_enc; label = "heading the wheels imply", color = C_ENC, base...)
plot!(p4; ylabel = "heading [deg]", xlabel = "time [s]")

p = plot(p1, p2, p3, p4; layout = (4, 1), size = (1300, 1000), dpi = 130, left_margin = 10Plots.mm)
out = joinpath(PLOTS, "trajectory_" * replace(basename(path), ".arrow" => "") * ".png")
savefig(p, out)
println("written to ", out)
