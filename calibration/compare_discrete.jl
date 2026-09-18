# Is the discrete deploy model worth fitting against, or is the continuous
# approximation good enough?
#
# Simulates both with the SAME parameters, the ones `calibrate_actuator.jl`
# fitted on the continuous model, and compares each against the robot log:
#
#   continuous  DyadBotComponents.NonidealReferenceDyadBot
#               CascadeController, ideal measurement plus a theta_bias term
#   discrete    BalansBotDeploy.IdentificationLoop
#               the deployed loop: 200 Hz DiscreteCascadeCore, the real
#               TiltEstimator, the device sensors, the limiter and the clock
#
# Both are driven with the reference the robot logged. The comparison is on the
# raw logged signals, so neither model is given a reconstructed tilt.
#
# Run:  julia +1.12 --project=. compare_discrete.jl data/run_clean.arrow 60 180

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))
using BalansBotDeploy, DiscreteComponents
import SynchToolkit          # the clock pass; DiscreteComponents loads it but does not export it
using DataInterpolations, Arrow, CSV, DataFrames, Statistics, Plots
using ModelingToolkit: t_nounits as t
mkpath(PLOTS)

path = length(ARGS) >= 1 ? ARGS[1] : joinpath(DATA, "run_clean.arrow")
t0 = length(ARGS) >= 3 ? parse(Float64, ARGS[2]) : 60.0
t1 = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 180.0

# ------------------------------------------------------------------ data ----
raw = read_log(path)
# convert_log gives timestamp, plant.theta (the deployed estimator replayed),
# plant.x and pos_ref in metres; the raw device channels come along for the
# comparison, because the model produces those itself.
conv = convert_log(raw)
data = hcat(conv, raw[!, [:accel_x, :accel_z, :gyro_y]])
win = data[(data.timestamp .>= t0) .& (data.timestamp .<= t1), :]
ts = win.timestamp .- win.timestamp[1]
println(nrow(win), " samples, ", round(ts[end], digits = 1), " s")

const REF = LinearInterpolation(win.pos_ref, ts; extrapolation = ExtrapolationType.Constant)
ref_fun(tt) = REF(tt)
@register_symbolic ref_fun(tt)

# The parameters the continuous fit found on this same window.
FIT = TOML.parsefile(joinpath(DATA, "calibrated_parameters_actuator.toml"))
println("parameters from the continuous fit: ", FIT)

# ----------------------------------------------------------------- models ----
"The continuous model, as calibrate_actuator.jl builds it."
function build_continuous()
    @named model = DyadBotComponents.NonidealReferenceDyadBot(;
        phi0 = win[1, "plant.theta"], theta_bias = FIT["theta_bias"],
        plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
        plant__Ib = FIT["Ib"], plant__Iw = FIT["Iw"], plant__d = FIT["d"],
        plant__k_scale = FIT["k_scale"], plant__db = FIT["db"], plant__tau_c = FIT["tau_c"],
        plant__motor__w_eps = 0.3,
        plant__wheelinertia__phi__initial = -win[1, "plant.x"] / MEAS.R, GAINS...)
    @named top = System([model.pos_reference ~ ref_fun(t)], t; systems = [model])
    return mtkcompile(top)
end

"The deployed loop with the reference as an input and the non-ideal plant."
function build_discrete()
    c = CFG["controller"]
    az = deg2rad(CFG["robot_log"]["angle_zero_deg"])
    @named model = BalansBotDeploy.IdentificationLoop(;
        phi0 = win[1, "plant.theta"],
        # theta_bias of the continuous model is the mismatch between where the IMU
        # sits and what the estimator calls zero, so it enters here as the mount.
        theta_mount = az + FIT["theta_bias"], angle_zero = az, alpha = CFG["robot_log"]["alpha"],
        k_angle = c["k_angle"], Ti_angle = c["Ti_angle"], Td_angle = c["Td_angle"],
        k_pos = c["k_pos"], Ti_pos = c["Ti_pos"], Td_pos = c["Td_pos"],
        controller__angle_controller__Nd = c["Nd_angle"], controller__angle_controller__wd = c["wd_angle"],
        controller__pos_controller__Nd = c["Nd_pos"],
        plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
        plant__Ib = FIT["Ib"], plant__Iw = FIT["Iw"], plant__d = FIT["d"],
        plant__k_scale = FIT["k_scale"], plant__db = FIT["db"], plant__tau_c = FIT["tau_c"],
        plant__motor__w_eps = 0.3,
        plant__wheelinertia__phi__initial = -win[1, "plant.x"] / MEAS.R)
    @named top = System([model.pos_reference ~ ref_fun(t)], t; systems = [model])
    return mtkcompile(top; additional_passes = [SynchToolkit.compile_lustre])
end

rms(a, b) = sqrt(sum(abs2, a .- b) / length(a))
results = Dict{String, Any}()

for (label, builder) in (("continuous", build_continuous), ("discrete", build_discrete))
    println("\n--- ", label, " ---")
    try
        @time "  build" sys = builder()
        @time "  solve" sol = solve(ODEProblem(sys, [], (0.0, ts[end])), Rodas5P();
                                    abstol = 1e-9, reltol = 1e-9)
        if !SciMLBase.successful_retcode(sol)
            println("  solve failed: ", sol.retcode); continue
        end
        x = sol(ts; idxs = sys.model.plant.x).u
        th = sol(ts; idxs = sys.model.plant.theta).u
        results[label] = (; x, th, steps = sol.stats.naccept)
        println("  retcode ", sol.retcode, ", ", sol.stats.naccept, " steps")
        println("  position error against the log: ", round(100rms(x, win[!, "plant.x"]), digits = 3), " cm")
    catch err
        println("  FAILED: ", first(sprint(showerror, err), 400))
        results[label] = nothing
    end
end

# ---------------------------------------------------------------- verdict ----
println()
if get(results, "discrete", nothing) === nothing
    println("The discrete model did not run; the failure is above.")
elseif get(results, "continuous", nothing) === nothing
    println("The continuous model did not run; the failure is above.")
else
    ec = 100rms(results["continuous"].x, win[!, "plant.x"])
    ed = 100rms(results["discrete"].x, win[!, "plant.x"])
    between = 100rms(results["continuous"].x, results["discrete"].x)
    println("position error against the log [cm]: continuous ", round(ec, digits = 3), "   discrete ", round(ed, digits = 3))
    println("difference between the 2 models [cm]: ", round(between, digits = 3))
    println(between < 0.3 * min(ec, ed) ?
        "The 2 models differ by much less than either differs from the robot: the continuous approximation is not what limits the fit." :
        "The 2 models differ by a useful fraction of the residual: fitting the discrete model is worth the cost.")
end

const C_DATA, C_CONT, C_DISC, C_REF = "#8a8985", "#eb6834", "#2a78d6", "#1baf7a"
p = plot(ts, 100 .* win.pos_ref; label = "reference", color = C_REF)
plot!(p, ts, 100 .* win[!, "plant.x"]; label = "robot", color = C_DATA)
haskey(results, "continuous") && results["continuous"] !== nothing &&
    plot!(p, ts, 100 .* results["continuous"].x; label = "continuous model", color = C_CONT)
haskey(results, "discrete") && results["discrete"] !== nothing &&
    plot!(p, ts, 100 .* results["discrete"].x; label = "discrete deploy model", color = C_DISC)
plot!(p; ylabel = "position [cm]", xlabel = "time [s]", title = "Same parameters, both models, against the robot",
      size = (1100, 500), lw = 1.5, framestyle = :box, gridalpha = 0.25,
      foreground_color_legend = nothing, background_color_legend = nothing, dpi = 130, left_margin = 8Plots.mm)
savefig(p, joinpath(PLOTS, "discrete_vs_continuous.png"))
println("plot written to ", joinpath(PLOTS, "discrete_vs_continuous.png"))
