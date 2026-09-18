# Fit the plant to a robot run with a position-reference square wave, and show
# how closely the model follows the robot.
#
# The robot runs BalanceDyadBot's loop with the filtered square-wave reference
# that the runtime generates (options ref_amplitude, ref_period, ref_start,
# ref_filter). The model is `ReferenceDyadBot`, driven with the reference the
# robot LOGGED, so both see the same excitation sample by sample. One long
# window is one `Experiment`: the loop is stable, so the whole run can be
# compared, unlike the open-loop plant.
#
# Measured and fixed:  M, m, R, L, controller gains (measured_parameters.toml)
# Fitted:              Ib, Iw, d  (log10 search)
#
# Run:
#   julia +1.12 --project=. calibrate_tracking.jl data/run.arrow [t_start t_end]
#   julia +1.12 --project=. calibrate_tracking.jl synthetic                # dry run on model-made data
#
# The window defaults to 1 s before the first reference step until the end of
# the log. Start it while the robot stands calm, because the model's controller
# starts with empty integrators.
# Outputs: data/calibrated_parameters_tracking.toml, plots/tracking_fit.png

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))   # read_log, convert_log (its main does not run)
using DyadModelOptimizer, DataInterpolations, CSV, DataFrames, Plots, Random
using DiffEqBase: BrownFullBasicInit
using ModelingToolkit: t_nounits as t
mkpath(PLOTS)

# ---------------------------------------------------------------- data ----
synthetic = !isempty(ARGS) && ARGS[1] == "synthetic"
if synthetic
    # A 40 s run of the model itself with the 'true' parameters of make_release_data.jl
    # and the robot's reference settings, plus measurement noise.
    const TRUTH = (Ib = 1.6 * NOMINAL.Ib, Iw = 2.5 * NOMINAL.Iw, d = 3.0e-4)
    ts = collect(0.0:0.005:40.0)
    sq = [tt <= 5.0 ? 0.0 : 0.05 * (4floor((tt - 5) / 10) - 2floor(2(tt - 5) / 10) + 1) for tt in ts]
    function filter2(sq; g = 1 - exp(-0.005 / 0.1))
        x1 = 0.0; x2 = 0.0; ref = similar(sq)
        for (k, u) in enumerate(sq); x1 += g * (u - x1); x2 += g * (x1 - x2); ref[k] = x2; end
        return ref
    end
    ref = filter2(sq)
    data = DataFrame("timestamp" => ts, "pos_ref" => ref, "plant.theta" => zeros(length(ts)), "plant.x" => zeros(length(ts)))
else
    isempty(ARGS) && error("usage: calibrate_tracking.jl <run.arrow|run.csv> [t_start t_end]  or  calibrate_tracking.jl synthetic")
    data = convert_log(read_log(ARGS[1]))
end

# ---------------------------------------------------------------- window ----
first_step = findfirst(abs.(data.pos_ref) .> 1e-4)
t_start = length(ARGS) >= 3 ? parse(Float64, ARGS[2]) : (first_step === nothing ? 0.0 : max(0.0, data.timestamp[first_step] - 1.0))
t_end = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : data.timestamp[end]
win = data[(data.timestamp .>= t_start) .& (data.timestamp .<= t_end), :]
win.timestamp .-= win.timestamp[1]
println("window ", t_start, " to ", t_end, " s: ", nrow(win), " samples")

# The logged reference drives the model through a registered interpolation.
const REF_INTERP = ConstantInterpolation(win.pos_ref, win.timestamp; extrapolation = ExtrapolationType.Constant)
ref_fun(tt) = REF_INTERP(tt)
@register_symbolic ref_fun(tt)

"ReferenceDyadBot with the measured plant, the robot's gains, the logged reference, and the logged initial state."
function build_tracking(; phi0, x0, Ib = NOMINAL.Ib, Iw = NOMINAL.Iw, d = NOMINAL.d)
    @named model = DyadBotComponents.ReferenceDyadBot(; phi0,
        plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
        plant__Ib = Ib, plant__Iw = Iw, plant__d = d,
        plant__wheelinertia__phi__initial = -x0 / MEAS.R, GAINS...)
    @named top = System([model.pos_reference ~ ref_fun(t)], t; systems = [model])
    return mtkcompile(top)
end

if synthetic
    sys_true = build_tracking(; phi0 = 0.0, x0 = 0.0, TRUTH...)
    sol = solve(ODEProblem(sys_true, [], (0.0, win.timestamp[end])), Rodas5P(); abstol = 1e-10, reltol = 1e-10)
    rng = Xoshiro(3)
    win[!, "plant.theta"] = sol(win.timestamp; idxs = sys_true.model.plant.theta).u .+ 0.002 .* randn(rng, nrow(win))
    win[!, "plant.x"] = sol(win.timestamp; idxs = sys_true.model.plant.x).u .+ 0.0005 .* randn(rng, nrow(win))
    println("synthetic truth: ", TRUTH)
end

sys = build_tracking(; phi0 = win[1, "plant.theta"], x0 = win[1, "plant.x"])
leaf = (Ib = sys.model.plant.body_mass.I, Iw = sys.model.plant.wheelinertia.I, d = sys.model.plant.simplemotor.damper.d)

# --------------------------------------------------------------- fit ----
exp1 = Experiment(win, sys; depvars = [sys.model.plant.theta => Symbol("plant.theta"), sys.model.plant.x => Symbol("plant.x")],
                  alg = Rodas5P(), abstol = 1e-9, reltol = 1e-9, initializealg = BrownFullBasicInit(), name = "tracking")
invprob = InverseProblem(exp1, [
    leaf.Ib => (NOMINAL.Ib, 0.2 * NOMINAL.Ib, 5 * NOMINAL.Ib, :log10),
    leaf.Iw => (NOMINAL.Iw, 0.2 * NOMINAL.Iw, 20 * NOMINAL.Iw, :log10),
    leaf.d => (NOMINAL.d, 1e-5, 5e-3, :log10),
])
@time "calibrate" result = calibrate(invprob, SingleShooting(maxiters = 200))
fitted = (Ib = result[1], Iw = result[2], d = result[3])
all(1e-7 .< collect(fitted) .< 1e-1) || error("calibration returned values outside the search bounds: $fitted")
println("\nparameter   start        fitted       fitted/start", synthetic ? "   true" : "")
for k in (:Ib, :Iw, :d)
    f, n = getproperty(fitted, k), getproperty(NOMINAL, k)
    print(rpad(k, 11), rpad(round(n, sigdigits = 4), 13), rpad(round(f, sigdigits = 4), 13), rpad(round(f / n, digits = 2), 15))
    synthetic ? println(round(getproperty(TRUTH, k), sigdigits = 4)) : println()
end
open(joinpath(DATA, "calibrated_parameters_tracking.toml"), "w") do io
    println(io, "# Fitted by calibrate_tracking.jl on ", synthetic ? "synthetic data" : ARGS[1], ", window $t_start to $t_end s.")
    println(io, "Ib = ", fitted.Ib); println(io, "Iw = ", fitted.Iw); println(io, "d = ", fitted.d)
end

# ------------------------------------------------------------- compare ----
run_model(pars) = solve(ODEProblem(sys, pars, (0.0, win.timestamp[end])), Rodas5P(); abstol = 1e-9, reltol = 1e-9)
s_start = run_model([])
s_fit = run_model([leaf.Ib => fitted.Ib, leaf.Iw => fitted.Iw, leaf.d => fitted.d])
ts = win.timestamp
rms(a, b) = sqrt(sum(abs2, a .- b) / length(a))
th_s, th_f = s_start(ts; idxs = sys.model.plant.theta).u, s_fit(ts; idxs = sys.model.plant.theta).u
x_s, x_f = s_start(ts; idxs = sys.model.plant.x).u, s_fit(ts; idxs = sys.model.plant.x).u
println("\ntilt RMS error [deg]:     start ", round(rad2deg(rms(th_s, win[!, "plant.theta"])), digits = 3), "   calibrated ", round(rad2deg(rms(th_f, win[!, "plant.theta"])), digits = 3))
println("position RMS error [cm]:  start ", round(100rms(x_s, win[!, "plant.x"]), digits = 3), "   calibrated ", round(100rms(x_f, win[!, "plant.x"]), digits = 3))
println("position RMS around the reference [cm]: robot ", round(100rms(win[!, "plant.x"], win.pos_ref), digits = 2), "   calibrated model ", round(100rms(x_f, win.pos_ref), digits = 2))

# Only the data and the calibrated model are plotted. The starting guess is far
# enough off to set the axis limits and hide the difference that matters; its
# error stays in the printed report above.
const C_DATA, C_FIT, C_REF = "#8a8985", "#2a78d6", "#1baf7a"
p1 = plot(ts, rad2deg.(win[!, "plant.theta"]); label = synthetic ? "synthetic robot" : "robot", color = C_DATA)
plot!(p1, ts, rad2deg.(th_f); label = "model, calibrated", color = C_FIT)
plot!(p1; ylabel = "tilt [deg]", title = "Square-wave tracking: robot against model")
p2 = plot(ts, 100 .* win.pos_ref; label = "reference", color = C_REF)
plot!(p2, ts, 100 .* win[!, "plant.x"]; label = synthetic ? "synthetic robot" : "robot", color = C_DATA)
plot!(p2, ts, 100 .* x_f; label = "model, calibrated", color = C_FIT)
plot!(p2; ylabel = "position [cm]", xlabel = "time [s]")
p = plot(p1, p2; layout = (2, 1), size = (1100, 750), lw = 1.5, framestyle = :box, gridalpha = 0.25,
         foreground_color_legend = nothing, background_color_legend = nothing, dpi = 130, left_margin = 8Plots.mm)
savefig(p, joinpath(PLOTS, synthetic ? "tracking_fit_synthetic.png" : "tracking_fit.png"))
println("plot written to ", joinpath(PLOTS, synthetic ? "tracking_fit_synthetic.png" : "tracking_fit.png"))
