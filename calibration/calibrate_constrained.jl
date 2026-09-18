# Refit the plant with the bench ratios held fixed.
#
# calibrate_actuator.jl fits Ib, Iw, d, k_scale, db and tau_c to a closed-loop
# run. That fit is badly conditioned: Iw and d trade against each other, and the
# result (k_tau/d = 78 rad/s per duty) disagrees with the bench by 3.2 times.
# The bench measures the drive with the wheels off the ground and no feedback,
# so it determines the RATIOS directly:
#
#     k_tau / d = 24.0 rad/s per duty      (steady-speed line, motor A 24.2, B 23.8)
#     Iw / d    = 0.0473 s                 (step time constant, A 47.3 ms, B 50 ms)
#
# ONLY Iw/d IS USED. Pinning k_tau/d as well was tried on 2026-09-18 and
# failed: `k_scale` ran to its upper bound on both passes and the position error
# rose from 1.71 cm to 5.24 cm, worse than no fit at all. The reason is that the
# bench measured the steady-speed line over duties 0.1 to 0.9, that is wheel
# speeds of 2.4 to 21 rad/s. A straight line through that range folds the
# Coulomb friction into an effective viscous `d`. The balance loop lives 1 to 2
# orders of magnitude slower, where that `d` over-damps badly. So k_tau/d = 24
# is a large-signal number and does not describe the robot.
#
# Iw/d survives the same test: the bench step gives 47.3 ms and the closed-loop
# fit gives 48.7 ms independently, and the chirp is a clean first-order roll-off
# at 47 ms up to 11 Hz. That is the constraint worth keeping, and it is exactly
# the Iw-against-d trade that made the unconstrained fit ill-conditioned.
#
# This script therefore fits
#
#     Ib, d, k_scale, db, tau_c, theta_bias
#
# and DERIVES Iw = 0.0473 * d.
#
# Run:  julia +1.12 --project=. calibrate_constrained.jl data/run_clean.arrow 60 180

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))
using DyadModelOptimizer, DataInterpolations, DataFrames, Plots, Printf
using DiffEqBase: BrownFullBasicInit
using ModelingToolkit: t_nounits as t
mkpath(PLOTS)

const IW_OVER_D = 0.0473       # s, bench step and chirp, and the closed-loop fit
const K_DUTY = 10.0            # the runtime's duty per N*m, as the robot ran it

path = length(ARGS) >= 1 ? ARGS[1] : joinpath(DATA, "run_clean.arrow")
t0 = length(ARGS) >= 3 ? parse(Float64, ARGS[2]) : 60.0
t1 = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 180.0

data = convert_log(read_log(path))
win = data[(data.timestamp .>= t0) .& (data.timestamp .<= t1), :]
win.timestamp .-= win.timestamp[1]
@printf("%s, %.0f to %.0f s: %d samples\n", basename(path), t0, t1, nrow(win))

const REF_INTERP = LinearInterpolation(win.pos_ref, win.timestamp; extrapolation = ExtrapolationType.Constant)
ref_fun(tt) = REF_INTERP(tt)
@register_symbolic ref_fun(tt)

"The tracking model with Iw pinned to d through the bench time constant."
function build_pinned(d; Ib, k_scale, db, tau_c, theta_bias)
    @named model = DyadBotComponents.NonidealReferenceDyadBot(;
        phi0 = win[1, "plant.theta"], theta_bias,
        plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
        plant__Ib = Ib, plant__Iw = IW_OVER_D * d, plant__d = d,
        plant__k_scale = k_scale, plant__db = db, plant__tau_c = tau_c,
        plant__motor__w_eps = 0.3,
        plant__wheelinertia__phi__initial = -win[1, "plant.x"] / MEAS.R, GAINS...)
    @named top = System([model.pos_reference ~ ref_fun(t)], t; systems = [model])
    return mtkcompile(top)
end

# Start from the unconstrained fit, which is the best guess of the level.
PREV = TOML.parsefile(joinpath(DATA, "calibrated_parameters_actuator.toml"))
state = (Ib = PREV["Ib"], d = PREV["d"], k_scale = PREV["k_scale"], db = PREV["db"],
         tau_c = PREV["tau_c"], theta_bias = PREV["theta_bias"])

"Warn when a fitted value sits on a search bound: that is a failed fit, not a result."
function at_bounds(name, v, lo, hi)
    (v <= lo * 1.001 || v >= hi * 0.999) &&
        @printf("  !! %s = %.4g sits on its bound [%.3g, %.3g]: the fit is constrained, not converged\n", name, v, lo, hi)
end

rms(a, b) = sqrt(sum(abs2, a .- b) / length(a))
local sys, leaf, fitted
const BOUNDS = (Ib = (0.3 * NOMINAL.Ib, 3 * NOMINAL.Ib), d = (1e-4, 5e-2), k_scale = (0.3, 3.0),
                db = (2e-5, 0.05), tau_c = (1e-6, 0.05))
for iter in 1:4
    @printf("\n--- pass %d: Iw pinned to d = %.3e  ->  Iw = %.3e ---\n", iter, state.d, IW_OVER_D * state.d)
    global sys = build_pinned(state.d; Ib = state.Ib, k_scale = state.k_scale,
                              db = state.db, tau_c = state.tau_c, theta_bias = state.theta_bias)
    global leaf = (Ib = sys.model.plant.body_mass.I,
                   d = sys.model.plant.motor.damper.d,
                   k_scale = sys.model.plant.motor.actuator.k_scale,
                   db = sys.model.plant.motor.actuator.db,
                   tau_c = sys.model.plant.motor.actuator.tau_c,
                   theta_bias = sys.model.bias_source.k)
    exp1 = Experiment(win, sys;
                      depvars = [sys.model.plant.theta => Symbol("plant.theta"),
                                 sys.model.plant.x => Symbol("plant.x")],
                      alg = Rodas5P(), abstol = 1e-9, reltol = 1e-9,
                      initializealg = BrownFullBasicInit(), name = "tracking")
    invprob = InverseProblem(exp1, [
        leaf.Ib => (state.Ib, BOUNDS.Ib..., :log10),
        leaf.d => (state.d, BOUNDS.d..., :log10),
        leaf.k_scale => (state.k_scale, BOUNDS.k_scale..., :log10),
        leaf.db => (state.db, BOUNDS.db..., :log10),
        leaf.tau_c => (state.tau_c, BOUNDS.tau_c..., :log10),
        leaf.theta_bias => (state.theta_bias, -0.03, 0.03),
    ])
    @time "  calibrate" r = calibrate(invprob, SingleShooting(maxiters = 200))
    global fitted = (Ib = r[1], d = r[2], k_scale = r[3], db = r[4], tau_c = r[5], theta_bias = r[6])
    @printf("  Ib %.4e  d %.4e  k_scale %.4f  db %.3e  tau_c %.3e  theta_bias %+.3f deg\n",
            fitted.Ib, fitted.d, fitted.k_scale, fitted.db, fitted.tau_c, rad2deg(fitted.theta_bias))
    for k in keys(BOUNDS); at_bounds(String(k), getproperty(fitted, k), getproperty(BOUNDS, k)...); end
    moved = abs(fitted.d - state.d) / state.d
    global state = fitted
    if moved < 0.02
        @printf("  d moved %.1f %%: the pinning has converged.\n", 100moved)
        break
    end
    @printf("  d moved %.1f %%: rebuild with the new Iw.\n", 100moved)
end

sys = build_pinned(fitted.d; Ib = fitted.Ib, k_scale = fitted.k_scale, db = fitted.db,
                   tau_c = fitted.tau_c, theta_bias = fitted.theta_bias)
sol = solve(ODEProblem(sys, [], (0.0, win.timestamp[end])), Rodas5P(); abstol = 1e-9, reltol = 1e-9)
ts = win.timestamp
x_f = sol(ts; idxs = sys.model.plant.x).u
th_f = sol(ts; idxs = sys.model.plant.theta).u

println("\n=== constrained fit against the unconstrained one ===")
@printf("%-12s %-14s %-14s %s\n", "parameter", "unconstrained", "constrained", "note")
for (k, note) in ((:Ib, "free"), (:d, "free"), (:k_scale, "free, sets the level"), (:db, "free"), (:tau_c, "free"))
    @printf("%-12s %-14.4g %-14.4g %s\n", k, PREV[String(k)], getproperty(fitted, k), note)
end
@printf("%-12s %-14.4g %-14.4g %s\n", "Iw", PREV["Iw"], IW_OVER_D * fitted.d, "derived: Iw/d = $IW_OVER_D s")
@printf("%-12s %-14.4g %-14.4g %s\n", "k_tau/d", 0.1 * PREV["k_scale"] / PREV["d"], 0.1 * fitted.k_scale / fitted.d, "bench says 24.0 at high speed")
@printf("%-12s %-14.3f %-14.3f %s\n", "theta_bias", rad2deg(PREV["theta_bias"]), rad2deg(fitted.theta_bias), "deg")
@printf("\nk_duty the runtime should use: %.2f  (it runs %.1f)\n", 1 / (fitted.k_scale / K_DUTY), K_DUTY)
@printf("tilt RMS %.3f deg, position RMS %.3f cm against the robot\n",
        rad2deg(rms(th_f, win[!, "plant.theta"])), 100rms(x_f, win[!, "plant.x"]))

open(joinpath(DATA, "calibrated_parameters_constrained.toml"), "w") do io
    println(io, "# calibrate_constrained.jl on ", basename(path), ", $t0 to $t1 s.")
    println(io, "# Iw is DERIVED from d by the bench time constant Iw/d = $IW_OVER_D s.")
    for (k, v) in (("Ib", fitted.Ib), ("Iw", IW_OVER_D * fitted.d), ("d", fitted.d), ("k_scale", fitted.k_scale),
                   ("db", fitted.db), ("tau_c", fitted.tau_c), ("theta_bias", fitted.theta_bias))
        println(io, k, " = ", v)
    end
end

const C_DATA, C_FIT, C_REF = "#8a8985", "#2a78d6", "#1baf7a"
p1 = plot(ts, rad2deg.(win[!, "plant.theta"]); label = "robot", color = C_DATA)
plot!(p1, ts, rad2deg.(th_f); label = "twin, bench-constrained", color = C_FIT)
plot!(p1; ylabel = "tilt [deg]", title = "Bench-constrained twin against the robot")
p2 = plot(ts, 100 .* win.pos_ref; label = "reference", color = C_REF)
plot!(p2, ts, 100 .* win[!, "plant.x"]; label = "robot", color = C_DATA)
plot!(p2, ts, 100 .* x_f; label = "twin, bench-constrained", color = C_FIT)
plot!(p2; ylabel = "position [cm]", xlabel = "time [s]")
p = plot(p1, p2; layout = (2, 1), size = (1100, 750), lw = 1.5, framestyle = :box, gridalpha = 0.25,
         foreground_color_legend = nothing, background_color_legend = nothing, dpi = 130, left_margin = 8Plots.mm)
savefig(p, joinpath(PLOTS, "constrained_fit.png"))
println("plot written to ", joinpath(PLOTS, "constrained_fit.png"))
