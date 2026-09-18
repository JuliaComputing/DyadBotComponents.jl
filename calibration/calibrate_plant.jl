# Calibrate the plant parameters of the balancing robot with DyadModelOptimizer.
#
# The model is `BalanceDyadBot`: the cascade-controlled planar robot with a
# position reference of 0, which is what the robot does on the floor. Each
# release test in data/ becomes one `Experiment`; all of them share one search
# space, so the three releases fit the same plant.
#
# Measured and fixed (measured_parameters.toml):  M, m, R, L and the controller gains
# Fitted:                                         Ib, Iw, d
#
# The three fitted parameters span decades, so the search runs on log10 of each.
#
# Run:
#   julia +1.12 --project=. calibrate_plant.jl                 # continuous data
#   julia +1.12 --project=. calibrate_plant.jl discrete        # 200 Hz sampled data
#   julia +1.12 --project=. calibrate_plant.jl robot           # releases cut from a robot log
#
# Inputs:  data/release_<tag>_*.csv from make_release_data.jl, or robot logs with
#          the columns timestamp, plant.theta, plant.x.
# Outputs: data/calibrated_parameters_<tag>.toml and plots/calibration_<tag>.png

include(joinpath(@__DIR__, "common.jl"))
using DyadModelOptimizer, CSV, DataFrames, Plots
# OrdinaryDiffEq v7 no longer exports the DAE initialization algorithms.
using DiffEqBase: BrownFullBasicInit
mkpath(PLOTS)

tag = length(ARGS) >= 1 && ARGS[1] in ("discrete", "robot") ? ARGS[1] : "continuous"

println("measured and fixed: ", MEAS)
println("controller gains:   ", GAINS)
println("starting guesses:   ", NOMINAL)

# ------------------------------------------------------------ experiments ----
paths = sort(filter(p -> startswith(basename(p), "release_$(tag)_"), readdir(DATA; join = true)))
isempty(paths) && error("no release data for tag '$tag' in $DATA; run make_release_data.jl first")

frames = [CSV.read(p, DataFrame) for p in paths]
# The release angle is read off the log, not fitted: the first tilt sample is the
# angle at which the robot was let go. It has to be baked into the model, so each
# release gets its own compiled system (see `build` in common.jl).
systems = map(frames) do df
    build(; phi0 = df[1, "plant.theta"], x0 = df[1, "plant.x"])
end
for (i, (p, df)) in enumerate(zip(paths, frames))
    println("experiment $i: ", basename(p), "  ", nrow(df), " rows, release ",
            round(df[1, "plant.theta"], digits = 4), " rad")
end

# `initializealg`: the model carries Bool parameters (the `render` flags of the
# multibody shapes). When the optimizer changes a parameter, the default
# initialization re-solve builds a parameter object of a different concrete type
# and the integrator rejects it (MTKParameters with Vector{Bool} vs BitVector).
# The differential states do not depend on the fitted parameters, so a plain
# algebraic-state initialization is enough and avoids that path.
experiments = map(enumerate(zip(frames, systems))) do (i, (df, sys))
    Experiment(df, sys; alg = Rodas5P(), abstol = 1e-10, reltol = 1e-10,
               initializealg = BrownFullBasicInit(), name = "release_$i")
end

# `mtkcompile` keeps only the leaf parameters, so the search space names those
# rather than plant.Ib, plant.Iw and plant.d. Every system is named `model`, so
# one set of symbols addresses all three experiments.
leaf = leaf_params(first(systems))
search_space = [
    leaf.Ib => (NOMINAL.Ib, 0.2 * NOMINAL.Ib, 5 * NOMINAL.Ib, :log10),
    leaf.Iw => (NOMINAL.Iw, 0.2 * NOMINAL.Iw, 20 * NOMINAL.Iw, :log10),
    leaf.d => (NOMINAL.d, 1e-5, 5e-3, :log10),
]
invprob = JointInverseProblem(experiments, search_space)

# -------------------------------------------------------------- calibrate ----
@time "calibrate" result = calibrate(invprob, SingleShooting(maxiters = 200))
println("\ncalibration result:\n", result)

fitted = (Ib = result[1], Iw = result[2], d = result[3])
# The search runs on log10 and the result comes back in the original units.
# Check that, so a change in the transform cannot make the numbers silently wrong.
all(1e-7 .< collect(fitted) .< 1e-1) ||
    error("calibration returned values outside the search bounds: $fitted; the log10 transform may not have been inverted")

open(joinpath(DATA, "calibrated_parameters_$(tag).toml"), "w") do io
    println(io, "# Fitted by calibrate_plant.jl on the '$tag' release data.")
    println(io, "Ib = ", fitted.Ib)
    println(io, "Iw = ", fitted.Iw)
    println(io, "d = ", fitted.d)
end

if tag == "robot"
    println("\nparameter   start        fitted       fitted / start")
    for k in (:Ib, :Iw, :d)
        f, n = getproperty(fitted, k), getproperty(NOMINAL, k)
        println(rpad(k, 11), rpad(round(n, sigdigits = 4), 13), rpad(round(f, sigdigits = 4), 13), round(f / n, digits = 2))
    end
else
    # The synthetic data has a known answer. TRUTH lives in make_release_data.jl,
    # repeated here so this script does not depend on it.
    truth = (Ib = 1.6 * NOMINAL.Ib, Iw = 2.5 * NOMINAL.Iw, d = 3.0e-4)
    println("\nparameter   start        fitted       true         error")
    for k in (:Ib, :Iw, :d)
        f, t, n = getproperty(fitted, k), getproperty(truth, k), getproperty(NOMINAL, k)
        println(rpad(k, 11), rpad(round(n, sigdigits = 4), 13), rpad(round(f, sigdigits = 4), 13),
                rpad(round(t, sigdigits = 4), 13), string(round(100 * (f - t) / t, digits = 1), " %"))
    end
end

# ------------------------------------------------------------------ plots ----
# Only the data and the calibrated model are plotted. The starting guess is far
# enough off to set the axis limits and hide the difference that matters; its
# error stays in the printed report above.
const C_DATA, C_FIT = "#8a8985", "#2a78d6"
rms(a, b) = sqrt(sum(abs2, a .- b) / length(a))
plts = []
for (i, (df, sys)) in enumerate(zip(frames, systems))
    ts = df.timestamp .- df.timestamp[1]
    s_start = run_release(sys, ts[end])
    s_fit = run_release(sys, ts[end];
                        overrides = [leaf.Ib => fitted.Ib, leaf.Iw => fitted.Iw, leaf.d => fitted.d])
    th_s, th_f = s_start(ts; idxs = sys.plant.theta).u, s_fit(ts; idxs = sys.plant.theta).u
    x_s, x_f = s_start(ts; idxs = sys.plant.x).u, s_fit(ts; idxs = sys.plant.x).u
    q1 = plot(ts, rad2deg.(df[!, "plant.theta"]); label = "data", color = C_DATA)
    plot!(q1, ts, rad2deg.(th_f); label = "calibrated", color = C_FIT)
    plot!(q1; ylabel = "tilt [deg]", title = "release $i")
    q2 = plot(ts, 100 .* df[!, "plant.x"]; label = "data", color = C_DATA)
    plot!(q2, ts, 100 .* x_f; label = "calibrated", color = C_FIT)
    plot!(q2; ylabel = "position [cm]", xlabel = "time [s]")
    push!(plts, plot(q1, q2; layout = (2, 1)))
    println("release $i  tilt RMS [deg]: start ", round(rad2deg(rms(th_s, df[!, "plant.theta"])), digits = 3),
            " calibrated ", round(rad2deg(rms(th_f, df[!, "plant.theta"])), digits = 3),
            "   position RMS [cm]: start ", round(100 * rms(x_s, df[!, "plant.x"]), digits = 3),
            " calibrated ", round(100 * rms(x_f, df[!, "plant.x"]), digits = 3))
end
p = plot(plts...; layout = (1, length(plts)), size = (450 * length(plts), 650), lw = 2,
         framestyle = :box, gridalpha = 0.25, foreground_color_legend = nothing,
         background_color_legend = nothing, dpi = 130, left_margin = 6Plots.mm)
savefig(p, joinpath(PLOTS, "calibration_$(tag).png"))
println("\nplot written to ", joinpath(PLOTS, "calibration_$(tag).png"))
