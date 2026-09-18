# Sweep the outer position loop on the calibrated twin, and report the duty
# ceiling as a first-class result.
#
# The robot found the ceiling the hard way on 2026-09-18: k_pos = 0.04 with
# Ti_pos = 3 spent 16 % of its ticks with the duty on the rail, and a saturated
# robot has NO steering, because `duty + y` and `duty - y` clamp to the same
# value. It turned 3079 counts in 29 s. So every candidate here is judged on 3
# things and a candidate that saturates is rejected whatever its tracking looks
# like:
#
#   saturation   the fraction of the run with |duty| >= 0.95 * max_duty
#   park         how far past the +-5 cm target the robot settles
#   overshoot    the peak past the settled value
#
# The twin is DyadBotComponents.NonidealReferenceDyadBot with the parameters
# from calibrate_constrained.jl, driven by the reference the robot logged. It
# reproduces run_clean to 1.71 cm RMS, which is the residual to keep in mind:
# differences between candidates smaller than that are not real.
#
# Run:  julia +1.12 --project=. sweep_kpos.jl [data/run_clean.arrow 60 180]

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))
using DataInterpolations, DataFrames, Plots, Printf, Statistics
using ModelingToolkit: t_nounits as t
mkpath(PLOTS)

const K_DUTY = 10.0            # the runtime's duty per N*m, as the robot runs it
const MAX_DUTY = 1.0

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

const FIT = TOML.parsefile(joinpath(DATA, "calibrated_parameters_constrained.toml"))
println("twin parameters: ", FIT)

"The calibrated twin with one pair of outer-loop gains."
function build(k_pos, Ti_pos)
    gains = merge(GAINS, (controller__k_pos = k_pos, controller__Ti_pos = Ti_pos))
    @named model = DyadBotComponents.NonidealReferenceDyadBot(;
        phi0 = win[1, "plant.theta"], theta_bias = FIT["theta_bias"],
        plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
        plant__Ib = FIT["Ib"], plant__Iw = FIT["Iw"], plant__d = FIT["d"],
        plant__k_scale = FIT["k_scale"], plant__db = FIT["db"], plant__tau_c = FIT["tau_c"],
        plant__motor__w_eps = 0.3,
        plant__wheelinertia__phi__initial = -win[1, "plant.x"] / MEAS.R, gains...)
    @named top = System([model.pos_reference ~ ref_fun(t)], t; systems = [model])
    return mtkcompile(top)
end

"The 3 numbers that decide a candidate, from one solved run."
function score(ts, x, tau)
    duty = clamp.(K_DUTY .* tau, -MAX_DUTY, MAX_DUTY)
    sat = 100 * count(>=(0.95 * MAX_DUTY), abs.(duty)) / length(duty)
    A = maximum(abs, win.pos_ref)
    lvl = sign.(round.(win.pos_ref ./ A, digits = 1))
    park = Float64[]; over = Float64[]
    for k in 2:length(lvl)
        (lvl[k] != lvl[k - 1] && lvl[k] != 0) || continue
        tk = win.timestamp[k]
        (tk < 15.0 || tk + 5.0 > ts[end]) && continue
        w = (ts .>= tk) .& (ts .<= tk + 4.0)
        s = (ts .>= tk + 3.5) .& (ts .<= tk + 5.0)
        (count(w) < 50 || count(s) < 20) && continue
        xf = mean(x[s])
        pk = lvl[k] > 0 ? maximum(x[w]) : minimum(x[w])
        push!(park, abs(xf - lvl[k] * A))
        push!(over, abs(pk - xf))
    end
    return (; sat, park = isempty(park) ? NaN : mean(park),
            over = isempty(over) ? NaN : mean(over),
            rms = sqrt(mean(abs2, x .- win.pos_ref)), peak_duty = maximum(abs, duty))
end

const K_POS = (0.033, 0.04, 0.05, 0.06, 0.07, 0.09, 0.12)
const TI_POS = (10.0, 6.0)

println("\n k_pos  Ti_pos | saturated | park  | overshoot | RMS   | peak duty | verdict")
println("-------------------------------------------------------------------------------")
results = Dict{Tuple{Float64, Float64}, Any}()
for Ti in TI_POS, kp in K_POS
    try
        sys = build(kp, Ti)
        sol = solve(ODEProblem(sys, [], (0.0, win.timestamp[end])), Rodas5P(); abstol = 1e-8, reltol = 1e-8)
        if !SciMLBase.successful_retcode(sol)
            @printf("%6.3f  %5.1f  | %s\n", kp, Ti, "solve failed: $(sol.retcode)")
            continue
        end
        ts = win.timestamp
        x = sol(ts; idxs = sys.model.plant.x).u
        tau = sol(ts; idxs = sys.model.controller.torque).u
        s = score(ts, x, tau)
        results[(kp, Ti)] = (; s, x)
        verdict = s.sat > 0.5 ? "REJECT: saturates" : s.sat > 0.05 ? "marginal" : "ok"
        @printf("%6.3f  %5.1f  |  %6.2f %%  | %4.1f cm |  %4.2f cm  | %4.2f cm|   %5.3f   | %s\n",
                kp, Ti, s.sat, 100s.park, 100s.over, 100s.rms, s.peak_duty, verdict)
    catch err
        @printf("%6.3f  %5.1f  | FAILED: %s\n", kp, Ti, first(sprint(showerror, err), 120))
    end
end

println("\nThe twin reproduces run_clean to 1.71 cm RMS. Differences smaller than that")
println("between 2 candidates are inside the model error and mean nothing.")

ok = [(k, v) for (k, v) in results if v.s.sat <= 0.05 && isfinite(v.s.park)]
if !isempty(ok)
    # argmin(f, collection) returns the ELEMENT, not an index into it.
    best = argmin(kv -> kv[2].s.park, ok)
    @printf("\nsmallest park with no saturation: k_pos = %.3f, Ti_pos = %.1f -> park %.1f cm, overshoot %.2f cm, peak duty %.3f\n",
            best[1][1], best[1][2], 100best[2].s.park, 100best[2].s.over, best[2].s.peak_duty)
end

C = ["#8a8985", "#2a78d6", "#eb6834", "#1baf7a", "#9b59b6", "#d4a017", "#c0392b"]
p = plot(win.timestamp, 100 .* win.pos_ref; label = "reference", color = "#000000", ls = :dash, lw = 1)
plot!(p, win.timestamp, 100 .* win[!, "plant.x"]; label = "robot (k_pos 0.033)", color = C[1], lw = 2)
for (i, kp) in enumerate(K_POS)
    haskey(results, (kp, 10.0)) || continue
    plot!(p, win.timestamp, 100 .* results[(kp, 10.0)].x; label = "twin k_pos $kp", color = C[min(i + 1, end)], lw = 1.2)
end
plot!(p; ylabel = "position [cm]", xlabel = "time [s]", title = "k_pos sweep on the calibrated twin, Ti_pos = 10",
      size = (1300, 600), framestyle = :box, gridalpha = 0.25, xlims = (0, 60),
      foreground_color_legend = nothing, background_color_legend = nothing, dpi = 130, left_margin = 8Plots.mm)
savefig(p, joinpath(PLOTS, "kpos_sweep.png"))
println("plot written to ", joinpath(PLOTS, "kpos_sweep.png"))
