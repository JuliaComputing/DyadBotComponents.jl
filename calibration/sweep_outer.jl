# Tune the OUTER loop on the calibrated twin to reduce the overshoot.
#
# On run_amp10 the robot reaches +-15 cm against a +-10 cm reference: about 5 cm
# past, every cycle, on both sides. The twin reproduces that (peak 16.0 cm), so
# it is the right place to ask what reduces it.
#
# The score is the PEAK EXCURSION PAST THE REFERENCE, not the RMS. RMS is
# dominated by the sweep between the levels, where most of the samples are, and
# the sweep is not what needs fixing.
#
# Rejected outright, whatever the overshoot:
#   saturation  |duty| >= 0.95 for more than 0.5 % of the run. A saturated robot
#               has NO steering, because `duty + y` and `duty - y` clamp to the
#               same value, and it turned 3079 counts in 29 s when that happened.
#   creep       the robot must still ARRIVE. A setting that overshoots less by
#               never getting there is not an improvement, so the settled
#               distance from the reference is reported next to the overshoot.
#
# Run:  julia +1.12 --project=. sweep_outer.jl [data/run_amp10.arrow 20 185]

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))
using DataInterpolations, DataFrames, Plots, Printf, Statistics
using ModelingToolkit: t_nounits as t
mkpath(PLOTS)

const K_DUTY, MAX_DUTY = 10.0, 1.0
path = length(ARGS) >= 1 ? ARGS[1] : joinpath(DATA, "run_amp10.arrow")
t0 = length(ARGS) >= 3 ? parse(Float64, ARGS[2]) : 20.0
t1 = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 185.0

data = convert_log(read_log(path))
win = data[(data.timestamp .>= t0) .& (data.timestamp .<= t1), :]; win.timestamp .-= win.timestamp[1]
const REF = LinearInterpolation(win.pos_ref, win.timestamp; extrapolation = ExtrapolationType.Constant)
ref_fun(tt) = REF(tt); @register_symbolic ref_fun(tt)
const F = TOML.parsefile(joinpath(DATA, "calibrated_parameters_actuator.toml"))
const A = maximum(abs, win.pos_ref)
const LVL = sign.(round.(win.pos_ref ./ A, digits = 1))

function build(; k_pos, Ti_pos, Td_pos)
    g = merge(GAINS, (controller__k_pos = k_pos, controller__Ti_pos = Ti_pos, controller__Td_pos = Td_pos))
    @named model = DyadBotComponents.NonidealReferenceDyadBot(; phi0 = win[1, "plant.theta"],
        theta_bias = F["theta_bias"], plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
        plant__Ib = F["Ib"], plant__Iw = F["Iw"], plant__d = F["d"], plant__k_scale = F["k_scale"],
        plant__db = F["db"], plant__tau_c = F["tau_c"], plant__tau_s = F["tau_s"], plant__w_s = F["w_s"],
        plant__motor__w_eps = 0.05,
        plant__wheelinertia__phi__initial = -win[1, "plant.x"] / MEAS.R, g...)
    @named top = System([model.pos_reference ~ ref_fun(t)], t; systems = [model])
    return mtkcompile(top)
end

"Overshoot past the reference level, settled distance from it, and saturation."
function score(ts, x, tau)
    duty = clamp.(K_DUTY .* tau, -MAX_DUTY, MAX_DUTY)
    sat = 100 * count(>=(0.95 * MAX_DUTY), abs.(duty)) / length(duty)
    over = Float64[]; settled = Float64[]
    for k in 2:length(LVL)
        (LVL[k] != LVL[k-1] && LVL[k] != 0) || continue
        tk = ts[k]; (tk < 10.0 || tk + 5.0 > ts[end]) && continue
        w = (ts .>= tk) .& (ts .<= tk + 4.0); s = (ts .>= tk + 3.5) .& (ts .<= tk + 5.0)
        (count(w) < 50 || count(s) < 20) && continue
        target = LVL[k] * A
        pk = LVL[k] > 0 ? maximum(x[w]) : minimum(x[w])
        push!(over, abs(pk) - abs(target))
        push!(settled, abs(mean(x[s])) - abs(target))
    end
    return (; sat, over = mean(over), settled = mean(settled), peak = maximum(abs, x))
end

println("robot: peak ", round(100maximum(abs, win[!, "plant.x"]), digits = 1), " cm against a reference of +-",
        round(100A, digits = 1), " cm\n")
println(" k_pos  Ti_pos  Td_pos | overshoot | settled | peak  | saturated | verdict")
println("-----------------------------------------------------------------------------")
results = Dict{Any, Any}()
for Td in (3.0, 6.0, 10.0), Ti in (6.0, 10.0), kp in (0.03, 0.05, 0.07)
    try
        sys = build(; k_pos = kp, Ti_pos = Ti, Td_pos = Td)
        sol = solve(ODEProblem(sys, [], (0.0, win.timestamp[end])), Rodas5P(); abstol = 1e-8, reltol = 1e-8)
        if !SciMLBase.successful_retcode(sol)
            @printf("%6.3f  %6.1f  %6.1f | solve failed: %s\n", kp, Ti, Td, sol.retcode); continue
        end
        ts = win.timestamp
        x = sol(ts; idxs = sys.model.plant.x).u
        tau = sol(ts; idxs = sys.model.controller.torque).u
        s = score(ts, x, tau)
        results[(kp, Ti, Td)] = (; s, x)
        v = s.sat > 0.5 ? "REJECT: saturates" : s.settled < -0.03 ? "REJECT: never arrives" : "ok"
        @printf("%6.3f  %6.1f  %6.1f | %+8.2f cm| %+6.2f cm| %5.1f | %6.2f %%  | %s\n",
                kp, Ti, Td, 100s.over, 100s.settled, 100s.peak, s.sat, v)
    catch err
        @printf("%6.3f  %6.1f  %6.1f | FAILED: %s\n", kp, Ti, Td, first(sprint(showerror, err), 100))
    end
end

ok = [(k, v) for (k, v) in results if v.s.sat <= 0.5 && v.s.settled >= -0.03]
if !isempty(ok)
    best = argmin(kv -> kv[2].s.over, ok)
    @printf("\nleast overshoot that still arrives and does not saturate:\n")
    @printf("  k_pos = %.3f, Ti_pos = %.1f, Td_pos = %.1f -> overshoot %+.2f cm, settled %+.2f cm, peak %.1f cm\n",
            best[1][1], best[1][2], best[1][3], 100best[2].s.over, 100best[2].s.settled, 100best[2].s.peak)
    cur = get(results, (0.05, 6.0, 6.0), nothing)
    cur === nothing || @printf("  the robot runs k_pos 0.05, Ti_pos 6, Td_pos 6 -> overshoot %+.2f cm, peak %.1f cm\n",
                               100cur.s.over, 100cur.s.peak)
end
println("\nThe twin under-sticks: it creeps 4.0 cm/s during the holds where the robot")
println("moves 2.7, so absolute overshoot is optimistic. Trust the ORDERING, not the values.")
