# One command to tune the controller against a robot run.
#
#     julia +1.12 --project=. tune.jl data/run_amp10.arrow
#     julia +1.12 --project=. tune.jl data/run_amp10.arrow --refit
#
# It does 3 things, in order, and stops at the first one that fails:
#
#   1 CHECK   that the gains in measured_parameters.toml match the run. A fit or
#             a sweep against the wrong controller is worthless: the plant
#             parameters absorb the mismatch and the twin then falls over at the
#             gains the robot runs. That happened on 2026-09-18 and cost an
#             invalid fit, an invalid sweep and a commit that had to be
#             corrected. THIS CHECK IS THE POINT OF THE SCRIPT.
#   2 FIT     only with --refit. 49 minutes. Needed when the ROBOT changes, not
#             when the gains change: the gains are inputs to the fit, not
#             outputs, so a gain change needs step 1 and step 3 only.
#   3 SWEEP   3 minutes. Scores candidates on overshoot, rejects anything that
#             saturates or fails to arrive, and prints what to set.
#
# The gains cannot be read from the log, because the runtime does not record
# them. Until it does, `tune.jl` compares against a file you must keep correct,
# and prints the robot's options.toml side by side if the robot is reachable.

include(joinpath(@__DIR__, "common.jl"))
using Printf

const HOST = get(ENV, "BALANSBOT_HOST", "balansbot@raspberrypi.local")

"Read the robot's [controller_params], or nothing if it cannot be reached."
function robot_gains()
    try
        out = read(`ssh -o ConnectTimeout=6 -o BatchMode=yes $HOST
                    "sed -n '/controller_params/,\$p' ~/BalansBot/BalansBotDyad/BalansRuntime/options.toml"`, String)
        g = Dict{String, Float64}()
        for line in split(out, '\n')
            m = match(r"^\"controller\.([A-Za-z_.]+)\"\s*=\s*([-0-9.eE]+)", strip(line))
            m === nothing || (g[m.captures[1]] = parse(Float64, m.captures[2]))
        end
        return isempty(g) ? nothing : g
    catch
        return nothing
    end
end

"""
    check_gains() -> Bool

Compare measured_parameters.toml with the robot. Returns false and explains if
they disagree. A disagreement is not a warning: every number downstream of it is
wrong, so the script stops.
"""
function check_gains()
    c = CFG["controller"]
    mine = Dict("k_angle" => c["k_angle"], "Ti_angle" => c["Ti_angle"], "Td_angle" => c["Td_angle"],
                "k_pos" => c["k_pos"], "Ti_pos" => c["Ti_pos"], "Td_pos" => c["Td_pos"],
                "pos_controller.Nd" => c["Nd_pos"], "angle_controller.Nd" => c["Nd_angle"],
                "pos_controller.y_max" => get(c, "y_max_pos", deg2rad(25.0)))
    theirs = robot_gains()
    if theirs === nothing
        println("The robot is not reachable, so the gains cannot be checked.")
        println("measured_parameters.toml says:")
        for (k, v) in sort(collect(mine)); @printf("    %-24s %.4g\n", k, v); end
        println("MAKE SURE these are the gains the run used. If they are not, everything below is wrong.")
        return true
    end
    bad = String[]
    println(" ", rpad("gain", 24), rpad("config", 12), rpad("robot", 12), "")
    for (k, v) in sort(collect(mine))
        r = get(theirs, k, nothing)
        mark = r === nothing ? "  (robot uses the model default)" :
               isapprox(v, r; rtol = 1e-6) ? "" : "   <-- DISAGREE"
        isempty(mark) || r === nothing || push!(bad, k)
        @printf(" %-24s %-12.4g %-12s%s\n", k, v, r === nothing ? "-" : string(round(r, sigdigits = 5)), mark)
    end
    if !isempty(bad)
        println("\n", length(bad), " gain(s) disagree: ", join(bad, ", "))
        println("The robot is the source of truth. Edit [controller] in measured_parameters.toml to match,")
        println("then rerun. Do NOT fit or sweep until they agree.")
        return false
    end
    println("\ngains agree.")
    return true
end

path = isempty(ARGS) ? joinpath(DATA, "run_amp10.arrow") : ARGS[1]
refit = "--refit" in ARGS
println("run: ", basename(path), refit ? "   (will refit, about 49 minutes)" : "")
println(repeat("-", 70))
flush(stdout)
check_gains() || exit(1)

if refit
    println("\n=== fitting the plant, about 49 minutes ===")
    flush(stdout)
    run(`julia +1.12 --project=$(HERE) $(joinpath(HERE, "calibrate_actuator.jl")) $path 20 185`)
end

println("\n=== sweeping the outer loop, about 3 minutes ===")
flush(stdout)
run(`julia +1.12 --project=$(HERE) $(joinpath(HERE, "sweep_outer.jl")) $path`)
