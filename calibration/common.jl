# Shared setup for the plant calibration: the measured parameters, the robot's
# controller gains, and the model builder. Included by make_release_data.jl and
# calibrate_plant.jl so that the data and the fit always use the same model.

using DyadBotComponents, ModelingToolkit, OrdinaryDiffEq, TOML

const HERE = @__DIR__
const DATA = joinpath(HERE, "data")
const PLOTS = joinpath(HERE, "plots")

const CFG = TOML.parsefile(joinpath(HERE, "measured_parameters.toml"))
const MEAS = (M = CFG["M"], m = CFG["m"], R = CFG["R"], L = CFG["L"])

# Controller gains, from the [controller] table. They are not fitted: the robot
# runs a controller whose gains are known, and the plant is what is unknown.
const GAINS = let c = CFG["controller"]
    (; controller__k_angle = c["k_angle"], controller__Ti_angle = c["Ti_angle"],
       controller__Td_angle = c["Td_angle"], controller__k_pos = c["k_pos"],
       controller__Ti_pos = c["Ti_pos"], controller__Td_pos = c["Td_pos"],
       controller__angle_controller__Nd = c["Nd_angle"],
       controller__angle_controller__wd = get(c, "wd_angle", 1.0),
       controller__pos_controller__Nd = c["Nd_pos"],
       # The model defaults to a 25 degree limit that the robot has never used.
       # Leaving it out let the twin command 3 times the lean the robot can, and
       # the sweep on 2026-09-18 diverged at the robot's own gains because of it.
       controller__pos_controller__y_max = get(c, "y_max_pos", deg2rad(25.0)),
       controller__pos_controller__y_min = -get(c, "y_max_pos", deg2rad(25.0)))
end

# Starting guesses for the three fitted parameters, from the measured geometry:
# the body as a uniform rod of length 2L about its centre, and the wheels as a
# solid disc. Both are underestimates of a real robot, which is the point of
# fitting them.
Ib_estimate() = MEAS.M * (2 * MEAS.L)^2 / 12
Iw_estimate() = MEAS.m * MEAS.R^2 / 2
const NOMINAL = (Ib = Ib_estimate(), Iw = Iw_estimate(), d = 1.0e-4)

"""
    build(; phi0, x0 = 0.0, Ib, Iw, d, discrete = false) -> compiled System

`BalanceDyadBot` with the measured parameters, the robot's gains, and the three
fitted parameters applied. `x0` is the wheel position at the release, measured
from the point the controller holds as zero: a robot log keeps the runtime's
zero, and the outer loop pulls toward it, not toward the release point.
The release angle `phi0` is baked in at construction:
the Dyad model sets the body's initial tilt through an initialization equation,
which wins over anything passed to the `ODEProblem`, so each release needs its
own model. Every instance is named `model`, so the same symbolic parameter names
address all of them and one search space covers every release.
"""
function build(; phi0, x0 = 0.0, Ib = NOMINAL.Ib, Iw = NOMINAL.Iw, d = NOMINAL.d, discrete = false)
    ctor = discrete ? DyadBotComponents.DiscreteBalanceDyadBot : DyadBotComponents.BalanceDyadBot
    # In PlanarDyadBot the wheel-axis position is x = -R * wheel angle.
    sys = ctor(; name = :model, phi0,
               plant__M = MEAS.M, plant__m = MEAS.m, plant__R = MEAS.R, plant__L = MEAS.L,
               plant__Ib = Ib, plant__Iw = Iw, plant__d = d,
               plant__wheelinertia__phi__initial = -x0 / MEAS.R, GAINS...)
    return mtkcompile(sys)
end

# The leaf parameters that the Dyad-level Ib, Iw and d are bound to. `mtkcompile`
# keeps only these; `plant.Ib` and friends do not survive it, so the search space
# and any run-time override must use these names.
leaf_params(sys) = (Ib = sys.plant.body_mass.I, Iw = sys.plant.wheelinertia.I,
                    d = sys.plant.simplemotor.damper.d)

"Solve one release to `t_end`."
function run_release(sys, t_end; abstol = 1e-10, reltol = 1e-10, overrides = [])
    return solve(ODEProblem(sys, overrides, (0.0, t_end)), Rodas5P(); abstol, reltol)
end
