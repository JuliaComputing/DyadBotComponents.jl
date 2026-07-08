# PID autotuning of the inner angle loop of the balancing robot.
#
# The loop from motor torque to body tilt angle is opened at the analysis
# points `u` (controller output) and `y` (angle measurement) of
# AngleControlledDyadBot, and a PID controller is optimized subject to
# constraints on the sensitivity functions.
#
# The resulting parameters (displayed as `optimized_params`) can be written
# back into the defaults of `k_angle`, `Ti_angle`, `Td_angle` (and `Nd`) in
# dyad/closed_loop.dyad.

cd(@__DIR__)
using Pkg
Pkg.activate(".")
using DyadBotComponents
using MultibodyComponents
using ModelingToolkit
using OrdinaryDiffEq
using DyadControlSystems, ControlSystemsBase, ControlSystemsMTK
using Plots
import DyadControlSystems as JSC

## Simulate the closed loop with the current gains
@named model = DyadBotComponents.AngleControlledDyadBot(phi0 = deg2rad(-10))
ssys = multibody(model)
prob = ODEProblem(ssys, [], (0.0, 1.0))
sol = solve(prob, Rodas5P())
plot(sol, idxs = [ssys.plant.theta, ssys.plant.x, ssys.plant.torque])
hline!([0], l = (:dash, :black), primary = false)

# Model instance for linearization and tuning: zero initial tilt, so that
# the model's initial equations are consistent with the upright trim.
@named tuning_model = DyadBotComponents.AngleControlledDyadBot(phi0 = 0)

# Operating point for the linearizations: the upright equilibrium. Opening a
# loop leaves the analysis-point inputs free, so they are given values here,
# together with the trim values of the mechanical states.
using ModelingToolkit: Symbolics
overrides = Dict{Symbolics.SymbolicT, Symbolics.SymbolicT}()
for (k, v) in [
    ssys.angle_controller.u_m => 0.0
    ssys.angle_controller.integrator.y => 0.0
    ssys.angle_controller.derivative.x => 0.0
    ssys.gain.u => 0.0
    ssys.plant.torque => 0.0
    ssys.plant.body_mass.body.phi => 0.0
    ssys.plant.body_mass.body.w => 0.0
    ssys.plant.wheelinertia.phi => 0.0
    ssys.plant.wheelinertia.w => 0.0
]
    overrides[Symbolics.unwrap(k)] = v
end

## Loop-shaping quantities at the current gains
S = get_named_sensitivity(tuning_model, tuning_model.y; op = overrides, MultibodyComponents.linsys...) |> sminreal
Ms, ws = hinfnorm2(S)
bodeplot(S, title = "\$S(s)\$ angle controlled", plotphase = false, legend = :bottomright)
hline!([Ms], l = (:dash, :black), label = "\$M_S = \$$(round(Ms, digits = 2))")

## PID autotuning
spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :AngleTuning,
    model = tuning_model,
    measurement = "y",
    control_input = "u",
    step_input = "u",
    step_output = "y",
    Ts = 0.01,           # Sample time
    duration = 3.0,      # Simulation duration
    Ms = 1.6,            # Sensitivity peak constraint
    Mt = 1.6,            # Complementary sensitivity peak constraint
    Mks = 2000.0,        # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    num_frequencies = 200,
    soft = true,
    homotopy = true,
    overrides,
)

asol = JSC.run_analysis(spec)

plot(asol.sol)

# Visualize results
Splot = JSC.artifacts(asol, :SensitivityFunctions)
response_plot = JSC.artifacts(asol, :OptimizedResponse)
nyquist_plot = JSC.artifacts(asol, :NyquistPlot)
optimized_params = JSC.artifacts(asol, :OptimizedParameters)

display(optimized_params)
