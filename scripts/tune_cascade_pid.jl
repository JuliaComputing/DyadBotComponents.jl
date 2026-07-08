# PID autotuning of the outer position loop of the cascade-controlled
# balancing robot, followed by robustness analysis of the full cascade.
#
# The inner angle loop is kept at its current gains while the outer loop is
# opened at the analysis points `u2` (position-controller output) and `y2`
# (position measurement) of CascadeControlledDyadBot.
#
# The resulting parameters (displayed as `cascade_optimized_params`) can be
# written back into the defaults of `k_pos`, `Ti_pos`, `Td_pos` (and `Nd`) in
# dyad/closed_loop.dyad.

cd(@__DIR__)
using Pkg
Pkg.activate(".")
using DyadBotComponents
using MultibodyComponents
using ModelingToolkit
using OrdinaryDiffEq
using DyadControlSystems, ControlSystemsBase, ControlSystemsMTK
using LinearAlgebra
using Plots
import DyadControlSystems as JSC

## Simulate the closed loop with the current gains
@named model = DyadBotComponents.CascadeControlledDyadBot()
ssys = multibody(model)
prob = ODEProblem(ssys, [], (0.0, 20.0))
sol = solve(prob, Rodas5P())
plot(sol, idxs = [ssys.plant.theta, ssys.plant.x, ssys.square.y, ssys.pos_controller.y])

# Operating point for the linearizations: the upright equilibrium with zero
# reference. Opening a loop leaves the analysis-point inputs free, so they
# are given values here, together with the trim values of the mechanical
# states.
using ModelingToolkit: Symbolics
overrides = Dict{Symbolics.SymbolicT, Symbolics.SymbolicT}()
for (k, v) in [
    ssys.square.amplitude => 0.0
    ssys.angle_controller.u_m => 0.0
    ssys.angle_controller.integrator.y => 0.0
    ssys.angle_controller.derivative.x => 0.0
    ssys.pos_controller.u_m => 0.0
    ssys.pos_controller.integrator.y => 0.0
    ssys.pos_controller.derivative.x => 0.0
    ssys.gain.u => 0.0
    ssys.gain1.u => 0.0
    ssys.plant.torque => 0.0
    ssys.plant.body_mass.body.phi => 0.0
    ssys.plant.body_mass.body.w => 0.0
    ssys.plant.wheelinertia.phi => 0.0
    ssys.plant.wheelinertia.w => 0.0
    ssys.firstorder.x => 0.0
    ssys.firstorder1.x => 0.0
]
    overrides[Symbolics.unwrap(k)] = v
end

## Autotuning of the outer position loop
# Model instance for linearization and tuning: zero initial tilt, so that
# the model's initial equations are consistent with the upright trim.
@named cascade_tuning_model = DyadBotComponents.CascadeControlledDyadBot(phi0 = 0)

cascade_spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :CascadePositionTuning,
    model = cascade_tuning_model,
    measurement = "y2",
    control_input = "u2",
    step_input = "u2",
    step_output = "y2",
    ref = 0.0,
    Ts = 0.01,
    duration = 15.0,
    Ms = 1.4,
    Mt = 1.9,
    Mks = 1000.0,
    wl = 1e-2,
    wu = 1e3,
    num_frequencies = 300,
    soft = false,
    exact_hessian = false,
    scale = true,
    auto_resolve = true,
    homotopy = true,
    overrides,
)

cascade_asol = JSC.run_analysis(cascade_spec)

plot(cascade_asol.sol)

cascade_Splot = JSC.artifacts(cascade_asol, :SensitivityFunctions)
cascade_nyquist_plot = JSC.artifacts(cascade_asol, :NyquistPlot)
cascade_optimized_params = JSC.artifacts(cascade_asol, :OptimizedParameters)

display(cascade_optimized_params)

## Robustness analysis of the full cascade at the current gains
# MIMO output sensitivity at both measurements
S2 = get_named_sensitivity(cascade_tuning_model,
    [cascade_tuning_model.y, cascade_tuning_model.y2]; op = overrides, MultibodyComponents.linsys...)
Ms2, ws2 = hinfnorm2(S2)
sigmaplot(S2)
hline!([Ms2], l = (:dash, :black), label = "\$M_S = \$$(round(Ms2, digits = 2))")

L2 = inv(S2) - ss(I(2))
bodeplot(L2)

# Input sensitivity at the torque input
Si2 = get_named_sensitivity(cascade_tuning_model, [cascade_tuning_model.u]; op = overrides, MultibodyComponents.linsys...)
w = exp10.(LinRange(-1, 3, 1000))
Msi2, _ = hinfnorm2(Si2)
bodeplot(Si2, w, plotphase = false)
hline!([Msi2], l = (:dash, :black), label = "\$M_S = \$$(round(Msi2, digits = 2))")

# Input loop transfer and diskmargin
Li2 = get_named_looptransfer(cascade_tuning_model, [cascade_tuning_model.u]; op = overrides, MultibodyComponents.linsys...) |> minreal
marginplot(Li2, w, adjust_phase_start = true)

dmi2 = diskmargin(Li2)
plot(dmi2)
