# LQG controller design for the balancing robot.
#
# An LQG controller with integral action on the position and reference
# feedforward is synthesized from a linearization of LQGTuningDyadBot with
# the loop opened at the analysis point `u`. The combined
# feedforward-plus-feedback controller has 6 inputs
# [r_x, r_xd, r_theta, r_thetad, x, theta] and one output (motor torque).
# Its state-space matrices are written to CSV files in ../data/, from where
# they are loaded by the component LQGControlledDyadBot.
#
# After running this script, recompile the Dyad package (`dyad compile .` in
# the package root) if the matrix dimensions changed, and simulate
# LQGControlledDyadBot to verify the design (final section).

cd(@__DIR__)
using Pkg
Pkg.activate(".")
using DyadBotComponents
using MultibodyComponents
using ModelingToolkit
using OrdinaryDiffEq
using DyadControlSystems, ControlSystemsBase, ControlSystemsMTK
using DelimitedFiles
using Plots
import DyadControlSystems as JSC

@named tuning_model = DyadBotComponents.LQGTuningDyadBot()
tssys = multibody(tuning_model)

# Operating point for the linearization: the upright equilibrium. Opening the
# loop at `u` leaves the analysis-point inputs free, so they are given values
# here, together with the trim values of the mechanical states (initial
# conditions are codegened to initial equations by dyad and those are ignored
# by linearization).
using ModelingToolkit: Symbolics
overrides = Dict{Symbolics.SymbolicT, Symbolics.SymbolicT}()
for (k, v) in [
    tssys.g_x.u => 0.0
    tssys.g_xd.u => 0.0
    tssys.g_theta.u => 0.0
    tssys.g_thetad.u => 0.0
    tssys.plant.torque => 0.0
    tssys.plant.body_mass.body.phi => 0.0
    tssys.plant.body_mass.body.w => 0.0
    tssys.plant.wheelinertia.phi => 0.0
    tssys.plant.wheelinertia.w => 0.0
]
    overrides[Symbolics.unwrap(k)] = v
end

# LQG Analysis Specification
# - 2 measurements: position (x) and tilt angle (theta)
# - 4 controlled outputs: position, velocity, angle and angular rate
# - 1 control input: motor torque
lqg_spec = JSC.LQGAnalysisSpec(;
    name = :DyadBotLQG,
    model = tuning_model,
    measurement = ["y_x", "y_theta"],
    controlled_output = ["y_x", "y_xd", "y_theta", "y_thetad"],
    control_input = ["u"],
    loop_openings = ["u"],
    q1_diag = [10.0, 0.1, 1, 0.1],  # Penalty on controlled outputs
    q2_diag = [0.0001],             # Penalty on control input
    r1_diag = [1.0],                # Dynamics noise covariance
    r2_diag = [0.001, 0.001],       # Measurement noise covariance (x, theta)
    wl = 1e-2,
    wu = 314,
    num_frequencies = 200,
    integrator_indices = [1],       # Integral action on the position
    integrator_r1_diag = [0.1],
    overrides,
)

lqg_asol = JSC.run_analysis(lqg_spec)

# Visualize results
step_response = JSC.artifacts(lqg_asol, :StepResponse)
gang_of_four = JSC.artifacts(lqg_asol, :GangOfFour)
bode_plot = JSC.artifacts(lqg_asol, :BodePlot)
margin_plot = JSC.artifacts(lqg_asol, :MarginPlot)
controller_gain = lqg_asol.L
observer_gain = lqg_asol.K

## Extract the state-space matrices of the combined feedforward + feedback
# controller `Cfffb = [Cff -Cfb]` and store them for LQGControlledDyadBot.
# The controller acts around the operating point (u0, y0); for the upright
# equilibrium both are zero so no compensation is required.
@show lqg_asol.u0 lqg_asol.y0
maximum(abs, [lqg_asol.u0; lqg_asol.y0]) < 1e-6 ||
    @warn "Nonzero operating point; LQGControlledDyadBot does not compensate for it"
A, B, C, D = ssdata(ss(lqg_asol.Cfffb))
@assert size(B, 2) == 6 "expected 6 controller inputs (4 references + 2 measurements)"
@assert size(C, 1) == 1 "expected 1 controller output (torque)"

dir = joinpath(@__DIR__, "..", "data")
mkpath(dir)
writedlm(joinpath(dir, "lqg_A.csv"), A, ',')
writedlm(joinpath(dir, "lqg_B.csv"), B, ',')
writedlm(joinpath(dir, "lqg_C.csv"), C, ',')
writedlm(joinpath(dir, "lqg_D.csv"), D, ',')
open(joinpath(dir, "lqg_nx.txt"), "w") do io
    println(io, size(A, 1))
end

## Fundamental limitations and robustness
rhp_pole, _ = findmax(real, poles(lqg_asol.P))
# The gain crossover must exceed twice the RHP pole

Li = lqg_asol.Cfb * system_mapping(lqg_asol.P_ext)
dmi = diskmargin(Li, offset = 0)
plot(dmi)
marginplot(Li, adjust_phase_start = false)
vline!([2 * rhp_pole], l = (:dash, :black), label = "Fundamental limitation")

## Verify the closed loop (requires the Dyad package to be recompiled after
# the CSV files were written, so that LQGControlledDyadBot loads the new
# controller)
@named cl_model = DyadBotComponents.LQGControlledDyadBot(phi0 = deg2rad(-10))
ssys = multibody(cl_model)
prob = ODEProblem(ssys, [], (0.0, 20.0))
sol = solve(prob, Rodas5P())
@show sol.retcode
plot(sol, idxs = [ssys.plant.theta, ssys.plant.x, ssys.plant.torque])
hline!([0 prob.ps[ssys.step_height]], l = (:dash, :black), primary = false)
