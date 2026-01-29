# cd(joinpath(@__DIR__, ".."))
# using Pkg
# Pkg.activate(".")
include("planar_multibody.jl")


using OrdinaryDiffEq
using DyadControlSystems, ControlSystemsBase, ControlSystemsMTK
using Plots
import DyadControlSystems as JSC
using LinearAlgebra
connect = ModelingToolkit.connect

# @named plant = FlatDyadBot()
# plant = complete(plant)
# inputs = [plant.control_input.u]
# outputs = [plant.theta_output.u]


# P0 = named_ss(plant, inputs, outputs; op = Dict([inputs .=> 0; plant.b_trans=>0; plant.b_rot=>0; plant.x => big(0.0)]), allow_input_derivatives=true)

# Pt = named_ss(plant, inputs, outputs; op = Dict([inputs .=> 0; plant.b_trans=>1; plant.b_rot=>0; plant.x => big(0.0)]), allow_input_derivatives=true)

# Pr = named_ss(plant, inputs, outputs; op = Dict([inputs .=> 0; plant.b_trans=>0; plant.b_rot=>1; plant.x => big(0.0)]), allow_input_derivatives=true)

# Ptr = named_ss(plant, inputs, outputs; op = Dict([inputs .=> 0; plant.b_trans=>1; plant.b_rot=>1; plant.x => big(0.0)]), allow_input_derivatives=true)

# Ps = [P0, Pt, Pr, Ptr]

# nyquistplot(8 .* Ps)


#=
The loop closed around the angle only has a zero in the origin corresponding to a constant wheel rotation and non-zero input exactly matching the friction torque. This zero attracts the RHP pole and the RHP pole will never move into the LHP unless we move the zero away from the origin by adding integral action. If the friction is zero, the zero cancels with a pole in the origin.
=#

@component function AngleControlledPlanarSegway(; name)
    pars = @parameters begin
        theta_ref = deg2rad(0)  # Reference angle (upright)
    end

    systems = @named begin
        plant = PlanarMultibodybot()
        controller = Blocks.LimPID(k=0.487401, Ti=0.0587352, Td=0.0420526, Nd=119.368)#, u_max=15)
        ref = Blocks.Constant(k=theta_ref)
    end

    eqs = [
        connect(ref.output, :r, controller.reference)
        connect(plant.theta_output, :y, controller.measurement)
        connect(controller.ctr_output, :u, plant.control_input)
    ]

    System(eqs, t, [], pars; systems, name)
end

@named segway_model = AngleControlledPlanarSegway()
ssys = multibody(segway_model)

x0 = [
    ssys.plant.body.phi => deg2rad(-10)  # Adjust for offset 
    ssys.plant.body.w => 0.0  # Adjust for offset 
]

prob = ODEProblem(ssys, x0, (0.0, 1.0))
sol = solve(prob, Rodas5P())
plot(sol, idxs=[ssys.plant.theta_output.u, ssys.plant.x_output.u, ssys.plant.tau]); hline!([0], l=(:dash, :black), primary=false)

import GLMakie
Multibody.render(segway_model, sol, 0.0, lookat=[0,0.1,0], x=0, y=0.1, z=-0.5)[1]


##

L = -get_named_looptransfer(segway_model, segway_model.y; Multibody.linsys...) |> sminreal
PS = named_ss(segway_model, segway_model.u, segway_model.y; loop_openings = [], Multibody.linsys...) |> sminreal
lsys = named_ss(segway_model, segway_model.u, segway_model.y; loop_openings = [segway_model.y, segway_model.u], Multibody.linsys...) |> sminreal

C = named_ss(segway_model, segway_model.r, segway_model.u; loop_openings = [segway_model.y], Multibody.linsys...) |> sminreal



"""
    zero_direction(lsys, z)

Compute the zero direction of the linear system `lsys` at the zero `z`. Returns the nullspace of the Rosenbrock matrix at `z`, where the first `nx` entries correspond to the state direction and the last `nu` entries correspond to the input direction. If `z` is not a zero, the nullspace will be empty.
"""
function zero_direction(lsys, z)
    (; A,B,C,D) = lsys
    R = [z*I-A -B; C D]
    nullspace(R)
end


S = get_named_sensitivity(segway_model, segway_model.y; Multibody.linsys...) |> sminreal
Ms, ws = hinfnorm2(S)
bodeplot(S, title="\$S(s)\$ angle controlled", plotphase=false, legend=:bottomright)
hline!([Ms], l=(:dash, :black), label="\$M_S = \$$(round(Ms, digits=2))")
##
# PID Autotuning Analysis
using DyadControlSystems
import DyadControlSystems as JSC

@named tuning_model = AngleControlledPlanarSegway()

spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :SegwayTuning,
    model = tuning_model,
    measurement = "y",
    control_input = "u",
    step_input = "u",
    step_output = "y",
    Ts = 0.01,           # Sample time
    duration = 3.0,      # Simulation duration
    Ms = 1.6,            # Sensitivity peak constraint
    Mt = 1.6,            # Complementary sensitivity peak constraint
    Mks = 2000.0,         # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    # ki_ub = 0.0,         # Tune PD controller
    num_frequencies = 200,
    soft = true,
    homotopy = true,
    # auto_resolve = false,
)

# Run the autotuning analysis
asol = JSC.run_analysis(spec)

plot(asol.sol)

# Visualize results
Splot = JSC.artifacts(asol, :SensitivityFunctions)
response_plot = JSC.artifacts(asol, :OptimizedResponse)
nyquist_plot = JSC.artifacts(asol, :NyquistPlot)
optimized_params = JSC.artifacts(asol, :OptimizedParameters)

display(optimized_params)

# ==============================================================================
##
# ==============================================================================

# Cascade control: outer velocity loop + inner angle loop
@component function CascadeControlledPlanarSegway(; name)
    pars = @parameters begin
        x_ref = 0.15  # Reference velocity
    end

    systems = @named begin
        plant = PlanarMultibodybot()
        # Inner loop: angle controller
        inner_controller = Blocks.LimPID(k=0.487401, Ti=0.0587352, Td=0.0420526, Nd=119.368)
        # Outer loop: velocity controller
        outer_controller = Blocks.LimPID(k=0.0666576, Ti=5.25024, Td=4.81393, Nd=4.76616, wd=1, wp=1, u_max=deg2rad(25.0))
        neg_gain = Blocks.Gain(k=1)
        ref = Blocks.Step(height=x_ref, start_time=10)
        # Add pi offset to inner loop reference
        pi_offset = Blocks.Constant(k=0)
        add_pi = Blocks.Add(k1=1, k2=1)
        ref_filter = Blocks.FirstOrder(T=1)
    end

    eqs = [
        # Outer loop: velocity reference -> angle reference
        connect(ref.output, ref_filter.input)
        connect(ref_filter.output, :r2, outer_controller.reference)
        connect(plant.x_output, neg_gain.input)
        connect(neg_gain.output, :y2, outer_controller.measurement)

        # Add pi to outer controller output for inner loop reference
        connect(outer_controller.ctr_output, :u2, add_pi.input1)
        connect(pi_offset.output, add_pi.input2)

        # Inner loop: angle reference -> torque
        connect(add_pi.output, inner_controller.reference)
        connect(plant.theta_output, :y, inner_controller.measurement)
        connect(inner_controller.ctr_output, :u, plant.control_input)
    ]

    System(eqs, t, [], pars; systems, name)
end

@named cascade_model = CascadeControlledPlanarSegway()
cascade_ssys = multibody(cascade_model)

x0 = [
    cascade_ssys.plant.body.phi => deg2rad(5)
    # cascade_ssys.plant.wheelJoint.color => [0,0,0,0.5]
    cascade_ssys.plant.wheelJoint.frame_a.render => true
    cascade_ssys.plant.wheelJoint.frame_a.length => 0.05
    cascade_ssys.plant.wheelJoint.frame_a.radius => 0.004
]

cascade_prob = ODEProblem(cascade_ssys, x0, (0.0, 20.0))
cascade_sol = solve(cascade_prob, Rodas5P(), dt=0.005, adaptive=false)
plot(cascade_sol, idxs=[cascade_ssys.plant.theta_output.u, cascade_ssys.plant.x_output.u, cascade_ssys.outer_controller.ctr_output.u]); hline!([0 0.15], l=(:dash, :black), primary=false, ylims=(-1, 3.5), size=(800,1600))

Multibody.render(cascade_model, cascade_sol, 0.0, lookat=[0,0.2,0], x=0, y=0.2, z=-0.5)[1]
# Multibody.render(cascade_model, cascade_sol, lookat=[0,0.2,0], x=0, y=0.2, z=-0.5, timescale=1)
# plot(cascade_sol, idxs=[cascade_ssys.pi_offset.output.u])

##

lsys = named_ss(cascade_model, cascade_model.u2, cascade_model.y2; loop_openings = [cascade_model.u2, cascade_model.y2], Multibody.linsys...) |> sminreal



##
# PID Autotuning for outer velocity loop
@named cascade_tuning_model = CascadeControlledPlanarSegway()

cascade_spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :CascadeVelocityTuning,
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
    # kd_ub = 0.0, # Tune PI controller
    # kp_guess = 0.55,
    # ki_guess = 0.55,
    # kd_guess = 1e-3,
    soft = false,
    exact_hessian = false,
    scale = true,
    auto_resolve = true,
    homotopy = true,
)

cascade_asol = JSC.run_analysis(cascade_spec)

plot(cascade_asol.sol)

cascade_Splot = JSC.artifacts(cascade_asol, :SensitivityFunctions)
# cascade_response_plot = JSC.artifacts(cascade_asol, :OptimizedResponse)
cascade_nyquist_plot = JSC.artifacts(cascade_asol, :NyquistPlot)
cascade_optimized_params = JSC.artifacts(cascade_asol, :OptimizedParameters)

display(cascade_optimized_params)

##

S2 = get_named_sensitivity(cascade_tuning_model, [cascade_tuning_model.y, cascade_tuning_model.y2]; Multibody.linsys...)
Ms2, ws2 = hinfnorm2(S2)
sigmaplot(S2); hline!([Ms2], l=(:dash, :black), label="\$M_S = \$$(round(Ms2, digits=2))")
# bodeplot(S2, plotphase=false)

L2 = inv(S2) - ss(I(2))
bodeplot(L2)

## Swms = freqresp(S2, ws2)

Si2 = get_named_sensitivity(cascade_tuning_model, [cascade_tuning_model.u]; Multibody.linsys...)
w = exp10.(LinRange(-1, 3, 1000))
Msi2, ws2 = hinfnorm2(Si2)
bodeplot(Si2, w, plotphase=false); hline!([Msi2], l=(:dash, :black), label="\$M_S = \$$(round(Msi2, digits=2))")

Li2 = minreal(inv(Si2) - 1, 1e-12)
Li2 = convert(StateSpace{Continuous, Float64}, Li2)
nyquistplot(Li2)
marginplot(Li2, w, adjust_phase_start=true)


Li22 = get_named_looptransfer(cascade_tuning_model, [cascade_tuning_model.u]; Multibody.linsys...) |> minreal

dmi2 = diskmargin(Li2)
plot(dmi2)

dmi22 = diskmargin(Li22)
plot(dmi22)




# ==============================================================================
## Feedforward generation
# ==============================================================================

T0 = ss(tf([1], [0.1, 1]))
Tr = T0^4

Ryur = DyadControlSystems.feedforward_generator(cascade_tuning_model; # requires master branch
    Tr,
    measurement = [cascade_tuning_model.y, cascade_tuning_model.y2],
    controlled_output = [cascade_tuning_model.y2],
    control_input = [cascade_tuning_model.u],
    Multibody.linsys...
)

Ryur = balance_statespace(Ryur)[1]

@component function FilteredCascadeControlledFlatDyadBot(; name)
    pars = @parameters begin
        x_ref = 0.15  # Reference velocity
    end

    systems = @named begin
        plant = PlanarMultibodybot()
        # Inner loop: angle controller
        inner_controller = Blocks.LimPID(k=22.3659, Ti=false, Td=0.0632787, Nd=54.6269)
        # Outer loop: velocity controller
        outer_controller = Blocks.LimPID(k=0.0666576, Ti=5.25024, Td=4.81393, Nd=4.76616, wd=1, wp=1, u_max=deg2rad(25.0))
        neg_gain = Blocks.Gain(k=1)
        ref = Blocks.Step(height=x_ref, start_time=5)
        # Add pi offset to inner loop reference
        pi_offset = Blocks.Constant(k=0)
        add_pi = Blocks.Add3(k1=1, k2=1)

        torque_input = Blocks.Add()
        R = Blocks.StateSpace(ssdata(Ryur)...)
    end

    eqs = [


        connect(ref.output, :r2, R.input)
        R.output.u[1] ~ add_pi.input3.u
        R.output.u[2] ~ outer_controller.reference.u
        R.output.u[3] ~ torque_input.input2.u


        # Outer loop: velocity reference -> angle reference
        connect(plant.x_output, neg_gain.input)
        connect(neg_gain.output, :y2, outer_controller.measurement)

        # Add pi to outer controller output for inner loop reference
        connect(outer_controller.ctr_output, :u2, add_pi.input1)
        connect(pi_offset.output, add_pi.input2)

        # Inner loop: angle reference -> torque
        connect(add_pi.output, inner_controller.reference)
        connect(plant.theta_output, :y, inner_controller.measurement)
        connect(inner_controller.ctr_output, :u, torque_input.input1)
        connect(torque_input.output, plant.control_input)
    ]

    System(eqs, t, [], pars; systems, name)
end

@named filtered_model = FilteredCascadeControlledFlatDyadBot()
filtered_ssys = multibody(filtered_model)

Pcl_angle = named_ss(filtered_model, [filtered_model.r2], [filtered_model.y]; allow_input_derivatives=true, Multibody.linsys...)
Pcl_pos = named_ss(filtered_model, [filtered_model.r2], [filtered_model.y2]; allow_input_derivatives=true, Multibody.linsys...)

Pcl_angle = minreal(Pcl_angle)
Pcl_pos = minreal(Pcl_pos)

w = exp10.(LinRange(-2.5, 1.5, 2000))
bodeplot([Pcl_angle, Pcl_pos, Tr], w, legend=:bottom, plotphase=false, legendfontsize=8, hz=true, label=["Gar" "Gpr" "Tr"])


x0 = [
    filtered_ssys.inner_controller.k  => optimized_params[1, :Kp_standard];
    filtered_ssys.inner_controller.Ti => false;#optimized_params[1, :Ti_standard];
    filtered_ssys.inner_controller.Td => optimized_params[1, :Td_standard];
    filtered_ssys.inner_controller.Nd => optimized_params[1, :Nd];
    filtered_ssys.outer_controller.k  => cascade_optimized_params[1, :Kp_standard];
    filtered_ssys.outer_controller.Ti => cascade_optimized_params[1, :Ti_standard];

    filtered_ssys.plant.body.phi => deg2rad(-10)
]


filtered_prob = ODEProblem(filtered_ssys, x0, (0.0, 10.0))
filtered_sol = solve(filtered_prob, Rodas5P(), dt=0.01, dtmin=0.01, adaptive=false)
plot(filtered_sol, idxs=[filtered_ssys.plant.theta_output.u, filtered_ssys.plant.x_output.u, filtered_ssys.outer_controller.ctr_output.u, filtered_ssys.plant.tau]); hline!([0 0.15], l=(:dash, :black), primary=false, ylims=(-1, 3.5), size=(800,1600), legend=:right)

plot(filtered_sol, idxs=[filtered_ssys.R.input.u; filtered_ssys.R.output.u; ])

# Multibody.render(filtered_model, filtered_sol, lookat=[0,0.2,0], x=0, y=0.2, z=-0.5, timescale=1)