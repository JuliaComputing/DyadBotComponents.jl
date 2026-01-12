# cd(joinpath(@__DIR__, ".."))
# using Pkg
# Pkg.activate(".")
include("planar_flat.jl")  # Get FlatDyadBot


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

@component function AngleControlledFlatDyadBot(; name)
    pars = @parameters begin
        theta_ref = deg2rad(180)  # Reference angle (upright)
    end

    systems = @named begin
        plant = FlatDyadBot()
        controller = Blocks.LimPID(k=15.6, Ti=Inf, Td=0.16, Nd=25, u_max=7)
        ref = Blocks.Constant(k=theta_ref)
    end

    eqs = [
        connect(ref.output, :r, controller.reference)
        connect(plant.theta_output, :y, controller.measurement)
        connect(controller.ctr_output, :u, plant.control_input)
    ]

    System(eqs, t, [], pars; systems, name)
end

@named flat_model = AngleControlledFlatDyadBot()
model = complete(flat_model)
ssys = structural_simplify(model)

x0 = [
    ssys.plant.theta => deg2rad(160)
]

prob = ODEProblem(ssys, x0, (0.0, 5.0))
sol = solve(prob, Rodas5P())
plot(sol, idxs=[ssys.plant.theta, ssys.plant.x, ssys.plant.tau]); hline!([π], l=(:dash, :black), primary=false)

##
S = get_named_sensitivity(flat_model, flat_model.y)
Ms, ws = hinfnorm2(S)
bodeplot(S, title="\$S(s)\$ angle controlled", plotphase=false, legend=:bottomright)
hline!([Ms], l=(:dash, :black), label="\$M_S = \$$(round(Ms, digits=2))")
##
# PID Autotuning Analysis
using DyadControlSystems
import DyadControlSystems as JSC

@named tuning_model = AngleControlledFlatDyadBot()

spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :SegwayTuning,
    model = tuning_model,
    measurement = "y",
    control_input = "u",
    step_input = "u",
    step_output = "y",
    Ts = 0.01,           # Sample time
    duration = 25.0,      # Simulation duration
    Ms = 1.5,            # Sensitivity peak constraint
    Mt = 1.5,            # Complementary sensitivity peak constraint
    Mks = 400.0,         # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    ki_ub = 0.0,         # Tune PD controller
    num_frequencies = 200,
    soft = true,
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

##
# Cascade control: outer velocity loop + inner angle loop
@component function CascadeControlledFlatDyadBot(; name)
    pars = @parameters begin
        x_ref = 0.15  # Reference velocity
    end

    systems = @named begin
        plant = FlatDyadBot()
        # Inner loop: angle controller
        inner_controller = Blocks.LimPID(k=15.6, Ti=Inf, Td=0.16, Nd=25, u_max=7)
        # Outer loop: velocity controller
        outer_controller = Blocks.LimPID(k=0.54, Ti=2.48, Td=0, Nd=600, wd=1, wp=0.5)
        neg_gain = Blocks.Gain(k=1)
        ref = Blocks.Step(height=x_ref, start_time=5)
        # Add pi offset to inner loop reference
        pi_offset = Blocks.Constant(k=pi)
        add_pi = Blocks.Add(k1=1, k2=1)
    end

    eqs = [
        # Outer loop: velocity reference -> angle reference
        connect(ref.output, :r2, outer_controller.reference)
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

@named cascade_model = CascadeControlledFlatDyadBot()
cascade_model = complete(cascade_model)
cascade_ssys = structural_simplify(cascade_model)

x0 = [
    cascade_ssys.plant.theta => deg2rad(170)
]

cascade_prob = ODEProblem(cascade_ssys, x0, (0.0, 20.0), dtmax=0.01)
cascade_sol = solve(cascade_prob, Rodas5P())
plot(cascade_sol, idxs=[cascade_ssys.plant.theta, cascade_ssys.plant.x_dot, cascade_ssys.plant.x, cascade_ssys.outer_controller.ctr_output.u, cascade_ssys.plant.tau]); hline!([π 0.15], l=(:dash, :black), primary=false, ylims=(-1, 3.5), size=(800,1600))

# plot(cascade_sol, idxs=[cascade_ssys.pi_offset.output.u])

##
# PID Autotuning for outer velocity loop
@named cascade_tuning_model = CascadeControlledFlatDyadBot()

cascade_spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :CascadeVelocityTuning,
    model = cascade_tuning_model,
    measurement = "y2",
    control_input = "u2",
    step_input = "u2",
    step_output = "y2",
    ref = 0.0,
    Ts = 0.01,
    duration = 20.0,
    Ms = 1.4,
    Mt = 1.9,
    Mks = 100.0,
    wl = 1e-2,
    wu = 1e3,
    num_frequencies = 300,
    kd_ub = 0.0, # Tune PI controller
    # kp_guess = 0.55,
    # ki_guess = 0.55,
    # kd_guess = 1e-3,
    soft = true,
    exact_hessian = false,
    homotopy = false,
)

cascade_asol = JSC.run_analysis(cascade_spec)

plot(cascade_asol.sol)

cascade_Splot = JSC.artifacts(cascade_asol, :SensitivityFunctions)
# cascade_response_plot = JSC.artifacts(cascade_asol, :OptimizedResponse)
cascade_nyquist_plot = JSC.artifacts(cascade_asol, :NyquistPlot)
cascade_optimized_params = JSC.artifacts(cascade_asol, :OptimizedParameters)

display(cascade_optimized_params)

##

S2 = get_named_sensitivity(cascade_tuning_model, [cascade_tuning_model.y, cascade_tuning_model.y2])
Ms2, ws2 = hinfnorm2(S2)
sigmaplot(S2); hline!([Ms2], l=(:dash, :black), label="\$M_S = \$$(round(Ms2, digits=2))")
# bodeplot(S2, plotphase=false)

L2 = inv(S2) - ss(I(2))
bodeplot(L2)

## Swms = freqresp(S2, ws2)

Si2 = get_named_sensitivity(cascade_tuning_model, [cascade_tuning_model.u])
w = exp10.(LinRange(-1, 3, 1000))
Msi2, ws2 = hinfnorm2(Si2)
bodeplot(Si2, w, plotphase=false); hline!([Msi2], l=(:dash, :black), label="\$M_S = \$$(round(Msi2, digits=2))")

Li2 = minreal(inv(Si2) - 1, 1e-12)
Li2 = convert(StateSpace{Continuous, Float64}, Li2)
nyquistplot(Li2)
marginplot(Li2, w, adjust_phase_start=true)


Li22 = get_named_looptransfer(cascade_tuning_model, [cascade_tuning_model.u]) |> minreal

dmi2 = diskmargin(Li2)
plot(dmi2)

dmi22 = diskmargin(Li22)
plot(dmi22)
