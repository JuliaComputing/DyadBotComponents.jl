using ModelingToolkit
import ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Mechanical.Rotational 
import ModelingToolkitStandardLibrary.Blocks 
using OrdinaryDiffEq
using ControlSystemsBase, ControlSystemsMTK


#= Multibody model, does not work due to MTK bugs
using Multibody
import Multibody.PlanarMechanics as Pl
# A simple model of a planar segway

@component function PlanarDyadBot(; name)
    pars = @parameters begin
        r_cm = 0.1
    end
    systems = @named begin
        body = Pl.Body(m = 0.1, I = 0.01, phi=0)
        translation_cm = Pl.FixedTranslation(r = [0, r_cm])
        wheelJoint = Pl.SlipBasedWheelJoint(
            radius = 0.25,
            r = [0, 1],
            mu_A = 1,
            mu_S = 0.7,
            N = 1000,
            sAdhesion = 0.04,
            sSlide = 0.12,
            vAdhesion_min = 0.05,
            vSlide_min = 0.15,
            phi_roll = 0)
        inertia = Rotational.Inertia(J = 0.01, phi = 0, w = 0)
    end

    vars = @variables begin
    end

    eqs = [
        connect(wheelJoint.flange_a, inertia.flange_b)
        connect(wheelJoint.frame_a, translation_cm.frame_a)
        connect(translation_cm.frame_b, body.frame_a)
        wheelJoint.v_lat ~ 0
    ]

    System(eqs, t, vars, pars; systems, name)
end


@named model = PlanarDyadBot()
model = complete(model)

ssys = mtkcompile(model)

guesses = [
    collect(ssys.wheelJoint.v) .=> 0;
]


prob = ODEProblem(ssys, [ssys.body.r => [0, 1.0]], (0.0, 10.0); guesses)

sol = solve(prob, Rodas5P())

=#
##


# ==============================================================================
## FlatDyadBot - Equation-based planar segway model
# ==============================================================================

@component function FlatDyadBot(; name)
    pars = @parameters begin
        M = 1.0     # Body mass
        m = 0.1     # Wheel mass
        R = 0.1     # Wheel radius
        L = 0.5     # Distance from wheel axis to body center of mass
        Ic = 0.1    # Body moment of inertia
        Iw = 0.01   # Wheel moment of inertia
        g = 9.81    # Gravity
        b = 1.0     # Damping coefficient
    end

    systems = @named begin
        control_input = Blocks.RealInput()
        x_output = Blocks.RealOutput()
        theta_output = Blocks.RealOutput()
        x_dot_output = Blocks.RealOutput()
        theta_dot_output = Blocks.RealOutput()
    end

    vars = @variables begin
        x(t) = 0.0          # Horizontal position
        theta(t) = deg2rad(180)      # Body angle (from vertical down)
        x_dot(t) = 0.0      # Horizontal velocity
        theta_dot(t) = 0.0  # Angular velocity
        x_ddot(t)           # Horizontal acceleration
        theta_ddot(t)       # Angular acceleration
        tau(t)              # Input torque
    end

    # Mass matrix elements
    # M11 = (M+m) + Iw/R^2
    # M12 = M21 = M*L*cos(theta)
    # M22 = Ic + M*L^2

    # RHS = G - C + B*tau where:
    # G = [0; -M*L*g*sin(theta)]
    # C = [-M*L*theta_dot^2*sin(theta) - (b/R^2)*x_dot; b*theta_dot]
    # B*tau = [tau/R; -tau]

    eqs = [
        # Connect input/outputs
        tau ~ -control_input.u
        x_output.u ~ x
        theta_output.u ~ theta
        x_dot_output.u ~ x_dot
        theta_dot_output.u ~ theta_dot

        # Kinematic equations
        D(x) ~ x_dot
        D(theta) ~ theta_dot
        D(x_dot) ~ x_ddot
        D(theta_dot) ~ theta_ddot

        # Mass matrix equation: M * [x_ddot; theta_ddot] = RHS
        # Row 1: ((M+m) + Iw/R^2)*x_ddot + M*L*cos(theta)*theta_ddot = RHS1
        # Row 2: M*L*cos(theta)*x_ddot + (Ic + M*L^2)*theta_ddot = RHS2

        ((M + m) + Iw/R^2) * x_ddot + M*L*cos(theta) * theta_ddot ~
            M*L*theta_dot^2*sin(theta) - (b/R^2)*x_dot + tau/R

        M*L*cos(theta) * x_ddot + (Ic + M*L^2) * theta_ddot ~
            -M*L*g*sin(theta) - 0.1*b*theta_dot - tau

        # tau ~ 0
    ]

    guesses = [
        x_ddot => -1
        # theta_ddot => 0
        # tau => 0
    ]

    System(eqs, t, vars, pars; systems, name, guesses)
end
##
@component function AngleControlledFlatDyadBot(; name)
    pars = @parameters begin
        theta_ref = deg2rad(180)  # Reference angle (upright)
    end

    systems = @named begin
        plant = FlatDyadBot()
        controller = Blocks.LimPID(k=10.5, Ti=17.6, Td=0.23, Nd=38)
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
    Mt = 1.9,            # Complementary sensitivity peak constraint
    Mks = 400.0,         # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    num_frequencies = 200,
    soft = true,
    # timeweight = true,
    # soft_penalty = 1e7,
    # exact_hessian = true,
    # kp_guess = 1.0,
    # ki_guess = 0.1,
    # kd_guess = 0.0,
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
        inner_controller = Blocks.LimPID(k=10.5, Ti=17.6, Td=0.23, Nd=38)
        # Outer loop: velocity controller
        outer_controller = Blocks.LimPID(k=0.82, Ti=1.3, Td=0.02, Nd=140, wd=0, wp=0.5)
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
    duration = 6.0,
    Ms = 1.6,
    Mt = 1.5,
    Mks = 100.0,
    wl = 1e-2,
    wu = 1e4,
    num_frequencies = 300,
    kp_guess = 0.55,
    ki_guess = 0.55,
    kd_guess = 1e-3,
    # kp_lb = -100.0,
    # kp_ub = 0.0,
    # ki_lb = -10.0,
    # ki_ub = 0.0,
    # kd_lb = -10.0,
    # kd_ub = 0.0,
    soft = false,
    exact_hessian = true,
)

cascade_asol = JSC.run_analysis(cascade_spec)

plot(cascade_asol.sol)

cascade_Splot = JSC.artifacts(cascade_asol, :SensitivityFunctions)
# cascade_response_plot = JSC.artifacts(cascade_asol, :OptimizedResponse)
cascade_nyquist_plot = JSC.artifacts(cascade_asol, :NyquistPlot)
cascade_optimized_params = JSC.artifacts(cascade_asol, :OptimizedParameters)

display(cascade_optimized_params)