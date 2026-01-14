# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
using ModelingToolkit
using Multibody
import ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Mechanical.Rotational

# A 3D model of a segway using RollingWheelSet (does not allow tilting of body)


##

# ==============================================================================
# DyadBot3DWheels - 3D segway using individual RollingWheel components
# This allows the body to tilt forward/backward (unlike RollingWheelSet which
# keeps wheels vertical)
# ==============================================================================

@component function DyadBot3DWheels(; name)
    pars = @parameters begin
        track = 0.13#, [description = "Distance between wheels"]
        body_height = 0.1#, [description = "Height of body above wheel axle"]
        wheel_radius = 0.04
        body_mass = 0.1
        wheel_mass = 0.05
        wheel_I_axis = 5e-5
        wheel_I_long = 1e-5
        d_wheel = 0.1, [description = "Wheel rotational damping coefficient"]
    end

    systems = @named begin
        world = World()

        # Two individual rolling wheels (allows tilting)
        wheel_left = SlippingWheel(
            angles = nothing,
            der_angles = nothing,
            radius = wheel_radius,
            m = wheel_mass,
            I_axis = wheel_I_axis,
            I_long = wheel_I_long,
            state = false,
        )
        wheel_right = SlippingWheel(
            angles = nothing,
            der_angles = nothing,
            radius = wheel_radius,
            m = wheel_mass,
            I_axis = wheel_I_axis,
            I_long = wheel_I_long,
            # iscut = true,  # Avoid over-constraining
            state = true,
        )

        # Revolute joints for wheel spin
        revolute_left = Revolute(n = [0, 0, 1], axisflange = true, phi0=nothing, w0=nothing)
        revolute_right = Revolute(n = [0, 0, 1], axisflange = true, phi0=nothing, w0=nothing)

        # Axis connecting wheels
        rod_left = FixedTranslation(r = [0, 0, track/2])
        rod_right = FixedTranslation(r = [0, 0, -track/2])

        # Axis body (small mass at center)
        axis_body = Body(m = 0.01, r_cm = [0, 0, 0])

        # Main body extending upward
        body = BodyShape(
            m = body_mass,
            I_22 = 0.01*0.03^2,
            I_11 = 0.01*0.05^2,
            I_33 = 0.01*0.05^2,
            r = [0, body_height, 0],
        )

        # Rotational damping for wheel joints
        damper_left = Rotational.Damper(d = d_wheel)
        damper_right = Rotational.Damper(d = d_wheel)
    end

    eqs = [
        # Connect wheels to revolute joints
        connect(wheel_left.frame_a, revolute_left.frame_b)
        connect(wheel_right.frame_a, revolute_right.frame_b)

        # Connect revolute joints to axis via rods
        connect(revolute_left.frame_a, rod_left.frame_b)
        connect(revolute_right.frame_a, rod_right.frame_b)

        # Connect rods to axis center
        connect(rod_left.frame_a, axis_body.frame_a)
        connect(rod_right.frame_a, axis_body.frame_a)

        # Connect body to axis center
        connect(axis_body.frame_a, body.frame_a)

        # Connect dampers to wheel revolute joints
        connect(damper_left.flange_a, revolute_left.axis)
        connect(damper_left.flange_b, revolute_left.support)
        connect(damper_right.flange_a, revolute_right.axis)
        connect(damper_right.flange_b, revolute_right.support)
    ]

    System(eqs, t, [], pars; systems, name)
end

##

@named model_wheels = DyadBot3DWheels()

ssys = multibody(model_wheels)


x0 = [
    # collect(ssys.axis_body.Q̂) .=> [1, 0, 0, 0];
]

guesses = Dict([
    unknowns(ssys) .=> 0.0;
    # (ssys.wheel_right.wheeljoint.vContact_0)[3]  => 0.0
    # (ssys.wheel_left.wheeljoint.vContact_0)[1]  => 0.0
    # (ssys.axis_body.a_0)[2]  => 0.0
    # D(D((ssys.wheel_right.wheeljoint.angles)[3]))  => 0.0
    # (ssys.wheel_left.wheeljoint.e_axis_0)[2]  => 0.0
    # (ssys.wheel_left.wheeljoint.e_axis_0)[1]  => 0.0
    # D(D((ssys.wheel_left.wheeljoint.delta_0)[3]))  => 0.0
    # D((ssys.wheel_left.wheeljoint.delta_0)[3])  => 0.0
    # (ssys.wheel_left.body.a_0)[3]  => 0.0
    # D(D((ssys.wheel_right.wheeljoint.angles)[2]))  => 0.0
    # D((ssys.wheel_right.wheeljoint.delta_0)[1])  => 0.0
    # D(D((ssys.wheel_right.wheeljoint.delta_0)[3]))  => 0.0
])

initsys_mtkcompile_kwargs = (allow_symbolic=true, )

prob = ODEProblem(ssys, x0, (0.0, 5.0); guesses, missing_guess_value = MissingGuessValue.Constant(0.0001));

using OrdinaryDiffEq
sol_wheels = solve(prob, Rodas5P())

##

