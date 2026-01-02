using ModelingToolkit
using Multibody
import ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Mechanical.Rotational
using OrdinaryDiffEq

# A 3D model of a segway using RollingWheelSet (does not allow tilting of body)

# @component function DyadBot3D(; name)
#     pars = @parameters begin
#         wheel_radius = 0.04
#         body_mass = 0.1
#         body_height = 0.1
#         body_width = 0.13
#         body_depth = 0.07
#     end

#     systems = @named begin
#         world = World()
#         wheels = RollingWheelSet(
#             radius = wheel_radius,
#             m_wheel = 0.05, # mass of one wheel
#             I_axis = 5e-5,  # moment of inertia of one wheel around the rotation axis
#             I_long = 1e-5,
#             track = body_width,   # distance between wheels
#         )
#         body = BodyShape(
#             m = body_mass,
#             I_22 = 0.01*0.03^2, # Inertia around vertical axis, a very rough approximation
#             I_11 = 0.01*0.05^2, # Total guesses
#             I_33 = 0.01*0.05^2, 
#             r = [0.0, body_height, 0.0], # Vector from `frame_a` to `frame_b` (head) resolved in `frame_a`
#             height = body_depth, # we use the word "height" for the length of r
#             width = body_width,
#         )  
#     end

#     eqs = [
#         connect(wheels.frame_middle, body_offset.frame_a)
#         connect(body_offset.frame_b, body.frame_a)
#     ]

#     System(eqs, t, [], pars; systems, name)
# end


# @named model = DyadBot3D()
# model = complete(model)

# ssys = structural_simplify(multibody(model))

# prob = ODEProblem(ssys, [], (0.0, 5.0))

# sol = solve(prob, Rodas5P())

# ##
# import GLMakie


##

# ==============================================================================
# DyadBot3DWheels - 3D segway using individual RollingWheel components
# This allows the body to tilt forward/backward (unlike RollingWheelSet which
# keeps wheels vertical)
# ==============================================================================

@component function DyadBot3DWheels(; name)
    body_height = 0.1#, [description = "Height of body above wheel axle"]
    track = 0.13#, [description = "Distance between wheels"]
    pars = @parameters begin
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
            radius = wheel_radius,
            m = wheel_mass,
            I_axis = wheel_I_axis,
            I_long = wheel_I_long,
            state = false,
        )
        wheel_right = SlippingWheel(
            radius = wheel_radius,
            m = wheel_mass,
            I_axis = wheel_I_axis,
            I_long = wheel_I_long,
            # iscut = true,  # Avoid over-constraining
            state = true,
        )

        # Revolute joints for wheel spin
        revolute_left = Revolute(n = [0, 0, 1], axisflange = true, phi0=0, w0=0)
        revolute_right = Revolute(n = [0, 0, 1], axisflange = true, phi0=0, w0=0)

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
model_wheels = complete(model_wheels)

linsys = (; allow_symbolic = false, inline_linear_sccs = true, analytical_linear_scc_limit = 5, reassemble_alg = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 5))
ssys = structural_simplify(multibody(model_wheels); linsys...)


x0 = [
    # collect(ssys.axis_body.Q̂) .=> [1, 0, 0, 0];
]

guesses = unknowns(ssys) .=> 0.0

prob = ODEProblem(ssys, x0, (0.0, 5.0); guesses)

sol_wheels = solve(prob, Rodas5P())

##

