# Multibody model, does not work due to MTK bugs
using Multibody, ModelingToolkit, Test
t = Multibody.t
import Multibody.PlanarMechanics as Pl
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Blocks
using Random
@component function PlanarMultibodybot(; name, g=9.81)
    pars = @parameters begin
        M = 1.0,     [description="Body mass"]
        m = 0.1,     [description="Wheel mass"]
        R = 0.1,     [description="Wheel radius"]
        L = 0.5,     [description="Distance from wheel axis to body center of mass"]
        Ic = 0.1,    [description="Body moment of inertia"]
        Iw = 0.01,   [description="Wheel moment of inertia"]
        g = g,    [description="Gravity"]
        b_rot = 0.01,   [description="Rotational damping coefficient"]
    end

    systems = @named begin
        body = Pl.Body(m = M, I = Ic, phi=0, w=0.0, radius=0.05, gy=-g)
        wheel_body = Pl.Body(m = m, I = Iw, radius=0.01, gy=-g, phi=0, w=0.00)
        translation_cm = Pl.FixedTranslation(r = [0, L], radius=0.01)
        wheel_rotation = Pl.Revolute(axisflange=true, render=false)
        torque = Rotational.Torque(use_support=true)
        # wheelJoint = Pl.OneDOFSlippingWheelJoint(
        #     radius = R,
        #     x = 0,
        #     v = 0,
        #     mu_A = 1,
        #     mu_S = 0.7,
        #     sAdhesion = 0.04,
        #     sSlide = 0.12,
        #     vAdhesion_min = 0.01,
        #     vSlide_min = 0.03,
        # )

        wheelJoint = Pl.OneDOFRollingWheelJoint(
            radius = R,
            x = nothing,
            v = nothing,
        )

        # Control input and state outputs
        control_input = Blocks.RealInput()
        x_output = Blocks.RealOutput()
        theta_output = Blocks.RealOutput()
        x_dot_output = Blocks.RealOutput()
        theta_dot_output = Blocks.RealOutput()

        # Rotational damping (body tilt)
        rot_damper = Rotational.Damper(d = b_rot)
        fixed = Rotational.Fixed()

    end

    vars = @variables begin
        tau(t), [description="Applied torque"]
    end

    eqs = [
        # Mechanical connections
        connect(wheelJoint.frame_a, wheel_body.frame_a, wheel_rotation.frame_a)
        connect(wheel_rotation.frame_b, translation_cm.frame_a)
        connect(translation_cm.frame_b, body.frame_a)

        # Torque from control input (negate as in flat model)
        tau ~ control_input.u
        torque.tau.u ~ tau
        # tau ~ 0
        connect(torque.flange, wheel_rotation.flange_a, rot_damper.flange_a)
        connect(torque.support, wheel_rotation.support, rot_damper.flange_b)

        x_output.u ~ -wheelJoint.x
        theta_output.u ~ body.phi
        x_dot_output.u ~ -wheelJoint.v
        theta_dot_output.u ~ wheel_rotation.w
    ]

    System(eqs, t, vars, pars; systems, name)
end


# @named model = PlanarMultibodybot()
# ssys = multibody(model)

# guesses = []


# prob = ODEProblem(ssys,
#     [
#         ssys.wheelJoint.frame_a.render => true, ssys.wheelJoint.frame_a.length => 0.1, ssys.wheelJoint.frame_a.radius => 0.002,
#         ssys.body.frame_a.render => true, ssys.body.frame_a.length => 0.1, ssys.body.frame_a.radius => 0.002,
#     ], (0.0, 2); 
#     missing_guess_value = MissingGuessValue.Constant(0.0),
#     # missing_guess_value = MissingGuessValue.Random(Random.GLOBAL_RNG),
# )


# using OrdinaryDiffEq
# sol = solve(prob, Rodas5P(autodiff=true))#, dt=0.01, dtmax=0.01, force_dtmin=true)


# plot(sol)

# Multibody.render(model, sol, 0.0, lookat=[0,0.1,0], x=0, y=0.1, z=-0.5)[1]


##

