using ModelingToolkit
using Multibody
import ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq

# A 3D model of a segway using RollingWheelSet

@component function DyadBot3D(; name)
    pars = @parameters begin
        wheel_radius = 0.04
        body_mass = 0.1
        body_height = 0.1
        body_width = 0.13
        body_depth = 0.07
    end

    systems = @named begin
        world = World()
        wheels = RollingWheelSet(
            radius = wheel_radius,
            m_wheel = 0.05, # mass of one wheel
            I_axis = 5e-5,  # moment of inertia of one wheel around the rotation axis
            I_long = 1e-5,
            track = body_width,   # distance between wheels
        )
        body = BodyShape(
            m = body_mass,
            I_22 = 0.01*0.03^2, # Inertia around vertical axis, a very rough approximation
            I_11 = 0.01*0.05^2, # Total guesses
            I_33 = 0.01*0.05^2, 
            r = [0.0, body_height, 0.0], # Vector from `frame_a` to `frame_b` (head) resolved in `frame_a`
            height = body_depth, # we use the word "height" for the length of r
            width = body_width,
        )  
    end

    eqs = [
        connect(wheels.frame_middle, body_offset.frame_a)
        connect(body_offset.frame_b, body.frame_a)
    ]

    System(eqs, t, [], pars; systems, name)
end


@named model = DyadBot3D()
model = complete(model)

ssys = structural_simplify(multibody(model))

prob = ODEProblem(ssys, [], (0.0, 5.0))

sol = solve(prob, Rodas5P())

##
import GLMakie


##

