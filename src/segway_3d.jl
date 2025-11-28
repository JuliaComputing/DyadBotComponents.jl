using ModelingToolkit
using Multibody
import ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq

# A 3D model of a segway using RollingWheelSet

@component function DyadBot3D(; name)
    pars = @parameters begin
        wheel_radius = 0.25
        body_mass = 10.0
        body_height = 0.8
    end

    systems = @named begin
        world = World()
        wheels = RollingWheelSet(
            radius = wheel_radius,
            m_wheel = 1.0,
            I_axis = 0.1,
            I_long = 0.1,
            track = 0.5,        # distance between wheels
        )
        body_offset = FixedTranslation(r = [0, body_height/2, 0])  # raise body above wheel axis
        body = BodyShape(m = body_mass, r = [0.1, body_height, 0.1])  # thin tall body
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

