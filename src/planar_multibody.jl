
# Multibody model, does not work due to MTK bugs
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

