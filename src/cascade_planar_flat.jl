using  ModelingToolkitParameters

Base.@kwdef mutable struct CascadeControlledFlatDyadBotParams <: Params
    # systems
    plant::FlatDyadBotParams = FlatDyadBotParams()
    inner_controller::LimPIDParams = LimPIDParams(k=15.6, Ti=Inf, Td=0.16, Nd=25, u_max=7)
    outer_controller::LimPIDParams = LimPIDParams(k=0.54, Ti=2.48, Td=0, Nd=600, wd=1, wp=0.5)
end


# Cascade control: outer velocity loop + inner angle loop
@component function CascadeControlledFlatDyadBot(; name)
    vars = @variables begin
        x_ref(t), [input=true]
    end

    systems = @named begin
        plant = FlatDyadBot()
        # Inner loop: angle controller
        inner_controller = Blocks.LimPID(k=15.6, Ti=Inf, Td=0.16, Nd=25, u_max=7)
        # Outer loop: velocity controller
        outer_controller = Blocks.LimPID(k=0.54, Ti=2.48, Td=0, Nd=600, wd=1, wp=0.5)
        # neg_gain = Blocks.Gain(k=1)
        # ref = Blocks.Step(height=x_ref, start_time=5)
        # ref = Blocks.Square(;  smooth = true)
        # Add pi offset to inner loop reference
        pi_offset = Blocks.Constant(k=pi)
        add_pi = Blocks.Add(k1=1, k2=1)
    end

    eqs = [
        # Outer loop: velocity reference -> angle reference
        # connect(ref.output, :r2, outer_controller.reference)
        outer_controller.reference.u ~ x_ref
        connect(plant.x_output, :y2, outer_controller.measurement)

        # Add pi to outer controller output for inner loop reference
        connect(outer_controller.ctr_output, :u2, add_pi.input1)
        connect(pi_offset.output, add_pi.input2)

        # Inner loop: angle reference -> torque
        connect(add_pi.output, inner_controller.reference)
        connect(plant.theta_output, :y, inner_controller.measurement)
        connect(inner_controller.ctr_output, :u, plant.control_input)
    ]

    System(eqs, t, vars, []; systems, name)
end
