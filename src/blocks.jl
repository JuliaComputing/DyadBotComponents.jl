import ModelingToolkitStandardLibrary: Blocks
import ModelingToolkitParameters: Params
using ModelingToolkit

Base.@kwdef mutable struct AddParams <: Params
    # parameters
    k1::Real = 1.0
    k2::Real = 1.0
end

Base.@kwdef mutable struct GainParams <: Params
    # parameters
    u_start::Real = 0.0
    y_start::Real = 0.0
    k::Real
end

Base.@kwdef mutable struct Add3Params <: Params
    # parameters
    k1::Real = 1.0
    k2::Real = 1.0
    k3::Real = 1.0
end

Base.@kwdef mutable struct LimiterParams <: Params
    # parameters
    u_start::Real = 0.5
    y_start::Real = 0.5
    y_max::Real = 1
    y_min::Real = 0
end

Base.@kwdef mutable struct ConstantParams <: Params
    # parameters
    k::Real = 0.0
end

Base.@kwdef mutable struct IntegratorParams <: Params
    # parameters
    u_start::Real = 0.0
    y_start::Real = 0.0
    k::Real = 1
end

Base.@kwdef mutable struct DerivativeParams <: Params
    # parameters
    u_start::Real = 0.0
    y_start::Real = 0.0
    T::Real = NaN
    k::Real = 1
end

Base.@kwdef mutable struct LimPIDParams <: Params
    # systems
    addP::AddParams = AddParams(k1 = 1, k2 = -1)
    gainPID::GainParams = GainParams()
    addPID::Add3Params = Add3Params()
    limiter::LimiterParams = LimiterParams()
    addSat::AddParams = AddParams(k1 = 1, k2 = -1)
    gainTrack::GainParams = GainParams()
    addI::Add3Params = Add3Params(k1 = 1, k2 = -1, k3 = 1)
    int::IntegratorParams = IntegratorParams()
    addD::AddParams = AddParams(k1 = 1, k2 = -1)
    der::DerivativeParams = DerivativeParams()

    # function LimPIDParams(;
    #     k = 1, Ti = false, Td = false, wp = 1, wd = 1,
    #     Ni = Ti == 0 ? Inf : √(max(Td / Ti, 1.0e-6)),
    #     Nd = 10,
    #     u_max = Inf,
    #     u_min = u_max > 0 ? -u_max : -Inf,
    #     gains = false,
    #     int__x = 0.0,
    #     der__x = 0.0
    # )
    #     P = new()
    #     P.gainPID.k = k
    #     P.limiter.y_max = u_max
    #     P.limiter.y_min = u_min
    #     P.gainTrack.k = 1/(k*Ni)
    #     P.int.k = 1/Ti
    #     P.int.x = int__x
    #     P.der.k = Td
    #     P.der.T = 1/Nd
    #     P.der.x = der__x
        
    # end
end



@component function LimPID(;
        name, k = 1, Ti = false, Td = false, wp = 1, wd = 1,
        Ni = Ti == 0 ? Inf : √(max(Td / Ti, 1.0e-6)),
        Nd = 10,
        u_max = Inf,
        u_min = u_max > 0 ? -u_max : -Inf,
        gains = false,
        int__x = 0.0,
        der__x = 0.0
    )
    with_I = true # !isequal(Ti, false)
    with_D = true # !isequal(Td, false)
    with_AWM = Ni != Inf
    if gains
        Ti = k / Ti
        Td = Td / k
    end
    # @symcheck Ti ≥ 0 ||
    #     throw(ArgumentError("Ti out of bounds, got $(Ti) but expected Ti ≥ 0"))
    # @symcheck Td ≥ 0 ||
    #     throw(ArgumentError("Td out of bounds, got $(Td) but expected Td ≥ 0"))
    # @symcheck u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))
    # @symcheck Nd > 0 ||
    #     throw(ArgumentError("Nd out of bounds, got $(Nd) but expected Nd > 0"))

    pars = []
    # @parameters begin
        # k = k, [description = "Proportional gain"]
        # Ti = Ti, [description = "Integrator time constant"]
        # Td = Td, [description = "Derivative time constant"]
        # wp = wp, [description = "Set-point weighting in the proportional part"]
        # wd = wd, [description = "Set-point weighting in the derivative part"]
        # Ni = Ni, [description = "Anti-windup tracking gain"]
        # Nd = Nd, [description = "Derivative limit"]
        # u_max = u_max, [description = "Upper saturation limit"]
        # u_min = u_min, [description = "Lower saturation limit"]
    # end
    @named reference = Blocks.RealInput()
    @named measurement = Blocks.RealInput()
    @named ctr_output = Blocks.RealOutput() # control signal
    @named addP = Blocks.Add(k1 = wp, k2 = -1)
    @named gainPID = Blocks.Gain(; k)
    @named addPID = Blocks.Add3()
    @named limiter = Blocks.Limiter(y_max = u_max, y_min = u_min)
    if with_I
        if with_AWM
            @named addI = Blocks.Add3(k1 = 1, k2 = -1, k3 = 1)
            @named addSat = Blocks.Add(k1 = 1, k2 = -1)
            @named gainTrack = Blocks.Gain(k = 1 / (k * Ni))
        else
            @named addI = Blocks.Add(k1 = 1, k2 = -1)
        end
        @named int = Blocks.Integrator(k = 1 / Ti, x = int__x)
    else
        @named Izero = Blocks.Constant(k = 0)
    end
    if with_D
        @named der = Blocks.Derivative(k = Td, T = 1 / Nd, x = der__x)
        @named addD = Blocks.Add(k1 = wd, k2 = -1)
    else
        @named Dzero = Blocks.Constant(k = 0)
    end

    sys = [reference, measurement, ctr_output, addP, gainPID, addPID, limiter]
    if with_I
        if with_AWM
            push!(sys, [addSat, gainTrack]...)
        end
        push!(sys, [addI, int]...)
    else
        push!(sys, Izero)
    end
    if with_D
        push!(sys, [addD, der]...)
    else
        push!(sys, Dzero)
    end

    eqs = [
        connect(reference, addP.input1),
        connect(measurement, addP.input2),
        connect(addP.output, addPID.input1),
        connect(addPID.output, gainPID.input),
        connect(gainPID.output, limiter.input),
        connect(limiter.output, ctr_output),
    ]
    if with_I
        push!(eqs, connect(reference, addI.input1))
        push!(eqs, connect(measurement, addI.input2))
        if with_AWM
            push!(eqs, connect(limiter.input, addSat.input2))
            push!(eqs, connect(limiter.output, addSat.input1))
            push!(eqs, connect(addSat.output, gainTrack.input))
            push!(eqs, connect(gainTrack.output, addI.input3))
        end
        push!(eqs, connect(addI.output, int.input))
        push!(eqs, connect(int.output, addPID.input3))
    else
        push!(eqs, connect(Izero.output, addPID.input3))
    end
    if with_D
        push!(eqs, connect(reference, addD.input1))
        push!(eqs, connect(measurement, addD.input2))
        push!(eqs, connect(addD.output, der.input))
        push!(eqs, connect(der.output, addPID.input2))
    else
        push!(eqs, connect(Dzero.output, addPID.input2))
    end

    #parameter_dependencies = [
    # push!(eqs, addP.k1 ~ wp)
    # push!(eqs, gainPID.k ~ k)
    # push!(eqs, gainTrack.k ~ 1 / (k * Ni))
    # push!(eqs, limiter.y_max ~ u_max)
    # push!(eqs, limiter.y_min ~ u_min)
    # push!(eqs, int.k ~ 1 / Ti)
    # if with_D
    #     push!(eqs, der.k ~ Td)
    #     push!(eqs, der.T ~ 1 / Nd)
    #     push!(eqs, addD.k1 ~ wd)
    # end
    #]

    System(eqs, t, [], pars; name = name, systems = sys)
end


@kwdef mutable struct ControllerParams <: Params
  #set_point::Float64 = 0
  kp::Float64 = 1.0
  ki::Float64 = 0.1
  kd::Float64 = 0.0
end

"""
  Controller\\_(set\\_point=0, kp=1, ki=0.1, kd=0, initial\\_output=0, output\\_limit=0)

PI controller
"""
@component function Controller(; name)
  params = @parameters begin
    kp
    ki
    kd 
    # initial_output
    # output_limit
  end
  systems = @named begin
    measurement = Blocks.RealInput()
    ctr_output = Blocks.RealOutput()
    reference = Blocks.RealInput()
  end
  vars = @variables begin
    x(t), [guess=0]
    dx(t), [guess=0]
    ddx(t), [guess=0]
    y(t)=0
    dy(t), [guess=0]
  end
  eqs = [

    x ~ reference.u - measurement.u
    ctr_output.u ~ y

    D(y) ~ dy
    D(measurement.u) ~ dx # Using measurement.u here because ERROR: ArgumentError: Differential(t)(Differential(t)(Differential(t)(x_ref(t)))) is present in the system but x_ref(t) is not an unknown.
    D(dx) ~ ddx

    dy ~ kp*(dx + ki*x + kd*ddx)

  ]
  return System(eqs, t, vars, params; systems, name)
end






