import ModelingToolkitParameters: Params

@kwdef mutable struct LimPIDParams <: Params
        k = 1
        Ti = 0.1
        Td = 0.01
        wp = 1
        wd = 1
        Ni = Ti == 0 ? Inf : √(max(Td / Ti, 1.0e-6))
        Nd = 10
        u_max = Inf
        u_min = u_max > 0 ? -u_max : -Inf
end

Base.@kwdef mutable struct GainParams <: Params
    # parameters
    u_start::Real = 0.0
    y_start::Real = 0.0
    k::Real
end

Base.@kwdef mutable struct SquareParams <: Params
    # parameters
    frequency::Real = 1.0
    amplitude::Real = 1.0
    offset::Real = 0.0
    start_time::Real = 0.0
end

Base.@kwdef mutable struct ConstantParams <: Params
    # parameters
    k::Real = 0.0
end

Base.@kwdef mutable struct AddParams <: Params
    # parameters
    k1::Real = 1.0
    k2::Real = 1.0
end

