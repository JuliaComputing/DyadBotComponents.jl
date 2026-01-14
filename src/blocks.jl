import ModelingToolkitStandardLibrary: Blocks
import ModelingToolkitParameters: Params
using ModelingToolkit

@kwdef mutable struct ControllerParams <: Params
  #set_point::Float64 = 0
  k::Float64 = 1.0
  Ti::Float64 = 10.0
  Td::Float64 = 1e-6
  N::Float64 = 10.0
end

"""
  Controller\\_(set\\_point=0, k, Ti, Td, initial\\_output=0, output\\_limit=0)

PID controller

``k(1 + 1/(T_i s) + T_d s)``
"""
@component function Controller(; name)
  params = @parameters begin
    k
    Ti
    Td
    N
  end
  systems = @named begin
    measurement = Blocks.RealInput()
    ctr_output = Blocks.RealOutput()
    reference = Blocks.RealInput()
  end
  vars = @variables begin
    e(t), [guess=0]
    x(t)=0, [guess=0]
    y(t), [guess=0]
    yf(t)=0, [guess=0]
    u(t), [guess=0]
  end
  Tf = Td/N
  eqs = [
    y ~ measurement.u
    e ~ reference.u - y
    ctr_output.u ~ u
    D(x) ~ e / Ti

    Tf*D(yf) + yf ~ -D(y) # D(e) # ERROR: ArgumentError: Differential(t)(bot₊d_u2(t)) is present in the system but bot₊d_u2(t) is not an unknown.
    u ~ k*(e + x + Td*yf)

    # Tf*D(yf) + yf ~ e
    # u ~ k*(e + x + Td*D(yf))
  ]
  return System(eqs, t, vars, params; systems, name)
end






