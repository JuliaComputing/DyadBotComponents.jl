import ModelingToolkitStandardLibrary: Blocks
import ModelingToolkitParameters: Params
using ModelingToolkit

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






