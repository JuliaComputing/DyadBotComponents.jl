using ModelingToolkit
import ModelingToolkit: t_nounits as t, D_nounits as D
# import ModelingToolkitStandardLibrary.Mechanical.Rotational 
import ModelingToolkitStandardLibrary.Blocks 
import ModelingToolkitParameters: Params
# using OrdinaryDiffEq
# using DyadControlSystems, ControlSystemsBase, ControlSystemsMTK
# # using Plots
# import DyadControlSystems as JSC
# # using LinearAlgebra


# # Plots.default(size=(1200,1200))
# connect = ModelingToolkit.connect


##


# ==============================================================================
## FlatDyadBot - Equation-based planar segway model
# ==============================================================================

Base.@kwdef mutable struct FlatDyadBotParams <: Params
    # parameters
    M::Real = 1.0
    m::Real = 0.1
    R::Real = 0.1
    L::Real = 0.5
    Ic::Real = 0.1
    Iw::Real = 0.01
    g::Real = 9.81
    b_trans::Real = 1.0
    b_rot::Real = 1.0
end


@component function FlatDyadBot(; name)
    pars = @parameters begin
        M = 1.0,     [description="Body mass"]
        m = 0.1,     [description="Wheel mass"]
        R = 0.1,     [description="Wheel radius"]
        L = 0.5,     [description="Distance from wheel axis to body center of mass"]
        Ic = 0.1,    [description="Body moment of inertia"]
        Iw = 0.01,   [description="Wheel moment of inertia"]
        g = 9.81,    [description="Gravity"]
        b_trans = 10.0,     [description="Translational Damping coefficient"]
        b_rot = 10.0,     [description="Rotational Damping coefficient"]
    end

    systems = @named begin
        control_input = Blocks.RealInput()
        x_output = Blocks.RealOutput()
        theta_output = Blocks.RealOutput()
        x_dot_output = Blocks.RealOutput()
        theta_dot_output = Blocks.RealOutput()
    end

    vars = @variables begin
        x(t) = 0.0,          [description="Horizontal position"]
        theta(t) = deg2rad(180),      [description="Body angle (from vertical down)"]
        x_dot(t) = 0.0,      [description="Horizontal velocity"]
        theta_dot(t) = 0.0,  [description="Angular velocity"]
        x_ddot(t),           [description="Horizontal acceleration"]
        theta_ddot(t),       [description="Angular acceleration"]
        tau(t),              [description="Input torque"]
    end

    # Mass matrix elements
    # M11 = (M+m) + Iw/R^2
    # M12 = M21 = M*L*cos(theta)
    # M22 = Ic + M*L^2

    # RHS = G - C + B*tau where:
    # G = [0; -M*L*g*sin(theta)]
    # C = [-M*L*theta_dot^2*sin(theta) - (b/R^2)*x_dot; b*theta_dot]
    # B*tau = [tau/R; -tau]

    eqs = [
        # Connect input/outputs
        tau ~ -control_input.u
        x_output.u ~ x
        theta_output.u ~ theta
        x_dot_output.u ~ x_dot
        theta_dot_output.u ~ theta_dot

        # Kinematic equations
        D(x) ~ x_dot
        D(theta) ~ theta_dot
        D(x_dot) ~ x_ddot
        D(theta_dot) ~ theta_ddot

        # Mass matrix equation: M * [x_ddot; theta_ddot] = RHS
        # Row 1: ((M+m) + Iw/R^2)*x_ddot + M*L*cos(theta)*theta_ddot = RHS1
        # Row 2: M*L*cos(theta)*x_ddot + (Ic + M*L^2)*theta_ddot = RHS2

        ((M + m) + Iw/R^2) * x_ddot + M*L*cos(theta) * theta_ddot ~
            M*L*theta_dot^2*sin(theta) - (b_trans/R^2)*x_dot + tau/R

        M*L*cos(theta) * x_ddot + (Ic + M*L^2) * theta_ddot ~
            -M*L*g*sin(theta) - b_rot*theta_dot - tau

        # tau ~ 0
    ]

    guesses = [
        x_ddot => -1
        # theta_ddot => 0
        # tau => 0
    ]

    System(eqs, t, vars, pars; systems, name, guesses)
end
##
