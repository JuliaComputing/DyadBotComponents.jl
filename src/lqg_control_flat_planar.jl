include("planar_flat.jl")  # Get FlatDyadBot

using DyadControlSystems
import DyadControlSystems as JSC
using LinearAlgebra

# Create the plant directly - no wrapper component needed

@component function Add4(; name, K)
  __vars = Any[]
  __systems = System[]
  __guesses = Dict()
  __defaults = Dict()
  __initialization_eqs = []
  __eqs = Equation[]

  r,c = size(K)
  ### Symbolic Parameters
  __params = @parameters begin
    K[1:r, 1:c] = K
  end
  append!(__vars, @variables (input1(t)::Real), [input = true])
  append!(__vars, @variables (input2(t)::Real), [input = true])
  append!(__vars, @variables (input3(t)::Real), [input = true])
  append!(__vars, @variables (input4(t)::Real), [input = true])
  append!(__vars, @variables (output(t)::Real), [output = true])
  __constants = Any[]
  push!(__eqs, output ~ dot(K, [input1, input2, input3, input4]))
  return System(__eqs, t, __vars, __params; systems=__systems, defaults=__defaults, guesses=__guesses, name, initialization_eqs=__initialization_eqs)
end

@component function LQGFlatDyadBot2(; name)
    systems = @named begin
        plant = FlatDyadBot()
        L = Add4(K = ones(1, 4))
    end

    eqs = [
        # Expose outputs with analysis point names for LQG
        connect(plant.x_output.u, :y_x, L.input1)           # Position measurement
        connect(plant.x_dot_output.u, :y_xd, L.input2)           
        connect(plant.theta_output.u, :y_theta, L.input3)   # Angle measurement
        connect(plant.theta_dot_output.u, :y_thetad, L.input4)  
        connect(L.output, :u, plant.control_input.u)        # Control input
    ]

    System(eqs, t, [], []; systems, name)
end



@named plant = LQGFlatDyadBot2()

# LQG Analysis Specification
# - 2 measurements: position (x) and angle (theta)
# - 2 controlled outputs: same as measurements (regulate both)
# - 1 control input: torque (tau)
lqg_spec = JSC.LQGAnalysisSpec(;
    name = :SegwayLQG,
    model = plant,
    measurement = ["y_x", "y_theta"],           # What we measure
    controlled_output = ["y_x", "y_xd", "y_theta", "y_thetad"],     # What we want to control
    control_input = ["u"],                       # Control input
    q1_diag = [10.0, 0.1, 1, 0.1],     # Penalty on controlled outputs (x, theta)
    q2_diag = [0.0001],          # Penalty on control input (tau)
    r1_diag = [1.0],          # Disturbance noise covariance
    r2_diag = [0.001, 0.001],     # Measurement noise covariance (x, theta)
    wl = 1e-2,
    wu = 1e3,
    num_frequencies = 200,
)

# Run the LQG analysis
lqg_asol = JSC.run_analysis(lqg_spec)

# Visualize results
step_response = JSC.artifacts(lqg_asol, :StepResponse)
gang_of_four = JSC.artifacts(lqg_asol, :GangOfFour)
bode_plot = JSC.artifacts(lqg_asol, :BodePlot)
margin_plot = JSC.artifacts(lqg_asol, :MarginPlot)
controller_gain = lqg_asol.L
observer_gain = lqg_asol.K

display(controller_gain)
display(observer_gain)

##

# cm = complete(plant)
ssys = mtkcompile(plant)
op = [
    ssys.L.K => -lqg_asol.L / lqg_asol.P_reduced.C
]


prob = ODEProblem(ssys, op, (0, 1))
sol = solve(prob, Rodas5P())
plot(sol, ylims=(-5, 5)); hline!([π], l=(:dash, :black), primary=false)


##

get_Cfffb(; name) = System(ss(lqg_asol.Cfffb); name)
get_Cff(; name) = System(ss(lqg_asol.Cff); name)
get_Cfb(; name) = System(ss(lqg_asol.Cfb); name)

@component function LQGFlatDyadBot4(; name)
    systems = @named begin
        plantfffb = FlatDyadBot()
        C = get_Cfffb()
        plant = FlatDyadBot()
        Cff = get_Cff()
        Cfb = get_Cfb()
        add = Blocks.Add(k1=1, k2=-1)
    end

    eqs = [
        C.input.u[1] ~ ifelse(t>5, 0.15, 0) # rx
        C.input.u[2] ~ 0 # rxd
        C.input.u[3] ~ 0 # rtheta
        C.input.u[4] ~ 0 # rthetad
        C.input.u[5] ~ plantfffb.x_output.u # x
        C.input.u[6] ~ plantfffb.theta_output.u-pi # theta
        connect(C.output, :ufffb, plantfffb.control_input)

        Cff.input.u[1] ~ ifelse(t>5, 0.15, 0) # rx
        Cff.input.u[2] ~ 0 # rxd
        Cff.input.u[3] ~ 0 # rtheta
        Cff.input.u[4] ~ 0 # rthetad
        Cfb.input.u[1] ~ plant.x_output.u # x
        Cfb.input.u[2] ~ plant.theta_output.u-pi # theta
        connect(Cff.output, add.input1)
        connect(Cfb.output, add.input2)
        connect(add.output, :u, plant.control_input)
    ]
    guesses = [
        plantfffb.theta_ddot => 0,
        plant.theta_ddot => 0
    ]
    System(eqs, t, [], []; systems, name, guesses)
end

@named lqg_cl = LQGFlatDyadBot4()

ssys = mtkcompile(lqg_cl)

op = [
    # collect(ssys.C.x) .=> nothing
    # ssys.C.output.u => 0
    # D(ssys.C.output.u) => 0
    ssys.plant.theta => deg2rad(160)
    ssys.plantfffb.theta => deg2rad(160)
]

prob = ODEProblem(ssys, op, (0, 20))
sol = solve(prob, Rodas5P())
plot(sol, idxs=[
    ssys.plant.theta
    ssys.plant.x
    ssys.plant.tau
    ssys.plantfffb.theta
    ssys.plantfffb.x
    ssys.plantfffb.tau
]); hline!([π 0.3], l=(:dash, :black), primary=false, ylims=(-1, 4))


##

LT = lqg_asol.Cfb*system_mapping(lqg_asol.P_ext)

plot(diskmargin(LT))

##
L2 = get_named_looptransfer(lqg_cl, [lqg_cl.u])
L2 = -L2#minreal(L2, 1e-12)

bodeplot([LT, L2]) # why not identical?


dmi2 = diskmargin(L2)
plot!(dmi2)

# TODO: test the Cff and Cfb separate contorllers, they do not appear to agree with Cfffb at all