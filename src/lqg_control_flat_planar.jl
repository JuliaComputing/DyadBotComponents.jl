include("planar_flat.jl")  # Get FlatDyadBot


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
        L = Add4(K = 1e-10*ones(1, 4))
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
    name = :DyadBotLQG,
    model = plant,
    measurement = ["y_x", "y_theta"],           # What we measure
    controlled_output = ["y_x", "y_xd", "y_theta", "y_thetad"],     # What we want to control
    control_input = ["u"],                       # Control input
    loop_openings = ["u"],
    q1_diag = [10.0, 0.1, 1, 0.1],     # Penalty on controlled outputs (x, theta)
    q2_diag = [0.0001],          # Penalty on control input (tau)
    r1_diag = [1.0],          # Disturbance noise covariance
    r2_diag = [0.001, 0.001],     # Measurement noise covariance (x, theta)
    wl = 1e-2,
    wu = 314,
    num_frequencies = 200,
    integrator_indices = [1],
    integrator_r1_diag = [0.1],
    # disc = "zoh",
    # Ts = 0.01,
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


##

@component function LQGFlatDyadBot4(; name)
    pars = @parameters begin
        step_time = 5
        step_height = 0.15
    end
    systems = @named begin
        plant = FlatDyadBot()
        C = DyadControlSystems.get_Cfffb(lqg_asol)
    end

    eqs = [
        C.input.u[1] ~ ifelse(t>step_time, step_height, 0) # rx
        C.input.u[2] ~ 0 # rxd
        C.input.u[3] ~ 0 # rtheta
        C.input.u[4] ~ 0 # rthetad
        C.input.u[5] ~ plant.x_output.u # x
        C.input.u[6] ~ plant.theta_output.u # theta (operating-point adjustment is already taken care of in C)
        connect(C.output, :u, plant.control_input)
    ]
    guesses = [
        plant.theta_ddot => 0,
    ]
    System(eqs, t, [], pars; systems, name, guesses)
end

@named lqg_cl = LQGFlatDyadBot4()

ssys = mtkcompile(lqg_cl)

op = [
    ssys.plant.theta => deg2rad(160)
    ssys.plant.theta => deg2rad(160)
]

prob = ODEProblem(ssys, op, (0, 20))
sol = solve(prob, Rodas5P())
plot(sol, idxs=[
    ssys.plant.theta
    ssys.plant.x
    ssys.plant.tau
]); hline!([π prob.ps[ssys.step_height]], l=(:dash, :black), primary=false, ylims=(-1, 4))


##

rhp_pole, _ = findmax(real, poles(lqg_asol.P))
# fundamental limitation due to RHP pole ω_gc > 2p = 2*

Lo = system_mapping(lqg_asol.P_ext)*lqg_asol.Cfb
Li = lqg_asol.Cfb*system_mapping(lqg_asol.P_ext)

# ##
L2 = get_named_looptransfer(lqg_cl, [lqg_cl.u])
L2 = -minreal(L2, 1e-8)

S2 = get_named_sensitivity(lqg_cl, [lqg_cl.u])



dmi = diskmargin(Li, offset=0)
dmi2 = diskmargin(L2, offset=0) # offset due to hard to cancel pole/zero pair in origin

plot(dmi)
plot!(dmi2)

marginplot([Li, L2], adjust_phase_start=false) # Verify equal
vline!([2*rhp_pole], l=(:dash, :black), label="Fundamental limitation")

# ##
# P_ext = system_mapping(lqg_asol.P_ext)
# @named model = FlatDyadBot()
# cm = complete(model)
# P2 = named_ss(cm, [cm.control_input.u], [cm.x, cm.theta], op=Dict(cm.control_input.u => 0), allow_input_derivatives=true)
# bodeplot([P_ext, P2], adjust_phase_start=false)