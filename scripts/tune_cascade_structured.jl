# Structured tuning of both cascade PID controllers of the balancing robot
# at once.
#
# Instead of tuning the loops one at a time (scripts tune_angle_pid.jl and
# tune_cascade_pid.jl), a single StructuredAutoTuningProblem holds the
# objectives and all six tunable gains (proportional gain, integral and
# derivative time constants of both controllers), exposed as top-level
# parameters of CascadeControlledDyadBot.
#
# Objectives:
#   * MaximumSensitivityObjective / MaximumComplementarySensitivityObjective
#     per loop (SISO, at y for the angle loop and y2 for the position loop).
#     The two measurements are a tilt angle and a position: their scales
#     differ and the loops are nested, so a joint MIMO singular-value bound
#     over [y, y2] is dominated by the unit mismatch (the maximum singular
#     value of the 2x2 sensitivity exceeds 14 even for a well-tuned design)
#     and is not a meaningful robustness measure here.
#   * PoleLocationObjective: the sensitivity bounds alone do not guarantee
#     stability over a finite frequency grid, so right-half-plane closed-loop
#     poles are forbidden explicitly.

cd(@__DIR__)
using Pkg
Pkg.activate(".")
using DyadBotComponents
using MultibodyComponents
using ModelingToolkit
using DyadControlSystems
using DyadControlSystems.MPC
using ControlSystemsBase
using OrdinaryDiffEq
using Plots

# Zero initial tilt, so that the model's initial equations are consistent
# with the upright trim used for linearization
@named model = DyadBotComponents.CascadeControlledDyadBot(phi0 = 0)
m = ModelingToolkit.toggle_namespacing(model, false)
ssys = multibody(model)

# Complementary-sensitivity weights: |W| ≈ Mt in the bandwidth with a
# first-order roll-off above ωc (both closed loops roll off at first order,
# so a steeper weight can never be satisfied at high frequency). The inner
# angle loop is fast, the outer position loop is slow (crossover around half
# a rad/s, limited by the unstable pendulum pole), so the two loops get
# different roll-off frequencies.
Mt = 2.0
ωc_angle = 200.0
WT_angle = tf(Mt * ωc_angle, [1.0, ωc_angle])
ωc_pos = 10.0
WT_pos = tf(Mt * ωc_pos, [1.0, ωc_pos])

# Sensitivity weights. The outer position loop owns the low-frequency range:
# its weight combines the peak bound Ms with disturbance rejection below the
# corner ωb (which must sit well below the loop crossover). The inner angle
# loop does not reject low-frequency output disturbances on its own (the
# outer loop commands the angle at DC), so it only gets a flat peak bound.
Ms = 2.0
WS_angle = tf(Ms)
ωb_pos = 2π * 0.005
WS_pos = tf([Ms, 0.0], [1.0, ωb_pos])

# Frequency and time grids (the time grid is unused by the frequency-domain
# objectives but required by the problem constructor)
w = exp10.(LinRange(-2, 3, 300))
t = 0:0.01:3

# Operating point: the upright configuration with all states pinned to zero.
# The square-wave amplitude is deliberately NOT zeroed here: the operating
# point is also applied to the simulation objective's problem, which must
# run with the reference active. The source value does not affect the
# linearizations (it enters through linear blocks only).
op = Dict([
    # ssys.controller.angle_controller.u_m => 0
    ssys.controller.angle_controller.integrator.y => 0
    ssys.controller.angle_controller.derivative.x => 0
    # ssys.controller.pos_controller.u_m => 0
    ssys.controller.pos_controller.integrator.y => 0
    ssys.controller.pos_controller.derivative.x => 0
    ssys.plant.body_mass.body.phi => 0
    ssys.plant.body_mass.body.w => 0
    ssys.plant.wheelinertia.phi => 0
    ssys.plant.wheelinertia.w => 0
    ssys.firstorder.x => 0
    ssys.firstorder1.x => 0
])
operating_points = [op]

# Simulation objective: track the filtered position reference. Without it,
# collapsing the position-controller gain toward zero satisfies all
# frequency-domain bounds trivially (an open position loop has benign
# sensitivity functions), so tracking performance must be rewarded
# explicitly.
simprob = ODEProblem(ssys, [], (0.0, 20.0))
terr = 0:0.1:20
function tracking_cost(sol)
    sol.retcode == ReturnCode.Success || return 1e6 # penalize gains for which the simulation fails
    e = sol(terr, idxs = ssys.plant.x).u .- sol(terr, idxs = ssys.firstorder1.y).u
    1e2 * sum(abs2, e) / length(terr)
end
simobj = SimulationObjective(; costfun = tracking_cost, prob = simprob,
    solve_args = (Rodas5P(),), solve_kwargs = (; reltol = 1e-8, abstol = 1e-8))

objectives = [
    simobj,
    MaximumSensitivityObjective(WS_angle, m.controller.y),
    MaximumSensitivityObjective(WS_pos, m.controller.y2),
    MaximumComplementarySensitivityObjective(WT_angle, m.controller.y),
    MaximumComplementarySensitivityObjective(WT_pos, m.controller.y2),
    # The sensitivity bounds alone do not guarantee stability over a finite
    # frequency grid; forbid right-half-plane closed-loop poles. The margin
    # (required decay rate of every pole) must be small: the position loop of
    # this robot is slow, with its dominant poles well below one rad/s.
    PoleLocationObjective(output = m.controller.y, margin = 0.01, reduce = true),
]

# All six tunable gains and their box constraints
tunable_parameters = [
    m.k_angle => (1e-2, 10.0), m.Ti_angle => (1e-3, 10.0), m.Td_angle => (1e-3, 1.0),
    m.k_pos => (1e-3, 1.0), m.Ti_pos => (1e-1, 100.0), m.Td_pos => (1e-2, 100.0),
]
# Initial guess: the currently tuned gains
x0 = [0.487401, 0.0587352, 0.0420526, 0.0666576, 5.25024, 4.81393]

# linearize_kwargs = MultibodyComponents.linsys enables the accurate SYMBOLIC
# state-selection / reassembly during every linearization rather than the
# numerical fallback.
prob = StructuredAutoTuningProblem(model, w, t, objectives, operating_points,
    tunable_parameters; linearize_kwargs = MultibodyComponents.linsys)
res = solve(prob, x0,
    MPC.IpoptSolver(verbose = true, exact_hessian = false, acceptable_iter = 4,
        tol = 1e-2, acceptable_tol = 1e-1, max_iter = 500, printerval = 2,
        mu_strategy = "monotone"))

# Tuned gains (res.sol.u is ordered as tunable_parameters; trailing entries
# are slacks)
display(res)

# For every objective and operating point, the closed-loop transfer function
# against its weight bound
plot(res, size = (1200, 800)) |> display

## Analyze the optimal parameters
using ControlSystemsMTK

# Operating point including the optimized parameters
opopt = merge(operating_points[1], Dict(DyadControlSystems.optmap(res)))

S_angle = get_named_sensitivity(model, [m.controller.y]; op = opopt, balance = true,
    MultibodyComponents.linsys...) |> sminreal
S_pos = get_named_sensitivity(model, [m.controller.y2]; op = opopt, balance = true,
    MultibodyComponents.linsys...) |> sminreal
bodeplot([S_angle, S_pos], plotphase = false, lab = ["So angle" "So position"])

Ms_angle = hinfnorm2(S_angle)[1]
Ms_pos = hinfnorm2(S_pos)[1]
@show Ms_angle Ms_pos
max(Ms_angle, Ms_pos) < Ms || @warn "Robustness objective not met"

## Simulate with the optimized parameters (only the gains are overridden;
# the initial conditions and reference of the model are kept)
oprob = ODEProblem(ssys, Dict(DyadControlSystems.optmap(res)), (0.0, 20.0))
sol = solve(oprob, Rodas5P())
@show sol.retcode
plot(sol, idxs = [ssys.plant.theta, ssys.plant.x, ssys.square.y])



#=
optimization status : Success
objective status    
        :SimulationObjective => 0.03544817713041285
        :MaximumSensitivityObjective => 0.0
        :MaximumSensitivityObjective => 0.0
        :MaximumComplementarySensitivityObjective => 0.0
        :MaximumComplementarySensitivityObjective => 0.0
        :PoleLocationObjective => 0.0
minimizer (obtain with `optmap(res)`
        k_angle => 0.25223177002412134
        Ti_angle => 4.492649470796655
        Td_angle => 0.07535177198242332
        k_pos => 0.7382833570482313
        Ti_pos => 41.62571931622867
        Td_pos => 0.537697918931463
objective value     : 0.06372773315974468
=#