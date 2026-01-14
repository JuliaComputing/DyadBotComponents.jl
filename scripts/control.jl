using DyadBotComponents
using DyadControlSystems
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t
import DyadControlSystems as JSC
using ModelingToolkitParameters
# using WGLMakie

# get parameter defaults
bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()

@component function CascadeControlledFlatDyadBotInput(; name, bot_params)
    systems = @named begin
        bot = DyadBotComponents.CascadeControlledFlatDyadBot()
        position = Blocks.Constant(k=0.0) #NOTE: setting to 0.1 causes ERROR: ArgumentError: invalid argument #4 to LAPACK call
    end
    eqs = [
        ModelingToolkit.connect(bot.ref, position.output)
    ]
    defaults=bot=>bot_params
    System(eqs, t, [], []; systems, name, defaults)
end

@named sys = CascadeControlledFlatDyadBotInput(; bot_params)
spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :sys,
    model = sys,
    measurement = "bot.y",
    control_input = "bot.u",
    step_input = "bot.u",
    step_output = "bot.y",
    Ts = 0.01,           # Sample time
    duration = 5.0,      # Simulation duration
    Ms = 1.5,            # Sensitivity peak constraint
    Mt = 1.5,            # Complementary sensitivity peak constraint
    Mks = 500.0,         # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    ki_ub = 0.0,         # Tune PD controller
    num_frequencies = 200,
    soft = true,
    loop_openings = ["bot.y2", "bot.u2"]
)

asol = JSC.run_analysis(spec)

# DyadControlSystems.launch_pid_autotuning_designer(spec)

arts = DyadControlSystems.artifacts(asol, :OptimizedParameters)

# update parameters from auto-tune
bot_params.inner_controller.k = arts[1, :Kp_standard]
bot_params.inner_controller.Ti = arts[1, :Ti_standard]
bot_params.inner_controller.Td = arts[1, :Td_standard]
bot_params.inner_controller.N = arts[1, :Nd]

bot_params.outer_controller.k = 0
bot_params.outer_controller.Ti = 1e3
bot_params.outer_controller.Td = 1

# Simulate Tilted Robot
bot_params.plant.theta_init = deg2rad(190)
ssys = mtkcompile(sys)



prob = ODEProblem(ssys, ssys.bot => bot_params, (0, 10))
sol = solve(prob; abstol=1e-5, reltol=1e-5)

using Plots
plot(sol; idxs=ssys.bot.plant.theta); hline!([pi])


# Implement Tuned Defaults
@named sys2 = CascadeControlledFlatDyadBotInput(; bot_params)

# Next Step - Solve Outer Loop ..........
spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :sys2,
    model = sys2,
    measurement = "bot.y2",
    control_input = "bot.u2",
    step_input = "bot.u2",
    step_output = "bot.y2",
    Ts = 0.01,           # Sample time
    duration = 5.0,      # Simulation duration
    Ms = 1.5,            # Sensitivity peak constraint
    Mt = 1.5,            # Complementary sensitivity peak constraint
    Mks = 400.0,         # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    kd_ub = 0.0,         # Tune PI controller
    num_frequencies = 200,
    soft = true
)

asol = JSC.run_analysis(spec)

arts = DyadControlSystems.artifacts(asol, :OptimizedParameters)

bot_params.outer_controller.k = arts[1, :Kp_standard]
bot_params.outer_controller.Ti = arts[1, :Ti_standard]
bot_params.outer_controller.Td = arts[1, :Td_standard]
bot_params.outer_controller.N = 1 # arts[1, :Nd] #NOTE: avoid NaN

# Simulate Off Position Robot
bot_params.plant.x_init = 0.1

ssys = mtkcompile(sys)
prob = ODEProblem(ssys, ssys.bot => bot_params, (0, 10))
sol = solve(prob; abstol=1e-5, reltol=1e-5)

using Plots
p1=plot(sol; idxs=ssys.bot.plant.theta); hline!([pi])
p2=plot(sol; idxs=ssys.bot.plant.x)
plot([p1, p2]; layout=(2,1))