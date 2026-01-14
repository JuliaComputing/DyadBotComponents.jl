using DyadBotComponents
using DyadControlSystems
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t
import DyadControlSystems as JSC
using ModelingToolkitParameters
using WGLMakie

# get parameter defaults
bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()

@component function CascadeControlledFlatDyadBotInput(; name)
    systems = @named begin
        bot = DyadBotComponents.CascadeControlledFlatDyadBot()
        position = Blocks.Constant(k=0)
    end
    eqs = [
        ModelingToolkit.connect(bot.ref, position.output)
    ]
    defaults=bot=>bot_params
    System(eqs, t, [], []; systems, name, defaults)
end

@named sys = CascadeControlledFlatDyadBotInput()
spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :sys,
    model = sys,
    measurement = "bot.y",
    control_input = "bot.u",
    step_input = "bot.u",
    step_output = "bot.y",
    Ts = 0.01,           # Sample time
    duration = 25.0,      # Simulation duration
    Ms = 1.5,            # Sensitivity peak constraint
    Mt = 1.5,            # Complementary sensitivity peak constraint
    Mks = 400.0,         # Control sensitivity constraint
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
bot_params.inner_controller.kp = arts[1, :Kp_standard]
bot_params.inner_controller.ki = 1/arts[1, :Ti_standard]
bot_params.inner_controller.kd = arts[1, :Td_standard]

# Simulate Tilted Robot
bot_params.plant.theta_init = 0.75*π
sys = mtkcompile(sys)
prob = ODEProblem(sys, sys.bot => bot_params, (0, 1))
sol = solve(prob)

plot(sol; idxs=sys.bot.plant.theta)


# Next Step - Solve Outer Loop ..........
# @named bot = DyadBotComponents.CascadeControlledFlatDyadBot()

# spec = JSC.PIDAutotuningAnalysisSpec(;
#     name = :bot,
#     model = bot,
#     measurement = "y2",
#     control_input = "u2",
#     step_input = "u2",
#     step_output = "y2",
#     Ts = 0.01,           # Sample time
#     duration = 25.0,      # Simulation duration
#     Ms = 1.5,            # Sensitivity peak constraint
#     Mt = 1.5,            # Complementary sensitivity peak constraint
#     Mks = 400.0,         # Control sensitivity constraint
#     wl = 1e-2,           # Lower frequency bound
#     wu = 1e3,            # Upper frequency bound
#     kd_ub = 0.0,         # Tune PI controller
#     num_frequencies = 200,
#     soft = true
# )