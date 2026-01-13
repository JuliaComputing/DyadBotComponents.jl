using DyadBotComponents
using DyadControlSystems
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t
import DyadControlSystems as JSC
using ModelingToolkitParameters
using WGLMakie

bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()

@component function CascadeControlledFlatDyadBotInput(; name)
    systems = @named begin
        plant = DyadBotComponents.CascadeControlledFlatDyadBot()
        position = Blocks.Constant(k=0)
    end
    eqs = [
        ModelingToolkit.connect(plant.ref, position.output)
    ]
    defaults=plant=>bot_params
    System(eqs, t, [], []; systems, name, defaults)
end

@named bot = CascadeControlledFlatDyadBotInput()
spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :bot,
    model = bot,
    measurement = "plant.y",
    control_input = "plant.u",
    step_input = "plant.u",
    step_output = "plant.y",
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
    loop_openings = ["plant.y2", "plant.u2"]
)

asol = JSC.run_analysis(spec)

# DyadControlSystems.launch_pid_autotuning_designer(spec)

arts = DyadControlSystems.artifacts(asol, :OptimizedParameters)

bot_params.inner_controller.kp = arts[1, :Kp_standard]
bot_params.inner_controller.ki = 1/arts[1, :Ti_standard]
bot_params.inner_controller.kd = arts[1, :Td_standard]

# Simulate Tilted Robot
bot_params.plant.theta_init = 0.75*π
sys = mtkcompile(bot)
prob = ODEProblem(sys, sys.plant => bot_params, (0, 1))
sol = solve(prob)

plot(sol; idxs=sys.plant.plant.theta)



@named bot = DyadBotComponents.CascadeControlledFlatDyadBot()

spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :bot,
    model = bot,
    measurement = "y2",
    control_input = "u2",
    step_input = "u2",
    step_output = "y2",
    Ts = 0.01,           # Sample time
    duration = 25.0,      # Simulation duration
    Ms = 1.5,            # Sensitivity peak constraint
    Mt = 1.5,            # Complementary sensitivity peak constraint
    Mks = 400.0,         # Control sensitivity constraint
    wl = 1e-2,           # Lower frequency bound
    wu = 1e3,            # Upper frequency bound
    kd_ub = 0.0,         # Tune PI controller
    num_frequencies = 200,
    soft = true
)