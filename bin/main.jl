using DyadBotComponents
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitInputs
using ModelingToolkitParameters
using OrdinaryDiffEq
using SciMLBase
using ShoelaceWidgets
using Bonito
using WGLMakie


# using ModelingToolkitStandardLibrary.Blocks

# params(Blocks.Add)
# params(Blocks.Gain)
# params(Blocks.Add3)
# params(Blocks.Limiter; defaults=(y_max=1, y_min=0))
# params(Blocks.Constant)
# params(Blocks.Integrator)
# params(Blocks.Derivative)
# params(DyadBotComponents.LimPID; defaults=(k=15.6, Ti=Inf, Td=0.16, Nd=25, u_max=7))
# params(DyadBotComponents.LimPID; defaults=(k=0.54, Ti=2.48, Td=0, Nd=600, wd=1, wp=0.5))

#=
[x] format the page
[x] add visual of the set point
[x] change the set point to a square cycle (target, cycle time, on/off)
[ ] add save functionality
[ ] publish the app on JuliaHub


[ ] get Fredrik's app up and running

[ ] add a scrolling plot
=#

struct FileDownload
    value::Observable{String}
end

function Bonito.jsrender(session::Session, x::FileDownload)

    dom = DOM.div()

    download = js"""
        function (value) {

            console.log('test');

            // Create a Blob with the content
            const blob = new Blob([value], { type: 'text/plain' });

            // Create a temporary download link
            const link = document.createElement('a');
            link.href = URL.createObjectURL(blob);
            link.download = `parameters.toml`;

            // Trigger the download
            link.click();

            // Clean up
            URL.revokeObjectURL(link.href);
        }
    """
    onjs(session, x.value, download)

    return Bonito.jsrender(session, dom)
end

# download = FileDownload("")
# app = App() do session
#     DOM.html(
#         DOM.head(
#             get_shoelace()...
#         ),
#         DOM.body(
#             DOM.h1("Hello World"),
#             download
#         )
#     )
# end
# download.value[]="Hello World"

struct FileInput
    value::Observable{String}
    click::Observable{Bool} 
    label::String
    id::String
end

FileInput(label::String, id::String) = FileInput(Observable(""), Observable(false), label, id)

function Bonito.jsrender(session::Session, x::FileInput)

    button = ShoelaceWidgets.sl_button(x.label)
    input = DOM.input(type="file", accept=".toml", id=x.id, style="display: none;")

    change = js"""
    function onload(element) {
        function onchange(e) {

            const file = element.files[0];

            if (!file) {
                alert('Please select a file first');
                return;
            }

            const reader = new FileReader();

            reader.onload = function(ee) {
                const contents = ee.target.result;
                const filename = file.name;
                // Update the observable with file contents
                $(x.value).notify(contents);
            };

            reader.onerror = function(ee) {
                alert('Error reading file');
            };

            reader.readAsText(file);

        }
        element.addEventListener(`change`, onchange);
    }
    """
    Bonito.onload(session, input, change)

    click = js"""
    function onload(element) {
      function click() {
        document.getElementById($(x.id)).click();
      }
      element.addEventListener('click', click);
    }
    """
    Bonito.onload(session, button, click)

    dom = DOM.div(
        input,
        button
    )
    return Bonito.jsrender(session, dom)
end

# --------------------------------------------------------------
# Precompiled Constants
# --------------------------------------------------------------
@named bot = DyadBotComponents.CascadeControlledFlatDyadBot()
bot_ns = ModelingToolkit.toggle_namespacing(bot, false)
inputs = [bot_ns.x_ref]
sys = mtkcompile(bot; inputs)
sys, input_functions = ModelingToolkitInputs.build_input_functions(sys, inputs)
bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()
bot_setters = ModelingToolkitParameters.cache(sys, DyadBotComponents.CascadeControlledFlatDyadBotParams)
prob = ODEProblem(sys, sys => bot_params, (0, 100))

#=
integrator = init(prob)
integrator.ps[getproperty(sys, (:plant, :b_rot))]
for i=1:5
    set_input!(input_functions, integrator, sys.x_ref, 1.0)
    step!(integrator, 20, true)
end
plot(integrator.sol; idxs=bot.plant.x)
=#




# --------------------------------------------------------------
# APP
# --------------------------------------------------------------
const STYLE_CSS = """
    * {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
    }

    body {
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    }

    table td {
        text-align: center;
        vertical-align: middle;
    }

    /* Top Banner */
    .top-banner {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 1.5rem 2rem;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
    flex-shrink: 0;
    }

    .top-banner h1 {
    margin: 0;
    font-size: 1.5rem;
    }

    /* Content Area (scrollable) */
    .content {
    flex: 1;
    overflow-y: auto;
    padding: 2rem;
    background: #f5f5f5;
    }

    /* Content Layout */
    .content-container {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
    height: 100%;
    }

    /* Row 1: Two cards side by side */
    .row-1 {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.5rem;
    }

    .row-1 sl-card {
    height: 100%;
    }

    /* Row 2: 3 spots */
    .row-2 {
    display: grid;
    grid-template-columns: 1fr 1fr 1fr;
    gap: 1.5rem;
    }

    .row-2 sl-card {
    height: 100%;
    }

    .row-3 {
    display: flex;
    }

    .row-3 sl-card {
    flex: 1;
    }

    /* Row 4: Fill remaining space */
    .row-4 {
    flex: 1;
    background: white;
    border: 2px #ccc;
    border-radius: 8px;
    padding: 2rem;
    display: flex;
    align-items: center;
    justify-content: center;
    color: #666;
    font-size: 1.2rem;
    }

    /* Bottom Footer */
    .bottom-footer {
    background: #2c3e50;
    color: white;
    padding: 1rem 2rem;
    box-shadow: 0 -2px 8px rgba(0, 0, 0, 0.1);
    flex-shrink: 0;
    }

    .bottom-footer p {
    margin: 0;
    }

    /* Card content styling */
    sl-card::part(body) {
    padding: 1.5rem;
    }

    sl-card::part(header) {
    font-weight: 600;
    font-size: 1.1rem;
    }

    sl-button,
    sl-select {
    margin: 2px;
    }

    /* input[type="file"]::file-selector-button {
        font-family: var(--sl-font-sans);
        font-size: var(--sl-button-font-size-medium);
        padding: var(--sl-spacing-small) var(--sl-spacing-medium);
        background-color: var(--sl-color-primary-600);
        color: var(--sl-color-neutral-0);
        border: none;
        border-radius: var(--sl-border-radius-medium);
        cursor: pointer;
        font-weight: var(--sl-font-weight-semibold);
        margin-right: var(--sl-spacing-small);
    }

    input[type="file"]::file-selector-button:hover {
        background-color: var(--sl-color-primary-500);
    } */
"""

@kwdef struct AppState
    integrator = Ref{SciMLBase.DEIntegrator}()

    x_ref = SLInput(0.5; label="target position [m]")
    cycle_time = SLInput(5.0; label="cycle time [s]")
    cycle_amplitude = SLInput(0.5; label="cycle amplitude [m]")
    cycle = SLCheckbox("enable cycle"; checked=true)

    run = SLButton("Run")
    stop = SLButton("Stop"; disabled=true)
    update = SLButton("update") #TODO: remove, not needed

    loop = Ref(true)
    loop_thread = Ref{Task}()

    # left1 = SLButton("<")
    # left2 = SLButton("<<")
    # left3 = SLButton("<<<")
    # pause = SLButton("||")
    # right1 = SLButton(">")
    # right2 = SLButton(">>")
    # right3 = SLButton(">>>")

    
    


    L = SLInput(0.0; label="length [m]")
    R = SLInput(0.0; label="wheel radius [m]")
    b_trans = SLInput(0.0; label="damping (translational) [N-s/m]")
    b_rot = SLInput(0.0; label="damping (rotational) [N-s/rad]")

    kpo = SLInput(0.0)
    kio = SLInput(0.0)
    kdo = SLInput(0.0)

    kpi = SLInput(0.0)
    kii = SLInput(0.0) #TODO: put to Inf
    kdi = SLInput(0.0)




    f1 = Figure()
    ax1 = Axis(f1[1,1])

    f2 = Figure()
    ax2 = Axis(f2[1,1])

    bot_params = Ref(DyadBotComponents.CascadeControlledFlatDyadBotParams())

    base_transform = Ref{Makie.Transformation}()
    spoke_transform = Ref{Makie.Transformation}()
    body_transform = Ref{Makie.Transformation}()


    time = Observable(0.0)
    footer_time = Observable("time: N/A")

    N = round(Int,30/0.1)
    past_position = Observable(zeros(N))
    past_target = Observable(zeros(N))

    fileinput = FileInput("Load Parameters","parameters_input")
    download = FileDownload("")
    filedownload = SLButton("Save Parameters")

    file = SLTextarea("no file loaded")

end

cycle(app::AppState, x::Bool) = app.x_ref.disabled[] = x
time(app::AppState, x::Float64) = app.footer_time[] =  "time: $(round(x; digits=2)) [s]"

function stop(app::AppState, x)

    app.loop[] = false
    wait(app.loop_thread[])
    app.run.disabled[] = false
    app.stop.disabled[] = true
    app.L.disabled[] = false
    app.R.disabled[] = false

end

function filedownload(app::AppState, x)
    params = ModelingToolkitParameters.parameters_to_string(app.bot_params[])
    app.download.value[] = params
end

function fileinput(app::AppState, x)
    app.file.value[] = x
    app.bot_params[] = ModelingToolkitParameters.string_to_parameters(x, DyadBotComponents.CascadeControlledFlatDyadBotParams)
    params_to_control(app)
end


function build_fig!(app::AppState)

    empty!(app.ax1)
    empty!(app.ax2)

    r = app.bot_params[].plant.R #wheel radius
    y1 = r

    L = app.bot_params[].plant.L # body length
    y2 = y1 + L

    w1 = r*0.5 
    w2 = r*2   

    # frame
    poly!(app.ax1, Rect(-3/2, 0 , 3, 3); color=:transparent)  

    # wheel
    app.base_transform[] = Transformation(origin = Vec3d(0.0, y1, 0))
    app.spoke_transform[] = Transformation(app.base_transform[]; origin = Vec3d(0.0, y1, 0))
    linesegments!(app.ax1, [-r,+r], [y1,y1]; transformation=app.spoke_transform[], color=:blue)
    linesegments!(app.ax1, [0.0,0.0], [y1-r,y1+r]; transformation=app.spoke_transform[], color=:blue)
    poly!(app.ax1, Circle(Point2f(0.0,y1), r); color=:transparent, strokecolor=:blue, strokewidth=2, transformation=app.spoke_transform[])

    # body
    app.body_transform[] = Transformation(app.base_transform[]; origin = Vec3d(0.0, y1, 0))
    poly!(app.ax1, Rect(0.0-w1/2, y1, w1, y2-y1); color=:transparent, strokecolor=:red, strokewidth=2, transformation=app.body_transform[])
    poly!(app.ax1, Rect(0.0-w2/2, y2, w2, w2); color=:transparent, strokecolor=:green, strokewidth=2, transformation=app.body_transform[])

    # set point
    vlines!(app.ax1, app.x_ref.value)




    time = 0:0.1:(app.N-1)*(0.1)

    lines!(app.ax2, time, app.past_position; label="position")
    lines!(app.ax2, time, app.past_target; label="target")
    axislegend(app.ax2)

  
end

parameter_control_map(app::AppState) = [
    (:plant, :R) => app.R
    (:plant, :L) => app.L
    (:plant, :b_rot) => app.b_rot
    (:plant, :b_trans) => app.b_trans
    (:outer_controller, :kp) => app.kpo
    (:outer_controller, :ki) => app.kio
    (:outer_controller, :kd) => app.kdo
    (:inner_controller, :kp) => app.kpi
    (:inner_controller, :ki) => app.kii
    (:inner_controller, :kd) => app.kdi
    
]

function Base.getproperty(value, names::NTuple{N, Symbol}) where N
    y = value
    for i=1:N-1
        y = getproperty(y, names[i])
    end

    return getproperty(y, names[N])
end

function Base.setproperty!(value, names::NTuple{N, Symbol}, x) where N
    y = value
    for i=1:N-1
        y = getproperty(y, names[i])
    end

    return setproperty!(y, names[N], x)
end

function control_to_params(app::AppState)
    # update bot_params from controls
    for (path, control) in parameter_control_map(app)
        setproperty!(app.bot_params[], path, control.value[])
    end
end

function params_to_control(app::AppState)
    for (path, control) in parameter_control_map(app)
        control.value[] = getproperty(app.bot_params[], path)
    end
end



function update_ps(app::AppState, path, x)

    if isdefined(app.integrator, 1)
        app.integrator[].ps[getproperty(sys, path)] = x
    end

end

function update_vec!(y::Observable{Vector{Float64}}, x::Float64)
    popat!(y[], 1)
    push!(y[], x)
    notify(y) 
end

function run(app::AppState, x)

    control_to_params(app) #updates app.bot_params
    
    build_fig!(app)

    # Update Model Parameters
    # -- precompiled constants: prob, bot_setters, sys
    prob_ = remake(prob, bot_setters, sys => app.bot_params[])

    # Initialize the integrator
    app.integrator[] = init(prob_, Rodas5P())

    step_size = 0.1
    app.loop[] = true

    r = app.bot_params[].plant.R

    cycle_current = 0

    app.loop_thread[] = Threads.@spawn while app.loop[]
        t0 = Base.time()

        # input 
        if app.cycle.value[]
            if cycle_current > app.cycle_time.value[]
                app.x_ref.value[] = -app.cycle_amplitude.value[] * sign(app.x_ref.value[])
                cycle_current = 0
            end
        end
        set_input!(input_functions, app.integrator[], sys.x_ref, app.x_ref.value[])

        # step
        step!(app.integrator[], 0.1, true)
        app.time[] = app.integrator[].t
        cycle_current += 0.1

        # update plot
        x = app.integrator[][sys.plant.x]
        theta = app.integrator[][sys.plant.theta]

        WGLMakie.translate!(app.base_transform[], x)
        WGLMakie.rotate!(app.body_transform[], -theta + π)
        WGLMakie.rotate!(app.spoke_transform[], -x/r) 

        update_vec!(app.past_position, x)
        update_vec!(app.past_target, app.x_ref.value[])
        
        # real time pause
        compute_time = Base.time() - t0
        sleep_time = max(step_size - compute_time, 0)*0.95 
        sleep(sleep_time)
    end

    app.run.disabled[] = true
    app.stop.disabled[] = false
    app.L.disabled[] = true
    app.R.disabled[] = true

end


function build_app!(app::AppState)

    on(Base.Fix1(run, app), app.run.value; weak=false)
    on(Base.Fix1(stop, app), app.stop.value; weak=false)

    on(Base.Fix1(time, app), app.time; weak=false)

    on(Base.Fix1(cycle, app), app.cycle.value; weak=false)

    for (path, control) in parameter_control_map(app)
       on(x->update_ps(app, path, x), control.value; weak=false)
    end

    on(Base.Fix1(fileinput, app), app.fileinput.value; weak=false)
    on(Base.Fix1(filedownload, app), app.filedownload.value; weak=false)
end

function initialize()
    app = AppState()

    #sync controls to the parameters
    params_to_control(app)

    build_app!(app)
    return app
end

function get_head()
    DOM.head(
        DOM.title("JuliaHub | DyadBot"),
        DOM.head(
            get_shoelace()...
        ),
        DOM.style(STYLE_CSS)
    )
end

function get_body(session, app::AppState)

    DOM.body(

        DOM.div(
            DOM.h1("DyadBot App"); 
            class="top-banner"
        ),

        DOM.div(
            DOM.div(

                DOM.div(
                    sl_card(
                        DOM.div("Bot Parameters"; slot="header"), 
                        DOM.div(app.R, app.L, app.b_rot, app.b_trans)
                    ), 
                    sl_card(
                        DOM.div("Controller Parameters"; slot="header"), 
                        
                        DOM.table(
                            DOM.tr(
                                DOM.td(""), DOM.td("P"), DOM.td("I"), DOM.td("D")
                            ),
                            DOM.tr(
                                DOM.td("outer"), DOM.td(app.kpo), DOM.td(app.kio), DOM.td(app.kdo)
                            ),
                            DOM.tr(
                                DOM.td("inner"), DOM.td(app.kpi), DOM.td(app.kii), DOM.td(app.kdi)
                            )
                        )
                        
                    )
                    ;class="row-1"
                ),
                DOM.div(
                    sl_card(
                        DOM.div("Control"; slot="header"), 
                        DOM.div(app.run, app.stop),
                        app.x_ref,
                        app.cycle_time,
                        app.cycle_amplitude,
                        app.cycle
                    ),
                    app.f1,
                    app.f2
                    ;class="row-2"
                ),
                DOM.div(
                    sl_card(
                        DOM.div("File"; slot="header"), 
                        DOM.div(app.fileinput, app.filedownload),
                        # app.file, #file display
                        app.download #hidden component
                    )
                    ;class="row-3"
                ),
                DOM.div(
                    
                    ;class="row-4"
                )

                ;class="content-container"
            )
            ;class="content"
        ),

        DOM.div(
            sl_tag(app.footer_time; pill=true, primary=true),
            ;class="bottom-footer"
        ),       

    )
end

# app = initialize()   
# get_html(session) = DOM.html(
#         get_head(),
#         get_body(session, app)
#     )

# # App(get_html)

function get_html_app_session(session)
    app = initialize()   
    DOM.html(
            get_head(),
            get_body(session, app)
        )
end


# ============================================================================
# publish
# ============================================================================

proxy = get(ENV, "JULIAHUB_APP_URL", "")
if isempty(proxy)
    @info "No Bonito proxy found in environment variable JULIAHUB_APP_URL"
else
    @info "Using Bonito proxy from JULIAHUB_APP_URL: $proxy"
end
port = get(ENV, "PORT", "8080") # it's guaranteed this exists on JuliaHub
@info "Constructing Bonito server on 0.0.0.0:$port $(isempty(proxy) ? "" : "with proxy $proxy")"
server = Bonito.Server("0.0.0.0", parse(Int, port); proxy_url=proxy, verbose=-1)
route!(server, "/" => App(get_html_app_session))
# Bonito.HTTPServer.openurl(Bonito.online_url(server, "/"));

@info "Starting Bonito server"
Bonito.HTTPServer.start(server)

@info "Server successfully started, waiting on connections"

# Wait for the server to exit, because if running in an app, the app will
# exit when the script is done.  This makes sure that the app is only closed
# if (a) the server closes, or (b) the app itself times out and is killed externally.
wait(server)

#=
using DyadControlSystems
import DyadControlSystems as JSC

spec = JSC.PIDAutotuningAnalysisSpec(;
    name = :CascadeTuning,
    model = bot,
    measurement = "y",
    control_input = "u",
    step_input = "u",
    step_output = "y",
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
)

DyadControlSystems.launch_pid_autotuning_designer(spec)
=#