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
prob = ODEProblem(sys, [sys.x_ref => 1.0], (0, 100))


#=
integrator = init(prob)
for i=1:5
    set_input!(input_functions, integrator, sys.x_ref, 1.0)
    step!(integrator, 20, true)
end
plot(integrator.sol; idxs=bot.plant.x)
=#






# --------------------------------------------------------------
# APP
# --------------------------------------------------------------
const STYLE_CSS = read(joinpath(@__DIR__, "style.css"), String)

@kwdef struct AppState
    integrator = Ref{SciMLBase.DEIntegrator}()

    x_ref = Ref(0.0)

    run = SLButton("run")
    stop = SLButton("stop"; disabled=true)

    loop = Ref(true)
    loop_thread = Ref{Task}()

    left1 = SLButton("<")
    left2 = SLButton("<<")
    left3 = SLButton("<<<")
    pause = SLButton("||")
    right1 = SLButton(">")
    right2 = SLButton(">>")
    right3 = SLButton(">>>")

    L = SLInput(0.5; label="length [m]")
    R = SLInput(0.1; label="wheel radius [m]")

    kpo = SLInput(0.54)
    kio = SLInput(2.48)
    kdo = SLInput(0.0)

    kpi = SLInput(15.6)
    kii = SLInput(1e12) #TODO: put to Inf
    kdi = SLInput(0.16)


    f = Figure()
    ax = Axis(f[1,1])

    bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()

    base_transform = Ref{Makie.Transformation}()
    spoke_transform = Ref{Makie.Transformation}()
    body_transform = Ref{Makie.Transformation}()


    time = Observable(0.0)
    footer_time = Observable("time: N/A")


end

left1(app::AppState, x) = app.x_ref[] = -0.5
left2(app::AppState, x) = app.x_ref[] = -1.0
left3(app::AppState, x) = app.x_ref[] = -1.5

pause(app::AppState, x) = app.x_ref[] = 0

right1(app::AppState, x) = app.x_ref[] = +0.5
right2(app::AppState, x) = app.x_ref[] = +1.0
right3(app::AppState, x) = app.x_ref[] = +1.5

time(app::AppState, x::Float64) = app.footer_time[] =  "time: $(round(x; digits=2)) [s]"

function stop(app::AppState, x)

    app.loop[] = false
    wait(app.loop_thread[])
    app.run.disabled[] = false
    app.stop.disabled[] = true

end


function build_fig!(app::AppState)

    empty!(app.ax)

    r = app.bot_params.plant.R #wheel radius
    y1 = r

    L = app.bot_params.plant.L # body length
    y2 = y1 + L

    w1 = r*0.5 
    w2 = r*2   

    # frame
    poly!(app.ax, Rect(-3/2, 0 , 3, 3); color=:transparent)  

    # wheel
    app.base_transform[] = Transformation(origin = Vec3d(0.0, y1, 0))
    app.spoke_transform[] = Transformation(app.base_transform[]; origin = Vec3d(0.0, y1, 0))
    linesegments!(app.ax, [-r,+r], [y1,y1]; transformation=app.spoke_transform[], color=:blue)
    linesegments!(app.ax, [0.0,0.0], [y1-r,y1+r]; transformation=app.spoke_transform[], color=:blue)
    poly!(app.ax, Circle(Point2f(0.0,y1), r); color=:transparent, strokecolor=:blue, strokewidth=2, transformation=app.spoke_transform[])

    # body
    app.body_transform[] = Transformation(app.base_transform[]; origin = Vec3d(0.0, y1, 0))
    poly!(app.ax, Rect(0.0-w1/2, y1, w1, y2-y1); color=:transparent, strokecolor=:red, strokewidth=2, transformation=app.body_transform[])
    poly!(app.ax, Rect(0.0-w2/2, y2, w2, w2); color=:transparent, strokecolor=:green, strokewidth=2, transformation=app.body_transform[])

  
end

parameter_control_map(app::AppState) = [
    (:plant, :R) => app.R
    (:plant, :L) => app.L
    (:outer_controller, :k) => app.kpo
    (:outer_controller, :Ti) => app.kio
    (:outer_controller, :Td) => app.kdo
    (:inner_controller, :k) => app.kpi
    (:inner_controller, :Ti) => app.kii
    (:inner_controller, :Td) => app.kdi
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

function run(app::AppState, x)

    # update bot_params from controls
    for (path, control) in parameter_control_map(app)
        setproperty!(app.bot_params, path, control.value[])
    end
    
    build_fig!(app)

    # Update Model Parameters
    # -- precompiled constants: prob, bot_setters, sys
    prob_ = remake(prob, bot_setters, sys => app.bot_params)

    # Initialize the integrator
    app.integrator[] = init(prob_, Rodas5P())

    step_size = 0.1
    app.loop[] = true

    r = app.bot_params.plant.R

    app.loop_thread[] = Threads.@spawn while app.loop[]
        t0 = Base.time()

        # input 
        set_input!(input_functions, app.integrator[], sys.x_ref, app.x_ref[])

        # step
        step!(app.integrator[], 0.1, true)
        app.time[] = app.integrator[].t

        # update plot
        x = app.integrator[][sys.plant.x]
        theta = app.integrator[][sys.plant.theta]

        WGLMakie.translate!(app.base_transform[], x)
        WGLMakie.rotate!(app.body_transform[], -theta + π)
        WGLMakie.rotate!(app.spoke_transform[], -x/r) 
        
        # real time pause
        compute_time = Base.time() - t0
        sleep_time = max(step_size - compute_time, 0)*0.95 
        sleep(sleep_time)
    end

    app.run.disabled[] = true
    app.stop.disabled[] = false

end


function build_app!(app::AppState)

    on(Base.Fix1(left1, app), app.left1.value; weak=false)
    on(Base.Fix1(left2, app), app.left2.value; weak=false)
    on(Base.Fix1(left3, app), app.left3.value; weak=false)
    on(Base.Fix1(pause, app), app.pause.value; weak=false)
    on(Base.Fix1(right1, app), app.right1.value; weak=false)
    on(Base.Fix1(right2, app), app.right2.value; weak=false)
    on(Base.Fix1(right3, app), app.right3.value; weak=false)

    on(Base.Fix1(run, app), app.run.value; weak=false)
    on(Base.Fix1(stop, app), app.stop.value; weak=false)

    on(Base.Fix1(time, app), app.time; weak=false)
    
end

function initialize()
    app = AppState()
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

function get_footer(app::AppState)
    DOM.div(
        sl_tag(app.footer_time; pill=true, primary=true),
        class="footer"
    )
end

function get_body(session, app::AppState)
    DOM.body(


        DOM.h1("DyadBot App"),
        DOM.div(
            sl_card(DOM.div(DOM.h3("Bot Parameters"); slot="header"), DOM.div(app.R, app.L)), 
            sl_card(
                DOM.div(DOM.h3("Controller Parameters"); slot="header"), 
                
                DOM.table(
                    DOM.tr(
                        DOM.td("-"), DOM.td("P"), DOM.td("I"), DOM.td("D")
                    ),
                    DOM.tr(
                        DOM.td("outer"), DOM.td(app.kpo), DOM.td(app.kio), DOM.td(app.kdo)
                    ),
                    DOM.tr(
                        DOM.td("inner"), DOM.td(app.kpi), DOM.td(app.kii), DOM.td(app.kdi)
                    )
                )
                
            )
        ),
        
        DOM.div(app.run, app.stop),
        DOM.hr(), #---------------
        DOM.div(app.left3, app.left2, app.left1, app.pause, app.right1, app.right2, app.right3),
        DOM.hr(), #---------------
        app.f,

        

        # footer --------------------------------------------
        get_footer(app)
        

    )
end

app = initialize()   
get_html(session) = DOM.html(
        get_head(),
        get_body(session, app)
    )

# App(get_html)


server = Bonito.Server("0.0.0.0", 9500)
route!(server, "/" => App(get_html))
Bonito.HTTPServer.openurl(Bonito.online_url(server, "/"));
