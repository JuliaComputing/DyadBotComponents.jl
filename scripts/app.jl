using DyadBotComponents
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitInputs
using ModelingToolkitParameters
using OrdinaryDiffEq
using SciMLBase

# --------------------------------------------------------------
# Precompiled Constants
# --------------------------------------------------------------
@named bot = DyadBotComponents.CascadeControlledFlatDyadBot()
bot_ns = ModelingToolkit.toggle_namespacing(bot, false)
inputs = [bot_ns.v_ref]
sys = mtkcompile(bot; inputs)
sys, input_functions = ModelingToolkitInputs.build_input_functions(sys, inputs)
bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()
prob = ODEProblem(sys, sys => bot_params, (0, Inf))
bot_setters = cache(sys, DyadBotComponents.CascadeControlledFlatDyadBotParams);



# --------------------------------------------------------------
# APP
# --------------------------------------------------------------
using ShoelaceWidgets
using Bonito
using WGLMakie

@kwdef struct AppState
    integrator = Ref{SciMLBase.DEIntegrator}()

    v_ref = Ref(0.0)

    setup = SLButton("setup")
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

    kpo = SLInput(0.54; label="outer P")
    kio = SLInput(2.48; label="outer i")
    kdo = SLInput(0.0; label="outer d")

    kpi = SLInput(15.6; label="inner P")
    kii = SLInput(1e12; label="inner i") #TODO: put to Inf
    kdi = SLInput(0.16; label="inner d")


    f = Figure()
    ax = Axis(f[1,1])

    bot_params = DyadBotComponents.CascadeControlledFlatDyadBotParams()

    base_transform = Ref{Makie.Transformation}()
    spoke_transform = Ref{Makie.Transformation}()
    body_transform = Ref{Makie.Transformation}()

end

left1(app::AppState, x) = app.v_ref[] = -0.5
left2(app::AppState, x) = app.v_ref[] = -1.0
left3(app::AppState, x) = app.v_ref[] = -1.5

pause(app::AppState, x) = app.v_ref[] = 0

right1(app::AppState, x) = app.v_ref[] = +0.5
right2(app::AppState, x) = app.v_ref[] = +1.0
right3(app::AppState, x) = app.v_ref[] = +1.5

function stop(app::AppState, x)

    app.loop[] = false
    wait(app.loop_thread[])
    app.setup.disabled[] = false
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

    # frame
    poly!(app.ax, Rect(-3/2, 0 , 3, 3); color=:transparent)    
end


function setup(app::AppState, x)

    app.bot_params.plant.R = app.R.value[]
    app.bot_params.plant.L = app.L.value[]
    
    app.bot_params.outer_controller.k = app.kpo.value[]
    app.bot_params.outer_controller.Ti = app.kio.value[]
    app.bot_params.outer_controller.Td = app.kdo.value[]

    app.bot_params.inner_controller.k = app.kpi.value[]
    app.bot_params.inner_controller.Ti = app.kii.value[]
    app.bot_params.inner_controller.Td = app.kdi.value[]

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
        set_input!(input_functions, app.integrator[], sys.v_ref, app.v_ref[])

        # step
        step!(app.integrator[], 0.1, true)

        # update plot
        x = app.integrator[][sys.plant.x]
        theta = app.integrator[][sys.plant.theta]

        translate!(app.base_transform[], x)
        rotate!(app.body_transform[], -theta + π)
        rotate!(app.spoke_transform[], -x/r) 
        
        # real time pause
        compute_time = Base.time() - t0
        sleep_time = max(step_size - compute_time, 0)
        sleep(sleep_time)
    end

    app.setup.disabled[] = true
    app.stop.disabled[] = false

end


function build_app!(app::AppState)

    on(Base.Fix1(left1, app), app.left1.value)
    on(Base.Fix1(left2, app), app.left2.value)
    on(Base.Fix1(left3, app), app.left3.value)
    on(Base.Fix1(pause, app), app.pause.value)
    on(Base.Fix1(right1, app), app.right1.value)
    on(Base.Fix1(right2, app), app.right2.value)
    on(Base.Fix1(right3, app), app.right3.value)

    on(Base.Fix1(setup, app), app.setup.value)
    on(Base.Fix1(stop, app), app.stop.value)
    
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
        )
    )
end



function get_body(session, app::AppState)
    DOM.body(

        DOM.div(
            DOM.h1("DyadBot App"),
            app.R,
            app.L,
            DOM.hr(), #---------------    
            app.kpo,
            app.kio,
            app.kdo,
            app.kpi,
            app.kii,
            app.kdi,
            DOM.hr(), #---------------
            DOM.div(app.setup, app.stop),
            DOM.hr(), #---------------
            DOM.div(app.left3, app.left2, app.left1, app.pause, app.right1, app.right2, app.right3),
            DOM.hr(), #---------------
            app.f
        )
        

    )
end

app = initialize()   
    

App() do session
    DOM.html(
        get_head(),
        get_body(session, app)
    )
end