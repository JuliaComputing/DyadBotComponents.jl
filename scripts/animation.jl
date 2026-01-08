using DyadBotComponents
using WGLMakie
using ModelingToolkit
using OrdinaryDiffEq

t_end = 100
@mtkcompile bot = DyadBotComponents.CascadeControlledFlatDyadBot()
prob = ODEProblem(bot, [bot.x_ref => 0.75], (0,t_end))
sol = solve(prob)

plot(sol; idxs=[bot.ref.output.u, bot.plant.x])
plot(sol; idxs=bot.plant.theta)


function get_figure()

    defs = ModelingToolkit.defaults(bot)

    r = defs[bot.plant.R] #wheel radius
    y1 = r

    L = defs[bot.plant.L] # body length
    y2 = y1 + L

    w1 = r*0.5 
    w2 = r*2


    f = Figure(size = (600, 600); aspect=DataAspect())
    ax = Axis(f[1, 1])
    

    # wheel
    base_transform = Transformation(origin = Vec3d(0.0, y1, 0))
    spoke_transform = Transformation(base_transform; origin = Vec3d(0.0, y1, 0))
    linesegments!(ax, [-r,+r], [y1,y1]; transformation=spoke_transform, color=:blue)
    linesegments!(ax, [0.0,0.0], [y1-r,y1+r]; transformation=spoke_transform, color=:blue)
    poly!(ax, Circle(Point2f(0.0,y1), r); color=:transparent, strokecolor=:blue, strokewidth=2, transformation=spoke_transform)

    # body
    body_transform = Transformation(base_transform; origin = Vec3d(0.0, y1, 0))
    poly!(ax, Rect(0.0-w1/2, y1, w1, y2-y1); color=:transparent, strokecolor=:red, strokewidth=2, transformation=body_transform)
    poly!(ax, Rect(0.0-w2/2, y2, w2, w2); color=:transparent, strokecolor=:green, strokewidth=2, transformation=body_transform)

    # frame
    poly!(ax, Rect(-3/2, 0 , 3, 3); color=:transparent)

    return f, base_transform, spoke_transform, body_transform
end

f, base_transform, spoke_transform, body_transform = get_figure()

function animate(;framerate = 30, filename="dyadbot.mp4")

    nframes = 500
    times = range(0, t_end, length=nframes)

    translate!(base_transform, 0.0)
    rotate!(body_transform, 0.0)
    rotate!(spoke_transform, 0.0)

    record(f, filename, times;
            framerate) do t
        x = sol(t; idxs=bot.plant.x)
        theta = sol(t; idxs=bot.plant.theta)

        # x = wheel*r
        
        translate!(base_transform, x)
        rotate!(body_transform, -theta + π)
        rotate!(spoke_transform, -x/r) 
    end

end

animate(; framerate=5, filename = "dyadbot_slow.mp4")
animate(; framerate=30, filename = "dyadbot_fast.mp4")

