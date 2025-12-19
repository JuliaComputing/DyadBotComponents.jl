using DyadBotComponents
using WGLMakie
using ModelingToolkit
using OrdinaryDiffEq

f = Figure(size = (600, 600); aspect=DataAspect())

ax = Axis(f[1, 1])

#TODO: connect to the bot parameters
y1 = 1.0
r = 1.0

y2 = 3.0
w = 2


base_transform = Transformation(origin = Vec3d(0.0, y1, 0))
spoke_transform = Transformation(base_transform; origin = Vec3d(0.0, y1, 0))
linesegments!(ax, [-r,+r], [y1,y1]; transformation=spoke_transform, color=:blue)
linesegments!(ax, [0.0,0.0], [y1-r,y1+r]; transformation=spoke_transform, color=:blue)
poly!(ax, Circle(Point2f(0.0,y1), r); color=:transparent, strokecolor=:blue, strokewidth=2, transformation=spoke_transform)

body_transform = Transformation(base_transform; origin = Vec3d(0.0, y1, 0))
poly!(ax, Rect(0.0-0.25/2, y1, 0.25, y2-y1); color=:transparent, strokecolor=:red, strokewidth=2, transformation=body_transform)
poly!(ax, Rect(0.0-w/2, y2, w, w); color=:transparent, strokecolor=:green, strokewidth=2, transformation=body_transform)

poly!(ax, Rect(-2.5, 0 , 5, 5); color=:transparent)

rotate!(spoke_transform, π*0.1)
rotate!(spoke_transform, π*0.2)
rotate!(spoke_transform, π*0.3)

translate!(base_transform, 0.5)

rotate!(body_transform, -π*0.1)



@mtkcompile bot = DyadBotComponents.CascadeControlledFlatDyadBot()
prob = ODEProblem(bot, [bot.x_ref => 0.75], (0,200))
sol = solve(prob)

plot(sol; idxs=bot.plant.x)
plot(sol; idxs=bot.plant.theta)
plot(sol; idxs=bot.ref.output.u)

nframes = 200
framerate = 30
times = range(0, 200, length=nframes)

translate!(base_transform, 0.0)
rotate!(body_transform, 0.0)
rotate!(spoke_transform, 0.0)

record(f, "dyadbot.mp4", times;
        framerate = 10) do t
    x = sol(t; idxs=bot.plant.x)
    theta = sol(t; idxs=bot.plant.theta)

    # x = wheel*r*π
    
    translate!(base_transform, x)
    rotate!(body_transform, -theta + π)
    rotate!(spoke_transform, -x/(r*π))
     
end