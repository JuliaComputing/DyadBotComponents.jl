# Compute a feedforward generator for the cascade-controlled balancing robot
# and store its state-space matrices as CSV files in ../data/.
#
# The feedforward generator drives the closed-loop system to follow a
# reference model Tr (a fourth-order low-pass filter). It has 1 input
# (position reference) and 3 outputs:
#   y[1]: angle feedforward
#   y[2]: filtered position reference
#   y[3]: torque feedforward
#
# The stored matrices are loaded by the component CascadeFFDyadBot through
# the loader functions ff_A(), ff_B(), ff_C(), ff_D(), ff_nx() defined in
# dyad/module.jl.

cd(@__DIR__)
using Pkg
Pkg.activate(".")
using DyadBotComponents
using MultibodyComponents
using ModelingToolkit
using OrdinaryDiffEq
using DyadControlSystems, ControlSystemsBase
using DelimitedFiles
using Plots

# Build model with analysis points. Zero initial tilt, so that the model's
# initial equations are consistent with the upright trim used for
# linearization.
@named model = DyadBotComponents.CascadeControlledDyadBot(phi0 = 0)
ssys = multibody(model)

# Reference model: 4th-order low-pass filter
T0 = ss(tf([1], [0.1, 1]))
Tr = T0^4

# Initial conditions are codegened to initial equations by dyad. Initial
# equations are ignored by linearization, so the upright equilibrium is
# specified explicitly here.
op = Dict([
    ssys.square.amplitude => 0
    ssys.angle_controller.u_m => 0
    ssys.pos_controller.u_m => 0
    ssys.plant.body_mass.body.phi => 0
    ssys.plant.body_mass.body.w => 0
    ssys.plant.wheelinertia.phi => 0
    ssys.plant.wheelinertia.w => 0
    ssys.firstorder.x => 0
    ssys.firstorder1.x => 0
])

Ryur = DyadControlSystems.feedforward_generator(model;
    Tr,
    measurement = [model.y, model.y2],
    controlled_output = [model.y2],
    control_input = [model.u],
    op,
    MultibodyComponents.linsys...
)

# feedforward_generator orders the outputs [control_input; measurement...];
# permute to the order assumed by CascadeFFDyadBot: [angle ff, position ref,
# torque ff], identified by the analysis-point names.
outnames = string.(Ryur.y)
perm = [findfirst(endswith("₊y"), outnames),
    findfirst(endswith("₊y2"), outnames),
    findfirst(endswith("₊u"), outnames)]
@assert !any(isnothing, perm) "unexpected generator output names: $outnames"
Ryur = Ryur[perm, :]

Ryur = balance_statespace(Ryur)[1]

bodeplot(Ryur, size = (800, 800))

# Save matrices to CSV files
dir = joinpath(@__DIR__, "..", "data")
mkpath(dir)
A, B, C, D = ssdata(Ryur)
writedlm(joinpath(dir, "ff_A.csv"), A, ',')
writedlm(joinpath(dir, "ff_B.csv"), B, ',')
writedlm(joinpath(dir, "ff_C.csv"), C, ',')
writedlm(joinpath(dir, "ff_D.csv"), D, ',')
open(joinpath(dir, "ff_nx.txt"), "w") do io
    println(io, size(A, 1))
end
