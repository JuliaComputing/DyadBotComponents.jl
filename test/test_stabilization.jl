# Closed-loop stabilization and setpoint-tracking tests for every DyadBot model.
#
# Each model is simulated through the same `DyadInterface.TransientAnalysis` path
# the generated case tests use (which applies the SynchToolkit clock pass needed
# by the discrete models). For every model we check that
#   * the solver reports success,
#   * the controller keeps the robot upright — the body tilt angle stays bounded
#     (never falls over) and settles back to ~0, and
#   * for the cascade models, the wheel position tracks the filtered position
#     reference by the end of the run (the outer loop reaches its setpoint).

using DyadBotComponents
using DiscreteComponents   # loads SynchToolkit so the clocked models get the Lustre pass
using DyadInterface
using ModelingToolkit
using Test

# Simulate a model to `stop` and return (model, solution).
function simulate_model(constructor; stop)
    model = constructor(; name = :model)
    result = TransientAnalysis(; model, alg = DyadInterface.ODEAlg.Auto(),
        start = 0.0, stop, abstol = 1e-8, reltol = 1e-8,
        automatic_discontinuity_detection = false)
    return model, DyadInterface.rebuild_sol(result)
end

succeeded(sol) = Symbol(sol.retcode) == :Success

# Maximum absolute value of a signal over the whole trajectory.
maxabs(sol, acc) = maximum(abs, sol[acc])

# Maximum absolute value of `acc` sampled over the time window [t0, t1].
function maxabs_window(sol, acc, t0, t1)
    ts = range(t0, t1; length = 50)
    maximum(abs, sol(ts; idxs = acc).u)
end

# Maximum absolute tracking error x - r sampled over the window [t0, t1].
function maxabs_error_window(sol, x, r, t0, t1)
    ts = range(t0, t1; length = 50)
    maximum(abs, sol(ts; idxs = x).u .- sol(ts; idxs = r).u)
end

# The body tilt angle must never exceed this — well below "fallen over".
const THETA_BOUND = 0.2
# Settled tolerances, evaluated over the final second of the run.
const THETA_SETTLE = 0.02
const POS_TRACK_TOL = 0.02

@testset "Stabilization and setpoint tracking" begin
    # Angle-only models: the sole setpoint is the upright pose (theta = 0). The
    # wheel position is intentionally uncontrolled and free to drift.
    @testset "$(nameof(ctor))" for ctor in (
        DyadBotComponents.AngleControlledDyadBot,
        DyadBotComponents.DiscreteAngleControlledDyadBot,
        DyadBotComponents.AngleControlledDyadBot3D,
    )
        stop = 5.0
        m, sol = simulate_model(ctor; stop)
        @test succeeded(sol)
        @test maxabs(sol, m.plant.theta) < THETA_BOUND          # never falls over
        @test maxabs_window(sol, m.plant.theta, stop - 1, stop) < THETA_SETTLE  # settles upright
    end

    # Cascade models: the inner loop holds the robot upright while the outer loop
    # drives the wheel position to the filtered square-wave reference.
    @testset "$(nameof(ctor))" for (ctor, reffun) in (
        (DyadBotComponents.CascadeControlledDyadBot, m -> m.firstorder1.y),
        (DyadBotComponents.DiscreteCascadeControlledDyadBot, m -> m.firstorder1.y),
        (DyadBotComponents.CascadeFFDyadBot, m -> m.pos_ref.y),
        (DyadBotComponents.DiscreteCascadeFFDyadBot, m -> m.pos_ref.y),
    )
        stop = 20.0
        m, sol = simulate_model(ctor; stop)
        @test succeeded(sol)
        @test maxabs(sol, m.plant.theta) < THETA_BOUND          # inner loop keeps it upright
        @test maxabs_window(sol, m.plant.theta, stop - 1, stop) < THETA_SETTLE
        # Outer loop reaches its setpoint: position tracks the reference at the
        # end of the run (the reference has been constant for the last 5 s).
        @test maxabs_error_window(sol, m.plant.x, reffun(m), stop - 1, stop) < POS_TRACK_TOL
    end

    # The ideal-rolling 3D robot moves in the vertical plane only, so its
    # closed-loop response must be identical to the planar model's up to solver
    # tolerance (observed max deviation ~2e-9 at the tolerances used here).
    @testset "Planar/3D equivalence" begin
        stop = 5.0
        m2d, sol2d = simulate_model(DyadBotComponents.AngleControlledDyadBot; stop)
        m3d, sol3d = simulate_model(DyadBotComponents.AngleControlledDyadBot3D; stop)
        @test succeeded(sol2d)
        @test succeeded(sol3d)
        ts = range(0, stop; length = 501)
        for (sig2d, sig3d) in ((m2d.plant.theta, m3d.plant.theta), (m2d.plant.x, m3d.plant.x))
            dev = maximum(abs, sol2d(ts; idxs = sig2d).u .- sol3d(ts; idxs = sig3d).u)
            @test dev < 1e-6
        end
    end
end
