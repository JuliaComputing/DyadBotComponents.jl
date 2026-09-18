# Make synthetic release data for the plant calibration.
#
# A release test is what the robot can do on the floor: hold it at a tilt, let
# go, and let the cascade controller balance it in place. The position reference
# stays at 0, so the release angle is the only excitation.
#
# This script plays the part of the robot until real logs exist. It uses the
# measured masses and lengths and the robot's gains from
# measured_parameters.toml, and differs from the calibration model only in the
# three fitted parameters: the true body inertia is 1.6 times the uniform-rod
# estimate, the true wheel inertia is 2.5 times the solid-disc estimate, and the
# true friction is 3 times the guess. A calibration that works recovers them.
#
# The CSV columns are the ones calibrate_plant.jl reads:
#   timestamp     time, s
#   plant.theta   body tilt, rad   (the robot gets this from the Kalman filter)
#   plant.x       wheel-axis position, m  (the robot gets this from the encoders)
#
# Run:
#   julia +1.12 --project=. make_release_data.jl              # continuous controller
#   julia +1.12 --project=. make_release_data.jl discrete     # 200 Hz sampled controller

include(joinpath(@__DIR__, "common.jl"))
using CSV, DataFrames, Random

const TRUTH = (Ib = 1.6 * NOMINAL.Ib, Iw = 2.5 * NOMINAL.Iw, d = 3.0e-4)

# Release angles, rad. Both signs and three sizes, so the fit sees more than one
# operating point. The torque limit is 0.1 N m and gravity asks for
# M * g * L * sin(theta), so a release much beyond 0.13 rad saturates the motor.
const RELEASES = (0.10, -0.06, 0.04)

const T_END = 3.0     # s, long enough for the wheel to come back to 0
const TS = 0.005      # s, the controller period and the log rate

# Measurement noise, matching the robot: the tilt estimate is good to a few
# milliradian and the encoder position to well under a millimetre.
const SIGMA_THETA = 0.002   # rad
const SIGMA_X = 0.0005      # m

function make_release(phi0; discrete = false, seed = 1, noise = true)
    sys = build(; phi0, TRUTH..., discrete)
    sol = run_release(sys, T_END)
    ts = collect(0.0:TS:T_END)
    theta = sol(ts; idxs = sys.plant.theta).u
    x = sol(ts; idxs = sys.plant.x).u
    if noise
        rng = Xoshiro(seed)
        theta = theta .+ SIGMA_THETA .* randn(rng, length(ts))
        x = x .+ SIGMA_X .* randn(rng, length(ts))
    end
    println("release ", rpad(phi0, 6), " retcode ", sol.retcode,
            "  max |theta| ", round(maximum(abs, theta), digits = 4),
            "  max |x| ", round(maximum(abs, x), digits = 4))
    return DataFrame("timestamp" => ts, "plant.theta" => theta, "plant.x" => x)
end

function make_all(; discrete = false, noise = true)
    mkpath(DATA)
    tag = discrete ? "discrete" : "continuous"
    for (i, phi0) in enumerate(RELEASES)
        df = make_release(phi0; discrete, seed = i, noise)
        CSV.write(joinpath(DATA, "release_$(tag)_$(i).csv"), df)
    end
    println("wrote ", length(RELEASES), " releases to ", DATA)
    println("measured and fixed: ", MEAS)
    println("starting guesses:   ", NOMINAL)
    println("true parameters:    ", TRUTH)
end

if abspath(PROGRAM_FILE) == @__FILE__
    make_all(; discrete = length(ARGS) >= 1 && ARGS[1] == "discrete")
end
