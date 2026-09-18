# Plant calibration

Fits the plant parameters of the balancing robot to release-test data with
[DyadModelOptimizer](https://help.juliahub.com/jsmo/stable/).

The model is `BalanceDyadBot` (`dyad/balance.dyad`): the cascade-controlled
planar robot with a **position reference of 0**, which is what the robot does on
the floor. `CascadeControlledDyadBot` tracks a square wave instead, which is a
tracking test and not how the robot is used.

## The split

| Parameter | Meaning | How it is obtained |
| --- | --- | --- |
| `M` | body mass | weighed, `measured_parameters.toml` |
| `m` | wheel mass, both wheels | weighed, `measured_parameters.toml` |
| `R` | effective rolling radius | rolled over a measured distance, `measured_parameters.toml` |
| `L` | wheel axis to body centre of mass | balanced on an edge, `measured_parameters.toml` |
| `Ib` | body moment of inertia | **fitted** |
| `Iw` | wheel + motor moment of inertia | **fitted** |
| `d` | motor viscous friction | **fitted** |

A scale and a ruler give the first four. The last three are the ones they cannot
give, and they set the natural frequency and the damping of the robot. Each spans
decades, so the search runs on `log10` of each.

## The experiment

A release test: hold the robot at a tilt, let go, and let the controller balance
it in place. The position reference stays at 0, so the release angle is the only
excitation. Three releases of different size and sign make three `Experiment`s
that share one search space, so all three fit the same plant.

## Running it

```
cd calibration
julia +1.12 --project=. make_release_data.jl      # synthetic robot -> data/
julia +1.12 --project=. calibrate_plant.jl        # fit -> data/calibrated_parameters_continuous.toml
```

`make_release_data.jl` stands in for the robot until real logs exist. It
simulates the model with plant parameters that differ from the defaults and adds
measurement noise, so the calibration can be checked against a known answer.

To use a real log instead, write one CSV per release into `data/` named
`release_continuous_<n>.csv` with these columns:

| Column | Meaning |
| --- | --- |
| `timestamp` | time, s |
| `plant.theta` | body tilt, rad, from the Kalman filter |
| `plant.x` | wheel-axis position, m, from the encoders |

Check that the release is hands-off before you use it. A stretch of a log where
the operator is holding the robot looks like good balancing and is not.

## The continuous approximation

The calibration model uses the continuous `CascadeController`. The robot runs the
200 Hz sampled `DiscreteCascadeController`, whose period is about 40 times shorter
than the closed-loop response, so the continuous model is a close stand-in and it
keeps the optimizer away from clocked solves.

To see what that approximation costs, make a discrete data set and fit it with
the same continuous model:

```
julia +1.12 --project=. make_release_data.jl discrete
julia +1.12 --project=. calibrate_plant.jl discrete
```

The gap between the two fitted parameter sets is the bias from the approximation.
