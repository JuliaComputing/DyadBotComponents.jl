# DyadBotComponents

Dyad models of a two-wheeled balancing robot, together with controller
architectures and controller-tuning scripts.

## Models (`dyad/`)

Plant models:

- `PlanarDyadBot`: planar model of the robot. A wheel rolling without slip,
  an inverted-pendulum body mounted on the wheel axis through a motor
  (`SimpleMotor`), and an IMU-like sensor. Inputs/outputs: motor torque in,
  position/velocity/tilt/tilt-rate out.
- `DyadBot3D`: three-dimensional model with two individually spinning wheels
  with slip-based ground contact, allowing the body to tilt.
- `RollingDyadBot3D`: three-dimensional model with an ideal no-slip rolling
  wheel axle (`MultibodyComponents.RollingWheelSet`) and individually driven
  wheels, so the robot can be steered by differential torque. Outputs
  odometric path position/velocity, tilt angle/rate and yaw (heading)
  angle/rate. With equal torque on both wheels it is dynamically equivalent
  to `PlanarDyadBot`.

Closed-loop models around the planar plant:

- `AngleControlledDyadBot`: single PID loop stabilizing the tilt angle.
- `CascadeControlledDyadBot`: cascade control with an inner tilt-angle loop
  and an outer position loop; all six controller gains are exposed as
  tunable top-level parameters.
- `CascadeFFDyadBot`: the cascade extended with a feedforward generator
  (state-space system loaded from `data/ff_*.csv`) providing filtered
  position reference, tilt-angle feedforward and torque feedforward.
- `LQGControlledDyadBot`: LQG controller with reference feedforward and
  integral action (state-space system loaded from `data/lqg_*.csv`).
- `LQGTuningDyadBot`: analysis model used by the LQG design script.

Closed-loop models around the 3D plant:

- `AngleControlledDyadBot3D`: `RollingDyadBot3D` stabilized by the same
  `AngleController` as `AngleControlledDyadBot`; with in-plane motion the
  closed-loop response is identical to the planar model
  (`test/test_stabilization.jl` verifies this).
- `YawControlledDyadBot3D`: `RollingDyadBot3D` with the balance/position
  cascade (`CascadeController`, position reference held at zero) and a
  separate `YawController` tracking a heading reference that steps at t = 5;
  a `ControlMixer` combines drive and yaw torque into the two wheel torque
  commands, so the robot spins in place to the new heading.

## Controller tuning scripts (`scripts/`)

Each script activates the `scripts/` environment and can be run directly:

- `tune_angle_pid.jl`: PID autotuning of the inner angle loop
  (`PIDAutotuningAnalysisSpec`).
- `tune_cascade_pid.jl`: PID autotuning of the outer position loop, plus
  robustness analysis (sensitivity functions, disk margins).
- `tune_cascade_structured.jl`: structured tuning of all six cascade gains
  simultaneously (`StructuredAutoTuningProblem`) with per-loop sensitivity
  bounds and a pole-location constraint.
- `compute_feedforward.jl`: computes the feedforward generator for
  `CascadeFFDyadBot` and stores it in `data/`.
- `tune_lqg.jl`: LQG design (`LQGAnalysisSpec`) producing the controller
  used by `LQGControlledDyadBot`, stored in `data/`.

## Julia utilities (`src/`)

- `IMUKalmanFilter`: discrete-time Kalman filter for tilt estimation from
  gyro and accelerometer measurements, suitable for embedded use.
- `balance_original.jl`: standalone Julia implementation of the original
  balance controller (not part of the package module).
