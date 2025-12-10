using StaticArrays

# System: x = [angle; bias], dynamics: ẋ = A_c*x + B_c*u, measurement: y = C*x
# Continuous: A_c = [0 -1; 0 0], B_c = [1; 0]
# Discretized with Ts = 0.005 (200 Hz)

const Ts = 0.005f0
const Ad = @SMatrix [1.0f0 -Ts; 0.0f0 1.0f0]   # exp(A_c * Ts), exact since A_c² = 0
const Bd = @SVector [Ts, 0.0f0]             # ∫₀ᵈᵗ exp(A_c*τ) B_c dτ
const Cd = @SMatrix [1.0f0 0.0f0]             # Measurement matrix (1×2)
const I2 = @SMatrix [1.0f0 0.0f0; 0.0f0 1.0f0]    # Identity for covariance update

"""
    DiscreteKalmanFilter

Discrete-time Kalman filter for angle estimation from gyro (input) and accelerometer (measurement).
State is `[angle, gyro_bias]`.
"""
mutable struct DiscreteKalmanFilter
    x::SVector{2, Float32}              # State [angle, bias]
    P::SMatrix{2, 2, Float32, 4}        # Covariance
    Q::SMatrix{2, 2, Float32, 4}        # Process noise covariance
    R::SMatrix{1, 1, Float32, 1}        # Measurement noise covariance
end

"""
    DiscreteKalmanFilter(; Q_angle=0.001f0, Q_bias=0.005f0, R_angle=0.5)

Create a discrete Kalman filter with specified noise covariances.
Default values from BalanceCar.h.
"""
function DiscreteKalmanFilter(; Q_angle::Float32=0.001f0, Q_bias::Float32=0.005f0, R_angle::Float32=0.5f0)
    x = @SVector zeros(Float32, 2)
    P = I2
    Q = @SMatrix [Q_angle 0.0f0; 0.0f0 Q_bias]
    R = @SMatrix [R_angle;;]
    DiscreteKalmanFilter(x, P, Q, R)
end

"""
    predict!(kf, u)

Prediction step: propagate state and covariance using gyro measurement `u` as input.

    x = A * x + B * u
    P = A * P * A' + Q
"""
function predict!(kf::DiscreteKalmanFilter, u::Float32)
    kf.x = Ad * kf.x + Bd * u
    kf.P = Ad * kf.P * Ad' + kf.Q
    nothing
end

"""
    correct!(kf, y)

Update step: correct state and covariance using accelerometer angle measurement `y`.

    S = C * P * C' + R
    K = P * C' / S
    x = x + K * (y - C * x)
    P = (I - K * C) * P
"""
function correct!(kf::DiscreteKalmanFilter, y::Float32)
    # Innovation covariance (1×1 matrix)
    S = Cd * kf.P * Cd' + kf.R

    # Kalman gain (2×1)
    K = kf.P * Cd' / S[1, 1]

    # Innovation (scalar wrapped as 1×1 for matrix multiply)
    innovation = SA[y] - Cd * kf.x

    # State update
    kf.x = kf.x + K * innovation

    # Covariance update
    kf.P = (I2 - K * Cd) * kf.P

    nothing
end

"""
    update!(kf, u, y)

Combined predict + correct step.

- `u`: gyro measurement (angular velocity)
- `y`: accelerometer angle measurement
"""
function update!(kf::DiscreteKalmanFilter, u::Float32, y::Float32)
    predict!(kf, u)
    correct!(kf, y)
    nothing
end

"""
    angle(kf)

Get the estimated angle from the filter state.
"""
angle(kf::DiscreteKalmanFilter) = kf.x[1]

"""
    bias(kf)

Get the estimated gyro bias from the filter state.
"""
bias(kf::DiscreteKalmanFilter) = kf.x[2]

# IMU calibration constants (from original C++ code)
const GYRO_OFFSET = 128.1f0
const GYRO_SCALE = 131.0f0

"""
    angle, gyro_x, gyro_z = compute_angles(kf, ax, ay, az, gx, gy, gz)

Apply calibration to raw IMU readings and compute
- `angle`: Tilt angle from accelerometer (radians)
- `gyro_x`: Calibrated gyro x (angular velocity) for tile derivative control
- `gyro_z`: Calibrated gyro z (angular velocity) for turn derivative control

Applies calibration to gyro readings and computes angle from accelerometer.

Use with Kalman filter like this:
```julia
angle, gyro_x, gyro_z = compute_angles(kf, ax, ay, az, gx, gy, gz)
update!(kf, gyro_x, angle)
```
That is, the gyro x reading is used as control input, and the angle computed from accelerometer as measurement.


# Arguments
- `kf`: DiscreteKalmanFilter instance
- `ax, ay, az`: Raw accelerometer readings (int16)
- `gx, gy, gz`: Raw gyroscope readings (int16)
"""
function compute_angles(kf::DiscreteKalmanFilter, ax::Integer, ay::Integer, az::Integer,
                       gx::Integer, gy::Integer, gz::Integer)
    # Calculate angle from accelerometer (radians to degrees)
    angle = atan(ay, az) * (180.0f0 / pi)

    # Apply gyro calibration offset and scale
    gyro_x = (gx - GYRO_OFFSET) / GYRO_SCALE


    # Store z-axis gyro (for turn control)
    gyro_z = -gz / GYRO_SCALE

    angle, gyro_x, gyro_z
end



# test
f = DiscreteKalmanFilter()

m_gyro_x = 0.1f0  # rad/s
m_angle = 0.05f0  # rad
update!(f, m_gyro_x, m_angle)
@show angle(f)
@show bias(f)