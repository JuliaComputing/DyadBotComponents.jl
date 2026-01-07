# Balance car controller translated from original Tumbller/BalanceCar.h
# Parallel PD/PI controller architecture with Kalman filter state estimation
include(joinpath(@__DIR__, "DiscreteKalmanFilter.jl"))
export BalanceController, compute_pwm!, handle_motion_mode!, balance_car!
export MotionMode, STOP, START, FORWARD, BACKWARD, LEFT, RIGHT
export car_stop!, apply_motor_output!

# =============================================================================
# Hardware Interface Stubs (override for actual hardware)
# =============================================================================

# Default no-op stubs - override these for actual hardware
function digitalWrite(pin, value)
    println("Stub - override for actual hardware")
end

function analogWrite(pin, value)
    println("Stub - override for actual hardware")
end

# Hardware pin constants (defaults - override for actual hardware)
const AIN1 = 0
const BIN1 = 0
const PWMA_LEFT = 0
const PWMB_RIGHT = 0
const STBY_PIN = 0

# =============================================================================
# Constants (matching original BalanceCar.h parameters)
# =============================================================================

# PID parameters
const KP_BALANCE = 55.0f0
const KD_BALANCE = 0.75f0
const KP_SPEED = 10.0f0
const KI_SPEED = 0.26f0
const KP_TURN = 2.5f0
const KD_TURN = 0.5f0

# Angle limits (degrees)
const BALANCE_ANGLE_MIN = -22.0f0
const BALANCE_ANGLE_MAX = 22.0f0

# PWM limits
const PWM_MAX = 255.0f0
const PWM_MIN = -255.0f0

# Integral anti-windup limit
const INTEGRAL_LIMIT = 3000.0f0

# Speed control decimation (runs every N cycles)
const SPEED_CONTROL_PERIOD = 8

# =============================================================================
# Motion Mode Enum
# =============================================================================

@enum MotionMode STOP START FORWARD BACKWARD LEFT RIGHT

# =============================================================================
# Controller State Struct
# =============================================================================

mutable struct BalanceController
    # Kalman filter for angle estimation
    kf::IMUKalmanFilter

    # Encoder pulse accumulators (signed based on PWM direction)
    encoder_left_pulse::Int
    encoder_right_pulse::Int

    # Speed control state (PI integrator)
    speed_filter::Float32
    speed_filter_old::Float32
    car_speed_integral::Float32
    speed_control_period_count::Int

    # Control outputs
    speed_control_output::Float32
    rotation_control_output::Float32

    # Setpoints
    setting_car_speed::Int
    setting_turn_speed::Int

    # PWM outputs
    pwm_left::Float32
    pwm_right::Float32

    # Motion mode
    motion_mode::MotionMode

    # Calibration offsets
    angle_zero::Float32
    angular_velocity_zero::Float32
end

"""
    BalanceController(; angle_zero=0.0f0, angular_velocity_zero=0.0f0)

Create a balance controller with default parameters matching the original Tumbller code.
"""
function BalanceController(; angle_zero::Float32=0.0f0, angular_velocity_zero::Float32=0.0f0)
    BalanceController(
        IMUKalmanFilter(),      # kf
        0,                      # encoder_left_pulse
        0,                      # encoder_right_pulse
        0.0f0,                  # speed_filter
        0.0f0,                  # speed_filter_old
        0.0f0,                  # car_speed_integral
        0,                      # speed_control_period_count
        0.0f0,                  # speed_control_output
        0.0f0,                  # rotation_control_output
        0,                      # setting_car_speed
        0,                      # setting_turn_speed
        0.0f0,                  # pwm_left
        0.0f0,                  # pwm_right
        STOP,                   # motion_mode
        angle_zero,
        angular_velocity_zero
    )
end

# =============================================================================
# Control Signal Computation (decoupled from motor output)
# =============================================================================

"""
    compute_pwm!(ctrl, encoder_count_left, encoder_count_right, ax, ay, az, gx, gy, gz)

Compute PWM control signals based on IMU and encoder data.
Updates `ctrl.pwm_left` and `ctrl.pwm_right` in place.

Returns the Kalman-filtered angle for use in motion mode handling.

# Control Architecture (parallel, not cascade):
- Balance control: PD on tilt angle
- Speed control: PI on wheel speed (runs every 8 cycles)
- Turn control: P + D on turn rate

Final: `pwm = balance - speed ± rotation`
"""
function compute_pwm!(ctrl::BalanceController,
                      encoder_count_left::Integer, encoder_count_right::Integer,
                      ax::Integer, ay::Integer, az::Integer,
                      gx::Integer, gy::Integer, gz::Integer)

    # Accumulate encoder pulses with sign based on current PWM direction
    ctrl.encoder_left_pulse += ctrl.pwm_left < 0 ? -encoder_count_left : encoder_count_left
    ctrl.encoder_right_pulse += ctrl.pwm_right < 0 ? -encoder_count_right : encoder_count_right

    # Get calibrated angles from Kalman filter
    accel_angle, gyro_x, gyro_z = compute_angles(ctrl.kf, ax, ay, az, gx, gy, gz)

    # Update Kalman filter
    update!(ctrl.kf, gyro_x, accel_angle)
    kalman_angle = angle(ctrl.kf)

    # Balance control (PD on tilt angle)
    balance_control_output = KP_BALANCE * (kalman_angle - ctrl.angle_zero) +
                             KD_BALANCE * (gyro_x - ctrl.angular_velocity_zero)

    # Speed control (PI, runs every 8 cycles = 40ms at 5ms sample time)
    ctrl.speed_control_period_count += 1
    if ctrl.speed_control_period_count >= SPEED_CONTROL_PERIOD
        ctrl.speed_control_period_count = 0

        # Average wheel speed from encoder pulses
        car_speed = (ctrl.encoder_left_pulse + ctrl.encoder_right_pulse) * 0.5f0
        ctrl.encoder_left_pulse = 0
        ctrl.encoder_right_pulse = 0

        # Low-pass filter
        ctrl.speed_filter = ctrl.speed_filter_old * 0.7f0 + car_speed * 0.3f0
        ctrl.speed_filter_old = ctrl.speed_filter

        # PI integrator with setpoint
        ctrl.car_speed_integral += ctrl.speed_filter - ctrl.setting_car_speed

        # Anti-windup
        ctrl.car_speed_integral = clamp(ctrl.car_speed_integral, -INTEGRAL_LIMIT, INTEGRAL_LIMIT)

        # Speed control output (negative feedback)
        ctrl.speed_control_output = -KP_SPEED * ctrl.speed_filter - KI_SPEED * ctrl.car_speed_integral

        # Turn control output (computed at same rate as speed control)
        ctrl.rotation_control_output = ctrl.setting_turn_speed + KD_TURN * gyro_z
    end

    # Combine control outputs
    ctrl.pwm_left = balance_control_output - ctrl.speed_control_output - ctrl.rotation_control_output
    ctrl.pwm_right = balance_control_output - ctrl.speed_control_output + ctrl.rotation_control_output

    # Clamp to PWM limits
    ctrl.pwm_left = clamp(ctrl.pwm_left, PWM_MIN, PWM_MAX)
    ctrl.pwm_right = clamp(ctrl.pwm_right, PWM_MIN, PWM_MAX)

    return kalman_angle
end

# =============================================================================
# Motion Mode Handling
# =============================================================================

"""
    handle_motion_mode!(ctrl, kalman_angle; key_flag='0')

Handle motion mode logic including angle limit detection and STOP mode behavior.
Modifies `ctrl.pwm_left`, `ctrl.pwm_right`, `ctrl.motion_mode`, and integrator state.

Returns `true` if motors should be stopped (car_stop! should be called).
"""
function handle_motion_mode!(ctrl::BalanceController, kalman_angle::Float32; key_flag::Char='0')
    should_stop_motors = false

    # Check angle limits - force STOP if exceeded (except during START or STOP modes)
    if ctrl.motion_mode != START && ctrl.motion_mode != STOP
        if kalman_angle < BALANCE_ANGLE_MIN || kalman_angle > BALANCE_ANGLE_MAX
            ctrl.motion_mode = STOP
            should_stop_motors = true
        end
    end

    # Handle STOP mode
    if ctrl.motion_mode == STOP
        if key_flag != '4'
            # Full stop - zero everything
            ctrl.car_speed_integral = 0.0f0
            ctrl.setting_car_speed = 0
            ctrl.pwm_left = 0.0f0
            ctrl.pwm_right = 0.0f0
            should_stop_motors = true
        else
            # key_flag == '4': Reset state but don't stop motors (allows restart)
            ctrl.car_speed_integral = 0.0f0
            ctrl.setting_car_speed = 0
            ctrl.pwm_left = 0.0f0
            ctrl.pwm_right = 0.0f0
        end
    end

    return should_stop_motors
end

# =============================================================================
# Motor Output Functions
# =============================================================================

"""
    car_stop!()

Stop both motors by setting PWM to 0 and appropriate direction pins.
"""
function car_stop!()
    digitalWrite(AIN1, 1)     # HIGH
    digitalWrite(BIN1, 0)     # LOW
    digitalWrite(STBY_PIN, 1) # HIGH (standby off = enabled, but PWM=0)
    analogWrite(PWMA_LEFT, 0)
    analogWrite(PWMB_RIGHT, 0)
end

"""
    apply_motor_output!(pwm_left, pwm_right)

Apply PWM signals to motors based on computed control values.
Handles direction pin setting based on PWM sign.
"""
function apply_motor_output!(pwm_left::Float32, pwm_right::Float32)
    # Left motor
    if pwm_left < 0
        digitalWrite(AIN1, 1)
        analogWrite(PWMA_LEFT, round(Int, -pwm_left))
    else
        digitalWrite(AIN1, 0)
        analogWrite(PWMA_LEFT, round(Int, pwm_left))
    end

    # Right motor
    if pwm_right < 0
        digitalWrite(BIN1, 1)
        analogWrite(PWMB_RIGHT, round(Int, -pwm_right))
    else
        digitalWrite(BIN1, 0)
        analogWrite(PWMB_RIGHT, round(Int, pwm_right))
    end
end

# =============================================================================
# Main Update Function
# =============================================================================

"""
    balance_car!(ctrl, encoder_count_left, encoder_count_right, ax, ay, az, gx, gy, gz; key_flag='0')

Main control loop update function. Combines:
1. `compute_pwm!` - Calculate control signals
2. `handle_motion_mode!` - Apply motion mode logic
3. `apply_motor_output!` or `car_stop!` - Output to motors

Call this at 200Hz (5ms interval) matching the original Tumbller timer.

# Arguments
- `ctrl`: BalanceController instance
- `encoder_count_left`, `encoder_count_right`: Encoder pulse counts since last call
- `ax, ay, az`: Raw accelerometer readings (Int16)
- `gx, gy, gz`: Raw gyroscope readings (Int16)
- `key_flag`: Remote control key flag (default '0')
"""
function balance_car!(ctrl::BalanceController,
                      encoder_count_left::Integer, encoder_count_right::Integer,
                      ax::Integer, ay::Integer, az::Integer,
                      gx::Integer, gy::Integer, gz::Integer;
                      key_flag::Char='0')

    # Compute control signals
    kalman_angle = compute_pwm!(ctrl, encoder_count_left, encoder_count_right,
                                ax, ay, az, gx, gy, gz)

    # Handle motion mode logic
    should_stop = handle_motion_mode!(ctrl, kalman_angle; key_flag)

    # Apply motor output
    if should_stop
        car_stop!()
    elseif ctrl.motion_mode != STOP
        apply_motor_output!(ctrl.pwm_left, ctrl.pwm_right)
    end

    return nothing
end


# =============================================================================
# Test no error
# =============================================================================
using Test
ctrl = BalanceController()
@test_nowarn balance_car!(ctrl, 1, 2, 3, 4, 5, 6, 7, 8)
@test ctrl.pwm_left == 0 # motion mode = stop

##

ctrl = BalanceController()
ctrl.motion_mode = FORWARD
@test_nowarn compute_pwm!(ctrl, 512, 512, 512, 512, 512, 512, 512, 512)
@test ctrl.pwm_left != 0