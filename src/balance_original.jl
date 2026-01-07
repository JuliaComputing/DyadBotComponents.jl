# Balance car controller translated from original Tumbller/BalanceCar.h
# Parallel PD/PI controller architecture with Kalman filter state estimation
include(joinpath(@__DIR__, "DiscreteKalmanFilter.jl"))
export BalanceController, compute_pwm!, handle_motion_mode!, balance_car!
export MotionMode, STOP, START, FORWARD, BACKWARD, LEFT, RIGHT
export car_stop!, apply_motor_output!
export ParallelPID

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
# ParallelPID Controller Parameters
# =============================================================================

"""
    ParallelPID

Controller parameters for parallel PD/PI architecture (original Tumbller).
Stores all tunable gains and limits for the balance car controller.
Default values match the original Tumbller code.

# Fields
- `kp_balance`, `kd_balance`: Balance control (PD on tilt angle)
- `kp_speed`, `ki_speed`: Speed control (PI on wheel speed)
- `kp_turn`, `kd_turn`: Turn control
- `balance_angle_min`, `balance_angle_max`: Angle limits (degrees)
- `pwm_min`, `pwm_max`: PWM output limits
- `integral_limit`: Anti-windup limit for speed integrator
- `speed_control_period`: Decimation period for speed control loop
"""
@kwdef mutable struct ParallelPID
    # Parameters
    # Balance control (PD)
    kp_balance::Float32 = 55.0f0
    kd_balance::Float32 = 0.75f0

    # Speed control (PI)
    kp_speed::Float32 = 10.0f0
    ki_speed::Float32 = 0.26f0

    # Turn control
    kp_turn::Float32 = 2.5f0
    kd_turn::Float32 = 0.5f0

    # Angle limits (degrees)
    balance_angle_min::Float32 = -22.0f0
    balance_angle_max::Float32 = 22.0f0

    # PWM limits
    pwm_min::Float32 = -255.0f0
    pwm_max::Float32 = 255.0f0

    # Integral anti-windup limit
    integral_limit::Float32 = 3000.0f0

    # Speed control decimation period
    speed_control_period::Int = 8

    # State (mutable, initialized to zero)
    speed_filter::Float32 = 0.0f0
    speed_filter_old::Float32 = 0.0f0
    car_speed_integral::Float32 = 0.0f0
    speed_control_period_count::Int = 0
    speed_control_output::Float32 = 0.0f0
    rotation_control_output::Float32 = 0.0f0
end

"""
    reset_state!(c::ParallelPID)

Reset controller state to initial values (zeros).
"""
function reset_state!(c::ParallelPID)
    c.speed_filter = 0.0f0
    c.speed_filter_old = 0.0f0
    c.car_speed_integral = 0.0f0
    c.speed_control_period_count = 0
    c.speed_control_output = 0.0f0
    c.rotation_control_output = 0.0f0
    nothing
end

# =============================================================================
# Motion Mode Enum
# =============================================================================

@enum MotionMode STOP START FORWARD BACKWARD LEFT RIGHT

# =============================================================================
# Controller State Struct
# =============================================================================

"""
    BalanceController{C}

Balance controller state with parametric controller type `C`.
The controller type determines the control algorithm (e.g., `ParallelPID`).
Controller-specific state is stored in the controller itself.
"""
mutable struct BalanceController{C}
    # Controller (contains parameters and algorithm-specific state)
    controller::C

    # Kalman filter for angle estimation
    kf::IMUKalmanFilter

    # Encoder pulse accumulators (signed based on PWM direction)
    encoder_left_pulse::Int
    encoder_right_pulse::Int

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
    BalanceController(; controller=ParallelPID(), angle_zero=0.0f0, angular_velocity_zero=0.0f0)

Create a balance controller with the specified controller type (defaults to ParallelPID).
"""
function BalanceController(;
    controller::C = ParallelPID(),
    angle_zero::Float32 = 0.0f0,
    angular_velocity_zero::Float32 = 0.0f0
) where C
    BalanceController{C}(
        controller,
        IMUKalmanFilter(),      # kf
        0,                      # encoder_left_pulse
        0,                      # encoder_right_pulse
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

# Control Architecture for ParallelPID (parallel, not cascade):
- Balance control: PD on tilt angle
- Speed control: PI on wheel speed (runs every N cycles)
- Turn control: P + D on turn rate

Final: `pwm = balance - speed ± rotation`
"""
function compute_pwm!(ctrl::BalanceController{ParallelPID},
                      encoder_count_left::Integer, encoder_count_right::Integer,
                      ax::Integer, ay::Integer, az::Integer,
                      gx::Integer, gy::Integer, gz::Integer)
    c = ctrl.controller

    # Accumulate encoder pulses with sign based on current PWM direction
    ctrl.encoder_left_pulse += ctrl.pwm_left < 0 ? -encoder_count_left : encoder_count_left
    ctrl.encoder_right_pulse += ctrl.pwm_right < 0 ? -encoder_count_right : encoder_count_right

    # Get calibrated angles from Kalman filter
    accel_angle, gyro_x, gyro_z = compute_angles(ax, ay, az, gx, gy, gz)

    # Update Kalman filter
    update!(ctrl.kf, gyro_x, accel_angle)
    kalman_angle = angle(ctrl.kf)

    # Balance control (PD on tilt angle)
    balance_control_output = c.kp_balance * (kalman_angle - ctrl.angle_zero) +
                             c.kd_balance * (gyro_x - ctrl.angular_velocity_zero)

    # Speed control (PI, runs every N cycles)
    c.speed_control_period_count += 1
    if c.speed_control_period_count >= c.speed_control_period
        c.speed_control_period_count = 0

        # Average wheel speed from encoder pulses
        car_speed = (ctrl.encoder_left_pulse + ctrl.encoder_right_pulse) * 0.5f0
        ctrl.encoder_left_pulse = 0
        ctrl.encoder_right_pulse = 0

        # Low-pass filter
        c.speed_filter = c.speed_filter_old * 0.7f0 + car_speed * 0.3f0
        c.speed_filter_old = c.speed_filter

        # PI integrator with setpoint
        c.car_speed_integral += c.speed_filter - ctrl.setting_car_speed

        # Anti-windup
        c.car_speed_integral = clamp(c.car_speed_integral, -c.integral_limit, c.integral_limit)

        # Speed control output (negative feedback)
        c.speed_control_output = -c.kp_speed * c.speed_filter - c.ki_speed * c.car_speed_integral

        # Turn control output (computed at same rate as speed control)
        c.rotation_control_output = ctrl.setting_turn_speed + c.kd_turn * gyro_z
    end

    # Combine control outputs
    ctrl.pwm_left = balance_control_output - c.speed_control_output - c.rotation_control_output
    ctrl.pwm_right = balance_control_output - c.speed_control_output + c.rotation_control_output

    # Clamp to PWM limits
    ctrl.pwm_left = clamp(ctrl.pwm_left, c.pwm_min, c.pwm_max)
    ctrl.pwm_right = clamp(ctrl.pwm_right, c.pwm_min, c.pwm_max)

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
    c = ctrl.controller
    should_stop_motors = false

    # Check angle limits - force STOP if exceeded (except during START or STOP modes)
    if ctrl.motion_mode != START && ctrl.motion_mode != STOP
        if kalman_angle < c.balance_angle_min || kalman_angle > c.balance_angle_max
            ctrl.motion_mode = STOP
            should_stop_motors = true
        end
    end

    # Handle STOP mode
    if ctrl.motion_mode == STOP
        if key_flag != '4'
            # Full stop - zero everything
            reset_state!(c)
            ctrl.setting_car_speed = 0
            ctrl.pwm_left = 0.0f0
            ctrl.pwm_right = 0.0f0
            should_stop_motors = true
        else
            # key_flag == '4': Reset state but don't stop motors (allows restart)
            reset_state!(c)
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