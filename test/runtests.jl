
using DyadBotComponents
using Test
    
# include("../generated/tests.jl")

# test
f = IMUKalmanFilter()

m_gyro_x = 0.0f0  # rad/s
m_angle = 0.05f0  # rad
for i = 1:100
    update!(f, m_gyro_x, m_angle)
end
@test angle(f) ≈ m_angle rtol = 1e-4
@test abs(bias(f)) < 1e-4


m_gyro_x = 0.1f0  # rad/s
m_angle = 0.05f0  # rad
for i = 1:100
    update!(f, m_gyro_x, m_angle)
end
@test angle(f) ≈ m_angle rtol = 1e-4
@test bias(f) ≈ 0.1 rtol = 1e-4