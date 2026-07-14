
using DyadBotComponents
using Test

@testset "DyadBotComponents" begin

@testset "Generated model tests" begin
    include("../generated/tests.jl")
end

@testset "IMUKalmanFilter" begin
    f = IMUKalmanFilter()

    m_gyro_x = 0.0f0  # rad/s
    m_angle = 0.05f0  # rad
    for i = 1:100
        update!(f, m_gyro_x, m_angle)
    end
    @test DyadBotComponents.angle(f) ≈ m_angle rtol = 1e-4
    @test abs(bias(f)) < 1e-4


    m_gyro_x = 0.1f0  # rad/s
    m_angle = 0.05f0  # rad
    for i = 1:100
        update!(f, m_gyro_x, m_angle)
    end
    @test DyadBotComponents.angle(f) ≈ m_angle rtol = 1e-4
    @test bias(f) ≈ 0.1 rtol = 1e-4
end

end

# Self-contained (defines its own `using`s and `@testset`), so include at top level.
include("test_stabilization.jl")
