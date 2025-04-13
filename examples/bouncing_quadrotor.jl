using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra
using HybridRobotDynamics

# Define bouncing quadrotor model
# Crazyflie 2.1 parameters: https://arxiv.org/pdf/1608.05786
e = 1.0
g = 9.81
m = 0.027
j = [1.436e-5, 1.395e-5,  2.173e-5]
b = 7.9379e-12 / 3.1582e-10 * ones(4)
c = 0.0283 * ones(4)
d = 0.0283 * ones(4)
system = bouncing_quadrotor(e, g, m, j, b, c, d)

# Define roll-out parameters
N = 50
Δt = 0.05

# Define control inputs and initial conditions
us = zeros(N * system.nu)
xic = [[0.0, 0.0, 5.0]; [1.0, 0.0, 0.0, 0.0]; [1.0, 0.0, 0.0]; [0.0, 0.0, 0.0]]
init_mode = :flight

# Roll out and visualize
rk4 = ExplicitIntegrator(:rk4)
xs = roll_out(system, rk4, N, Δt, us, xic, init_mode)
plot_2d_states(N, system.nx, (1,3), xs; title="Bouncing Quadrotor Roll-Out")
nothing
