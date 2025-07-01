using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using HybridRobotDynamics

# Define elastic bouncing ball model
system = hopper()

# Define roll-out parameters
N = 50
Δt = 0.05

# Define control inputs and initial conditions
us = zeros((N-1) * system.nu)
xic = [0.0; 5.0; 0.0; 4.0; ones(4)]
init_mode = :flight

# Roll out and visualize
rk4 = ExplicitIntegrator(:rk4)
xs = roll_out(system, rk4, N, Δt, us, xic, init_mode)
plot_2d_states(N, system.nx, (3,4), xs; title="Hopper Roll-Out", ylim=(-5,5))

nothing
