using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using HybridRobotDynamics

# Define elastic bouncing ball model
system = bouncing_ball()

# Define roll-out parameters
N = 50
Δt = 0.05

# Define control inputs and initial conditions
us = zeros(N * system.nu)
xic = [0.0; 5.0; 1.0; 0.0]
init_mode = :down

# Roll out and visualize
rk4 = ExplicitIntegrator(:rk4)
xs = roll_out(system, rk4, N, Δt, us, xic, init_mode)
plot_2d_states(N, system.nx, (1,2), xs; title="Bouncing Ball Roll-Out")
nothing
