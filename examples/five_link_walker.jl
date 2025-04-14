using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra
using HybridRobotDynamics

# Define RABBIT (Five-Link Walker) model parameters
g = 9.81  # gravity

# Link lengths (meters)
ls = [0.5, 0.5, 0.5, 0.5, 0.5]  

# Link masses (kg)
ms = [5.0, 5.0, 5.0, 5.0, 5.0]  

# Link inertias (kg*m^2) - planar scalar Izz only
Is = [0.1, 0.1, 0.1, 0.1, 0.1] 

# Generate the hybrid system
system = five_link_walker(Is, ls, ms, g)

# Roll-out simulation parameters
N = 100               # number of timesteps
Δt = 0.05             # timestep duration (s)

# Control inputs (no actuation for now)
us = zeros(N * system.nu)

# Initial conditions (position + velocity)
# q1...q5, dq1...dq5
xic = [0.0, -0.2, 0.2, -0.1, 0.0, 
       0.0, 0.0, 0.0, 0.0, 0.0]  

init_mode = :stance

# Integrator choice
rk4 = ExplicitIntegrator(:rk4)

# Roll out the trajectory
xs = roll_out(system, rk4, N, Δt, us, xic, init_mode)

# Plotting configuration
plot_2d_states(N, system.nx, (1,3), xs; title="Five-Link Walker Roll-Out")
nothing