using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra
using HybridRobotDynamics

<<<<<<< Updated upstream
# Define RABBIT (Five-Link Walker) model parameters
g = 9.81  # gravity

# Link lengths (meters)
ls = [0.5, 0.5, 0.5, 0.5, 0.5]  

# Link masses (kg)
ms = [5.0, 5.0, 5.0, 5.0, 5.0]  

# Link inertias (kg*m^2) - planar scalar Izz only
Is = [0.1, 0.1, 0.1, 0.1, 0.1] 
=======
using Plots
gr()
rm("walker.gif")

function compute_joint_positions(q::AbstractVector, ls::AbstractVector)
  x_com, y_com = q[1:2]
  θ₁, θ₂, θ₃, θ₄, θ₅ = q[3:7]
  l1, l2, l3, l4, l5 = ls

  hip = [x_com, y_com]

  θ_torso      = θ₅
  θ_stance_th  = θ_torso + θ₂
  θ_stance_leg = θ_stance_th + θ₁
  θ_swing_th   = θ_torso + θ₄
  θ_swing_leg  = θ_swing_th + θ₃

  torso_tip     = hip .+ l5 * [cos(θ_torso), sin(θ_torso)]

  stance_thigh  = hip .- l2 * [cos(θ_stance_th), sin(θ_stance_th)]
  stance_foot   = stance_thigh .- l1 * [cos(θ_stance_leg), sin(θ_stance_leg)]

  swing_thigh   = hip .- l4 * [cos(θ_swing_th), sin(θ_swing_th)]
  swing_foot    = swing_thigh .- l3 * [cos(θ_swing_leg), sin(θ_swing_leg)]

  return [
    hip stance_thigh stance_foot swing_thigh swing_foot torso_tip
    ]
end

function animate_walker(xs, ls, N, nx; title="Five-Link Walker", fps=20, gifname="walker.gif")
  anim = @animate for k in 1:N
      xk = xs[(k-1)*nx+1 : k*nx]
      joints = compute_joint_positions(xk, ls)

      # Extract x and y for each joint
      xpts = joints[1, :]
      ypts = joints[2, :]

      # Base scatter plot of all joints
      p = scatter(
          xpts, ypts,
          label = false,
          markersize = 5,
          markercolor = :blue,
          xlim = (-2, 3),
          ylim = (-2, 3),
          title = title,
          xlabel = "x",
          ylabel = "y",
          aspect_ratio = :equal,
      )

      # Add physical links (manually between joint indices)
      # Format: (joint1, joint2)
      links = [
          (1, 2),  # hip to stance thigh
          (2, 3),  # stance thigh to stance foot
          (1, 4),  # hip to swing thigh
          (4, 5),  # swing thigh to swing foot
          (1, 6)   # hip to torso tip
      ]

      for (i, j) in links

          if ((i,j) == (1, 2)) 
            _color = :green; 
          elseif ((i,j) == (2, 3)) 
            _color = :blue; 
          elseif ((i,j) == (1, 4))
             _color = :yellow; 
          elseif ((i,j) == (4, 5)) 
            _color = :red; 
          elseif ((i,j) == (1, 6)) 
            _color = :black; 
          else
            _color = :white;
          end
          
          plot!([xpts[i], xpts[j]], [ypts[i], ypts[j]], lw=2, color=_color, label=false)
      end
  end

  gif(anim, gifname, fps = fps)
end

function draw_walker(x, ls) 
  # Unpack state
  x_com, y_com = x[1:2]
  θ₁, θ₂, θ₃, θ₄, θ₅ = x[3:7]
  l1, l2, l3, l4, l5 = ls

  # Hip is the floating base
  hip = [x_com; y_com]

  # Absolute angles for each link
  θ_torso       = θ₅
  θ_stance_th   = θ_torso + θ₂
  θ_stance_leg  = θ_stance_th + θ₁
  θ_swing_th    = θ_torso + θ₄
  θ_swing_leg   = θ_swing_th + θ₃

  # Compute joint positions
  torso_tip     = hip .+ l5 * [cos(θ_torso), sin(θ_torso)]
  stance_thigh  = hip .- l2 * [cos(θ_stance_th), sin(θ_stance_th)]
  stance_foot   = stance_thigh .- l1 * [cos(θ_stance_leg), sin(θ_stance_leg)]
  swing_thigh   = hip .- l4 * [cos(θ_swing_th), sin(θ_swing_th)]
  swing_foot    = swing_thigh .- l3 * [cos(θ_swing_leg), sin(θ_swing_leg)]

  joints =  Dict(
    "torso_tip"     => torso_tip,
    "hip"           => hip,
    "stance_thigh"  => stance_thigh,
    "stance_foot"   => stance_foot,
    "swing_thigh"   => swing_thigh,
    "swing_foot"    => swing_foot
  )

  # Plot setup
  plt = plot(; aspect_ratio=1, legend=false, title="walker model", xlabel="x", ylabel="y")

  # Plot links
  plot!([torso_tip[1], hip[1]], [torso_tip[2], hip[2]], lw=3, color=:black, label="torso")
  plot!([hip[1], stance_thigh[1]], [hip[2], stance_thigh[2]], lw=3, color=:green, label="stance thigh")
  plot!([stance_thigh[1], stance_foot[1]], [stance_thigh[2], stance_foot[2]], lw=3, color=:blue, label="stance leg")
  plot!([hip[1], swing_thigh[1]], [hip[2], swing_thigh[2]], lw=3, color=:orange, label="swing thigh")
  plot!([swing_thigh[1], swing_foot[1]], [swing_thigh[2], swing_foot[2]], lw=3, color=:red, label="swing leg")

  # Plot joints
  for (label, pos) in joints
      scatter!([pos[1]], [pos[2]]; color=:gray, markersize=6)
      annotate!(pos[1], pos[2] + 0.05, text(label, :black, 8))
  end

  # CoM marker
  scatter!([x_com], [y_com]; marker=:x, color=:magenta, label="CoM")

  display(plt)
  savefig("walker.png")
end


# Define RABBIT (Five-Link Walker) model parameters
g = -9.81  # gravity

# Link lengths (meters)
ls = [
    0.4,  # l1: stance leg
    0.5,  # l2: stance thigh
    0.4,  # l3: swing leg
    0.5,  # l4: swing thigh
    0.6   # l5: torso
]

# Link masses (kg)
ms = [
    5.0,   # m1: stance leg
    7.0,   # m2: stance thigh
    5.0,   # m3: swing leg
    7.0,   # m4: swing thigh
    10.0   # m5: torso
]

# Link inertias (kg·m²)
Is = [
    (1/12) * ms[1] * ls[1]^2,  # I1: stance leg
    (1/12) * ms[2] * ls[2]^2,  # I2: stance thigh
    (1/12) * ms[3] * ls[3]^2,  # I3: swing leg
    (1/12) * ms[4] * ls[4]^2,  # I4: swing thigh
    (1/12) * ms[5] * ls[5]^2   # I5: torso
]

# # Initial conditions (position + velocity)
# # q1...q7, dq1...dq7
xic = [
    0.0, 0.75,                  # x_com, y_com
    π/5, -π/4, -π/5, π/4, π/2,  # θ₁ to θ₅
    0.0, 0.0,                   # ẋ_com, ẏ_com
    0.0, 0.0, 0.0, 0.0, 0.0     # joint velocities
]

function generate_stand_us(N::Int, dt::Float64; amp1=10.0, amp2=10.0, freq=1.0)
  # Create sinusoidal torque profiles for two hip actuators
  us = zeros(N * 2)
  for k in 1:2:2*N-1
      us[k] =  34 # Hip 1 (stance) // +:CCW    -:CW
      us[k+1] =  -23 # Hip 2 (swing) //  +:CW     -:CCW
  end
  return us
end

draw_walker(xic, ls)
>>>>>>> Stashed changes

# Generate the hybrid system
system = five_link_walker(Is, ls, ms, g)

# Roll-out simulation parameters
<<<<<<< Updated upstream
N = 100               # number of timesteps
Δt = 0.05             # timestep duration (s)

# Control inputs (no actuation for now)
us = zeros(N * system.nu)

# Initial conditions (position + velocity)
# q1...q5, dq1...dq5
xic = [0.0, -0.2, 0.2, -0.1, 0.0, 
       0.0, 0.0, 0.0, 0.0, 0.0]  

init_mode = :stance
=======
N =  100               # number of timesteps
Δt = 0.05             # timestep duration (s)

# Control inputs
# us = zeros(N * system.nu)
us = generate_stand_us(N-1, Δt)

init_mode = :flight
>>>>>>> Stashed changes

# Integrator choice
rk4 = ExplicitIntegrator(:rk4)

# Roll out the trajectory
xs = roll_out(system, rk4, N, Δt, us, xic, init_mode)

# Plotting configuration
<<<<<<< Updated upstream
plot_2d_states(N, system.nx, (1,3), xs; title="Five-Link Walker Roll-Out")
nothing
=======
animate_walker(xs, ls, N, system.nx; title="Walker Rollout", fps=20)
nothing
>>>>>>> Stashed changes
