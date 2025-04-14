"""
    bouncing_ball(e=1.0, g=9.81)

Returns the hybrid system model containing the modes and transitions of a planar (in)elastic bouncing ball.
"""
function bouncing_ball(
    e::Real = 1.0,
    g::Real = 9.81
)::HybridSystem
    # Model ballistic dynamics with thrust
    # State space: x, y, xdot, ydot
    nx = 4
    nu = 1
    flight_flow = (
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    ) -> [x[3:4]; 0.0; u[1] - g]

    # Define hybrid mode
    flight_mode = HybridMode(flight_flow)

    # Define impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[2]::DiffFloat
    R_impact = x::Vector{<:DiffFloat} -> [x[1:3]; -e*x[4]]::Vector{<:DiffFloat}
    impact = Transition(flight_flow, flight_flow, g_impact, R_impact)
    add_transition!(flight_mode, flight_mode, impact)

    # Create hybrid system dicts
    transitions = Dict(:impact => impact)
    modes = Dict(:flight => flight_mode)
    return HybridSystem(nx, nu, transitions, modes)
end

"""
"""
function bouncing_quadrotor(
    e::Float64,
    gravity::Float64,
    mass::Float64,
    principal_moments::Vector{Float64},
    torque_consts::Vector{Float64},
    rotor_xs::Vector{Float64},
    rotor_ys::Vector{Float64}
)::HybridSystem
    @assert size(principal_moments) == (3,)
    @assert size(torque_consts) == (4,)
    @assert size(rotor_xs) == (4,)
    @assert size(rotor_ys) == (4,)

    # Constants
    g = [0.0, 0.0, -gravity]
    J = diagm(principal_moments)
    K = [zeros(2,4); ones(1,4)]
    B = [
        rotor_ys';
        -rotor_xs';
        (torque_consts .* [-1.0, 1.0, -1.0, 1.0])'
    ]

    # Smooth dynamics
    function quadrotor_flow(
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    )::Vector{<:DiffFloat}
        q = x[4:7]
        v = x[8:10]
        ω = x[11:13]
        return [
            Qquat(q) * v;
            0.5 * Gquat(q) * ω;
            Qquat(q)'*g + K*u/mass - cross(ω, v);
            J \ (B*u - cross(ω, J*ω))
        ]
    end

    # Hybrid modes
    flight_mode = HybridMode(quadrotor_flow)

    # Impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[3]::DiffFloat
    R_impact = x::Vector{<:DiffFloat} -> [
        x[1:2]; 1e-3; x[4:9]; -e*x[10]; x[11:13]]::Vector{<:DiffFloat}
    impact = Transition(quadrotor_flow, quadrotor_flow, g_impact, R_impact)
    add_transition!(flight_mode, flight_mode, impact)

    # Create hybrid system
    nx = 13
    nu = 4
    transitions = Dict(:impact => impact)
    modes = Dict(:flight => flight_mode)
    return HybridSystem(nx, nu, transitions, modes)
end

"""
"""
function hopper(
    m1::Real = 5.0,    # body mass
    m2::Real = 1.0,    # foot mass
    e::Real = 0.0,     # coefficient of restitution of foot
    g::Real = 9.81,    # acceleration due to gravity
    Llb::Real = 0.5,
    Lub::Real = 1.5
)::HybridSystem
    # State space:
    #   body x, body y, foot x, foot y,
    #   body xdot, body ydot, foot xdot, foot ydot
    nx = 8
    nu = 2
    M = Diagonal([m1; m1; m2; m2])

    function get_length_vector(
        x::Vector{<:DiffFloat}
    )::Vector{<:DiffFloat}
        return x[1:2] - x[3:4]
    end

    function get_unit_length(
        x::Vector{<:DiffFloat}
    )::Vector{<:DiffFloat}
        L = get_length_vector(x)
        return L / norm(L)
    end

    function B_flight(
        x::Vector{<:DiffFloat}
    )::Matrix{<:DiffFloat}
        L1, L2 = get_unit_length(x)
        return [L1  L2; L2 -L1; -L1 -L2; -L2  L1]
    end

    function B_stance(
        x::Vector{<:DiffFloat}
    )::Matrix{<:DiffFloat}
        L1, L2 = get_unit_length(x)
        return [L1  L2; L2 -L1; zeros(2,2)]
    end

    function flow(
        control_allocation::Function,
        gravity::Vector{<:Real},
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    )::Vector{<:DiffFloat}
        B = control_allocation(x)
        return [x[5:8]; gravity + M\(B*u)]
    end

    # Define flight mode
    grav_flight = [0; -g; 0; -g]
    flight_flow = (
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    ) -> flow(B_flight, grav_flight, x, u)
    flight_mode = HybridMode(flight_flow)

    # Define stance mode
    grav_stance =  [0; -g; 0; 0]
    stance_flow = (
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    ) -> flow(B_stance, grav_stance, x, u)
    stance_mode = HybridMode(stance_flow)

    # Define impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[4]
    R_impact = x::Vector{<:DiffFloat} -> [x[1:6]; e * x[7:8]]
    impact = Transition(flight_flow, stance_flow, g_impact, R_impact)
    add_transition!(flight_mode, stance_mode, impact)

    # Define liftoff transition
    g_liftoff = x::Vector{<:DiffFloat} -> -x[4]
    R_liftoff = x::Vector{<:DiffFloat} -> x
    liftoff = Transition(stance_flow, flight_flow, g_liftoff, R_liftoff)
    add_transition!(stance_mode, flight_mode, liftoff)

    # Create hybrid system dicts
    transitions = Dict(
        :liftoff => liftoff,
        :impact => impact
    )
    modes = Dict(
        :flight => flight_mode,
        :stance => stance_mode
    )

    # Define length inequality constraint functions
    glb = x::Vector{<:DiffFloat} -> Llb - norm(get_length_vector(x))
    gub = x::Vector{<:DiffFloat} -> -Lub + norm(get_length_vector(x))
    g = x::Vector{<:DiffFloat} -> [glb(x); gub(x)]
    g_stage = (x::Vector{<:DiffFloat}, u::Vector{<:DiffFloat}) -> g(x)
    g_term = x::Vector{<:DiffFloat} -> g(x)

    return HybridSystem(
        nx, nu, transitions, modes;
        stage_ineq_constr=g_stage, term_ineq_constr=g_term
    )
end

"""
## Inputs:
Is = [I1 I2 I3 I4 I5] :: link inertias (rotational :: Izz)
ms = [m1 m2 m3 m4 m5] :: link masses

## Summary:
Returns the hybrid system model containing the modes and 
transitions of a planar 5-link biped walker.

# Notes:
q :: link angles
q1 : stance foot angle
q2 : stance leg angle
q3 : swing leg angle
q4 : swing thigh angle
q5 : torso angle
"""
function five_link_walker(
    Is::Vector{<:Real},        # Link inertias
    ls::Vector{<:Real},        # Link lengths
    ms::Vector{<:Real},        # Link masses
    g::Real = 9.81             # Gravitational constant
)::HybridSystem

  nx = 10  # 5 positions + 5 velocities
  nu = 2   # 2 actuators (hips)

  # Unpack link parameters
  I1, I2, I3, I4, I5 = Is
  l1, l2, l3, l4, l5 = ls
  m1, m2, m3, m4, m5 = ms
  
  # Actuation mapping (only joints actuated)
  B = [ zeros(3,2); I(2) ]

  function com_position(q, i)
    if i == 1 # stance leg
        x = -(l1/2)*sin(q[1])
        y =  (l1/2)*cos(q[1])
    elseif i == 2 # swing leg
        x = -l1*sin(q[1]) + (l2/2)*sin(q[2])
        y =  l1*cos(q[1]) - (l2/2)*cos(q[2])
    elseif i == 3 # stance thigh
        x = -(l1 + l3/2)*sin(q[1])
        y =  (l1 + l3/2)*cos(q[1])
    elseif i == 4 # swing thigh
        x = -l1*sin(q[1]) + (l2 + l4/2)*sin(q[2])
        y =  l1*cos(q[1]) - (l2 + l4/2)*cos(q[2])
    elseif i == 5 # torso
        x = -(l1 + l3 + l5/2)*sin(q[1])
        y =  (l1 + l3 + l5/2)*cos(q[1])
    end
    return [x; y]
  end

  function com_height(q, i)
      return com_position(q, i)[2]
  end

  # Swing foot position and Jacobian
  function p_swing(q)
    # Stance leg: l1, Swing leg: l2
    x = -l1*sin(q[1]) + l2*sin(q[2])
    y =  l1*cos(q[1]) - l2*cos(q[2])
    return [x; y]
  end

  function swing_foot_height(q)
      return p_swing(q)[2]
  end

  function swing_foot_jacobian(q)
    Symbolics.jacobian(p_swing(q), q)
  end

  # Mass Matrix M(q)
  function M(q)
    n = length(q)

    @variables dq[1:n]
    dq = collect(dq)

    # Total kinetic energy = translational + rotational
    Tke = 0 
    for i = 1:n
      p_com = com_position(q, i)
      J_i = Symbolics.jacobian(p_com, q)
      v_com = J_i * dq
      ω_i = dq[i]
      Tke += 0.5 * ms[i] * Symbolics.dot(v_com, v_com) + 0.5 * Is[i] * ω_i^2
    end

    # Mass matrix from Hessian of T wrt dq
    return Symbolics.simplify(Symbolics.hessian(Tke, dq))
  end

  # Coriolis Matrix C(q, dq) from M(q)
  function C(q, dq)
    Mq = M(q)
    n = length(q)
    Cq = zeros(Symbolics.Num, n, n)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                c_ijk = 0.5 * (
                    Symbolics.derivative(Mq[i,j], q[k]) +
                    Symbolics.derivative(Mq[i,k], q[j]) -
                    Symbolics.derivative(Mq[j,k], q[i])
                )
                Cq[i,j] += c_ijk * dq[k]
            end
        end
    end
    return Cq
  end

  # Gravity Vector G(q)
  function G(q)
    n = length(q)
    Gq = zeros(Symbolics.Num, n)
    for i = 1:n
        Gq[i] = -ms[i] * g * Symbolics.derivative(com_height(q, i), q[i])
    end
    return Gq
  end

  ## Symbolic Variables
  @variables q[1:5] dq[1:5]

  ## Pre-compute dynamics matrices sym.
  M_func = M(q)
  C_func = C(q, dq)
  G_func = G(q)
  J_func = swing_foot_jacobian(q)

  # Build fast eval functions
  eval_M = (Symbolics.build_function(M_func, q, expression=Val{false}) |> eval)[1]
  eval_C = (Symbolics.build_function(C_func, [q; dq], expression=Val{false}) |> eval)[1]
  eval_G = (Symbolics.build_function(G_func, q, expression=Val{false}) |> eval)[1]
  eval_J = (Symbolics.build_function(J_func, q, expression=Val{false}) |> eval)[1]

  # Impact Reset Map
  function reset_velocity_after_impact(x)
    q = x[1:5]
    dq_pre = x[6:10]
    Mq = eval_M(q)
    J = eval_J(q)
    K = [Mq -J'; J zeros(2, 2)]
    rhs = [Mq * dq_pre; zeros(2)]
    sol = K \ rhs
    dq_post = sol[1:5]
    return [q; dq_post]
  end

  # walker dynamics
  function walker_flow(
      M_f::Function,
      C_f::Function,
      G_f::Function,
      B::Matrix{Float64},
      x::Vector{<:DiffFloat},
      u::Vector{<:DiffFloat}
  )::Vector{<:DiffFloat}
      q = x[1:5]
      q̇ = x[6:10]
      Mq = M_f(q)
      Cq = C_f(q, q̇)
      Gq = G_f(q)
      return [q̇; Mq \ (B * u - Cq * q̇ - Gq)]
  end

  flow = (
      x::Vector{<:DiffFloat},
      u::Vector{<:DiffFloat}
  ) -> walker_flow(eval_M, eval_C, eval_G, B, x, u)

  # Hybrid mode: single support
  stance_mode = HybridMode(flow)

  # Impact: swing foot strikes ground
  g_impact = x::Vector{<:DiffFloat} -> swing_foot_height(x[1:5])
  R_impact = x::Vector{<:DiffFloat} -> reset_velocity_after_impact(x)

  impact = Transition(flow, flow, g_impact, R_impact)
  add_transition!(stance_mode, stance_mode, impact)

  transitions = Dict(:impact => impact)
  modes = Dict(:stance => stance_mode)

  return HybridSystem(nx, nu, transitions, modes)
end
