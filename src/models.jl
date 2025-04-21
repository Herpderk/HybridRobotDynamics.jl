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
"""
function five_link_walker(
  Is::Vector{<:Real},        # Link inertias
  ls::Vector{<:Real},        # Link lengths
  ms::Vector{<:Real},        # Link masses
  g::Real = 9.81             # Gravitational constant
)::HybridSystem

    nx = 14  # 5 positions + 5 velocities + 2 COM + 2 COM velocities
    nu = 2   # 2 actuators (hips)

    # Unpack link parameters
    l1, l2, l3, l4, l5 = ls

    # Actuation mapping (only joints actuated)
    B = [   0.0  0.0;  # x_com
            0.0  0.0;  # y_com
            0.0  0.0;  # θ₁ = left leg
            1.0  0.0;  # θ₂ = left thigh (u₁)
            0.0  0.0;  # θ₃ = right leg
            0.0  1.0;  # θ₄ = right thigh (u₂)
            0.0  0.0;  # θ₅ = torso
    ]

    function com_position(q::AbstractVector, i::Int)
        hip = q[1:2]
        θ₁, θ₂, θ₃, θ₄, θ₅ = q[3:7]    
    
        if i == 1  # left leg
            θ_th = θ₅ + θ₂
            θ_leg = θ_th + θ₁
            left_thigh = hip .- (l2/2) * [cos(θ_th), sin(θ_th)]
            return left_thigh .- (l1/2) * [cos(θ_leg), sin(θ_leg)]
    
        elseif i == 2  # left thigh
            θ_th = θ₅ + θ₂
            return hip .- (l2/2) * [cos(θ_th), sin(θ_th)]
    
        elseif i == 3  # right leg
            θ_th = θ₅ + θ₄
            θ_leg = θ_th + θ₃
            right_thigh = hip .- (l4/2) * [cos(θ_th), sin(θ_th)]
            return right_thigh .- (l3/2) * [cos(θ_leg), sin(θ_leg)]
    
        elseif i == 4  # right thigh
            θ_th = θ₅ + θ₄
            return hip .- (l4/2) * [cos(θ_th), sin(θ_th)]
    
        elseif i == 5  # torso
            return hip .+ (l5/2) * [cos(θ₅), sin(θ₅)]
    
        else
            error("Invalid link index")
        end
    end

     # left foot position and Jacobian
    function p_left(q::AbstractVector)
        x_com, y_com = q[1:2]
        θ₁, θ₂, _, _, θ₅ = q[3:7]
        l1, l2 = ls[1:2]
    
        θ_thigh = θ₅ + θ₂
        θ_leg   = θ_thigh + θ₁
    
        hip = [x_com, y_com]
        thigh = hip .- l2 * [cos(θ_thigh), sin(θ_thigh)]
        foot  = thigh .- l1 * [cos(θ_leg), sin(θ_leg)]
        return foot
    end

    function left_foot_height(x::AbstractVector)
        q = x[1:7]
        return p_left(q)[2]
    end

    # right foot position and Jacobian
    function p_right(q::AbstractVector)
        x_com, y_com = q[1:2]
        _, _, θ₃, θ₄, θ₅ = q[3:7]
        l3, l4 = ls[3:4]
    
        θ_thigh = θ₅ + θ₄
        θ_leg   = θ_thigh + θ₃
    
        hip = [x_com, y_com]
        thigh = hip .- l4 * [cos(θ_thigh), sin(θ_thigh)]
        foot  = thigh .- l3 * [cos(θ_leg), sin(θ_leg)]
        return foot
    end

    function right_foot_height(x::AbstractVector)
        q = x[1:7]
        return p_right(q)[2]
    end
    
    function total_com(q::AbstractVector)
        M_total = sum(ms)
        com = zeros(2)
        for i in 1:5
            com += ms[i] * com_position(q, i)
        end
        return com / M_total
    end

    function com_position_from_theta(θ::AbstractVector{<:Real}, i::Int)
        q = vcat(zeros(eltype(θ), 2), θ)
        return com_position(q, i)
    end

    function M(q::AbstractVector{<:Real})
        θ = q[3:7]
        θ̇₀ = zeros(eltype(q), 5)
    
        T_of_θ̇ = θ̇ -> begin
            T = 0.0
            for i in 1:5
                J_i = ForwardDiff.jacobian(θ -> com_position_from_theta(θ, i), θ)
                v_i = J_i * θ̇
                T += 0.5 * ms[i] * dot(v_i, v_i) + 0.5 * Is[i] * θ̇[i]^2
            end
            return T
        end
    
        M_joint = ForwardDiff.hessian(T_of_θ̇, θ̇₀)
    
        m_total = sum(ms)
        M_full = zeros(eltype(q), 7, 7)
        M_full[1:2, 1:2] .= m_total * I(2)     # ẋ_com and ẏ_com
        M_full[3:7, 3:7] .= M_joint            # joint-space
        return M_full
    end

    function C(q::AbstractVector{<:Real}, q̇::AbstractVector{<:Real})
        θ = q[3:7]
        θ̇ = q̇[3:7]
        n = length(θ)
    
        Cθ = zeros(eltype(q), n, n)
    
        # Mθ: M as a function of θ only
        function Mθ(θ₁::AbstractVector)
            q_full = vcat([0.0, 0.0], θ₁)
            return M(q_full)[3:7, 3:7]     # extract joint part
        end
    
        # Loop over Christoffel symbols to construct joint-space C matrix
        for i in 1:n, j in 1:n
            for k in 1:n
                ∂M_ij_∂θk = ForwardDiff.gradient(θ₁ -> Mθ(θ₁)[i, j], θ)[k]
                ∂M_ik_∂θj = ForwardDiff.gradient(θ₁ -> Mθ(θ₁)[i, k], θ)[j]
                ∂M_jk_∂θi = ForwardDiff.gradient(θ₁ -> Mθ(θ₁)[j, k], θ)[i]
    
                Cθ[i, j] += 0.5 * (∂M_ij_∂θk + ∂M_ik_∂θj - ∂M_jk_∂θi) * θ̇[k]
            end
        end
    
        # Embed into full 7×7 matrix
        C_full = zeros(eltype(q), 7, 7)
        C_full[3:7, 3:7] .= Cθ
    
        return C_full
    end
    
    function G(q::AbstractVector{<:Real})
        θ = q[3:7]
        Gθ = zeros(eltype(q), 5)
    
        for i in 1:5
            y_i = θ -> com_position(vcat(zeros(2), θ), i)[2]
            Gθ[i] = -ms[i] * g * ForwardDiff.gradient(y_i, θ)[i]
        end
    
        G_full = zeros(eltype(q), 7)
        G_full[3:7] .= Gθ
        G_full[2] = -sum(ms) * g   # gravity acts downward on total mass
        return G_full
    end

    function c(q, q̇)
      return C(q, q̇) * q̇ + G(q)
    end
  
    function reset_velocity_after_impact(x, p_foot_fn)
        q = x[1:7]
        q̇ = x[8:14]
    
        Mq = M(q)                        # 7x7

        Jc = ForwardDiff.jacobian(p_foot_fn, q)  # 2×7
    
        A = [Mq -Jc'; Jc zeros(2, 2)]
        rhs = [Mq * q̇; zeros(2)]
    
        sol = A \ rhs
        q̇_post = sol[1:7]
    
        return vcat(q, q̇_post)
    end

    function foot_jacobian_and_jacdot(p_foot, q, q̇)
      J = ForwardDiff.jacobian(p_foot, q)
      J̇ = ForwardDiff.jacobian(q_ -> ForwardDiff.jacobian(p_foot, q_) * q̇, q)
      return J, J̇
  end

    function walker_flow(
        M_f::Function,
        C_f::Function,
        G_f::Function,
        B::Matrix{Float64},  # 5×2 actuation mapping
        x::Vector{<:DiffFloat},  # full state: [q; q̇], 14×1
        u::Vector{<:DiffFloat},  # input torque: 2×1
        mode::Symbol             # symbol to discern current mode
    )::Vector{<:DiffFloat}
    
        # Split state
        q = x[1:7]
        q̇ = x[8:14]

        # Evaluate dynamics
        Mq = M_f(q)       # 7×7
        Cq = C_f(q, q̇)    # 7×7
        Gq = G_f(q)       # 7×1
        Bu = B * u        # 7×1

        if mode == :flight
            q̈ = Mq \ (Bu - Cq * q̇ - Gq)
            return vcat(q̇, q̈)
        end
        
        # Select contact Jacobian
        Jc, J̇c = (mode == :left)  ? foot_jacobian_and_jacdot(p_left, q, q̇) :
                (mode == :right) ? foot_jacobian_and_jacdot(p_right, q, q̇) :
                error("Unknown mode: $mode")

        # KKT system
        A = [Mq  -Jc'; 
            Jc   zeros(eltype(q), size(Jc, 1), size(Jc, 1))]

        rhs = [Bu - Cq * q̇ - Gq;
            -J̇c * q̇]

        sol = A \ rhs
        q̈ = sol[1:7]

        return vcat(q̇, q̈)
    end

    left_flow = (
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    ) -> walker_flow(M, C, G, B, x, u, :left)

    right_flow = (
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    ) -> walker_flow(M, C, G, B, x, u, :right)

    flight_flow = (
        x::Vector{<:DiffFloat},
        u::Vector{<:DiffFloat}
    ) -> walker_flow(M, C, G, B, x, u, :flight)

    modes = Dict(
        :flight  => HybridMode(flight_flow),
        :left    => HybridMode(left_flow),
        :right   => HybridMode(right_flow),
    )

    g_FL = x -> left_foot_height(x) 
    R_FL = x -> reset_velocity_after_impact(x, p_left)
    impact_FL = Transition(flight_flow, left_flow, g_FL, R_FL)

    g_RL = x -> left_foot_height(x) 
    R_RL = x -> reset_velocity_after_impact(x, p_left)
    impact_RL = Transition(right_flow, left_flow, g_RL, R_RL)
    
    g_FR = x -> right_foot_height(x) 
    R_FR = x -> reset_velocity_after_impact(x, p_right)
    impact_FR = Transition(flight_flow, right_flow, g_FR, R_FR)

    g_LR = x -> right_foot_height(x) 
    R_LR = x -> reset_velocity_after_impact(x, p_right)
    impact_LR = Transition(left_flow, right_flow, g_LR, R_LR)

    add_transition!(modes[:flight], modes[:left], impact_FL)
    add_transition!(modes[:flight], modes[:right], impact_FR)
    
    add_transition!(modes[:left], modes[:right], impact_LR)
    add_transition!(modes[:right], modes[:left], impact_RL)
    
    transitions = Dict(
        :impact_FL => impact_FL,
        :impact_FR => impact_FR,
        :impact_RL => impact_RL,
        :impact_LR => impact_LR,
    )

    return HybridSystem(nx, nu, transitions, modes)
end
