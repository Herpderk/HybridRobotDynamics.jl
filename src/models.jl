"""
    bouncing_ball(e=1.0, g=9.81)

Returns the hybrid system model containing the modes and transitions of a planar (in)elastic bouncing ball.
"""
function bouncing_ball(
    m::Real = 1.0,
    e::Real = 1.0,
    g::Real = 9.81
)::HybridSystem
    # Model ballistic dynamics with thrust
    # State space: x, y, xdot, ydot
    nx = 4
    nu = 1
    nq = 2

    # Manipulator equation form
    B = (q, q̇) -> reshape([0.0, 1.0], 2, 1)::Matrix{<:DiffFloat}
    M = q -> Matrix{Float64}(m * I(2))::Matrix{<:DiffFloat}
    c = (q, q̇) -> m * [0.0, g]::Vector{<:DiffFloat}

    # We can define a ControlAffineFlow like this:
    #   manip = ManipulatorEquation(manip_args...)
    #   flow = ControlAffineFlow(manip)
    # Or like this:
    #   flow = ControlAffineFlow(manip_args...)

    # Hybrid modes
    qidx = 1:2
    q̇idx = 3:4
    flight_flow = ControlAffineFlow(qidx, q̇idx, B, M, c)
    up_mode = HybridMode(flight_flow)
    down_mode = HybridMode(flight_flow)

    # Apex transition
    g_apex = x -> x[4]::DiffFloat
    R_apex = x -> x::Vector{<:DiffFloat}
    apex = Transition(flight_flow, flight_flow, g_apex, R_apex)
    add_transition!(up_mode, down_mode, apex)

    # Define impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[2]::DiffFloat
    R_impact = x::Vector{<:DiffFloat} -> (
        [x[1]; 1e-9; x[3]; abs(e*x[4])]::Vector{<:DiffFloat}
    )
    impact = Transition(flight_flow, flight_flow, g_impact, R_impact)
    add_transition!(down_mode, up_mode, impact)

    # Create hybrid system dicts
    transitions = Dict(:apex => apex, :impact => impact)
    modes = Dict(:up => up_mode, :down => down_mode)
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
    function actuated(x::Vector{<:DiffFloat})::Matrix{<:DiffFloat}
        return [zeros(7,4); mass\K; J\B]
    end
    function unactuated(x::Vector{<:DiffFloat})::Vector{<:DiffFloat}
        q = x[4:7]
        v = x[8:10]
        ω = x[11:13]
        return [
            Qquat(q) * v;
            0.5 * Gquat(q) * ω;
            Qquat(q)'*g - cross(ω, v);
            -J \ cross(ω, J*ω)
        ]
    end
    flight_flow = ControlAffineFlow(actuated, unactuated)

    # Hybrid modes
    up_mode = HybridMode(flight_flow)
    down_mode = HybridMode(flight_flow)

    # Apex transition
    g_apex = x::Vector{<:DiffFloat} -> x[10]::DiffFloat
    R_apex = x::Vector{<:DiffFloat} -> x::Vector{<:DiffFloat}
    apex = Transition(flight_flow, flight_flow, g_apex, R_apex)
    add_transition!(up_mode, down_mode, apex)

    # Impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[3]::DiffFloat
    R_impact = x::Vector{<:DiffFloat} -> (
        [x[1:9]; abs(e*x[10]); x[11:13]]::Vector{<:DiffFloat}
    )
    impact = Transition(flight_flow, flight_flow, g_impact, R_impact)
    add_transition!(down_mode, up_mode, impact)

    # Create hybrid system
    nx = 13
    nu = 4
    transitions = Dict(:apex => apex, :impact => impact)
    modes = Dict(:up => up_mode, :down => down_mode)
    return HybridSystem(nx, nu, transitions, modes)
end

"""
"""
function hopper(
    m1::Real = 5.0,    # body mass
    m2::Real = 1.0,    # foot mass
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

    # Define flight mode
    grav_flight = [0; -g; 0; -g]
    flight_unactuated = x::Vector{<:DiffFloat} -> (
        [x[5:8]; grav_flight]::Vector{<:DiffFloat}
    )
    flight_actuated = x::Vector{<:DiffFloat} -> (
        (M \ B_flight(x))::Vector{<:DiffFloat}
    )
    flight_flow = ControlAffineFlow(flight_actuated, flight_unactuated)
    flight_mode = HybridMode(flight_flow)

    # Define stance mode
    grav_stance =  [0; -g; 0; 0]
    stance_unactuated = x::Vector{<:DiffFloat} -> (
        [x[5:8]; grav_stance]::Vector{<:DiffFloat}
    )
    stance_actuated = x::Vector{<:DiffFloat} -> (
        (M \ B_stance(x))::Vector{<:DiffFloat}
    )
    stance_flow = ControlAffineFlow(stance_actuated, stance_unactuated)
    stance_mode = HybridMode(stance_flow)

    # Define impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[4]
    R_impact = x::Vector{<:DiffFloat} -> [x[1:6]; zeros(2)]
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
