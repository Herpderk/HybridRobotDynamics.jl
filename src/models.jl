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
    qidx = 1:2
    q̇idx = 3:4

    # Manipulator equation form
    B = q -> reshape([0.0, 1.0], 2, 1)::Matrix{<:DiffFloat}
    M = q -> Matrix{Float64}(m * I(2))::Matrix{<:DiffFloat}
    c = (q, q̇) -> m * [0.0, g]::Vector{<:DiffFloat}

    # We can define a ControlAffineFlow like this:
    #   manip = ManipulatorEquation(manip_args...)
    #   flow = ControlAffineFlow(manip)
    # Or like this:
    #   flow = ControlAffineFlow(manip_args...)

    # Flow and hybrid mode
    flight_flow = ControlAffineFlow(qidx, q̇idx, B, M, c)
    flight_mode = HybridMode(flight_flow)

    # Impact transition
    g_impact = x -> x[2]::DiffFloat
    R_impact = x -> [x[1]; 1e-9; x[3]; abs(e*x[4])]::Vector{<:DiffFloat}
    impact = Transition(flight_flow, flight_flow, g_impact, R_impact)
    add_transition!(flight_mode, flight_mode, impact)

    # Hybrid system dicts
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

    # Hybrid mode
    flight_flow = ControlAffineFlow(actuated, unactuated)
    flight_mode = HybridMode(flight_flow)

    # Impact transition
    g_impact = x -> x[3]::DiffFloat
    R_impact = x -> (
        [x[1:2]; 1e-9; x[4:9]; abs(e*x[10]); x[11:13]]::Vector{<:DiffFloat}
    )
    impact = Transition(flight_flow, flight_flow, g_impact, R_impact)
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

    # q or x can be inputted into these auxiliary functions
    function get_length_vector(q::Vector{<:DiffFloat})::Vector{<:DiffFloat}
        return q[1:2] - q[3:4]
    end
    function get_unit_length(q::Vector{<:DiffFloat})::Vector{<:DiffFloat}
        L = get_length_vector(q)
        return L / norm(L)
    end

    # Manipulator equation form
    function B_flight(q::Vector{<:DiffFloat})::Matrix{<:DiffFloat}
        L1, L2 = get_unit_length(q)
        return [L1  L2; L2 -L1; -L1 -L2; -L2  L1]
    end

    function B_stance(x::Vector{<:DiffFloat})::Matrix{<:DiffFloat}
        L1, L2 = get_unit_length(x)
        return [L1  L2; L2 -L1; zeros(2,2)]
    end

    # Define flight mode
    grav_flight = [0; -g; 0; -g]
    flight_unactuated = x -> [x[5:8]; grav_flight]::Vector{<:DiffFloat}
    flight_actuated = x -> [zeros(4,2); M\B_flight(x)]::Matrix{<:DiffFloat}
    flight_flow = ControlAffineFlow(flight_actuated, flight_unactuated)
    flight_mode = HybridMode(flight_flow)

    # Define stance mode
    grav_stance =  [0; -g; 0; 0]
    stance_unactuated = x -> [x[5:8]; grav_stance]::Vector{<:DiffFloat}
    stance_actuated = x -> [zeros(4,2); M\B_stance(x)]::Matrix{<:DiffFloat}
    stance_flow = ControlAffineFlow(stance_actuated, stance_unactuated)
    stance_mode = HybridMode(stance_flow)

    # Impact transition
    g_impact = x -> x[4]::DiffFloat
    R_impact = x -> [x[1:3]; 1e-9; x[5:6]; zeros(2)]::Vector{<:DiffFloat}
    impact = Transition(flight_flow, stance_flow, g_impact, R_impact)
    add_transition!(flight_mode, stance_mode, impact)

    # Define liftoff transition
    g_liftoff = x -> -x[4]::DiffFloat
    R_liftoff = x -> x::Vector{<:DiffFloat}
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
    glb = x -> Llb - norm(get_length_vector(x))
    gub = x -> -Lub + norm(get_length_vector(x))
    glen = x -> [glb(x); gub(x)]
    g_stage = (x, u) -> glen(x)::Vector{<:DiffFloat}
    g_term = x -> glen(x)::Vector{<:DiffFloat}

    return HybridSystem(
        nx, nu, transitions, modes;
        stage_ineq_constr=g_stage, term_ineq_constr=g_term
    )
end
