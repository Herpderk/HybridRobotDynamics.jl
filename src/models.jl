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

    # Define hybrid modes
    up_mode = HybridMode(flight_flow)
    down_mode = HybridMode(flight_flow)

    # Define apex transition
    g_apex = x::Vector{<:DiffFloat} -> x[4]
    R_apex = x::Vector{<:DiffFloat} -> x
    apex = Transition(flight_flow, flight_flow, g_apex, R_apex)
    add_transition!(up_mode, down_mode, apex)

    # Define impact transition
    g_impact = x::Vector{<:DiffFloat} -> x[2]
    R_impact = x::Vector{<:DiffFloat} -> [x[1:3]; -e * x[4]]
    impact = Transition(flight_flow, flight_flow, g_impact, R_impact)
    add_transition!(down_mode, up_mode, impact)

    # Create hybrid system dicts
    transitions = Dict(
        :apex => apex,
        :impact => impact
    )
    modes = Dict(
        :upwards => up_mode,
        :downwards => down_mode
    )
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
