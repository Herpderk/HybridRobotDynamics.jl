"""
    SaltationMatrix(flow_I, flow_J, guard, reset)

Derives the saltation matrix function for a given hybrid transition.
"""
struct SaltationMatrix
    flow_I::Function
    flow_J::Function
    guard::Function
    reset::Function
end

"""
    saltation(x, u)

Computes the saltation matrix at the corresponding hybrid transition given a state and input.
"""
function (saltation::SaltationMatrix)(
    x::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat}
)::Matrix{<:DiffFloat}
    xJ = saltation.reset(x)
    g_grad = ForwardDiff.gradient(saltation.guard, x)
    R_jac = ForwardDiff.jacobian(saltation.reset, x)
    return (
        R_jac + (saltation.flow_J(xJ, u) - R_jac * saltation.flow_I(x, u))
              * g_grad' / (g_grad' * saltation.flow_I(x,u))
    )
end

"""
    Transition(flow_I, flow_J, guard, reset)

Contains all hybrid system objects pertaining to a hybrid transition. Automatically computes the corresponding saltation matrix expresison.
"""
struct Transition
    flow_I::Function
    flow_J::Function
    guard::Function
    reset::Function
    saltation::SaltationMatrix
    function Transition(
        flow_I::Function,
        flow_J::Function,
        guard::Function,
        reset::Function,
    )::Transition
        saltation = SaltationMatrix(flow_I, flow_J, guard, reset)
        return new(flow_I, flow_J, guard, reset, saltation)
    end
end

"""
    HybridMode(flow, transitions=Dict())

Contains the flow and feasible transitions of the corresponding hybrid mode. Transitions and adjacent modes are stored as key-value pairs.
"""
mutable struct HybridMode
    flow::Function
    transitions::Dict{Transition, HybridMode}
    function HybridMode(
        flow::Function,
        transitions::Dict{Transition, HybridMode} = Dict{Transition, HybridMode}(Dict())
    )::HybridMode
        return new(flow, transitions)
    end
end

"""
    add_transition!(mode_I, mode_J, transition)

Adds a given transition and adjacent mode_J to the feasible transition dictionary of mode_I.
"""
function add_transition!(
    mode_I::HybridMode,
    mode_J::HybridMode,
    transition::Transition
)::Nothing
    mode_I.transitions[transition] = mode_J
    return nothing
end

"""
    add_transition!(mode_I, mode_J, guard, reset)

Constructs a transition using the given guard, reset, and modes and adds it to the feasible transition dictionary of mode_I.
"""
function add_transition!(
    mode_I::HybridMode,
    mode_J::HybridMode,
    guard::Function,
    reset::Function
)::Nothing
    transition = Transition(mode_I.flow, mode_J.flow, guard, reset)
    add_transition!(mode_I, mode_J, transition)
    return nothing
end


"""
    HybridSystem(
        nx, nu, transitions;
        stage_ineq_constr=nothing, stage_eq_constr=nothing,
        term_ineq_constr=nothing, term_eq_constr=nothing
    )

Contains all hybrid system objects in addition to the system's state and input dimensions. Assumes the following forms for stage and terminal constraints:

stage_constr(x,u) (<)= 0

term_constr(x) (<)= 0
"""
mutable struct HybridSystem
    nx::Int
    nu::Int
    stage_ng::Int
    stage_nh::Int
    term_ng::Int
    term_nh::Int
    stage_ineq_constr::Union{Nothing, Function}
    stage_eq_constr::Union{Nothing, Function}
    term_ineq_constr::Union{Nothing, Function}
    term_eq_constr::Union{Nothing, Function}
    transitions::Dict{Symbol, Transition}
    modes::Dict{Symbol, HybridMode}
    function HybridSystem(
        nx::Int,
        nu::Int,
        transitions::Dict{Symbol, Transition},
        modes::Dict{Symbol, HybridMode};
        stage_ineq_constr::Union{Nothing, Function} = nothing,
        stage_eq_constr::Union{Nothing, Function} = nothing,
        term_ineq_constr::Union{Nothing, Function} = nothing,
        term_eq_constr::Union{Nothing, Function} = nothing
    )::HybridSystem
        # Pack stage and terminal constraints together
        stage_constrs = (stage_ineq_constr, stage_eq_constr)
        term_constrs = (term_ineq_constr, term_eq_constr)

        # Get dimensions of constraints
        xtest = zeros(nx)
        utest = zeros(nu)
        stage_dims = [
            isnothing(constr) ? 0 : length(constr(xtest, utest))
            for constr in stage_constrs
        ]
        term_dims = [
            isnothing(constr) ? 0 : length(constr(xtest))
            for constr in term_constrs
        ]
        return new(
            nx, nu, stage_dims..., term_dims...,
            stage_constrs..., term_constrs...,
            transitions, modes,
        )
    end
end

"""
    roll_out(system, integrator, N, Δt, us, x0, init_transition)

Simulates a given system forward in time given an explicit integrator, horizon parameters, control sequence, and initial conditions. Returns the rolled out state trajectory.
"""
function roll_out(
    system::HybridSystem,
    integrator::ExplicitIntegrator,
    N::Int,
    Δt::AbstractFloat,
    us::Vector{<:AbstractFloat},
    x0::Vector{<:AbstractFloat},
    init_mode::Symbol
)::Vector{<:AbstractFloat}
    # Init loop variables
    u_idx = [(1:system.nu) .+ (k-1)*system.nu for k = 1:N-1]
    xs = [zeros(system.nx) for k = 1:N]
    xs[1] = x0
    mode_I = system.modes[init_mode]

    # Roll out over time horizon
    for k = 1:N-1
        xk = xs[k]

        # Reset and update mode if guard is hit
        for (trans, mode_J) in mode_I.transitions
            if trans.guard(xk) <= 0.0
                xk = trans.reset(xk)
                mode_I = mode_J
                break
            end
        end

        # Integrate smooth dynamics
        xs[k+1] = integrator(mode_I.flow, xk, us[u_idx[k]], Δt)
    end
    return vcat(xs...)
end
