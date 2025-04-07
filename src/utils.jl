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
    g_grad = FD.gradient(saltation.guard, x)
    R_jac = FD.jacobian(saltation.reset, x)
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
    transitions::Dict{Transition, HybridMode} = Dict()
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
    transitions::Dict{Symbol, Transition}
    modes::Dict{Symbol, HybridMode}
    stage_ineq_constr::Union{Nothing, Function}
    stage_eq_constr::Union{Nothing, Function}
    term_ineq_constr::Union{Nothing, Function}
    term_eq_constr::Union{Nothing, Function}
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
            nx, nu, transitions, modes,
            stage_dims..., term_dims..., stage_constrs..., term_constrs...,
        )
    end
end
