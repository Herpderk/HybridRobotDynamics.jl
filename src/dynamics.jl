"""
"""
struct BlockIndices
    nΛ::Int
    M::Tuple{UnitRange{Int}, UnitRange{Int}}
    J::Tuple{UnitRange{Int}, UnitRange{Int}}
    Λ::Tuple{UnitRange{Int}, UnitRange{Int}}
end

function BlockIndices(
    nq::Int,
    nJ::Tuple{Int, Int}
)::BlockIndices
    Midx = (1:nq, 1:nq)
    Jidx = (nq .+ (1:nJ[1]), 1:nJ[2])
    Λidx = (nq .+ (1:nJ[1]), nq .+ (1:nJ[1]))
    nΛ = nJ[1]
    return BlockIndices(nΛ, Midx, Jidx, Λidx)
end

function BlockIndices(
    nq::Int,
    J::Function
)::BlockIndices
    qtest = zeros(nq)
    nJ = size(J(qtest))
    return BlockIndices(nq, nJ)
end


"""
"""
struct ManipulatorEquation
    nq::Int
    qidx::UnitRange{Int}
    q̇idx::UnitRange{Int}
    blockidx::BlockIndices
    B::Function     # B(q, q̇)
    M::Function     # M(q)
    c::Function     # c(q, q̇)
    J::Function     # J(q)
    J̇::Function     # J̇(q, q̇)
end

function ManipulatorEquation(
    qidx::UnitRange{Int},
    q̇idx::UnitRange{Int},
    B::Function,
    M::Function,
    c::Function,
    J::Function,
    J̇::Function
)::ManipulatorEquation
    # Assert position and velocity sizes
    nq = length(qidx)
    nq̇ = length(q̇idx)
    @assert nq == nq̇

    # Get sample function outputs
    qtest = zeros(nq)
    Btest = B(qtest, qtest)
    Mtest = M(qtest)
    ctest = c(qtest, qtest)
    Jtest = J(qtest)
    J̇test = J̇(qtest, qtest)

    # Assert function output sizes
    @assert size(Btest)[1] == nq
    @assert size(Mtest) == (nq, nq)
    @assert size(ctest) == (nq,)
    @assert size(Jtest)[2] == size(J̇test)[2] == nq
    @assert size(Jtest) == size(J̇test)

    # Assert function output types
    @assert typeof(Btest) <: Matrix{<:DiffFloat}
    @assert typeof(Mtest) <: Matrix{<:DiffFloat}
    @assert typeof(ctest) <: Vector{<:DiffFloat}
    @assert typeof(Jtest) <: Matrix{<:DiffFloat}
    @assert typeof(J̇test) <: Matrix{<:DiffFloat}

    blockidx = BlockIndices(nq, J)
    return ManipulatorEquation(nq, qidx, q̇idx, blockidx, B, M, c, J, J̇)
end

function ManipulatorEquation(
    qidx::UnitRange{Int},
    q̇idx::UnitRange{Int},
    B::Function,
    M::Function,
    c::Function
)::ManipulatorEquation
    nq = length(qidx)
    J = q::Vector{<:DiffFloat} -> zeros(0, nq)::Matrix{<:DiffFloat}
    J̇ = (q::Vector{<:DiffFloat}, q̇::Vector{<:DiffFloat}) -> (
        zeros(0, nq)::Matrix{<:DiffFloat}
    )
    return ManipulatorEquation(qidx, q̇idx, B, M, c, J, J̇)
end


"""
"""
function manipulator_block(
    manip::ManipulatorEquation,
    x::Vector{<:DiffFloat}
)::Matrix{<:DiffFloat}
    q = x[manip.qidx]
    M = manip.M(q)
    J = manip.J(q)
    Λ = zeros(manip.blockidx.nΛ, manip.blockidx.nΛ)
    return [M J'; J Λ]
end

"""
"""
function manipulator_inverses(
    manip::ManipulatorEquation,
    x::Vector{<:DiffFloat}
)::Tuple{Matrix{<:DiffFloat}, Matrix{<:DiffFloat}, Matrix{<:DiffFloat}}
    block = manipulator_block(manip, x)
    blockinv = block \ I
    Minv = blockinv[manip.blockidx.M...]
    Jinv = blockinv[manip.blockidx.J...]
    Λinv = blockinv[manip.blockidx.Λ...]
    return Minv, Jinv, Λinv
end

"""
"""
function unactuated_acceleration(
    manip::ManipulatorEquation,
    Minv::Matrix{<:DiffFloat},
    Jinv::Matrix{<:DiffFloat},
    x::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    q = x[manip.qidx]
    q̇ = x[manip.q̇idx]
    J̇ = manip.J̇(q, q̇)
    c = manip.c(q, q̇)
    q̈u = -Minv*c - Jinv'*J̇*q̇
    return q̈u
end

"""
"""
function unactuated_acceleration(
    manip::ManipulatorEquation,
    x::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    Minv, Jinv, Λinv = manipulator_inverses(manip, x)
    return unactuated_acceleration(manip, Minv, Jinv, x)
end

"""
"""
function actuation_mapping(
    manip::ManipulatorEquation,
    Minv::Matrix{<:DiffFloat},
    x::Vector{<:DiffFloat}
)::Matrix{<:DiffFloat}
    q = x[manip.qidx]
    q̇ = x[manip.q̇idx]
    return Minv * manip.B(q, q̇)
end

"""
"""
function actuation_mapping(
    manip::ManipulatorEquation,
    x::Vector{<:DiffFloat}
)::Matrix{<:DiffFloat}
    Minv, Jinv, Λinv = manipulator_inverses(manip, x)
    return actuation_mapping(manip, Minv, x)
end

"""
"""
function acceleration(
    manip::ManipulatorEquation,
    x::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    Minv, Jinv, Λinv = manipulator_inverses(manip, x)
    q̈u = unactuated_acceleration(manip, Minv, Jinv, x)
    q̈f = actuation_mapping(manip, Minv, x) * u
    return q̈u + q̈f
end


"""
"""
struct ControlAffineFlow
    manip::Union{ManipulatorEquation, Nothing}
    actuated::Function
    unactuated::Function
end

function ControlAffineFlow(
    actuated::Function,
    unactuated::Function
)::ControlAffineFlow
    return ControlAffineFlow(nothing, actuated, unactuated)
end

function ControlAffineFlow(
    manip::ManipulatorEquation
)::ControlAffineFlow
    qtest = zeros(manip.nq)
    nu = size(manip.B(qtest, qtest))[2]

    actuated = x::Vector{<:DiffFloat} -> (
        [zeros(manip.nq, nu); actuation_mapping(manip, x)]
    )::Matrix{<:DiffFloat}

    unactuated = x::Vector{<:DiffFloat} -> (
        [x[manip.q̇idx]; unactuated_acceleration(manip, x)]
    )::Vector{<:DiffFloat}

    return ControlAffineFlow(manip, actuated, unactuated)
end

function ControlAffineFlow(
    qidx::UnitRange{Int},
    q̇idx::UnitRange{Int},
    manipfuncs::Vararg{Function}
)::ControlAffineFlow
    return ControlAffineFlow(ManipulatorEquation(qidx, q̇idx, manipfuncs...))
end


"""
"""
function (flow::ControlAffineFlow)(
    x::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    if typeof(flow.manip) === ManipulatorEquation
        q̇ = x[flow.manip.q̇idx]
        q̈ = acceleration(flow.manip, x, u)
        return [q̇; q̈]
    else
        return flow.unactuated(x) + flow.actuated(x)*u
    end
end

const Flow = Union{ControlAffineFlow, Function}
