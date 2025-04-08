"""
"""
struct BlockIndices
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
    return new(Midx, Jidx, Λidx)
end

function BlockIndices(
    nq::Int,
    J::Function
)::BlockIndices
    qtest = zeros(nq)
    nJ = size(J(qtest, qtest))
    return BlockIndices(nq, nJ)
end

"""
"""
struct LagrangianDynamics
    nq::Int
    idx::BlockIndices
    M::Function     # M(q)
    c::Function     # c(q, q̇)
    J::Function     # J(q)
    J̇::Function     # J̇(q, q̇)
    B::Function     # B(q)
end

function LagrangianDynamics(
    nq::Int,
    M::Function,
    c::Function,
    J::Function,
    J̇::Function,
    B::Function
)::LagrangianDynamics
    idx = BlockIndices(nq, J)
    return new(nq, idx, M, c, J, J̇, B)
end

"""
"""
function get_lagrangian_block(
    dynamics::LagrangianDynamics,
    q::Vector{<:DiffFloat}
)::Matrix{<:DiffFloat}
    M = dynamics.M(q)
    J = dynamics.J(q)
    nΛ = size(J)[1]
    return [M J'; J zeros(nΛ, nΛ)]
end

"""
"""
function get_unforced_acceleration(
    dynamics::LagrangianDynamics,
    q::Vector{<:DiffFloat},
    q̇::Vector{<:DiffFloat}
)::Tuple{Vector{<:DiffFloat}, Matrix{<:DiffFloat}}
    block = get_lagrangian_block(dynamics, q)
    block_inv = inv(block)
    Minv = block_inv[dynamics.idx.M]
    Jinv = block_inv[dynamics.idx.J]
    J̇ = dynamics.J̇(q, q̇)
    c = dynamics.c(q, q̇)
    q̈u = -Minv*c - Jinv'*J̇*q̇
    return q̈u, Minv
end

"""
"""
function get_input_mapping(
    dynamics::LagrangianDynamics,
    Minv::Union{Nothing, Matrix{<:DiffFloat}},
    q::Vector{<:DiffFloat}
)::Matrix{<:DiffFloat}
    return Minv * dynamics.B(q)
end

"""
"""
function get_forced_acceleration(
    dynamics::LagrangianDynamics,
    Minv::Union{Nothing, Matrix{<:DiffFloat}},
    q::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    return get_input_mapping(dynamics, Minv, q) * u
end

"""
"""
function get_acceleration(
    dynamics::LagrangianDynamics,
    q::Vector{<:DiffFloat},
    q̇::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    q̈u, Minv = get_unforced_acceleration(dynamics, q, q̇)
    q̈f = get_forced_acceleration(dynamics, Minv, q, u)
    return q̈u + q̈f
end
