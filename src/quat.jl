const Hquat = [zeros(1, 3); I(3)]

"""
"""
function skew(a::Vector{<:DiffFloat})
    @assert size(a) == (3,)
    return [
        0.0 -a[3] a[2];
        a[3] 0.0 -a[1];
        -a[2] a[1] 0.0
    ]
end

"""
"""
function L_or_R(q::Vector, is_L::Bool)::Matrix
    @assert size(q) == (4,)
    s = q[1]
    v = q[2:4]
    sign = is_L ? 1 : -1
    return [s -v'; v s*I(3) + sign*skew(v)]
end

"""
"""
function Lquat(q::Vector{<:DiffFloat})::Matrix
    return L_or_R(q, true)
end

"""
"""
function Rquat(q::Vector{<:DiffFloat})::Matrix
    return L_or_R(q, false)
end

"""
"""
function Gquat(q::Vector{<:DiffFloat})::Matrix
    return Lquat(q) * Hquat
end

"""
"""
function Qquat(q::Vector{<:DiffFloat})::Matrix
    return Hquat' * Rquat(q)' * Lquat(q) * Hquat
end

"""
"""
function random_unit_quat()::Vector{Float64}
    unscaled_quat = randn(4)
    return unscaled_quat / norm(unscaled_quat)
end
