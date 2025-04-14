"""
"""
struct ControlAffineFlow
    actuated::Function
    unactuated::Function
end

"""
"""
function (flow::ControlAffineFlow)(
    x::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat}
)::Vector{<:DiffFloat}
    return flow.unactuated(x) + flow.actuated(x)*u
end

const Flow = Union{ControlAffineFlow, Function}
