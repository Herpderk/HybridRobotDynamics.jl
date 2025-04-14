module Explicit
    using ..HybridRobotDynamics: DiffFloat

    """
        rk4(dynamics, x0, u0, x1, Δt)

    Explicit integrator using the fourth order Runge-Kutta method. Written as an equality constraint.
    """
    function rk4(
        dynamics::Function,
        x0::Vector{<:DiffFloat},
        u0::Vector{<:DiffFloat},
        Δt::DiffFloat
    )::Vector{<:DiffFloat}
        k1 = dynamics(x0, u0)
        k2 = dynamics(x0 + Δt/2*k1, u0)
        k3 = dynamics(x0 + Δt/2*k2, u0)
        k4 = dynamics(x0 + Δt*k3, u0)
        return x0 + Δt/6*(k1 + 2*k2 + 2*k3 + k4)
    end
end

"""
    ExplicitIntegrator(method_name)

Callable struct that instantiates its integration method based on the corresponding integrator in the `Explicit` module. When called on, returns the forward-simulated state given a set of primal variables.
"""
struct ExplicitIntegrator
    method::Function
end

function ExplicitIntegrator(
    method_name::Symbol
)::ExplicitIntegrator
    method = get_module_function(Explicit, method_name)
    return ExplicitIntegrator(method)
end

"""
    integrator(dynamics, primals)

Callable struct method for the `ExplicitIntegrator` and `ImplicitIntegrator` structs. Returns either the forward-simulated state or defect residuals between states at adjacent time steps, respectively.
"""
function (integrator::ExplicitIntegrator)(
    dynamics::Function,
    x::Vector{<:DiffFloat},
    u::Vector{<:DiffFloat},
    Δt::DiffFloat
)::Vector{<:DiffFloat}
    return integrator.method(dynamics, x, u, Δt)
end
