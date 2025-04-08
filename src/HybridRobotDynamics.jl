module HybridRobotDynamics

using LinearAlgebra
using ForwardDiff
using Plots

export
        add_transition!,
        Transition,
        HybridMode,
        HybridSystem,
        ExplicitIntegrator,
        roll_out,
        plot_2d_states,
        bouncing_ball,
        hopper

include("utils.jl")
include("integrators.jl")
include("dynamics.jl")
include("models.jl")
include("plotting.jl")

end # module HybridRobotDynamics
