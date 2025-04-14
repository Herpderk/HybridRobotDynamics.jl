module HybridRobotDynamics

using LinearAlgebra
using ForwardDiff
using ModelingToolkit
using Symbolics
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
        bouncing_quadrotor,
        five_link_walker,
        hopper

include("utils.jl")
include("integrators.jl")
include("lagrangian.jl")
include("hybrid.jl")
include("quat.jl")
include("models.jl")
include("plotting.jl")

end # module HybridRobotDynamics
