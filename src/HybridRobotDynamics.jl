module HybridRobotDynamics

using LinearAlgebra
using ForwardDiff
using Plots

export
        ManipulatorEquation,
        manipulator_inverses,
        unactuated_acceleration,
        actuation_mapping,
        ControlAffineFlow,
        Transition,
        HybridMode,
        add_transition!,
        HybridSystem,
        ExplicitIntegrator,
        roll_out,
        plot_2d_states,
        bouncing_ball,
        bouncing_quadrotor,
        hopper

include("utils.jl")
include("dynamics.jl")
include("integrators.jl")
include("hybrid.jl")
include("quat.jl")
include("models.jl")
include("plotting.jl")

end # module HybridRobotDynamics
