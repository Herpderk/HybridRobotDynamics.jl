module HybridRobotDynamics

using ForwardDiff
using RigidBodyDynamics

export
        add_transition!,
        Transition,
        HybridMode,
        HybridSystem,
        ExplicitIntegrator,
        roll_out,
        bouncing_ball,
        hopper

include("utils.jl")
include("integrators.jl")
include("dynamics.jl")
include("models.jl")

end # module HybridRobotDynamics
