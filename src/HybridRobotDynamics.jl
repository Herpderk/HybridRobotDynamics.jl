module HybridRobotDynamics

using RigidBodyDynamics

export
        add_transition!,
        Transition,
        HybridMode,
        HybridSystem,
        bouncing_ball,
        hopper

include("utils.jl")
include("models.jl")

end # module HybridRobotDynamics
