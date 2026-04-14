using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT

using LinearAlgebra
using SparseArrays
using StaticArrays
using ElasticArrays
using TypeSortedCollections
using Rotations
using FiniteDiff

using RecursiveArrayTools

using EponymTuples
import GeometryBasics as GB
using Test
using StructArrays
using LoggingExtras
using EponymTuples
using ElasticArrays
using FileIO, MeshIO
import Rible.Meshes

try
    @eval using AbbreviatedStackTraces;
    @eval ENV["JULIA_STACKTRACE_ABBREVIATED"] = true;
    @eval ENV["JULIA_STACKTRACE_MINIMAL"] = true;
catch e
    @warn "error while importing AbbreviatedStackTraces" e
end

include("../../examples/robots/define_prototype.jl");

function actuate_bezier!(P3; t1)
    function bezier(t)
        P0 = 0    # Start position is 0 mm  
        P1 = 3   # The first control point can be set to 20 mm  
        P2 = P3 - 3

        if 0 <= t <= 1
            return 0.0
        elseif 1 < t <= t1
            return (1 - t / t1)^3 * P0 +
                   3 * (1 - t / t1)^2 * t / t1 * P1 +
                   3 * (1 - t / t1) * (t / t1)^2 * P2 +
                   (t / t1)^3 * P3
        else
            return P3
        end
        # Cubic Bézier curve formula   
    end

    RB.TimeFunctionPolicy(
        (
            # (t) -> [0.0, 0.1t, 0.1t]
            (t) -> [
                # 0,0,0,0,
                0.5bezier(t),
                0.5bezier(t),
                -0.25bezier(t),
                -0.25bezier(t),
                1.2bezier(t),
                1.2bezier(t),
                -0.8bezier(t),
                -0.8bezier(t),
                # 0
                sin(t/50)
            ]
        )
    )
end

global_logger(ConsoleLogger(stdout, Logging.Info; show_limited=false))

@testset "Cluster and RheonomicJointForce" begin

    dt = 5e-2; T = 40.0; T1 = 10.0

    bot_2seg = build_2_seg(2, [60.0, 60.0], [0, 90.0])
    policy = actuate_bezier!(40; t1=T)
    # policy = RB.NoPolicy()
    
    RB.solve!(
        RB.DynamicsProblem(bot_2seg;
            policy=policy,
            env=RB.GravityEnv()
            # RB.EmptyEnv()

        ),
        RB.DynamicsSolver(
            RB.Zhong06(),
            RB.MonolithicApparatusSolver(
                RB.SmoothedFischerBurmeister()
            ),
        );
        dt,
        tspan=(0.0, T + T1),
        ftol=1e-7, maxiters=300, verbose=false, exception=false
    )

    #note mid point [40.0103253138607,-59.175115505833126,-422.2716898268551]
    @test bot_2seg.traj[end].q[end-11:end-9] ≈ [39.98299975420927, -59.15317921373849, -422.2706219961591]  rtol = 1e-6
    #note midpoint [-0.6202082426844637,1.6935314796792662,-0.4810901968221382] 
    @test bot_2seg.traj[end].q̇[end-11:end-9] ≈ [-0.5991303756744005, 1.6575590163912663, -0.4736832879461835] rtol = 1e-4
end
