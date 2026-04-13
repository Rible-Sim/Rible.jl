
using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using Test
using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using ElasticArrays
using TypeSortedCollections
using TypedTables
using Rotations
using FiniteDiff
using CoordinateTransformations
using OffsetArrays
using RecursiveArrayTools
using CircularArrays
using Interpolations
using DataStructures
using EponymTuples
import GeometryBasics as GB
using Match
using CSV, Tables
using StructArrays
using LoggingExtras
using Unitful
using FileIO, MeshIO
import Rible.Meshes

@testset "ES Examples" begin

    # --- Shared Includes ---
    # These files are used by multiple examples
    include("../../examples/bodies/make_3d_bar.jl")
    include("../../examples/bodies/make_3d_plate.jl")
    include("../../examples/bodies/make_3d_tri.jl")
    include("../../examples/bodies/new_deck.jl")
    
    # --- Simple Example ---
    @testset "Simple" begin
        # Referenced in main.jl line 182: ../../../examples/robots/simple.jl
        include("../../examples/robots/simple.jl")
        
        c = 0.0
        z0 = 0.2
        ωz = 5.0
        mbar = 0.05
        
        # Case 1: Fixed
        newsim_fixed = simple(; c, z0, ωz, mbar)
       
        dt = 1e-3
        tend = 100.
        
        RB.solve!(
            RB.DynamicsProblem(newsim_fixed; env=RB.EmptyEnv()),
            RB.DynamicsSolver(RB.Zhong06());
            dt, tspan=(0.0, tend), maxiters=5, ftol=1e-14
        )
        @test true 
        
        # Case 2: Freed
        newsim_freed = simple(; c, z0, ωz, mbar, free=true)
        RB.solve!(
            RB.DynamicsProblem(newsim_freed; env=RB.EmptyEnv()),
            RB.DynamicsSolver(RB.Zhong06());
            dt, tspan=(0.0, tend), ftol=1e-14
        )
        @test true
    end

    # --- Bridge Example ---
    @testset "Bridge" begin
        # Referenced in main.jl line 510: ../../examples/robots/bridge3d.jl (relative to RibleTensegrity)
        include("../examples/robots/bridge3d.jl")

        n = 2
        newbridge = bridge3d(;n)
        
        # Static Equilibrium check logic
        μ = RT.inverse_for_restlength(newbridge, newbridge, RB.NoField(); fmin=1e4, eps_rel=1e-10)
        
        newbridge_modes = deepcopy(newbridge)
        RT.set_restlen!(newbridge_modes, μ)
        RB.update!(newbridge_modes.structure, RB.NoField())
        tension_min,tension_max = (RT.get_cables_tension(newbridge_modes.structure) |> extrema)
        @test tension_min ≈ 9999.999995840772
        @test tension_max ≈ 153902.36778065306
        RB.GDR!(newbridge_modes, RB.Gravity();β=1e-4, res=1e-9)
        q = newbridge_modes.traj.q[end]
        RB.set_initial!(newbridge_modes,q,zero(q))
        isequilibrium,_ = RT.check_static_equilibrium_output_multipliers(newbridge_modes.structure,RB.Gravity())
        @test isequilibrium
        # Check eigen
        eigenvals,_ = RB.undamped_eigen(newbridge_modes.structure,RB.Gravity())
        @test eigenvals ≈ [5.172897797597702, 12.080272827065247, 13.131973093116175, 22.311189506645164, 32.63682490310705, 71.27579981027554, 115.2437686464143, 117.2935023522459, 189.76928775686056, 264.6279658216149, 401.51753119866044, 408.8536915133325, 430.7450418234224, 453.5236028125474]
        
        # Dynamics with Load
        tend = 10.0
        dt = 1e-2
        ν = 0.453
        policy = RB.TimeFunctionPolicy((t) -> [2e4*(sin(ν*2π*t)+1)])
        newbridge_dyn = deepcopy(newbridge_modes)
        RB.solve!(
            RB.DynamicsProblem(newbridge_dyn; policy=policy),
            RB.DynamicsSolver(RB.Zhong06());
            dt, tspan=(0.0, tend), ftol=1e-6, exception=false
        )
        newbridge_dyn_b5p11 = RB.get_trajectory!(newbridge_dyn,5,11)
        @test newbridge_dyn_b5p11.u[end] ≈ [ -0.19196885402154784, -0.1189074548176583, 5.3543902307093205] 
    end

    # --- Embed Example ---
    @testset "Embed" begin
        # Referenced in main.jl line 779: ../../examples/robots/twotre3d.jl
        include("../examples/robots/twotre3d.jl")
        # Referenced in main.jl line 781: ../../../examples/robots/embed3d.jl
        include("../../examples/robots/embed3d.jl")
        
        m = 3
        α = 2π/m
        θ = 1.25α
        n = 4
        b = 0.14
        r = 0.04*sqrt(2)
        
        # newembed1
        newembed1 = embed3d(;
            r1= 0.03*sqrt(2),
            r,b,m,α,θ,n,
            outer = false,
        )
        μ1 = RT.inverse_for_restlength(newembed1, newembed1, RB.NoField(); fmin=100.0, eps_rel=1e-10)
        RT.set_restlen!(newembed1, μ1)
        RB.update!(newembed1.structure)
        κ = RT.get_cables_stiffness(newembed1.structure)
        f = RT.get_cables_tension(newembed1.structure) 
        @test minimum(f) ≈ 99.99999999999395 
        @test maximum(f) ≈ 106.49567191511044 
        isequilibrium,_ = RB.check_static_equilibrium_output_multipliers(newembed1.structure, RB.NoField())
        @test isequilibrium
        
        # Folded / Outer logic
        newembed0 = embed3d(;
             r1= 0.06*sqrt(2),
             r,b,m,α,θ,n
        )
        μ0 = RT.inverse_for_restlength(newembed0, newembed0, RB.NoField(); fmin=100.0, eps_rel=1e-10)
        
        newembed0_outer = embed3d(;
            r1= 0.06*sqrt(2),
            r,b,m,α,θ,n, outer = true
        )
        μ0o = RT.get_cables_restlen(newembed0_outer)
        μ0o[begin:3n*m] .= μ0
        μ1o = deepcopy(μ0o)
        μ1o[begin:3n*m] .= μ1
        RT.set_restlen!(newembed0_outer,μ0o)
        RB.GDR!(newembed0_outer, RB.Gravity())
        RB.set_initial!(newembed0_outer,newembed0_outer.traj.q[end],zero(newembed0_outer.traj.q[end]))
        RB.update!(newembed0_outer.structure,RB.Gravity())
        RB.check_static_equilibrium_output_multipliers(newembed0_outer.structure,RB.Gravity())
        eigenvals,_ = RB.undamped_eigen(newembed0_outer.structure,RB.Gravity())
        @test eigenvals[1:6] ≈ [714.1679254973039, 714.1679255067368, 1017.100607757013, 1542.0743114696083, 8962.462710104302, 12944.409367987175]
        # Policy construction logic
        function make_policy_test(μs, start, stop, len)
            nμ = length(μs[begin])
            RB.TimeFunctionPolicy(
                (
                    (t) -> begin
                         scaled_itps = extrapolate(
                            Interpolations.scale(
                                interpolate(
                                    reduce(hcat,μs),
                                    (NoInterp(),BSpline(Linear()))
                                ),
                                1:nμ, range(;start, stop, length = len,)
                            ),
                            (Throw(),Flat())
                        )
                        [scaled_itps(j,t) for j in 1:nμ]
                    end
                )
            )
        end
        
        μs_arr = [ begin μ = deepcopy(μ0o); μ[begin:3m*j] .= μ1[begin:3m*j]; μ end for j=0:n ]
        
        Td = 60.0
        tend = 60.0
        policy = make_policy_test(μs_arr, 0.0, Td, length(μs_arr))
        
        newembed0_outer_deploy = deepcopy(newembed0_outer)
        RB.solve!(
            RB.DynamicsProblem(newembed0_outer_deploy; policy=policy),
            RB.DynamicsSolver(RB.Zhong06());
            dt=4e-3, tspan=(0.0, tend), ftol=1e-10, exception=false
        )
        embed_b29p1 = RB.get_trajectory!(newembed0_outer_deploy,29,1)
        @test embed_b29p1.u[end] ≈ [-0.19486177830886528, -0.04729773686364182, 0.6920707228687331]
    end

end
