#mark 
#TestEnv.activate(); using Revise;
using Rible
using RibleTensegrity
import Rible as RB
import RibleTensegrity as RT
using Test
using LinearAlgebra
using SparseArrays
using StaticArrays
using CircularArrays
import CircularArrays as CA

using TypeSortedCollections
using LoggingExtras
using Rotations
using Unitful
using ElasticArrays
# IO
using FileIO, MeshIO
import GeometryBasics as GB
using EponymTuples
using StructArrays


# Include dependencies
# Adjust paths relative to RibleTensegrity/test/


# Include bodies/robots
include("../../examples/bodies/rigidbar.jl")
include("../../examples/robots/superball.jl")

include("../../examples/bodies/make_3d_bar.jl")
include("../../examples/bodies/make_3d_plate.jl")
include("../examples/robots/prism.jl")

include("../../examples/bodies/make_3d_tri.jl")
include("../examples/robots/tower.jl")

include("../../examples/bodies/build_2d_tri.jl")
include("../examples/robots/two_tri.jl")

include("../examples/robots/Tbars.jl")

const nofield = RB.NoField()
const gravity = RB.Gravity()

@testset "Stiffness Tests" begin

    @testset "Superball" begin
        # mark superball
        function build_nullspace_on_free(st)
            (;sys_full_coords_idx,bodyid2sys_dof_idx) = st.connectivity
            # Get intrinsic nullspace, excluding the first body
            Nin = RB.intrinsic_nullspace(st,st.state.system)[
                sys_full_coords_idx,
                reduce(vcat,bodyid2sys_dof_idx[2:end])
            ]
        end

        ballbot = superball(;
            θ = 0.0,
            l = 2.0/2,
            d = 2.0/2/2,
            z0 = 2.0/2,
            visible=true,
            loadmesh = false,
            constrained = false
        )
        bot = ballbot

        # Check static equilibrium
        RT.check_static_equilibrium_output_multipliers(bot.structure, nofield)
        
        RB.update!(bot.structure)
        f = RT.get_cables_tension(bot)

        function verify_lambda(st)
            T = RB.get_numbertype(st)
            λs = zeros(T,st.connectivity.num_of_bodies)
            foreach(st.bodies) do body
                (;prop,state) = body
                (;loci_states,origin_frame) = state
                for locus_state in loci_states
                    λs[prop.id] += -1/2*(locus_state.frame.position-origin_frame.position)'*locus_state.force
                end
            end
            λs
        end
        
        verify_lambda(bot.structure)
        q = RB.get_coords(bot.structure)
        q̌ = RB.get_free_coords(bot.structure)
        Ǎ = RB.cstr_jacobian(bot.structure,bot.structure.state.system)
        Ň = build_nullspace_on_free(bot.structure)
        
        Q̃ = RT.build_Q(bot.structure)
        L̂ = RT.build_L̂(bot.structure)

        @test rank(Ň) > 0 # Ensure non-trivial nullspace
        @test norm(Ǎ*Ň) < 1e-9

        Q̃L̂ = Q̃*L̂
        Bᵀ = -Q̃L̂
        ℬᵀ = transpose(Ň)*Bᵀ

        S,D = RT.static_kinematic_determine(ℬᵀ)
        ns = size(S,2)
        nk = size(D,2)
        
        @test ns == 1
        @test nk == 2
        k = RT.get_cables_stiffness(bot.structure)
        l = RT.get_cables_len(bot.structure)
        
        λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
        @test all(λ .≈ 2241.232864751818)
        Ǩa = -RB.cstr_forces_jacobian(bot.structure,q,λ)
        𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
        
        vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
        @test all(vals_𝒦a[1:10] .≈ -4482.465729503636)
        
        Ǩm = RT.build_material_stiffness_matrix!(bot.structure,q,k)
        Ǩg = RT.build_geometric_stiffness_matrix!(bot.structure,q,f)

        vec𝒦ps = [
            begin
                si = S[:,i]
                λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
                Ǩai = -RB.cstr_forces_jacobian(bot.structure,q,λi)
                Ǩgi = RT.build_geometric_stiffness_matrix!(bot.structure,q,si)
                𝒦pi = transpose(Ň)*(Ǩgi.+Ǩai)*Ň |> Symmetric 
                vec𝒦pi = vec(𝒦pi)
            end
            for i = 1:ns
        ]
        mat𝒦ps = reduce(hcat,vec𝒦ps)
        
        𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric 
        𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 
        𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
        𝒦 = 𝒦m.+ 𝒦p |> Symmetric
        
        vals_𝒦,vecs_𝒦 = eigen(𝒦)
        @test all(vals_𝒦 .> 0) # Should be mostly positive definite or zero
        
        # Optimization tests - reduced nullspace
        v = vecs_𝒦[:,1]
        Ňv = Ň*nullspace(v')
        
        r𝒦m = transpose(Ňv)*(Ǩm)*Ňv |> Symmetric
        vecr𝒦m = vec(r𝒦m)
        vecI = vec(Matrix(1.0I,size(r𝒦m)))
        
        r𝒦g = transpose(Ňv)*(Ǩg)*Ňv |> Symmetric
        r𝒦a = transpose(Ňv)*(Ǩa)*Ňv |> Symmetric
        r𝒦p = r𝒦g .+ r𝒦a
        
        vecr𝒦ps = [
            begin
                si = S[:,i]
                λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
                Ǩai = -RB.cstr_forces_jacobian(bot.structure,q,λi)
                Ǩgi = RT.build_geometric_stiffness_matrix!(bot.structure,q,si)
                r𝒦pi = transpose(Ňv)*(Ǩgi.+Ǩai)*Ňv |> Symmetric
                vec(r𝒦pi)
            end
            for i = 1:ns
        ]
        matr𝒦ps = reduce(hcat,vecr𝒦ps)
        
        # Maximum stiffness optimization
        ᾱ = [1.0]
        A = hcat(
            -Matrix(1.0I,ns,ns),
            ᾱ,
            zero(ᾱ)
        )
        b = [0.0]
        nx = ns+2
        result_max = RT.optimize_maximum_stiffness(matr𝒦ps,vecr𝒦m,vecI,A,b,nx)
        σ_max = result_max.x[end-1]
        ρ_max = result_max.x[end]
        
        @test σ_max ≈ 2025.3931308646224
        @test ρ_max ≈ 464.78761689658154
        
        r𝒦_max = r𝒦m + σ_max*reshape(matr𝒦ps*ᾱ,size(r𝒦m))
        vals_r𝒦_max, vecs_r𝒦_max = eigen(r𝒦_max)
        @test minimum(vals_r𝒦_max) ≈ ρ_max rtol=1e-6
        
        # Zero stiffness optimization
        result_zero = RT.optimize_zero_stiffness(matr𝒦ps,vecr𝒦m,vecI,
            hcat(
                -Matrix(1.0I,ns,ns),
                ᾱ,
            ),
            [0.0],
            ns+1,
            result_max.x[1:end-1]
        )
        σ_zero = result_zero.x[end]
        
        r𝒦_zero = r𝒦m + σ_zero*reshape(matr𝒦ps*ᾱ,size(r𝒦m))
        vals_r𝒦_zero, vecs_r𝒦_zero = eigen(r𝒦_zero)
        ρ_zero = vals_r𝒦_zero[1]
        
        @test abs(ρ_zero) ≈ 0.0 atol=1e-3
        @test σ_zero ≈ 4898.977685798693
    end

    @testset "Prism" begin
        # mark prism
        m = 3
        α = 2π/m
        θ = 1.25α
        n = 4
        b = 0.14
        r = 0.04*sqrt(2)
        prism1 = prisms(;
            r1= 0.03*sqrt(2),
            r,b,m,α,θ,n = 2,
            coordsType=RB.NCF.NC3D1P1V,
            barcolor = :slategrey,
            hasplate = true,
        )
        bot = prism1
        
        RT.check_static_equilibrium_output_multipliers(bot.structure, nofield)
        RB.update!(bot.structure)
        f = RT.get_cables_tension(bot)

        function build_nullspace_on_free(st)
            (;sys_free_coords_idx,bodyid2sys_full_coords,bodyid2sys_dof_idx) = st.connectivity
            Nin_unfixed = RB.intrinsic_nullspace(st,st.state.system)
            Nin = Nin_unfixed[
                sys_free_coords_idx,
                reduce(vcat,bodyid2sys_dof_idx[begin:end-1])
            ]
            # Build external constraint matrix
            Nex = zeros(eltype(q),30,12)
            for i = 1:6
                is = (i-1)*5
                js = (i-1)*2
                Nex[is+4:is+5,js+1:js+2] .= I(2)
            end
            cm = CircularArray(collect(1:3))
            for i = 1:3
                is = (3+cm[i+2]-1)*5
                js = (i-1)*2
                q_I = q[bodyid2sys_full_coords[i]]
                ri = @view q_I[1:3]
                u = @view q_I[4:6]
                v,w = RB.HouseholderOrthogonalization(u)
                Nex[is+1:is+3,js+1:js+2] = -RB.skew(0.14u)*[v w;]
            end
            Nin*Nex
        end

        q = RB.get_coords(bot.structure)
        q̌ = RB.get_free_coords(bot.structure)
        Ǎ = RB.cstr_jacobian(bot.structure,bot.structure.state.system)
        Ň = build_nullspace_on_free(bot.structure)
        
        Q̃ = RT.build_Q(bot.structure)
        L̂ = RT.build_L̂(bot.structure)
        
        @test rank(Ň) > 0
        @test norm(Ǎ*Ň) < 1e-9

        Q̃L̂ = Q̃*L̂
        Bᵀ = -Q̃L̂
        ℬᵀ = transpose(Ň)*Bᵀ
        
        S,D = RT.static_kinematic_determine(ℬᵀ)
        ns = size(S,2)
        nk = size(D,2)
        @test ns == 2
        @test nk == 2

        k = RT.get_cables_stiffness(bot.structure)
        Ǩm = RT.build_material_stiffness_matrix!(bot.structure,q,k)
        𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric
        vec𝒦m  = vec(𝒦m)
        vecI = vec(Matrix(1.0I,size(𝒦m)))
        
        ᾱ = [1.0,1.0]
        f = S*ᾱ 
        λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
        
        Ǩa = -RB.cstr_forces_jacobian(bot.structure,q,λ)
        𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
        
        Ǩg = RT.build_geometric_stiffness_matrix!(bot.structure,q,f)
        𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric
        𝒦p = 𝒦g.+ 𝒦a |> Symmetric
        𝒦 = 𝒦m.+ 𝒦p |> Symmetric
        
        vals_𝒦,vecs_𝒦 = eigen(𝒦)
        @test all(vals_𝒦 .> 0.0)
        
        # Optimization tests
        vec𝒦ps = [
            begin
                si = S[:,i]
                λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
                Ǩai = -RB.cstr_forces_jacobian(bot.structure,q,λi)
                Ǩgi = RT.build_geometric_stiffness_matrix!(bot.structure,q,si)
                𝒦pi = transpose(Ň)*(Ǩgi.+Ǩai)*Ň |> Symmetric
                vec(𝒦pi)
            end
            for i = 1:ns
        ]
        mat𝒦ps = reduce(hcat,vec𝒦ps)
        
        # Maximum stiffness optimization
        ᾱ = [1.0,1.0]
        A = hcat(
            -Matrix(1.0I,ns,ns),
            ᾱ,
            zero(ᾱ)
        )
        b = [0.0,0.0]
        nx = ns+2
        result_max = RT.optimize_maximum_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx)
        σ_max = result_max.x[end-1]
        ρ_max = result_max.x[end]
        
        @test σ_max ≈ 1431.125948398539
        @test ρ_max ≈ 41.44038830646161
        
        𝒦_max = 𝒦m + σ_max*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦_max, vecs_𝒦_max = eigen(𝒦_max)
        @test minimum(vals_𝒦_max) ≈ ρ_max rtol=1e-6
        
        # Zero stiffness optimization
        result_zero = RT.optimize_zero_stiffness(mat𝒦ps,vec𝒦m,vecI,
            hcat(
                -Matrix(1.0I,ns,ns),
                ᾱ,
            ),
            [0.0,0.0],
            ns+1,
            result_max.x[begin:end-1]
        )
        σ_zero = result_zero.x[end]
        
        𝒦_zero = 𝒦m + σ_zero*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦_zero, vecs_𝒦_zero = eigen(𝒦_zero)
        ρ_zero = vals_𝒦_zero[1]
        
        @test abs(ρ_zero) ≈ 0.001284515883509317 atol=1e-5
        @test σ_zero ≈ 1549.6228995399074
        
        # Second optimization set with different constraints
        σ̄ = [0,1]
        b = [500.0,0]
        A_α = hcat(
            -Matrix(1.0I,ns,ns),
            σ̄,
            zero(ᾱ)
        )
        
        result_max_α = RT.optimize_maximum_stiffness(mat𝒦ps,vec𝒦m,vecI,A_α,b,nx)
        σ_max_α = result_max_α.x[end-1]
        ρ_max_α = result_max_α.x[end]
        
        @test σ_max_α ≈ 891.393980846986
        @test ρ_max_α ≈ 26.881113750141196
        
        𝒦_max_α = 𝒦m + reshape(mat𝒦ps*(σ_max_α*σ̄+b),size(𝒦m))
        vals_𝒦_max_α, vecs_𝒦_max_α = eigen(𝒦_max_α)
        @test minimum(vals_𝒦_max_α) ≈ ρ_max_α rtol=1e-6
        
        result_zero_α = RT.optimize_zero_stiffness(mat𝒦ps,vec𝒦m,vecI,
            hcat(
                -Matrix(1.0I,ns,ns),
                σ̄,
            ),
            b,
            ns+1,
            result_max_α.x[begin:end-1]
        )
        σ_zero_α = result_zero_α.x[end]
        
        @test σ_zero_α ≈ 1658.5 atol = 1e-1
        
        𝒦_zero_α = 𝒦m + reshape(mat𝒦ps*(σ_zero_α*σ̄ + b),size(𝒦m))
        vals_𝒦_zero_α, vecs_𝒦_zero_α = eigen(𝒦_zero_α)
        ρ_zero_α = vals_𝒦_zero_α[1]
        
        @test ρ_zero_α ≈ 0 atol = 1e-1
    end

    @testset "Two Triangles" begin
        two = two_tri()
        bot = two
        RT.check_static_equilibrium_output_multipliers(bot.structure)

        function make_nullspace_on_free(st)    
            (;sys_free_coords_idx,bodyid2sys_dof_idx) = st.connectivity
            Nin = RB.intrinsic_nullspace(st,st.state.system)
            Nin[
                sys_free_coords_idx,
                reduce(vcat,bodyid2sys_dof_idx[2:end])
            ][:,end]
        end

        Ǎ = RB.cstr_jacobian(bot.structure,bot.structure.state.system)
        Ň = make_nullspace_on_free(bot.structure)
        
        @test rank(Ň) > 0
        @test norm(Ǎ*Ň) < 1e-9
        
        Q̃ = RT.build_Q(bot.structure)
        L̂ = RT.build_L̂(bot.structure)
        Q̃L̂ = Q̃*L̂
        Bᵀ = -Q̃L̂
        ℬᵀ = transpose(Ň)*Bᵀ
        
        S,D = RT.static_kinematic_determine(ℬᵀ;atol=1e-14)
        ns = size(S,2)
        nk = size(D,2)
        @test ns == 4
        @test nk == 0
    end

    @testset "Tower" begin
        towerbot = tower()
        bot = towerbot
        function build_nullspace_on_free(st)
            (;sys_free_coords_idx,bodyid2sys_full_coords,bodyid2sys_dof_idx) = st.connectivity
            Nin = RB.intrinsic_nullspace(st,st.state.system)[
                sys_free_coords_idx,
                reduce(vcat,bodyid2sys_dof_idx[begin+1:end])
            ]
            q = st.state.system.q
            Nex = zeros(eltype(q),11,5)
            is = 0; js = 0
            Nex[is+4:is+5,js+1:js+2] .= Matrix(1I,2,2)
            is = 5; js = 2
            Nex[is+4:is+6,js+1:js+3] .= Matrix(1I,3,3)

            q_I = q[bodyid2sys_full_coords[2]]
            u = @view q_I[4:6]
            v,w = RB.NCF.HouseholderOrthogonalization(u)
            is = 5; js = 0
            Nex[is+1:is+3,js+1:js+2] = -RB.skew(u)*[v w;]
            Nin*Nex
        end
        
        q = RB.get_coords(bot.structure)
        Ǎ = RB.cstr_jacobian(bot.structure,bot.structure.state.system)
        Ň = build_nullspace_on_free(bot.structure)
        
        @test rank(Ň) > 0
        @test norm(Ǎ*Ň) < 1e-9
        
        Q̃ = RT.build_Q(bot.structure)
        L̂ = RT.build_L̂(bot.structure)
        Q̃L̂ = Q̃*L̂
        Bᵀ = -Q̃L̂
        ℬᵀ = transpose(Ň)*Bᵀ
        
        S,D = RT.static_kinematic_determine(ℬᵀ)
        ns = size(S,2)
        nk = size(D,2)
        @test ns == 24
        @test nk == 0
        
        # Optimization tests
        k = RT.get_cables_stiffness(bot.structure)
        Ǩm = RT.build_material_stiffness_matrix!(bot.structure,q,k)
        𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric
        vec𝒦m = vec(𝒦m)
        vecI = vec(Matrix(1.0I,size(𝒦m)))
        
        # Build prestress stiffness matrices for each self-stress state
        vec𝒦ps = [
            begin
                si = S[:,i]
                λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
                Ǩai = -RB.cstr_forces_jacobian(bot.structure,q,λi)
                Ǩgi = RT.build_geometric_stiffness_matrix!(bot.structure,q,si)
                𝒦pi = transpose(Ň)*(Ǩgi.+Ǩai)*Ň |> Symmetric
                vec(𝒦pi)
            end
            for i = 1:ns
        ]
        mat𝒦ps = reduce(hcat,vec𝒦ps)
        
        # Use selected self-stress states (indices 8, 14, 24)
        isis = [8,14,24]
        ᾱ = zeros(ns)
        ᾱ[isis] .= [1,1,1]
        
        # Maximum stiffness optimization
        A = hcat(
            -Matrix(1.0I,ns,ns),
            ᾱ,
            zero(ᾱ)
        )
        b = zeros(ns)
        nx = ns+2
        result_max = RT.optimize_maximum_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx)
        σ_max = result_max.x[end-1]
        ρ_max = result_max.x[end]
        
        @test ρ_max ≈ 91.85359516039843
        
        𝒦_max = 𝒦m + σ_max*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦_max, vecs_𝒦_max = eigen(𝒦_max)
        @test minimum(vals_𝒦_max) ≈ ρ_max rtol=1e-6
        
        # Zero stiffness optimization
        result_zero = RT.optimize_zero_stiffness(mat𝒦ps,vec𝒦m,vecI,
            hcat(
                -Matrix(1.0I,ns,ns),
                ᾱ,
            ),
            zeros(ns),
            ns+1,
            result_max.x[1:end-1]
        )
        σ_zero = result_zero.x[end]
        
        @test σ_zero ≈ 95.88760950482862
        𝒦_zero = 𝒦m + σ_zero*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦_zero, vecs_𝒦_zero = eigen(𝒦_zero)
        ρ_zero = vals_𝒦_zero[1]
        
        @test abs(ρ_zero) < 1e-6  # Should be near zero
        @test σ_zero > σ_max  # Zero stiffness requires higher prestress
    end

    @testset "Tbars" begin
        tbbot = Tbars(;θ=0.0)
        bot = tbbot
        RT.check_static_equilibrium_output_multipliers(bot.structure)

        function nullspace_on_free(st)    
            (;sys_free_coords_idx,bodyid2sys_full_coords,bodyid2sys_dof_idx) = st.connectivity
            Nin = RB.intrinsic_nullspace(st,st.state.system)[
                sys_free_coords_idx,
                reduce(vcat,bodyid2sys_dof_idx[2:end])
            ]
            q = RB.get_coords(st)
            I3 = I(3)
            O3 = zero(I3)
            o3 = O3[:,1]
            qbar = q[bodyid2sys_full_coords[4]]
            ri = @view qbar[1:3]
            u = @view qbar[4:6]
            v,w = RB.NCF.HouseholderOrthogonalization(u)
            x,_,_ = ri
            Nslider1 = [
                o3;
                o3;;
            ] |> sparse
            Nslider2 = [
                 0;
                -x;
                 0;
                o3;;
            ] |> sparse
            Nbar = [
                o3;
                0;
                1;;
            ] |> sparse
            Nex = vcat(Nslider1,Nslider2,Nbar)
            Nin*Nex
        end
        q = RB.get_coords(bot.structure)
        Ǎ = RB.cstr_jacobian(bot.structure,bot.structure.state.system)
        Ň = nullspace_on_free(bot.structure)
        
        @test rank(Ň) > 0
        @test norm(Ǎ*Ň) < 1e-9
        
        Q̃ = RT.build_Q(bot.structure)
        L̂ = RT.build_L̂(bot.structure)
        Q̃L̂ = Q̃*L̂
        Bᵀ = -Q̃L̂
        ℬᵀ = transpose(Ň)*Bᵀ
        
        S,D = RT.static_kinematic_determine(ℬᵀ;atol=1e-14)
        ns = size(S,2)
        nk = size(D,2)
        @test ns == 4
        @test nk == 1
                
        # Get cable stiffness and lengths
        k = RT.get_cables_stiffness(bot.structure)

        l = RT.get_cables_len(bot.structure)

        # Build stiffness matrices for each self-stress state
        struct𝒦 = [
            begin
                s = S[:,i]        
                Ǩm = RT.build_material_stiffness_matrix!(bot.structure,q,s)
                𝒦m = transpose(Ň)*Ǩm*Ň 
                # s = S\f
                # @show s
                λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*s
                @show λ
                Ǩa = -RB.cstr_forces_jacobian(bot.structure,q,λ)

                Ǩg = RT.build_geometric_stiffness_matrix!(bot.structure,q,s)

                𝒦g = transpose(Ň)*(Ǩg)*Ň

                𝒦a = transpose(Ň)*(Ǩa)*Ň

                𝒦p = 𝒦g .+ 𝒦a
                (𝒦m = 𝒦m, 𝒦g = 𝒦g, 𝒦a = 𝒦a, 𝒦p = 𝒦p)
            end
            for i = 1:ns
        ] |> StructArray

        # Combine the stiffness matrices
        mat𝒦ms = reduce(hcat,struct𝒦.𝒦m)
        mat𝒦gs = reduce(hcat,struct𝒦.𝒦g)
        mat𝒦as = reduce(hcat,struct𝒦.𝒦a)
        mat𝒦ps = reduce(hcat,struct𝒦.𝒦p)
        @test mat𝒦ps ≈ [0.0  0.0  1.0  -1.0]
        
    end
end
