using Test
using Rible
import Rible as RB
using RibleQCF
import RibleQCF as QCF
using RibleExtraIntegrators
using LinearAlgebra
using StaticArrays
using TypeSortedCollections
using StructArrays
using ForwardDiff
using Statistics
using Rotations
using CoordinateTransformations
using EponymTuples
using FileIO

# Load the robot definition
# Path relative to RibleExtraIntegrators/test/slider_crank.jl -> ../../examples/robots/slider_crank.jl
include("../../examples/robots/slider_crank.jl")

@testset "Slider Crank Moreau vs Zhong06" begin
    # Setup problem
    # Natural coordinates
    sc_ref = slider_crank(;coordsType= :NCF)
    sc_moreau = slider_crank(;coordsType= :NCF)
    
    # Simulation parameters
    # dt = 5e-4
    # tspan = (0.0, 1.0)
    dt = 1e-3
    tspan = (0.0, 0.1) 

    prob_ref = RB.DynamicsProblem(sc_ref; env=RB.GravityEnv())
    prob_moreau = RB.DynamicsProblem(sc_moreau; env=RB.GravityEnv())

    # Solve with Zhong06 (Reference)
    RB.solve!(
        prob_ref,
        RB.DynamicsSolver(RB.Zhong06());
        dt, tspan, ftol=1e-12, maxiters=50, verbose=false, exception=true, progress=false,
    )

    # Solve with Moreau
    RB.solve!(
        prob_moreau,
        RB.DynamicsSolver(RibleExtraIntegrators.Moreau(0.5));
        dt, tspan, ftol=1e-12, maxiters=50, verbose=false, exception=true, progress=false,
    )

    # Compare trajectories
    # Get positions of a specific point on the system, e.g., link 2 midpoint
    ref_traj = RB.get_trajectory!(sc_ref, 4, 2)
    moreau_traj = RB.get_trajectory!(sc_moreau, 4, 2)

    # Check that they produce similar results
    # Ideally, they should be relatively close for this smooth-ish problem (no contact enabled in this run)
    diff = norm(ref_traj - moreau_traj)
    @info "Trajectory difference (Norm): $diff"
    @test diff < 1e-4

    # Energy Check
    me_ref = RB.mechanical_energy!(sc_ref, RB.Gravity())
    me_moreau = RB.mechanical_energy!(sc_moreau, RB.Gravity())
    
    # Check energy conservation/stability
    # For Zhong06 (symplectic), energy should be well conserved
    energy_var_ref = std(me_ref.E)
    @info "Zhong06 Energy Variation (Std): $energy_var_ref"
    @test energy_var_ref < 1e-6

    # For Moreau, energy behavior might be different but shouldn't explode
    energy_var_moreau = std(me_moreau.E)
    @info "Moreau Energy Variation (Std): $energy_var_moreau"
    @test energy_var_moreau < 1e-4
    
    # Compare final energy
    final_energy_diff = abs(me_ref.E[end] - me_moreau.E[end])
    @info "Final Energy Difference: $final_energy_diff"
    @test final_energy_diff < 1e-4

end
