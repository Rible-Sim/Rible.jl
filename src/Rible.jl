module Rible

# std
using LinearAlgebra
using SparseArrays
using Random
# arrays
using StaticArrays
using StructArrays
using RecursiveArrayTools
using SparseMatricesCSR
using BlockArrays
using BlockDiagonals
using SymmetricFormats
using TypeSortedCollections
using ComponentArrays
using Rotations
using CoordinateTransformations
using Tullio
# data structure
using Interpolations
using EponymTuples
using Unitful
import IterTools
using Match
using Tables
using TypedTables
# NN/diff
using ForwardDiff
using FiniteDiff
# solvers
using RowEchelon
using NLsolve
import LineSearches as LS
# print
using ProgressMeter
using Printf
using Logging
using LoggingExtras, Dates
# visualize/plot
import Makie
using LaTeXStrings
using Observables
import GeometryBasics as GB
import ColorTypes as CT
import MarchingCubes
import Meshes
const lift = map
# docs
using DocStringExtensions
# code
using PrecompileTools
using ConcreteStructs

using Pkg.Artifacts

const ASSETS_DIR = joinpath(@__DIR__, "..", "examples/assets")
const ARTIFACTS_TOML = joinpath(@__DIR__, "..", "Artifacts.toml")

function assetpath(files...)
    local_path = normpath(joinpath(ASSETS_DIR, files...))

    # Tier 1: file exists locally (test-required assets in git)
    if isfile(local_path)
        return local_path
    end

    # Tier 2: lazy-download via artifacts
    hash = artifact_hash("rible_extra_meshes", ARTIFACTS_TOML)
    if hash !== nothing
        ensure_artifact_installed("rible_extra_meshes", ARTIFACTS_TOML)
        return normpath(joinpath(artifact_path(hash), files...))
    end

    error("Asset not found: $(join(files, "/"))")
end

# Export main types and interfaces
export PresFreeCoordinates, CoordinatesState, TransLocalNonminimalCoordinatesState, InertiaCache

# Export rigid body types
export RigidBody, RigidBodyProperty, RigidBodyState, RigidBodyCache

# Export flexible body types  
export FlexibleBody, FlexibleBodyProperty, FlexibleBodyState, FlexibleBodyCoordinatesCache

# Export natural coordinate types
export NC, LNCData
export NC2D, NC3D
export NC2D2C, NC2D4C, NC2D6C, NC3D3C, NC3D6C, NC3D12C
export NC2D1P, NC2D1P1V, NC2D2P, NC1P2V, NC2P1V, NC3P
export NC3D1P, NC3D1P1V, NC3D2P, NC1P3V, NC2P2V, NC3P1V, NC4P

# Export interface methods (for extension)
export get_num_of_dims, get_num_of_coords, get_num_of_dof, get_num_of_cstr, get_num_of_intrinsic_cstr, get_num_of_local_dims
export to_local_coords, to_position, to_position_jacobian, to_velocity_jacobian
export find_rotation, find_angular_velocity, find_local_angular_velocity
export cstr_function, cstr_jacobian, cstr_hessians, cstr_forces_jacobian, cstr_velocity_jacobian
export cartesian_frame2coords, get_deform
export build_joint_cache, get_joint_violations!, get_joint_jacobian!
export add_joint_forces_jacobian!, get_joint_velocity_jacobian!
export kinetic_energy_coords, find_independent_free_idx, nullspace_mat
export make_cstr_function, make_cstr_jacobian
export has_constant_mass_matrix

# Export body interface methods
export get_id, get_mass, get_inertia
export update_inertia_cache!, lazy_update_state!, update_state!
export update_loci_states!, stretch_loci!, update_transformations!
export kinetic_energy, potential_energy_field, potential_energy_strain
export apply_field!, field_jacobian!, clear_forces!, reset!
export body_state2coords_state, strain!
export BodyJoint

# Export structure update pipeline
export update!, lazy_update!
export stretch!, update_bodies!, lazy_update_bodies!, update_apparatuses!
export assemble_forces!, assemble_M, assemble_M!

# Export constraint interface (structure-level)
export cstr_function!, cstr_jacobian!

# Export post-processing
export mechanical_energy!, get_mid_velocity!, get_trajectory!, get_velocity!, get_bodies

# Export visualization
export plot_traj!, savefig

# Export robot methods
export set_initial!, goto_step!, set_robot_state!

# Export locus, connection, and joint types
export Locus, Anchor, Hen2Egg, ProtoJoint
export NCF

# Export onboarding/public example API
export StaticContactSurfaces, HalfSpace, GravityEnv
export DynamicsProblem, DynamicsSolver, Simulator, solve!
export RestitutionFrictionCombined, NewtonRestitution, CoulombFriction
export Zhong06, InnerLayerContactSolver, InteriorPointMethod
export Robot, Structure, Connectivity, ControlHub, Coalition
export Locus, Signifier, CaptumGauge, ErrorGauge
export PositionCaptum, PosVelCaptum, ExternalForceActuator, NaiveOperator
export measure, measure!, measure_jacobian!, measure_gradient!

# Export objective and cost
export Objective, cost!, get_trajectory_cost_weights, get_num_of_actions
export cost_gradient!, cost_hessian!

# Export sensitivity solvers
export DiscreteAdjointDynamicsSolver, AdjointDynamicsSensitivitySolver
export DirectDynamicsSensitivitySolver

# Export integrators
export RKIntegrator

include("utils/utils.jl")

# Include abstract types and interfaces first
include("abstract_types.jl")

include("environments/env.jl")

include("environments/geometries/surfaces.jl")

include("environments/geometries/rigids.jl")


# Import interfaces
include("coordinates/coordinate_interfaces.jl")
# Include coordinate systems
include("coordinates/nonminimal_coordinates.jl")
include("coordinates/contact_state.jl")
include("coordinates/natural_coordinates/NCF.jl")
using .NCF

# Include body implementations
include("body/types.jl")
include("body/body_interfaces.jl")
include("body/implementations.jl")
include("body/coordinate_specific.jl")

include("force/spring_dampers.jl")
include("force/rheonomic_joint.jl")

include("joint/joints.jl")
include("joint/interfaces.jl")

include("apparatus/apparatus.jl")

include("structures/connectivity.jl")
include("structures/structure.jl")
include("structures/interfaces.jl")
include("structures/mass_matrices.jl")
include("structures/constraints.jl")
include("structures/auxiliary.jl")
include("structures/mutate.jl")

include("structures/linearization.jl")
include("robots/robot.jl")

include("coordinates/presfree_coordinates.jl")
include("structures/pres_free_state.jl")


include("mechanics/materials.jl")
include("mechanics/linearized_mechanics.jl")
include("mechanics/dynamic_relax.jl")
include("mechanics/generalized_forces.jl")
include("mechanics/contact.jl")

include("dynamics_solvers/solvers.jl")

include("policy/basics.jl")
include("policy/basis_openloop.jl")
include("policy/legacy_static_linear.jl")
include("policy/legacy_trajectory_policy.jl")
include("policy/feedback_static_linear.jl")
include("policy/feedback_proportional_derivative.jl")
include("policy/discrete_linear.jl")
include("policy/discrete_pd.jl")


include("control/control.jl")
include("control/cost.jl")
include("control/cost_helpers.jl")

include("postprocess/analysis.jl")

include("postprocess/vis/mesh.jl")
include("postprocess/vis/recipe.jl")
include("postprocess/vis/plot_traj.jl")

include("precompile.jl")

end
