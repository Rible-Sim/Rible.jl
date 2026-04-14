module RibleTensegrity

using LinearAlgebra
using SparseArrays
using StaticArrays
using FiniteDiff
using DocStringExtensions
using Rible
using Makie
using EponymTuples
using TypeSortedCollections
using Printf
import JuMP, COSMO, Clarabel
using Polyhedra

import CDDLib 
lib = CDDLib.Library()

import GeometryBasics as GB

# Import specific Rible components for clarity
import Rible: StructureState, StructureCache, AbstractStructure, Robot, Apparatus, ApparatusCache, PrototypeJoint, AbstractField
import Rible: check_jacobian_singularity, check_constraints_consistency, PresFreeConnectivity
import Rible: Hen2Egg, Signifier, Anchor, update_bodies!, update_apparatuses!, apply_force_to_bodies!, DistanceSpringDamperState
import Rible: get_num_of_dims, get_num_of_coords, get_numbertype, get_coords, get_free_coords
import Rible: get_params,get_num_of_params, get_ids, check_id_sanity
import Rible: update!, switch!, lazy_update_bodies!,zero!
import Rible: clear_forces!, execute!, ExternalForceActuator, NaiveOperator
import Rible: assemble_forces!, mass_center_force!, centrifugal_force!, concentrated_force!,apply_field!, get_params!, get_auxiliary
import Rible: get_initial, Gravity, NoField, DistanceSpringDamper
import Rible: AbstractForce, RegisterActuator, FunctionOperator
import Rible: prepare_cache!, get_appar_idx
import Rible: make_cstr_jacobian, make_cstr_function, cstr_forces_jacobian, cstr_forces_jacobian!, add_cstr_forces_jacobian!
import Rible: auxi_function!, auxi_jacobian!, gen_force_auxi_jacobian!
import Rible: update_apparatus!, AbstractJoint
import Rible: Connectivity, AbstractConnectivity
import Rible: add_tangent_stiffness_matrix!,add_tangent_damping_matrix!
import Rible: gen_force_para_jacobian!,stretch!
import Rible: gen_force_actu_jacobian!, check_static_equilibrium_output_multipliers

include("structs.jl")

include("cable_joint.jl")
include("cluster_joint.jl")

include("connectivity.jl")

include("structure.jl")

include("inverse_statics.jl")

include("stiffness.jl")

include("analysis.jl")

include("plot.jl")

end
