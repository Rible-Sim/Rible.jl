module RibleQCF

using LinearAlgebra
using Rotations
using Quaternions
using StaticArrays
using ForwardDiff
using DocStringExtensions
using EponymTuples
using PreallocationTools
using SparseArrays

import Rible: HouseholderOrthogonalization, skew
import Rible: get_num_of_cstr, get_num_of_coords
import Rible: get_num_of_dof, get_num_of_local_dims
import Rible: get_num_of_intrinsic_cstr
import Rible: to_local_coords, to_position, to_position_jacobian, to_velocity_jacobian

import Rible: get_deform
import Rible: nullspace_mat
import Rible: cstr_function!
import Rible: cstr_jacobian!
import Rible: cstr_forces_jacobian!
import Rible: add_cstr_forces_jacobian!
import Rible: cartesian_frame2coords
import Rible: find_rotation, find_local_angular_velocity
import Rible: build_joint_cache, get_joint_violations!
import Rible: get_joint_jacobian!
import Rible: add_joint_forces_jacobian!
import Rible: get_joint_velocity_jacobian!
import Rible: find_independent_free_idx

import Rible: RigidBodyProperty, RigidBodyState, AbstractBodyCache, RigidBodyCache, AbstractCoordinates
import Rible: RigidBodyCache, update_inertia_cache!, lazy_update_state!, update_state!,
    stretch_loci!, update_transformations!, update_loci_states!, get_cstr_idx,
    has_constant_mass_matrix, AbstractCoordinatesState, InertiaCache
import Rible: Axes

include("quaternion_coordinates/QC.jl")

include("quaternion_coordinates/utils.jl")

include("quaternion_coordinates/constraints.jl")

include("quaternion_coordinates/mass_matrix.jl")

include("quaternion_coordinates/functions.jl")
 
include("quaternion_coordinates/joints.jl")

include("quaternion_coordinates.jl")

end
