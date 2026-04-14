module NCF
using LinearAlgebra
using StaticArrays
using SparseArrays
using LuxurySparse
using SymmetricFormats
using BlockDiagonals
using Rotations
using EponymTuples
using ForwardDiff
using DocStringExtensions

# Import utilities from parent module
import ..Rible: HouseholderOrthogonalization, skew, Axes, GECP, CartesianFrame, rotation_matrix

# Import abstract types for inheritance
import ..Rible: AbstractNonminimalCoordinates

# Import interface methods to extend them
import ..Rible: get_num_of_dims, get_num_of_coords, get_num_of_dof, get_num_of_cstr, get_num_of_intrinsic_cstr, get_num_of_local_dims
import ..Rible: to_local_coords, to_position, to_position_jacobian, to_velocity_jacobian
import ..Rible: find_rotation, find_angular_velocity, find_local_angular_velocity
import ..Rible: cstr_function, cstr_function!, cstr_jacobian, cstr_jacobian!, cstr_hessians
import ..Rible: cstr_forces_jacobian, add_cstr_forces_jacobian!, cstr_velocity_jacobian!
import ..Rible: cartesian_frame2coords, get_deform
import ..Rible: build_joint_cache, get_joint_violations!, get_joint_jacobian!
import ..Rible: add_joint_forces_jacobian!, get_joint_velocity_jacobian!
import ..Rible: kinetic_energy_coords, find_independent_free_idx, nullspace_mat
import ..Rible: make_cstr_function, make_cstr_jacobian
import ..Rible: has_constant_mass_matrix, get_cstr_idx

# Include type definitions and interface implementations
include("./types.jl")
include("./interface_implementations.jl")
include("./constructors.jl")

# Include specific implementations
include("./constraints.jl")
include("./mass_matrix.jl")
include("./functions.jl")
include("./joints.jl")

# Export the main types and constructors
export NC, LNCData
export NC2D, NC3D
export NC2D2C, NC2D4C, NC2D6C, NC3D3C, NC3D6C, NC3D12C
export NC2D1P, NC2D1P1V, NC2D2P, NC1P2V, NC2P1V, NC3P
export NC3D1P, NC3D1P1V, NC3D2P, NC1P3V, NC2P2V, NC3P1V, NC4P

# Note: Interface methods are exported from the main Rible module
# This allows for proper method extension without conflicts

end


