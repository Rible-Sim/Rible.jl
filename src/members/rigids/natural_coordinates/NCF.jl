module NCF
using LinearAlgebra
using StaticArrays
using SparseArrays
using LuxurySparse
using ForwardDiff
using DocStringExtensions
using SymmetricFormats
using EponymTuples
using BlockDiagonals

import ..Rible: HouseholderOrthogonalization, skew, Axes, GECP, CartesianFrame

import ..Rible: get_num_of_cstr, get_num_of_coords
import ..Rible: get_num_of_dof, get_num_of_local_dims
import ..Rible: to_local_coords, to_transformation
import ..Rible: find_rotation, find_angular_velocity
import ..Rible: find_local_angular_velocity
import ..Rible: get_deform
import ..Rible: cstr_function
import ..Rible: cstr_jacobian
import ..Rible: cstr_forces_jacobian
import ..Rible: cartesian_frame2coords
import ..Rible: build_joint_cache, get_joint_violations!
import ..Rible: get_joint_jacobian!, get_joint_forces_jacobian!
import ..Rible: find_independent_free_idx
import ..Rible: nullspace_mat

include("./LNC.jl")

include("./constraints.jl")

include("./mass_matrix.jl")

include("./functions.jl")

include("./joints.jl")
end


