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

import ..Rible: HouseholderOrthogonalization, skew, Axes

export get_num_of_cstr, get_num_of_coords
export get_num_of_dof, get_num_of_local_dims
export to_local_coords, to_transformation
export find_rotation, find_angular_velocity

export make_cstr_function
export make_cstr_jacobian
export make_cstr_hessians
export make_cstr_forces_jacobian
export cartesian_frame2coords
export build_joint_cache, get_joint_violations!, get_joint_jacobian!

include("./LNC.jl")

include("./constraints.jl")

include("./mass_matrix.jl")

include("./functions.jl")

include("./joints.jl")
end


