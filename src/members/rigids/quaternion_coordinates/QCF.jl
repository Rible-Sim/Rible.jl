module QCF
using LinearAlgebra
using Rotations
using Quaternions
using StaticArrays
using ForwardDiff
using DocStringExtensions
using EponymTuples
using PreallocationTools

import ..Rible: HouseholderOrthogonalization, skew
import ..Rible: get_num_of_cstr, get_num_of_coords
import ..Rible: get_num_of_dof, get_num_of_local_dims
import ..Rible: to_local_coords, to_transformation

import ..Rible: make_cstr_function
import ..Rible: make_cstr_jacobian
import ..Rible: make_cstr_hessians
import ..Rible: make_cstr_forces_jacobian
import ..Rible: cartesian_frame2coords
import ..Rible: find_rotation
import ..Rible: build_joint_cache, get_joint_violations!, get_joint_jacobian!
import ..Rible: find_independent_idx

include("QC.jl")

include("utils.jl")

include("constraints.jl")

include("mass_matrix.jl")

include("functions.jl")
 
include("joints.jl")

end