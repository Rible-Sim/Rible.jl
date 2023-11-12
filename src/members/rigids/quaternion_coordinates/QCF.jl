module QCF
using LinearAlgebra
using Rotations
using Quaternions
using StaticArrays
using ForwardDiff
using DocStringExtensions

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


export get_num_of_cstr, get_num_of_coords
export get_num_of_dof, get_num_of_local_dims
export to_local_coords, to_transformation

export make_cstr_function
export make_cstr_jacobian
export make_cstr_hessians
export make_cstr_forces_jacobian
export cartesian_frame2coords
export find_rotation
export find_local_angular_velocity

include("QC.jl")

include("utils.jl")

include("constraints.jl")

include("mass_matrix.jl")

include("functions.jl")

end