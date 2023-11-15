module Rible


# arrays
using LinearAlgebra
using SparseArrays
using PreallocationTools
using StaticArrays
using StructArrays
using RecursiveArrayTools
using SparseMatricesCSR
using BlockArrays
using BlockDiagonals
using SymmetricFormats
using RowEchelon
using TypeSortedCollections
using TypedTables
# other
using Quaternions
using Rotations
using CoordinateTransformations
using Tullio
using EponymTuples
using Unitful
# solvers
using NLsolve
import JuMP, COSMO, Clarabel
using Polyhedra
import CDDLib 
lib = CDDLib.Library()
# diff
using ForwardDiff
using FiniteDiff
# poly
using HomotopyContinuation
using DynamicPolynomials
# graph
using Graphs
using MetaGraphsNext
# IO
using EzXML
# print
using ProgressMeter
using Printf
using Logging
# visualize/plot
import MakieCore
using Observables
import GeometryBasics as GB
import ColorTypes as CT
const lift = map
# docs
using DocStringExtensions
# code
using PrecompileTools
using Test

const ASSETS_DIR = joinpath(@__DIR__, "..", "assets")
assetpath(files...) = normpath(joinpath(ASSETS_DIR, files...))


include("utils.jl")

include("members/loci.jl")

include("members/rigids/natural_coordinates/NCF.jl")
using .NCF

import .NCF: make_cstr_function
import .NCF: make_cstr_jacobian
import .NCF: make_cstr_hessians
import .NCF: make_cstr_forces_jacobian
import .NCF: get_num_of_cstr, get_num_of_coords, get_num_of_dof
import .NCF: get_num_of_dims, get_num_of_local_dims
import .NCF: cartesian_frame2coords
import .NCF: find_rotation, find_angular_velocity
import .NCF: build_joint_cache, get_joint_violations!, get_joint_jacobian!

include("members/rigids/quaternion_coordinates/QCF.jl")
using .QCF

include("members/abstract_body.jl")

include("members/rigids/rigid_body.jl")

include("members/flexibles/cables.jl")

include("members/flexibles/finite_elements/ANCF.jl")
import .ANCF

include("members/flexibles/flexible_body.jl")

include("structures/connectivity.jl")

include("structures/structure.jl")

include("members/joints.jl")

include("structures/mutate.jl")

include("robots/robot.jl")

include("control/control.jl")

include("./get.jl")

include("mechanics/inverse_statics.jl")

include("mechanics/forward_statics.jl")

include("mechanics/dynamic_relax.jl")

include("mechanics/contact.jl")

include("mechanics/linearization.jl")

include("mechanics/stiffness.jl")

include("mechanics/solvers.jl")

include("misc/miscellaneous.jl")

include("mechanics/dynamics.jl")

include("postprocess/analysis.jl")

include("preprocess/URDF.jl")
import .URDF

include("postprocess/vis.jl")

include("precompile.jl")

end