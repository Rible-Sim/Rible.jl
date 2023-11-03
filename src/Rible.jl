module Rible

using Logging
using LinearAlgebra
using SparseArrays
using PreallocationTools
using Parameters
using StaticArrays
using StructArrays
using RecursiveArrayTools
using SparseMatricesCSR
using NLsolve
import JuMP, COSMO, Clarabel
using BlockArrays
using BlockDiagonals
using SymmetricFormats
using Polyhedra
import CDDLib 
lib = CDDLib.Library()
using Tullio
using Quaternions
using ForwardDiff
using FiniteDiff
using RowEchelon
using HomotopyContinuation
using DynamicPolynomials
using Graphs
using MetaGraphsNext
using TypeSortedCollections
using EzXML
using ProgressMeter, Printf
using EponymTuples
using DocStringExtensions
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

include("members/rigids/QBF.jl")
using .QBF

include("members/abstract_body.jl")

include("members/rigids/rigid_body.jl")

include("members/flexibles/cables.jl")

include("members/flexibles/finite_elements/ANCF.jl")
import .ANCF

include("members/flexibles/flexible_body.jl")

include("members/joints.jl")

include("structures/connectivity.jl")

include("structures/structure.jl")

include("structures/mutate.jl")

include("robots/robot.jl")

include("control/control.jl")

include("./get.jl")

include("mechanics/Constraints.jl")

include("mechanics/inverse_statics.jl")

include("mechanics/forward_statics.jl")

include("mechanics/dynamic_relax.jl")

include("mechanics/contact.jl")

include("mechanics/linearization.jl")

include("mechanics/stiffness.jl")

include("mechanics/solvers.jl")

include("misc/miscellaneous.jl")

include("postprocess/analysis.jl")

include("preprocess/URDF.jl")
import .URDF

include("precompile.jl")

end