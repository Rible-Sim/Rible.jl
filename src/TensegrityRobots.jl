module TensegrityRobots

using Logging
using LinearAlgebra
using SparseArrays
using PreallocationTools
using Parameters
using StaticArrays
using StructArrays
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
using TypeSortedCollections
using EzXML
using ProgressMeter, Printf
using EponymTuples
using DocStringExtensions

const ASSETS_DIR = joinpath(@__DIR__, "..", "assets")
assetpath(files...) = normpath(joinpath(ASSETS_DIR, files...))

rotation_matrix(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

include("coordinates/QBF.jl")
using .QBF
include("coordinates/NCF.jl")
using .NCF
include("coordinates/ANCF.jl")
import .ANCF

include("members/rigidbody.jl")

include("members/flexiblebody.jl")

include("members/cable.jl")

include("members/joints.jl")

include("robots/tensegrity.jl")

include("robots/clustertensegrity.jl")

include("control/control.jl")

include("mechanics/inverse_statics.jl")

include("mechanics/forward_statics.jl")

include("mechanics/dynamic_relax.jl")

include("mechanics/contact.jl")

include("mechanics/linearization.jl")

include("mechanics/stiffness.jl")

include("mechanics/solvers.jl")

include("misc/miscellaneous.jl")

include("misc/URDF.jl")
import .URDF

end
