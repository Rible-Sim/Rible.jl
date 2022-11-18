module TensegrityRobots

using Logging
using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using PreallocationTools
using Parameters
using StaticArrays
using StructArrays
using NLsolve
import JuMP
import COSMO
using BlockArrays
using BlockDiagonals
using Tullio
using Quaternions
using ForwardDiff
using FiniteDiff
using HomotopyContinuation
using DynamicPolynomials
using TypeSortedCollections
using EzXML
using ProgressMeter, Printf
using EponymTuples
using DocStringExtensions
rotation_matrix(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

include("quaternions.jl")
using .QuaternionCoordinates
include("naturalcoordinates.jl")
using .NaturalCoordinates

include("rigidbody.jl")

include("cable.jl")

include("tensegrity.jl")

include("inverse_statics.jl")

include("forward_statics.jl")

include("dynamic_relax.jl")

include("control.jl")


# include("plotting.jl")
include("contact.jl")
include("linearization.jl")
include("tangent.jl")
include("miscellaneous.jl")
include("solvers.jl")
include("planning.jl")
include("clustercables.jl")

end
