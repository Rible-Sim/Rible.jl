module TensegrityRobots

using Logging
using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using Parameters
using StaticArrays
using StructArrays
using NLsolve
import JuMP
import COSMO
using BlockArrays
using BlockDiagonals
using GeometryTypes
using ForwardDiff
using FiniteDiff
using HomotopyContinuation
using DynamicPolynomials
using TypeSortedCollections
using EzXML
using Cubature
using GaussQuadrature
using ProgressMeter, Printf
using EponymTuples
rotation_matrix(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

include("naturalcoordinates.jl")
using .NaturalCoordinates

include("rigidbody.jl")

include("cable.jl")

include("tensegrity.jl")

include("inverse_statics.jl")

include("forward_statics.jl")

include("control.jl")


# Write your package code here.
include("plotting.jl")
include("contact.jl")
include("linearization.jl")
include("tangent.jl")
include("miscellaneous.jl")
include("solvers.jl")
include("planning.jl")

end
