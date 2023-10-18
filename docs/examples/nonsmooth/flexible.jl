using LinearAlgebra, Statistics
using StaticArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using StructArrays
using BenchmarkTools
using TypeSortedCollections
using CoordinateTransformations
using Rotations
using ForwardDiff
using Interpolations
import GeometryBasics as GB
using OffsetArrays
using Printf
using Unitful
using Match
using EponymTuples
using Integrals
using Revise
import Meshes
using Meshing
import TensegrityRobots as TR
cd("examples/nonsmooth")
include("../analysis.jl"); includet("../analysis.jl")
include("../vis.jl"); includet("../vis.jl")
include("../dyn.jl"); includet("../dyn.jl")
figdir::String = raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ns"
