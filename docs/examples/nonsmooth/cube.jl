using LinearAlgebra
using StaticArrays
# using Parameters
using TypeSortedCollections
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using CoordinateTransformations
using BenchmarkTools
using Rotations
import GeometryBasics as GB
using Unitful, Match, Printf
import Meshes
# using MeshViz
# using Meshing
using EponymTuples
using OffsetArrays
using Revise
import TensegrityRobots as TR
# GLMakie.enable_SSAO[] = true
cd("examples/nonsmooth")
include("../analysis.jl"); includet("../analysis.jl")
include("../vis.jl"); includet("../vis.jl")
include("../dyn.jl"); includet("../dyn.jl")
