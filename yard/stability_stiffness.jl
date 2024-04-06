#note -- preamble
using Statistics
# arrays
using LinearAlgebra
using SparseArrays
using StaticArrays
using TypeSortedCollections
using ElasticArrays
using OffsetArrays
using RecursiveArrayTools
using CircularArrays
const CA = CircularArray
using TypedTables
using Interpolations
using StructArrays
using SymmetricFormats
# data
using DataStructures
# using CubicSplines
# import FLOWMath
using CoordinateTransformations
using Rotations
using Unitful
using EponymTuples
# visualize/plot
using Makie
import GLMakie as GM
import CairoMakie as CM
import WGLMakie as WM
GM.activate!()
Makie.inline!(false)
import GeometryBasics as GB
import Meshes
using Meshing
using Match
# IO
using FileIO, MeshIO
using JLD2
using CSV, Tables
using EzXML
# print
using LaTeXStrings
using Latexify
using TexTables
using PrettyTables
using Printf
auto_display(false)
# solver
import DifferentialEquations as DE
using NLsolve
using Arpack
using COSMO
import Clarabel
using Polyhedra
import CDDLib 
lib = CDDLib.Library()
# code
using BenchmarkTools
using Cthulhu