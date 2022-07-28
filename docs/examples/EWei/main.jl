using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
using HomotopyContinuation
using Printf
using GLMakie
using Unitful
using NLsolve
using Revise
using StructArrays
using EponymTuples
using LaTeXStrings
using StaticArrays
using TypeSortedCollections
using Match
import TensegrityRobots as TR
includet("define.jl")
includet("../vis.jl")

bot1 = BuildTail(4; β=.9)
ω², δq̌ = TR.undamped_eigen(bot1.tg)
@show ω = [sqrt(ωi) for ωi in ω²[1:end]]
