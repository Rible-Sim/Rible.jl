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
# includet("plotting.jl")
# includet("../analysis.jl")
# includet("../plot_helpers.jl")
# includet("dyfun.jl")
bot1 = BuildTail()


function dynfuncs(tr)
    @unpack tg = tr
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)

    Q̃=TR.build_Q̃(tg)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.update_clustercables_apply_forces!(tg)
        # F .= Q̃*TR.fvector(tgstruct)
        TR.assemble_forces!(F,tg)
    end

    M,Φ,A,F!,nothing
end


plot_traj!(bot)