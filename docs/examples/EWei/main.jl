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

bot1 = BuildTail(1; β=.9)
ω², δq̌ = TR.undamped_eigen(bot1.tg)
@show ω = [sqrt(ωi) for ωi in ω²[1:end]]

function dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,s,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.distribute_s̄!(tg,s)
        TR.update_cables!(tg)
        ## TR.apply_gravity!(tg)
        TR.generate_forces!(tg)
        TR.get_force!(F,tg)
        ## F .= 0
    end
    function apply_acu!(tg, t; dt=1e-2)
        function inner(t;dt=dt)
            a = 10.0
            if 0<t<1
                return -a*t*dt
            elseif 1<=t<5
                return -a*dt
            elseif 5<=t<6
                return a*t*dt - 6a*dt
            else
                return 0
            end
        end
        tg.clustercables[1].segs[1].state.restlen += inner(t)
    end
    # apply_acu! = nothing
    Jac_F! = true
    @eponymtuple(F!, Jac_F!, apply_acu!)
end
prob = TR.SimProblem(bot1,dynfuncs)
TR.solve!(prob,TR.FBZhong06();dt=0.01,tspan=(0.0,5.0),ftol=1e-7,verbose=true)