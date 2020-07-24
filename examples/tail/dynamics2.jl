using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
# using NLsolve
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
include("tail_define.jl")

n = 1
tail = make_tail(n)
# @code_warntype make_tail(n)
q0,q̇0,λ0 = TR.get_initial(tail)
q̇0[end-3:end] .= 5*[0.1,0.0,0.1,0.0]

# TR.get_nbodyconstraint(tail)
# TR.get_nbodydof(tail)
# TR.get_nbodycoords(tail)


function dynfuncs(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    #Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,A,F!,Jacs = dynfuncs(tail,q0)
# Φ(q0)
# @code_warntype Φ(q0)
# A(q0)
# @code_warntype A(q0)
F = similar(q0)
# @code_warntype F!(F,q0,q̇0,0.0)
dt = 0.01
prob = TS.DyProblem(dynfuncs(tail,q0),q0,q̇0,λ0,(0.0,20.0))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
# sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
@code_warntype TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)

M̂,Ĉ,K̂ = TR.linearize(tail,sol.qs[end],sol.q̇s[end],sol.λs[end])
M̄,C̄,K̄ = TR.frequencyshift(M̂,Ĉ,K̂,0.1)
M̃,K̃ = TR.enlarge(M̄,C̄,K̄)
eigenvalues = eigen(K̃,M̃).values
TR.∂Aᵀλ∂q(tail,sol.λs[end])

@code_warntype A(q0)

TR.∂Aᵀλ∂q(tail,sol.λs[end])
function ∂Aᵀλ∂q_finite(A,q,λ)
    function Aᵀλ(q)
        transpose(A(q))*λ
    end
    FiniteDiff.finite_difference_jacobian(Aᵀλ,q,relstep=1e-10)
end
fd = ∂Aᵀλ∂q_finite(A,sol.qs[end],sol.λs[end])
ad = TR.∂Aᵀλ∂q(tail,sol.λs[end])
tail.rigidbodies[1].state.cache.cfuncs.∂Aᵀλ∂q(4)

fd - ad
d = my∂𝐟∂q - ∂𝐟∂q
ḋ = my∂𝐟∂q̇ - ∂𝐟∂q̇
my∂𝐟∂q̇ = TR.build_tangent(tail,q0,q̇0)

my∂𝐟∂q = TR.build_tangent(tail,q0,q̇0)

my∂l̂∂q = TR.build_tangent(tail,q0,q̇0)

my∂f∂q = TR.build_tangent(tail,q0,q̇0)
my∂f∂q̇ = TR.build_tangent(tail,q0,q̇0)
@code_warntype TR.build_tangent(tail,q0,q̇0)
𝐟func = TR.make_testtangent(tail)
𝐟func(vcat(q0,q̇0))
using FiniteDiff
∂l̂∂q∂l̂∂q̇ = FiniteDiff.finite_difference_jacobian(𝐟func,vcat(q0,q̇0),relstep=1e-10)
∂l̂∂q = ∂l̂∂q∂l̂∂q̇[:,1:tail.ncoords]
∂l̂∂q̇ = ∂l̂∂q∂l̂∂q̇[:,tail.ncoords+1:end]

∂𝐟∂q∂𝐟∂q̇ = FiniteDiff.finite_difference_jacobian(𝐟func,vcat(q0,q̇0),relstep=1e-11)
∂𝐟∂q = ∂𝐟∂q∂𝐟∂q̇[:,1:tail.ncoords]
∂𝐟∂q̇ = ∂𝐟∂q∂𝐟∂q̇[:,tail.ncoords+1:end]

∂f∂q∂f∂q̇ = FiniteDiff.finite_difference_jacobian(𝐟func,vcat(q0,q̇0),relstep=1e-11)
∂f∂q = ∂f∂q∂f∂q̇[:,1:tail.ncoords]
∂f∂q̇ = ∂f∂q∂f∂q̇[:,tail.ncoords+1:end]




∂l∂q∂l∂q̇ = FiniteDiff.finite_difference_jacobian(𝐟func,vcat(q0,q̇0),relstep=1e-13)
∂l∂q = ∂l∂q∂l∂q̇[:,1:tail.ncoords]
∂l∂q̇ = ∂l∂q∂l∂q̇[:,tail.ncoords+1:end]

using ForwardDiff

function make_ATλ(λ)
    bps,q = TR.NaturalCoordinates.BP1P3V(zeros(3))
    cf = TR.NaturalCoordinates.CoordinateFunctions(bps)
    Φq = cf.Φq
    function inner_ATλ(q)
        transpose(Φq(q))*λ
    end
end

ATλ = make_ATλ([1,2,3,4,5,6])
ForwardDiff.jacobian(ATλ,rand(12))
