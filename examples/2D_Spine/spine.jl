using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
using LaTeXStrings
using Cthulhu
# using NLsolve
using Convex
using MosekTools
using Revise
# using TS
# using Robot2D
# const TR = Robot2D
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot


function make_spine(n)
    a = 0.1
    b = 0.1*√2
    α = 3/4*π
    d = 0.15

    function vert(i,r,θ,a,b,α)
        if i == 1
            movable = false
        else
            movable = true
        end
        CoM = zeros(2)
        ap1 = [a,0.0]
        ap2 = b*[cos( α),sin( α)]
        ap3 = b*[cos(-α),sin(-α)]
        aps = [ap1,ap2,ap3]
        m = 0.495
        inertia = 1.0
        prop = TR.RigidBodyProperty(i,movable,m,inertia,CoM,aps)
        state = TR.RigidBodyState(prop,r,θ)
        rb = TR.RigidBody(prop,state)
    end

    rs = [[i*d,0.0] for i = 0:n-1]
    θs = zeros(n)
    rbs = [vert(i,rs[i],θs[i],a,b,α) for i = 1:n]

    nstrings = 4*(n-1)
    k = 840.0
    c = 1000.0
    original_restlen = 0.0
    ss = [TR.SString2D(original_restlen,k,c) for i = 1:nstrings]

    acs = [TR.Actuator(SVector(ss[i])) for i = 1:nstrings]
    body2q_raw = [collect(4(i-1)+1:4i) for i = 1:n]
    body2q = TR.filter_body2q(body2q_raw,rbs)
    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 2:n
        push!(string2ap,(TR.ID(i-1,2),TR.ID(i,2)))
        push!(string2ap,(TR.ID(i-1,1),TR.ID(i,2)))
        push!(string2ap,(TR.ID(i-1,1),TR.ID(i,3)))
        push!(string2ap,(TR.ID(i-1,3),TR.ID(i,3)))
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end

spine = make_spine(2)
W = TR.build_W(spine)

Q̃ = TR.build_Q̃(spine)
# G = TR.build_G(spine)
nλ = TR.get_nconstraint(spine)

K = TR.build_K(spine)
RHS = TR.build_RHS(spine)

K̂ = TR.build_K̂(spine)
ℓ = TR.build_ℓ(spine)
k = [s.k for s in spine.strings]
Diagℓ = Diagonal(ℓ)
Diagk = Diagonal(k)
ns = length(k)
# Exℓ = vcat(zeros(nλ),ℓ)
# Exk = vcat(zeros(nλ),k)
#
# x = Variable(ns+nλ)
# u = x[nλ+1:end]
# obj = quadform(u,Diagk)-2sum(Exℓ.*Exk.*x)
#
# problem = minimize(obj,
#             K*x==RHS,0.00<=u,u<=ℓ)

u = Variable(ns,Positive())
obj = quadform(u,Diagk)-2sum(ℓ.*k.*u)
obj = sumsquares(Diagk*ℓ-u)
problem = minimize(obj,
        (I-W)*Q̃*K̂*u==(I-W)*RHS,u<=0.99ℓ)

solve!(problem,Mosek.Optimizer)

problem.status
problem.optval
u.value
K*x.value - RHS
f = k.*(ℓ-x.value[nλ+1:end])
TR.set_restlen!(spine,x.value[nλ+1:end])
q0,q̇0,λ0 = TR.get_initial(spine)
function dynfuncs(tg,q0)
    M = TR.build_massmatrix(tg)
    A = TR.build_A(tg)
    Φ = TR.build_Φ(tg,q0)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_strings_apply_forces!(tg)
        TR.apply_gravity!(tg)
        TR.assemble_forces!(F,tg)
    end
    M,Φ,A,F!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(spine,q0)
Φ(q0)
size(A(q0))[2]
dt = 0.01
prob = TS.DyProblem(dynfuncs(spine,q0),q0,q̇0,λ0,(0.0,10.0))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-13,verbose=true)
sol = TS
