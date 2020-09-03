using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
# using NLsolve
# using DifferentialEquations
# using Sundials
# using DASKR

using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")
include("../analysis.jl")

function simulate_freefall(ndof,k=1.e3,c=0.0)
    manipulator = man_ndof(ndof,k=k,c=c)
    q0,q̇0,λ0 = TR.get_initial(manipulator)

    dt = 0.01 # Same dt used for PID AND Dynamics solver

    function freefall(tgstruct,q0)

        M = TR.build_massmatrix(tgstruct)
        #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
        Φ = TR.build_Φ(tgstruct,q0)
        A = TR.build_A(tgstruct)

        Q̃ = TR.build_Q̃(tgstruct)
        G = TR.build_G(tgstruct)
        function F!(F,q,q̇,t)
            # TR.reset_forces!(tgstruct)
            # TR.distribute_q_to_rbs!(tgstruct,q,q̇)
            # TR.update_strings_apply_forces!(tgstruct)
            # F .= Q̃*TR.fvector(tgstruct)
            # TR.apply_gravity!(tgstruct)
            F .= G
            # TR.assemble_forces!(F,tgstruct)
            # @show isapprox(F,Q̃*TR.fvector(tgstruct))
        end

        M,Φ,A,F!,nothing

        #A,Φ,∂T∂q̇!,F!,M!,nothing
    end
    M,Φ,A,F!,Jacs = freefall(manipulator,q0)

    prob = TS.DyProblem(freefall(manipulator,q0),q0,q̇0,λ0,(0.0,100.0))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-13,verbose=true)
    
    manipulator, sol
end
man_freefall, sol_freefall = simulate_freefall(2)

_,_,_,_,es_err = analyse_energy(man_freefall,sol_freefall,gravity=true)
fig,ax = plt.subplots(1,1,num = 1,figsize=(8,6))
ax.plot(sol_freefall.ts,es_err)
ax.set_xlabel("Time(s)")
ax.set_ylabel("Rel. Energy Error")
ax.set_ylim(-10e-5,5e-5)
ax.set_xlim(sol_freefall.ts[1],sol_freefall.ts[end])
fig.savefig("man_freefall_energy.png",dpi=300,bbox_inches="tight")

bscene = plotstructure(man_freefall,sol_freefall,sliderplot)

function simulate_freevibra(ndof;k=3.e2,c=1.e2)
    manipulator = man_ndof(ndof,k=k,c=c)
    q0,q̇0,λ0 = TR.get_initial(manipulator)

    dt = 0.01 # Same dt used for PID AND Dynamics solver

    function freevibra(tgstruct,q0)

        M = TR.build_massmatrix(tgstruct)
        #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
        Φ = TR.build_Φ(tgstruct,q0)
        A = TR.build_A(tgstruct)
    
        Q̃ = TR.build_Q̃(tgstruct)
        G = TR.build_G(tgstruct)
        function F!(F,q,q̇,t)
            TR.reset_forces!(tgstruct)
            TR.distribute_q_to_rbs!(tgstruct,q,q̇)
            TR.update_strings_apply_forces!(tgstruct)
            F .= Q̃*TR.fvector(tgstruct)
            TR.apply_gravity!(tgstruct)
            # F .= G
            TR.assemble_forces!(F,tgstruct)
            # @show isapprox(F,Q̃*TR.fvector(tgstruct))
        end
    
        M,Φ,A,F!,nothing
    
        #A,Φ,∂T∂q̇!,F!,M!,nothing
    end
 
    M,Φ,A,F!,Jacs = freevibra(manipulator,q0)

    prob = TS.DyProblem(freevibra(manipulator,q0),q0,q̇0,λ0,(0.0,100.0))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-13,verbose=true)
    
    manipulator, sol
end

man_freevibra, sol_freevibra = simulate_freevibra(2,c=0.0)

_,_,_,_,es_err = analyse_energy(man_freevibra,sol_freevibra,gravity=true,elasticity=true)

fig,ax = plt.subplots(1,1,num = 1,figsize=(8,6))
ax.plot(sol_freevibra.ts,es_err)
ax.set_xlabel("Time(s)")
ax.set_ylabel("Rel. Energy Error")
ax.set_ylim(-5e-4,10e-4)
ax.set_xlim(sol_freevibra.ts[1],sol_freevibra.ts[end])
fig.savefig("man_freevibra_energy.png",dpi=300,bbox_inches="tight")


fig,ax = plt.subplots(1,1,num = 1,figsize=(8,6))
ax.plot(sol.ts,es)
ax.set_xlabel("Time(s)")
ax.set_ylabel("Energy")
ax.set_ylim(es[1]-0.1abs(es[1]),es[1]+0.1abs(es[1]))








q5 = [q[5] for q in sol.qs]
plt.plot(sol.ts,q5)

function DAE(M,Φ,A,F!,q0,λ0)
    nq = length(q0)
    nλ = length(λ0)
    nu = 2nq + nλ
    mr = norm(M,Inf)
    function f(resid,du,u,p,t)
        q = u[1:nq]
        q̇ = u[nq+1:2nq]
        λ = u[2nq+1:end]
        # q̇ = du[1:nq]
        q̈ = du[nq+1:2nq]
        F = similar(q)
        F!(F,q,q̇,t)
        fr = norm(F,Inf)
        scaling = mr + fr
        resid[1:nq] = du[1:nq] - q̇
        resid[nq+1:2nq] = M*q̈ + scaling*transpose(A(q))*λ - F
        resid[2nq+1:end] =  scaling*Φ(q)
    end
end
G = TR.build_G(manipulator)
q̈0 = inv(M)*G
u0 = vcat(q0,q̇0,zero(λ0))
du0 = vcat(q̇0,q̈0,zero(λ0))
dvars = vcat(ones(Bool,2length(q0)),zeros(Bool,length(λ0)))
dae = DAEProblem(DAE(M,Φ,A,F!,q0,λ0),du0,u0,(0.0,50.0);differential_vars=dvars)
daesol = solve(dae,IDA())
daesol[5,:]

plt.plot(daesol.t,daesol[5,:])


function freevibra(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        # F .= Q̃*TR.fvector(tgstruct)
        TR.apply_gravity!(tgstruct,factor=0.01)
        # F .= G
        TR.assemble_forces!(F,tgstruct)
        # @show isapprox(F,Q̃*TR.fvector(tgstruct))
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end

prob = TS.DyProblem(freevibra(manipulator,q0),q0,q̇0,λ0,(0.0,50.0))

sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-13,verbose=true)

kes = [TR.kinetic_energy_coords(manipulator,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]

gpes = 0.01*[TR.gravity_potential_energy(manipulator,q) for q in sol.qs]
epes = [TR.potential_energy(manipulator,q) for q in sol.qs]
es = kes + gpes + epes

es_err = (es.-es[1])./es

plt.plot(sol.ts,es_err)
plt.plot(sol.ts,es)
q5 = [q[5] for q in sol.qs]
plt.plot(sol.ts,q5)
