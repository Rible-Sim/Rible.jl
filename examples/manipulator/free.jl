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
using Makie
AbstractPlotting.__init__()
using CSV
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")
include("../analysis.jl")

k = 2.e3; c = 0.0
"""
    function simulate_freefall(ndof;k,c,unit="mks")
        manipulator = man_ndof(ndof;k,c,unit)
        q0,q̇0,λ0 = TR.get_initial(manipulator)

        dt = 0.01 # Same dt used for PID AND Dynamics solver

        function freefall(tgstruct,q0)

            M = TR.build_massmatrix(tgstruct)
            #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
            Φ = TR.build_Φ(tgstruct,q0)
            A = TR.build_A(tgstruct)
            G = TR.build_G(tgstruct)
            function F!(F,q,q̇,t)
                # TR.reset_forces!(tgstruct)
                # TR.distribute_q_to_rbs!(tgstruct,q,q̇)
                # TR.update_strings_apply_forces!(tgstruct)
                # TR.apply_gravity!(tgstruct,factor=0.1)
                # TR.assemble_forces!(F,tgstruct)
                F .= G
            end

            M,Φ,A,F!,nothing

            #A,Φ,∂T∂q̇!,F!,M!,nothing
        end
        M,Φ,A,F!,Jacs = freefall(manipulator,q0)

        prob = TS.DyProblem(freefall(manipulator,q0),q0,q̇0,λ0,(0.0,100.0))

        sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-9,verbose=true)

        manipulator, sol
    end
    man_freefall, sol_freefall = simulate_freefall(2;k,c)
    manmks = man_ndof(2;k,c,unit="mks")
    mancgs = man_ndof(2;k,c,unit="cgs")
    Mmks = TR.build_massmatrix(manmks)
    Mcgs = TR.build_massmatrix(mancgs)
    Gmks = TR.build_G(manmks,factor=0.001)
    Gcgs = TR.build_G(mancgs,factor=100.0)

    kes,epes,gpes,es,es_err = analyse_energy(man_freefall,sol_freefall,gravity=true,factor=1.0)
    fig,ax = plt.subplots(1,1,num = 1,figsize=(4,3))
    # ax.plot(sol_freefall.ts,es)
    ax.plot(sol_freefall.ts,es_err)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rel. Energy Error")
    ax.set_ylim(-2e-4,2e-4)
    ax.set_xlim(sol_freefall.ts[1],sol_freefall.ts[end])
    ax.ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    ax.grid(true)
    fig.savefig("man_freefall_energy.png",dpi=300,bbox_inches="tight")

    bscene = plotstructure(man_freefall,sol_freefall,sliderplot)
"""
function simulate_freevibra(ndof;k,c,unit="mks",dt=0.01)
    manipulator = man_ndof(ndof;k,c,unit)
    q0,q̇0,λ0 = TR.get_initial(manipulator)

    # dt = 0.01 # Same dt used for PID AND Dynamics solver

    function freevibra(tgstruct,q0)

        M = TR.build_massmatrix(tgstruct)
        Φ = TR.build_Φ(tgstruct,q0)
        A = TR.build_A(tgstruct)

        function F!(F,q,q̇,t)
            TR.reset_forces!(tgstruct)
            TR.distribute_q_to_rbs!(tgstruct,q,q̇)
            TR.update_strings_apply_forces!(tgstruct)
            TR.apply_gravity!(tgstruct)
            TR.assemble_forces!(F,tgstruct,factor=1)
            # @show isapprox(F,Q̃*TR.fvector(tgstruct))
        end

        M,Φ,A,F!,nothing
    end

    M,Φ,A,F!,Jacs = freevibra(manipulator,q0)

    prob = TS.DyProblem(freevibra(manipulator,q0),q0,q̇0,λ0,(0.0,10.0))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14)

    manipulator, sol
end

# man_freevibra, sol_freevibra = simulate_freevibra(2;k,c,dt=0.01)
man_freevibra001, sol_freevibra001 = simulate_freevibra(2;k,c,dt=0.001)
# man_freevibra0001, sol_freevibra0001 = simulate_freevibra(2;k,c,dt=0.0001)

plotstructure(man_freevibra001,sol_freevibra001,sliderplot)
adams = TR.AdamsResults("triangle_man2e3_energy.xml")
# adams001 = TR.AdamsResults("triangle_man2e3_dt001.xml")
# adams001_9 = TR.AdamsResults("triangle_man2e3_dt001_err-9.xml")
t = adams("time")
M4Y = adams("M4Y")
# t001 = adams001("time")
# M4Y001 = adams001("M4Y")
# t001_9 = adams001_9("time")
# M4Y001_9 = adams001_9("M4Y")
p3x = [q[5] for q in sol_freevibra.qs]
p3y = [q[6] for q in sol_freevibra.qs]
p4x = [q[7] for q in sol_freevibra.qs]
plt.plot(sol_freevibra.ts,p3x)
plt.plot(sol_freevibra.ts,p3y)
plt.plot(sol_freevibra.ts,p4x)
# p4y01 = [q[8] for q in sol_freevibra.qs]
# plt.plot(sol_freevibra.ts,p4y01)
p4y001 = [q[8] for q in sol_freevibra001.qs]
plt.plot(sol_freevibra001.ts,p4y001)
# p4y0001 = [q[8] for q in sol_freevibra0001.qs]
# plt.plot(sol_freevibra0001.ts,p4y0001)
# plt.plot(t001_9,M4Y001_9)
# plt.plot(t001,M4Y001)
plt.plot(t,M4Y)
plt.legend()

MB1G = adams("MB1G")
MA2G = adams("MA2G")
MB1E = adams("MB1E")
MA2E = adams("MA2E")
SA1E = adams("SA1E")
SA2E = adams("SA2E")
SB1E = adams("SB1E")
SB2E = adams("SB2E")
cable_energy = [a+b+c+d for (a,b,c,d) in zip(SA1E,SA2E,SB1E,SB2E)]
body_energy = [a+b+c+d for (a,b,c,d) in zip(MB1G,MA2G,MB1E,MA2E)]
adams_es = body_energy .+ cable_energy
MB1E
ps = string_potential(man_freevibra001,sol_freevibra001)

kes,epes,gpes,es,es_err = analyse_energy(man_freevibra001,sol_freevibra001,gravity=true,elasticity=true,factor=1)
fig,ax = plt.subplots(3,1,num = 1,figsize=(4,3))
ax[1].plot(sol_freevibra001.ts,kes)
ax[1].plot(t,MB1E+MA2E)
ax.plot(sol_freevibra001.ts,sum(ps))
ax.plot(t,cable_energy)

SA1E = adams("SA1E")
SA2E = adams("SA2E")
SB1E = adams("SB1E")
SB2E = adams("SB2E")

ax.plot(sol_freevibra001.ts,ps[1])
ax.plot(t,SA1E)
ax.plot(sol_freevibra001.ts,es.-es[1])
ax.plot(t,adams_es.-adams_es[1])

ax.plot(sol_freevibra001.ts,es_err)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Rel. Energy Error")
ax.set_ylim(-8e-4,8e-4)


ax.set_xlim(sol_freevibra.ts[1],sol_freevibra.ts[end])
ax.grid(true)
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
