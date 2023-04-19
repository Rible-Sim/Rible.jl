using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
scene = plot(rand(10))
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
using CSV
using Revise
using TensegrityRobots; const TR = TensegrityRobots
cd("examples/LC")
includet("mydef.jl")
includet("man_plotting.jl")
includet("../analysis.jl")


k = 2.e3; c = 0

manipulator = man_nd1(k,c)

plotstructure(manipulator)

function simulate_freevibra(ndof;k,c,unit="mks",dt=0.01)
    manipulator = man_nd1(k,c)
    q0,q̇0,λ0 = TR.get_initial(manipulator.tg)

    # q̇0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1]
    @show q̇0
    # dt = 0.01 # Same dt used for PID AND Dynamics solver

    function freevibra(tr)
        @unpack tg = tr
        M = TR.build_massmatrix(tg)
        Φ = TR.build_Φ(tg)
        A = TR.build_A(tg)

        function F!(F,q,q̇,t)
            TR.reset_forces!(tg)
            TR.distribute_q_to_rbs!(tg,q,q̇)
            TR.update_cables_apply_forces!(tg)
            # TR.apply_gravity!(tg)
            TR.assemble_forces!(F,tg)
            # @show isapprox(F,Q̃*TR.fvector(tg))
        end

        M,Φ,A,F!,nothing
    end

    prob = TR.DyProblem(freevibra(manipulator),q0,q̇0,λ0,(0.0,5.0))

    TR.solve!(manipulator,prob,TR.Zhong06(),dt=dt,ftol=1e-14)

    manipulator
end

# man_freevibra, sol_freevibra = simulate_freevibra(2;k,c,dt=0.01)
man_freevibra001 = simulate_freevibra(2;k,c,dt=0.01)

# man_freevibra0001, sol_freevibra0001 = simulate_freevibra(2;k,c,dt=0.0001)

plotstructure(manipulator)


manipulator = man_nd1(k,c;ratio=0.85)
q0,_ = TR.get_q(manipulator.tg)
A = TR.build_A(manipulator.tg)
Q̃ = TR.build_Q̃(manipulator.tg)
Γ = TR.build_Γ(manipulator.tg)
Aq = A(q0)
λ = transpose(Aq)\(Q̃*Γ)
transpose(Aq)*λ-Q̃*Γ
M̂,Ĉ,K̂ = TR.linearize!(manipulator.tg,q0,λ)
d = eigen(K̂,M̂)
