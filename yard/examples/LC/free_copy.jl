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
using Rible; const TR = Rible
cd("examples/LC")
includet("mydef.jl")
includet("man_plotting.jl")
includet("../analysis.jl")


k = 2.e3; c = 0

manipulator = man_nd1(k,c)

plotstructure(manipulator)

function simulate_freevibra(num_of_dof;k,c,unit="mks",dt=0.01)
    manipulator = man_nd1(k,c)
    q0,q̇0,λ0 = RB.get_initial(manipulator.st)

    # q̇0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1]
    @show q̇0
    # dt = 0.01 # Same dt used for PID AND Dynamics solver

    function freevibra(tr)
        @unpack st = tr
        M = RB.build_massmatrix(st)
        Φ = RB.build_Φ(st)
        A = RB.build_A(st)

        function F!(F,q,q̇,t)
            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st,q,q̇)
            RB.update_cables_apply_forces!(st)
            # RB.apply_gravity!(st)
            RB.assemble_forces!(F,st)
            # @show isapprox(F,Q̃*RB.fvector(st))
        end

        M,Φ,A,F!,nothing
    end

    prob = RB.DyProblem(freevibra(manipulator),q0,q̇0,λ0,(0.0,5.0))

    RB.solve!(manipulator,prob,RB.Zhong06(),dt=dt,ftol=1e-14)

    manipulator
end

# man_freevibra, sol_freevibra = simulate_freevibra(2;k,c,dt=0.01)
man_freevibra001 = simulate_freevibra(2;k,c,dt=0.01)

# man_freevibra0001, sol_freevibra0001 = simulate_freevibra(2;k,c,dt=0.0001)

plotstructure(manipulator)


manipulator = man_nd1(k,c;ratio=0.85)
q0,_ = RB.get_coords(manipulator.st)
A = RB.build_A(manipulator.st)
Q̃ = RB.build_Q̃(manipulator.st)
Γ = RB.build_Γ(manipulator.st)
Aq = A(q0)
λ = transpose(Aq)\(Q̃*Γ)
transpose(Aq)*λ-Q̃*Γ
M̂,Ĉ,K̂ = RB.linearize!(manipulator.st,q0,λ)
d = eigen(K̂,M̂)
