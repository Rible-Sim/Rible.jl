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
using Revise

using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

include("man_define.jl")
function eigenanalysis(n,α=π/6,k=3.e2)
    ndof = 6
    manipulator = man_ndof(ndof,k=k,c=0.0)

    Y = build_Y(manipulator)

    q0,_,_ = TR.get_initial(manipulator)

    λ0,Δu,_= TR.inverse(manipulator,deepcopy(manipulator),Y)
    TR.actuate!(manipulator,Δu)
    TR.reset_forces!(manipulator)
    TR.update_strings_apply_forces!(manipulator)
    @show transpose(TR.build_A(manipulator)(q0))*λ0 ≈ TR.build_Q̃(manipulator)*TR.fvector(manipulator)

    ωs = Vector{Vector{Float64}}()
    Zs = Vector{Matrix{Float64}}()

    ω0,Z0 = TR.undamped_eigen(manipulator,q0,λ0)

    push!(ωs,ω0)
    push!(Zs,Z0)

    for i = 1:n
        iθ = -i*α/n
        refman = man_ndof(ndof,θ=iθ,k=k,c=0.0) # reference
        refqi,_,_ = TR.get_initial(refman)
        refλi,Δu,a= TR.inverse(manipulator,refman,Y)

        actmani = deepcopy(manipulator)
        TR.actuate!(actmani,a)
        TR.reset_forces!(actmani)
        TR.distribute_q_to_rbs!(actmani,refqi)
        TR.update_strings_apply_forces!(actmani)

        @show transpose(TR.build_A(actmani)(refqi))*refλi ≈ TR.build_Q̃(actmani)*TR.fvector(actmani)

        ωi,Zi = TR.undamped_eigen(actmani,refqi,refλi)
        push!(ωs,ωi)
        push!(Zs,Zi)
    end
    ωs,Zs
end

n = 5
ωs,Zs = eigenanalysis(n,π/4)


for i = 1:n+1
    plt.plot(ωs[i],label="$i,k=300")
end
plt.legend()


ωs6,Zs6 = eigenanalysis(n,π/4,6.e2)

plt.plot(ωs6)
plt.legend()





function freevibra(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    function F!(F,q,q̇,t)

        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        # F .= Q̃*TR.fvector(tgstruct)
        # TR.apply_gravity!(tgstruct,factor=0.01)
        # F .= G
        TR.assemble_forces!(F,tgstruct)
        # @show isapprox(F,Q̃*TR.fvector(tgstruct))
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = freevibra(manipulator,q0)

prob = TS.DyProblem(freevibra(manipulator,q0),q0,q̇0,λ,(0.0,50.0))


function perturbation!(q0,q̇0)
    q̇0[end] = 0.0001
end
perturbation!(q0,q̇0)

sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-12,verbose=true)
δq = [q-q0 for q in sol.qs]

plt.plot(transpose(hcat(δq...)))

q = TR.undamped_modal_solve!(manipulator,q0,q̇0,λ,50.0,dt)
plt.plot(transpose(q))
