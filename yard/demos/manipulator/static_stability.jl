using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Observables
using Makie
import GLMakie
import CairoMakie
GLMakie.activate!()
# import PyPlot as plt
# plt.pygui(true)
using LaTeXStrings
using RecursiveArrayTools
using Cthulhu
using XLSX, DataFrames
using Unitful
# using NLsolve
# using DifferentialEquations
using ForwardDiff
using HomotopyContinuation
using Printf
using Revise
import Rible as RB
cd("examples/manipulator")
includet("man_define.jl")
includet("man_plotting.jl")
includet("../plot_helpers.jl")
set_theme!(font = "Nimbus Rom No9 L", fontsize = fontsize)


k = 292.14; c = 0.0;
man_inv = man_ndof(1;θ=0.0,k,c)
rb1prop = man_inv.st.rigidbodies[1].prop
rb2prop = man_inv.st.rigidbodies[2].prop
a1 = norm(rb1prop.aps[6] - rb1prop.aps[2])
a2 = norm(rb2prop.aps[5] - rb2prop.aps[2])
β1 = get_angle(rb1prop.aps[6] - rb1prop.aps[2],rb1prop.aps[7] - rb1prop.aps[2])
β2 = get_angle(rb2prop.aps[4] - rb2prop.aps[1],rb2prop.aps[5] - rb2prop.aps[1])
fig,ax,plt = scatter(rb2prop.aps, markersize=20)
text!(
    "p" .* string.(1:7),
    position = rb2prop.aps
)
ax.aspect = DataAspect()

θ0 = π-(β1/2+β2/2)

u = RB.get_original_restlen(man_inv)
L0 = u[1]

function theta2L1(θ,θ0,a1,a2)
    sqrt(a1^2 + a2^2 - 2a1*a2*cos(θ0 + θ))
end
function theta2L2(θ,θ0,a1,a2)
    sqrt(a1^2 + a2^2 - 2a1*a2*cos(θ0 - θ))
end

function make_f_static(θ0,a1,a2,k,L0,ΔL)
    function inner_f_static(θ)
        L1 = theta2L1(θ,θ0,a1,a2)
        L2 = theta2L2(θ,θ0,a1,a2)
        ret = k*a1*a2*((L1-L0+ΔL)*sin(θ0+θ)/L1 - (L2-L0-ΔL)*sin(θ0-θ)/L2)
    end
end
L1 = theta2L1(0.0,θ0,a1,a2)
L2 = theta2L2(0.0,θ0,a1,a2)

ΔL = 0.0
k = 800.0
f_static = make_f_static(θ0,a1,a2,k,L0,ΔL)
f_static(-1)


function grid_check(θ0,a1,a2,k,ΔL)
    L0_range = 0.0:0.001:0.2
    ζ = Vector{eltype(L0_range)}()
    for L0 in L0_range
        f_static = make_f_static(θ0,a1,a2,k,L0,ΔL)
        push!(ζ,ForwardDiff.derivative(f_static,0.0))
    end
    L0_range,ζ
end
L0_range,ζ = grid_check(θ0,a1,a2,k,ΔL)
plot(L0_range,ζ)
0.1355
0.174-3.9061/453.09

function free(bot)
    st = bot.st
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)
    Q̃ = RB.build_Q̃(st)
    function F!(F,q,q̇,t)
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        # RB.apply_gravity!(st)
        RB.assemble_forces!(F,st)
        # @show isapprox(F,Q̃*RB.fvector(st))
    end
    Jac_Γ = RB.build_Jac_Γ(st)
    function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
        ∂Γ∂q,∂Γ∂q̇ = Jac_Γ(q,q̇)
        Q̃*∂Γ∂q,Q̃*∂Γ∂q̇
    end
    M,Φ,A,F!,nothing
end

man_inv = man_ndof(1,[1.0,0.0];θ=deg2rad(1),k,c,restlen=0.1480,isvirtual=false)
dt = 0.01
prob = RB.DynamicsProblem(man_inv,free,(0.0,10.0))

RB.solve!(prob,RB.Zhong06();dt,ftol=1e-14)
plotstructure(man_inv)

A = [0.1 0; 0 0.1]
A^0
