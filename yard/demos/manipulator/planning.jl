using LinearAlgebra
using Parameters
using StaticArrays
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
using Cthulhu
using NLsolve
using RecursiveArrayTools
using DynamicPolynomials
using Distributions
using BenchmarkTools
using NearestNeighbors
using LightGraphs, MetaGraphs
using GraphPlot, Compose, Cairo
using HomotopyContinuation
using Revise
import Rible as RB

cd("examples/manipulator")
includet("man_define.jl")
includet("man_plotting.jl")
includet("man_compare.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
includet("planning_dev.jl")

num_of_dof = 6
k = 1.e3
c = 0.0
man = man_ndof(num_of_dof;θ=0.0,k,c)
function initialize_RRT(tr)
    @unpack st = tr
    q0,q̇0 = RB.get_coords(st)
    chart0 = RB.Chart(st,q0,q̇0)
    x0 = chart0.xc

    states = ElasticArray{eltype(x0)}(undef, length(x0), 0)
    append!(states,x0)
    atlas = [chart0]

    rrt = MetaDiGraph(path_digraph(1))
    set_prop!(rrt, :name, "RRT")
    set_prop!(rrt, 1, :state, states[1])
    set_prop!(rrt, 1, :chart, 1)

    rrt, atlas, states
end
rrt, atlas, states = initialize_RRT(man)
x_rand = sample_atlas(atlas)
function new_chart_condition(chart)
    @unpack xc, Uc, ε, cosα, ρ = chart
    xₖ
    xₖ₊₁
    yₖ
    xₖ₊₁
    c1 = norm(xₖ₊₁-(xc+Uc*yₖ₊₁)) - ε
    c2 = cosα - norm(yₖ₊₁-yₖ)/norm(xₖ₊₁-xₖ)
    c3 = norm(yₖ₊₁) - ρ
end


function dynamics(tr,action)
    @unpack st = tr
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)
    Q̃ = RB.build_Q̃(st)
    function F!(F,q,q̇,t)
        RB.actuate!(tr,action(t))
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        RB.apply_gravity!(st;factor=0.1)
        RB.assemble_forces!(F,st)
    end
    Jac_Γ = RB.build_Jac_Γ(st)
    function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
        ∂Γ∂q,∂Γ∂q̇ = Jac_Γ(q,q̇)
        Q̃*∂Γ∂q,Q̃*∂Γ∂q̇
    end
    M,Φ,A,F!,Jac_F!
end

function simulate_action(tr,altas,mg,x_k,x_g,a)
    ix = find_nearest_state(states,x0)
    x_near = states[ix]
    ic = get_prop(mg, ix, :chart)
    chart = atlas[ic]
    t = 0.0
    tm = 1.0
    δ = 0.1
    while norm(x_k-x_g) > δ && abs(t) <= abs(tm)
        nextstate(tr,x_k,a)
    end
end
function nextstate(tr,x_k,a,tm;dt=0.01,ftol=1e-14)
    @unpack st = tr
    @unpack ncoords = st
    q_k = x_k[1:ncoords]
    q̇_k = x_k[ncoords+1:2ncoords]
    RB.set_new_initial!(tr,q_k,q̇_k)
    prob = RB.DyProblem(dynamics(tr,a),tr,(0.0,tm))
    intor = RB.solve!(tr,prob,RB.Zhong06();dt,ftol,exception=false,progress=false)
    q_next = tr.traj.qs[end]
    q̇_next = tr.traj.q̇s[end]
    x_next = vcat(q_next,q̇_next)
    x_next, intor.convergence
end
actions = [sample_action(man,0.001) for i = 1:100]
function simulate_action(tr,x_k,x_rand,actions,tm)
    divergence_count = 0
    dis = Inf
    i_next = 0
    x_next = copy(x_k)
    for (i,action) in enumerate(actions)
        x_kplus1, convergence = nextstate(tr,x_k,action,tm)
        if !convergence
            divergence_count += 1
            continue
        end
        new_dis = norm(x_kplus1-x_rand)
        if new_dis < dis
            dis = new_dis
            i_next = i
            x_next .= x_kplus1
        end
    end
    divergence_count,dis
end
simulate_action(man,atlas[1].xc,x_rand,actions,1.0)
norm(atlas[1].xc-x_rand)
@time simulate_action(man,x0,1.0)
@descend_code_warntype nextstate(man,x0,sample_action(man),1.0)
