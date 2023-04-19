using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
#using ForwardDiff
# using DynamicPolynomials
# using Printf
using Cthulhu
using RecursiveArrayTools
using XLSX
using CSV
using FileIO
using DataFrames
using EasyFit
using CoordinateTransformations
using HomotopyContinuation
using TypeSortedCollections
using Interpolations
using GeometryBasics
using FileIO
using NLsolve
using Optim
using Revise
import TensegrityRobots as TR
cd("examples/ming")
includet("links_define.jl")
includet("links_plotting.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
set_theme!(font = "Times New Roman")

function fit_spring(sheetname,rl,dmin,dmax,do_plot=true)
    spring_df = DataFrame(
        XLSX.readtable(
        "弹簧数据.xlsx", sheetname,infer_eltypes=true)...)
    deformations = spring_df[:,"位移/mm"]
    forces = spring_df[:,"力/N"]
    # plt.plot(deformations,forces,ls="",marker="o")
    # plt.xlim(0,deformations[end])
    # plt.ylim(0,forces[end])

    s_indices = findall((x)->dmin<x<dmax,deformations)
    ds = deformations[s_indices]
    fs = forces[s_indices]
    s_fit = fitlinear(ds,fs)
    r = s_fit.b
    k = s_fit.a
    plot_fit_fs = [0,fs[end]]
    plot_fit_ds = (plot_fit_fs .-r)./k
    effective_rl = rl - r/k
    if do_plot
        fig, ax = plt.subplots(1,1)
        ax.plot(deformations,forces,ls="",marker="o")
        ax.plot(ds,fs,ls="",marker="o")
        ax.plot(plot_fit_ds,plot_fit_fs)
        # ax.set_xlim(0,deformations[end])
        ax.set_ylim(0,forces[end])
        fig
    end
    effective_rl, k
end

er1, k1 = fit_spring("Sheet1",40,1,65)
er2, k2 = fit_spring("Sheet2",40,5,80)
er3, k3 = fit_spring("Sheet3",30,1,15)

pyramid = read_rb_mesh()

load1df = CSV.File("527load1.ts",header=1) |> DataFrame
pp = [1e-3.*load1df[!, ["X$i","Y$i","Z$i"]] for i = 1:8]

spine2 = spine_true(2,-0.1,RotY(0.0);o1=zeros(3),X1=Matrix(RotY(0.0)),dir=[1,0,0])
rs = [[pp[j][1,i] for i = 1:3] for j = 1:4]
r̄s = spine2.tg.rigidbodies[1].prop.aps[[13,17,18,19]]

function optim_npoints(r̄s,rs)
    function npoints(x)
        ro = x[1:3]
        α,β,γ = x[4:6]
        R = RotXYZ(α,β,γ)
        one(r̄,r) = ro + R*r̄ - r
        weights = ones(length(rs))
        # weights[1] =
        sum([w*sum(one(r̄,r).^2) for (r̄,r,w) in zip(r̄s,rs,weights)])
    end
    Optim.optimize(npoints,zeros(6))
end
optim_res = optim_npoints(r̄s,rs)
optim_minimizer = Optim.minimizer(optim_res)
ro = optim_minimizer[1:3]
R = RotXYZ(optim_minimizer[4:6]...)
rs_opt = [ro + R*r̄  for r̄ in r̄s]
rs

spine2 = spine_true(2,-0.105,RotY(0.0);o1=ro,X1=Matrix(R),dir=[1,0,0])

plotstructure(spine2,pyramid)

function dynfuncs(tgbot)
    @unpack tg = tgbot
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)
    Q̃ = TR.build_Q̃(tg)
    function F!(F,q,q̇,t)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_strings_apply_forces!(tg)
        # TR.apply_gravity!(tg)
        TR.assemble_forces!(F,tg)
    end
    M,Φ,A,F!,nothing
end
function get_static_λ(bot,dyfuncs)
    @unpack tg = bot
    M,Φ,A,F!,_ = dyfuncs(bot)
    q = bot.traj.qs[end]
    F = similar(q)
    F!(F,q,zero(q),0.0)
    λ = transpose(A(q))\F
    TR.check_static_equilibrium(tg,q,λ)
    # @show transpose(A(q))*λ-F
    λ
end
function get_forward_start(tgbot,funcs)
    @unpack tg = tgbot
    q = tgbot.traj.qs[end]
    TR.reset_forces!(tg)
    TR.distribute_q_to_rbs!(tg,q)
    TR.update_strings_apply_forces!(tg)
    ℓ = [s.state.length for s in tg.strings]
    s = 1 ./ℓ
    λ = get_static_λ(tgbot,funcs)
    q,s,λ
end
function find_initial_position(bot)
    @unpack tg = bot
    q0,s0,λ0 = get_forward_start(bot,dynfuncs)
    u0 = TR.get_strings_restlen(bot)
    g0 = [0.0]
    Psys,vars = TR.forward_system(tg,reshape(TR.build_G!(tg),:,1))
    @unpack q,s,λ,u,g = vars
    P = subs(Psys.expressions,u=>u0,g=>g0)
    function f(x)
        q1,s1,λ1 = TR.split_by_lengths(x,[length(q0),length(s0),length(λ0)])
        to_number.(subs(P,q=>q1,s=>s1,λ=>λ1))
    end
    nlres = nlsolve(f,vcat(q0,s0,λ0),ftol=1e-14)
    @assert converged(nlres)
    qzero,szero,λzero = TR.split_by_lengths(nlres.zero,[length(q0),length(s0),length(λ0)])
end

qzero,szero,λzero = find_initial_position(spine2)

function dynamic_relax(tgbot,dynfuncs;dt=0.001,tend=10.0,verbose=false)
    prob = TR.SimProblem(tgbot,dynfuncs,(0.0,tend))
    TR.solve!(prob,TR.Zhong06();dt,ftol=1e-14,verbose)
    tgbot,dynfuncs
end
_,dynfuncs = dynamic_relax(spine2,dynfuncs)
qend,_ = TR.get_q(spine2.tg)

TR.set_new_initial!(spine2,qzero)
function plot_super(pp,bot,pyramid;linked=false)
    @unpack tg,traj = bot
    fig = Figure(resolution=(1920,1080))
    ax1 = fig[1,1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    pyramids,strings,point_spheres = plotstructure!(ax1,tg,pyramid)
    ax2 = fig[1,2] = Axis(fig)
    colsize!(fig.layout, 1, Relative(2/3))
    scatter!(ax2,pp[5][:,3],markersize=4)
    rp13 = get_trajectory(bot,2,13)
    scatter!(ax2,rp13[3,:],markersize=4,color=:blue)

    r = 0.002
    it = 1
    points = [Point3f0(pp[i][it,1],pp[i][it,2],pp[i][it,3]) for i = 1:8]
    annotations!(ax1,
            ["P$i" for i = 1:8],
            points,
            rotation = qrotation(Vec3f(0.0,0.0,1.0), pi)*qrotation(Vec3f(1.0,0.0,0.0), pi/2),
            align = (:left, :baseline),
            textsize = 0.02,
            space = :data
    )
    spheres = [Node(Makie.Sphere(points[i], r)) for i = 1:8]
    for i = 1:8
        mesh!(ax1,spheres[i],color=:lightgrey)
    end
    ls_step = labelslider!(fig,"step",1:size(pp[1],1))
    fig[2,2] = ls_step.layout
    on(ls_step.slider.value) do this_step
        # analyse_slackness(tg,sol.qs[this_step])
        for i = 1:8
            point = Point3f0(pp[i][this_step,1],pp[i][this_step,2],pp[i][this_step,3])
            spheres[i][] = Makie.Sphere(point, r)
        end
        if linked
            update_scene!(tg,ax1,strings,point_spheres,traj.qs[this_step])
        end
        # @show TR.get_strings_len(tg)
        @show TR.get_strings_deform(tg)
    end
    vlines!(ax2, ls_step.slider.value, color = :red)
    fig
end
plot_super(pp,spine2,pyramid)

# spine2 = spine_true(2,-0.1,RotY(0.0);o1=[0.134074,0.119507,0.234401],X1=Matrix(RotYZ(-π/2,π/3-0.18)),dir=[1,0,0])
range_index = [
                [106,207],
                [106,207],
                [106,207],
                [106,207],
                [106,207],
                [299,369],
                [689,769],
                [1065,1128],
                [1312,1404],
                [1554,1671],
                [1790,1959]
                ]
function averaging(pp,range_index)
        [begin
            df = DataFrame()
            names = ["X$j","Y$j","Z$j"]
            for i = 1:3
                df[!,names[i]] = [mean(pdata[is:ie,i]) for (is,ie) in range_index]
            end
            df
        end
            for (j,pdata) in enumerate(pp)]
end
pp_avg = averaging(pp,range_index)

plot_super(pp_avg,spine2,pyramid)

FileIO.save("qend.jld2", "qend", qend)
plotstructure(spine2,pyramid)

function static_loading(bot_input,dynfuncs)
    bot = deepcopy(bot_input)
    startsols = get_forward_start(bot,dynfuncs)
    @unpack tg, traj = bot
    C = tg.rigidbodies[2].state.cache.Cp[1]
    T = TR.build_Ti(tg,2)
    F1 = reshape(transpose(C*T)*[0.0,0.0,-1.0],:,1)
    G = TR.build_G!(tg)
    F = hcat(G,F1)
    u0 = TR.get_strings_restlen(bot)
    gs = [[0.0,0.0]]
    for i = 0:10
        push!(gs,[1.0,i*9.8*0.1])
    end
    parameters = [(u=u0,g=gi) for gi in gs]
    TR.reset!(bot)
    q0,s0,λ0 = startsols
    TR.check_static_equilibrium(tg,q0,λ0)

    TR.check_static_equilibrium(tg,qzero,λzero)
    seqs = TR.forward_multi_sequence(tg,(q=qzero,s=szero,λ=λzero),parameters,F;n=1)
    # append!(traj.qs,seq.q)
    # append!(traj.ts,[g1[1] for g1 in seq.g])
    bot,seqs
end

spine2_loaded, seqs = static_loading(spine2,dynfuncs)
TR.set_new_initial!(spine2_loaded,seqs[1][end].q)
for seq in seqs[begin+1:end]
    push!(spine2_loaded.traj.ts,seq[end].g[2])
    push!(spine2_loaded.traj.qs,seq[end].q)
end
plot_super(pp,spine2_loaded,pyramid)
plot_super(pp_avg,spine2_loaded,pyramid;linked=true)
@descend_code_warntype static_loading(spine2,dynfuncs)
plotstructure(spine2_loaded,pyramid)
