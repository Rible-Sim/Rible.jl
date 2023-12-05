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
using HomotopyContinuation
using Printf
using Revise
import Rible as RB
cd("examples/DualTriangles")
includet("man_define.jl")
includet("man_plotting.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_theme!(font = "Nimbus Rom No9 L", fontsize = fontsize)
# set_theme!(font = "Informal Roman", fontsize = fontsize)
# set_theme!(font = "Ink Free", fontsize = fontsize)
k = 453.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(1,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=false)
plotstructure(man_inv)
function get_actuate_seqs(bot,g=0.0)
    # plotstructure(man_inv)
    Y = RB.build_Y(bot)
    # q,_ = RB.get_coords(bot.structure)
    # λ,u,a = RB.inverse(bot,deepcopy(bot),Y)
    start_sol,parameters0 = RB.get_start_system(bot,Y)
    # start_sol
    # parameters0.u
    # parameters0.g
    # RB.actuate!(bot,a)
    # RB.check_static_equilibrium(bot.structure,q,λ)
    # 0.182-5.2024/292.14
    u1 = copy(parameters0.u)
    a1 = 0.01
    u1[[1]] .+= a1
    u1[[2]] .-= a1
    u2 = copy(parameters0.u)
    a2 = -0.01
    u2[[1]] .+= a2
    u2[[2]] .-= a2
    # u1,u2

    g1 = g
    g2 = g

    # result = RB.forward_once(bot.structure,(q,s,λ),(u,0.0),(u1,g1))

    parameters1 = (u=u1, g=[g1])
    parameters2 = (u=u2, g=[g2])
    parameters_seq = [parameters0,parameters1,parameters2]

    man_sys,_ = RB.forward_system(bot.structure)

    # seq = RB.forward_sequence(man_sys,start_sol,start_parameters,target_parameters)
    seqs = RB.forward_multi_sequence(man_sys,start_sol,parameters_seq;n=20)

    # RB.set_new_initial!(bot,result.q)
    # plotstructure(bot)
    # get_angles(bot)
    # plotstructure(bot,seqs[1])
    # plotstructure(bot,seqs[2])
    seqs[2]
end
seq = get_actuate_seqs(man_inv)
seq_g = get_actuate_seqs(man_inv,1.0)
plotstructure(man_inv,seq_g)
function plot_angles(bot,seq,seq_g)
    @unpack traj, st = bot

    df = DataFrame(XLSX.readtable("弹簧2验证水平0924.xlsx", "Sheet2";infer_eltypes=true)...)
    Δu_read = 100df.a
    θ_read = df.theta
    θ_g_read = df.theta_G

    # seq 1
    u = [u[1] for u in seq.u]
    RB.set_new_initial!(bot,seq.q[1])
    append!(traj.qs, seq.q[2:end])

    θ = [
        begin
            RB.distribute_q_to_rbs!(st,q)
            get_angles(bot)[1]
        end
        for q in seq.q
    ]

    # seq 2
    u_g = [u[1] for u in seq_g.u]
    RB.set_new_initial!(bot,seq_g.q[1])
    append!(traj.qs, seq_g.q[2:end])

    θ_g = [
        begin
            RB.distribute_q_to_rbs!(st,q)
            get_angles(bot)[1]
        end
        for q in seq_g.q
    ]


    fig_width = 245.71811u"pt"
    fig_height = 110.89u"pt"
    font_size = 10u"pt"
    dpi=300/(1u"inch")
    font_pixel = to_resolution(dpi,font_size)
    size = to_resolution.(dpi,(fig_width,fig_height))
    # @show resolution
    fig = Figure(;resolution)
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[1,2])
    scatter!(ax1,Δu_read,θ,markersize=20)
    scatter!(ax1,Δu_read,θ_read,markersize=20,marker='▲')
    scatter!(ax2,Δu_read,θ_g,markersize=20,label="Sim.")
    scatter!(ax2,Δu_read,θ_g_read,markersize=20,marker='▲',label="Exp.")
    axislegend(ax2,position = :rt)
    xlims!(ax1,-1.2,1.2)
    xlims!(ax2,-1.2,1.2)
    ylims!(ax1,-36,36)
    ylims!(ax2,-36,36)
    ax2.yticksvisible = false
    ax2.yticksvisible = false
    ax2.yticklabelsvisible = false
    ax1.xlabel = "Δu (cm)"
    ax1.ylabel = "θ (°)"
    ax2.xlabel = "Δu (cm)"
    fig
end

fig_angles = plot_angles(man_inv,seq,seq_g)
GLMakie.activate!()
CairoMakie.activate!(type="svg")
CairoMakie.save(raw"C:\Users\luo22\OneDrive - 中山大学\Papers\ForwardStatics\ieeeconf\images\fig_angles.svg", fig_angles)

k = 1203.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(4,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
plotstructure(man_inv)

function get_actuate_seqs(bot,g=0.0)
    # plotstructure(man_inv)
    Y = RB.build_Y(bot)
    # q,_ = RB.get_coords(bot.structure)
    # λ,u,a = RB.inverse(bot,deepcopy(bot),Y)
    start_sol,parameters0 = RB.get_start_system(bot,Y)
    # start_sol
    # parameters0.u
    # parameters0.g
    # RB.actuate!(bot,a)
    # RB.check_static_equilibrium(bot.structure,q,λ)
    # 0.182-5.2024/292.14
    u1 = copy(parameters0.u)
    u2 = copy(parameters0.u)
    a2 = 0.01
    u2[[1,3,5]] .+= a2
    u2[[2,4,6]] .-= a2
    u3 = copy(parameters0.u)
    u4 = copy(parameters0.u)
    a4 = -0.01
    u4[[1,3,5]] .+= a4
    u4[[2,4,6]] .-= a4
    # u1,u2

    g1 = 1.0g
    g2 = 1.0g
    g3 = 1.0g
    g4 = 1.0g

    # result = RB.forward_once(bot.structure,(q,s,λ),(u,0.0),(u1,g1))

    parameters1 = (u=u1, g=[g1])
    parameters2 = (u=u2, g=[g2])
    parameters3 = (u=u3, g=[g3])
    parameters4 = (u=u4, g=[g4])
    parameters_seq = [parameters0,parameters1,parameters2,parameters3,parameters4]

    man_sys,_ = RB.forward_system(bot.structure)

    # seq = RB.forward_sequence(man_sys,start_sol,start_parameters,target_parameters)
    seqs = RB.forward_multi_sequence(man_sys,start_sol,parameters_seq;n=10)

    # RB.set_new_initial!(bot,result.q)
    # plotstructure(bot)
    # get_angles(bot)
    # plotstructure(bot,seqs[1])
    # plotstructure(bot,seqs[2])
    seqs
end

k = 403.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs1 = get_actuate_seqs(man_inv,0.0)

man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs2 = get_actuate_seqs(man_inv,1.0)

k = 1203.09; c = 100.0;
man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs3 = get_actuate_seqs(man_inv,1.0)

plotstructure(man_inv,seqs[4])

function plot_sweeping!(ax,bot,seq1,seq2)
    @unpack st, traj = bot

    for i in 1:length(seq1.q)
        RB.distribute_q_to_rbs!(st,seq1.q[i])
        lower_bars = get_lower_bars(st)
        for bar in lower_bars
            linesegments!(ax, bar, color = :grey, linewidth = 3)
        end
    end

    RB.set_new_initial!(bot,seq1.q[begin])
    append!(traj.qs,seq1.q[begin+1:end])
    endpoint_line = get_trajectory(bot,4,2)
    # @show endpoint_line
    lines!(ax, endpoint_line, color = :grey, linewidth = 3)

    for i in 1:length(seq2.q)
        RB.distribute_q_to_rbs!(st,seq2.q[i])
        upper_bars = get_upper_bars(st)
        for bar in upper_bars
            linesegments!(ax, bar, color = :grey, linewidth = 3)
        end
    end

    RB.set_new_initial!(bot,seq2.q[begin])
    append!(traj.qs,seq2.q[begin+1:end])
    endpoint_line = get_trajectory(bot,4,2)
    # @show endpoint_line
    lines!(ax, endpoint_line, color = :grey, linewidth = 3)

    bars,strings = plotstructure!(ax,st)
    update_scene!(st,bars,strings,seq1.q[begin])
    ax.aspect = DataAspect()
end

function plot_sweeping_compare(bot,seqs1,seqs2,seqs3)
    fig_width = 1.5*245.71811u"pt"
    fig_height = 150.89u"pt"
    font_size = 10u"pt"
    dpi=300/(1u"inch")
    font_pixel = to_resolution(dpi,font_size)
    size = to_resolution.(dpi,(fig_width,fig_height))
    fig = Figure(;resolution)
    ax1 = fig[1, 1] = Axis(fig)
    ax2 = fig[1, 2] = Axis(fig)
    ax3 = fig[1, 3] = Axis(fig)
    plot_sweeping!(ax1,bot,seqs1[2],seqs1[4])
    plot_sweeping!(ax2,bot,seqs2[2],seqs2[4])
    plot_sweeping!(ax3,bot,seqs3[2],seqs3[4])
    for ax in [ax1,ax2,ax3]
        xlims!(ax,-0.1,0.7)
        ylims!(ax,-0.5,0.4)
        ax.xticks = [0.0,0.2,0.4,0.6]
        ax.yticks = [-0.4,-0.2,0.0,0.2]

    end
    ax2.yticksvisible = false
    ax2.yticksvisible = false
    ax2.yticklabelsvisible = false
    ax3.yticksvisible = false
    ax3.yticksvisible = false
    ax3.yticklabelsvisible = false
    ax1.ylabel = "y (m)"
    ax1.xlabel = "x (m)"
    ax2.xlabel = "x (m)"
    ax3.xlabel = "x (m)"
    for (ilabel,label) in enumerate(["(a)", "(b)", "(c)"])
        Label(fig.layout[1, ilabel, TopLeft()], label,
            textsize = font_pixel,
            # font = noto_sans_bold,
            padding = (0, 0, 0, 0),
            halign = :right)
    end
    fig
end

fig_sweeping_compare = plot_sweeping_compare(man_inv,seqs1,seqs2,seqs3)
GLMakie.activate!()
CairoMakie.activate!(type="svg")
CairoMakie.save(raw"D:\OneDrive - 中山大学\Papers\ForwardStatics\ieeeconf\images\fig_sweeping_compare.svg", fig_sweeping_compare)

plot_sweeping(man_inv,seqs[2],seqs[4])

d0 = RB.get_d(man_inv)
d1 = copy(d0)
d1[4] = 0.1
d1[5] = 0.1
d1
result = RB.deform_forward(man_inv,(q,s,λ),(d0,u,0.0),(d1,copy(u),0.0))
# result = RB.forward(man_inv,(q,s,λ,u,0.0),(copy(u),g1))
@code_warntype RB.defrom_forward(man_inv,(q,s,λ),(d0,u,0.0),(d1,copy(u),0.0))

rsols = real_solutions(result)
rsol1 = rsols[1]
qsol = rsol1[1:length(q)]
RB.distribute_q_to_rbs!(man_inv,qsol)
scene,_ = plotstructure(man_inv)
scene

fig = Figure(size = (1600, 1000))
ax1 = fig[1, 1] = Axis(fig,size = (800, 800))
bars,cables,_ = plotstructure!(ax1,man_inv)
ls_g = labelslider!(fig, "Gravity:", 0:0.0001:0.01;
            format = x -> "$(x)×9.8 m/s²")
set_startindex(ls_g.slider,0.005)

fig[2, 1] = ls_g.layout

nu = length(u)
ulabels = ["u$i:" for i in 1:nu]
uranges = [0:(l/100):l for l in ℓ]
formats = Ref(x -> begin @sprintf "%6fm" x end)
ls_u = labelslidergrid!(fig,ulabels,uranges;formats)
set_startindex.(ls_u.sliders,u)
fig[:,2] = ls_u.layout

targets_ob = [ls_g.slider.value]
for sl in ls_u.sliders
    push!(targets_ob,sl.value)
end

onany(targets_ob...) do targets...
    g1 = targets[1]
    u1 = [u for u in targets[2:end]]
    result =  RB.forward(man_inv,(q,s,λ),(u,0.0),(u1,g1))
    rsols = real_solutions(result)
    rsol1 = rsols[1]
    qsol = rsol1[1:man_inv.ncoords]
    update_scene!(man_inv,bars,cables,qsol)
end
colsize!(fig.layout, 1, Aspect(1,1))
colsize!(fig.layout, 2, Relative(2/5))
# rowsize!(fig.layout, 1, Relative(10/10))
fig
