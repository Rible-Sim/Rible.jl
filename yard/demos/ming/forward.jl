using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using GeometryBasics
using RecursiveArrayTools
using FileIO
using Meshing
using DynamicPolynomials
using Printf
using CoordinateTransformations
using ColorTypes
using Makie
using Statistics
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using TypeSortedCollections
using EponymTuples
using Cthulhu
using BenchmarkTools
# using Observables
using Unitful, Match
using HomotopyContinuation

using Revise
import Rible as RB

cd("examples/ming")
includet("links_define.jl")
# includet("../links/link_u_plotting.jl")
includet("../analysis.jl")
includet("../vis.jl")

const figdir=raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ff"

spine2 = spine_true(2,0.105,RotY(0.0))
# @code_warntype spine_true(2,0.105,RotY(0.0))
# Φ = RB.build_Φ(spine2.st)
# q, _ = RB.get_coords(spine2.st)
# Φ(q)
# @code_warntype Φ(q)
# A = RB.build_A(spine2.st)
# A(q)
# @code_warntype A(q)
plotstructure(spine2)

bodyid = 2
pid = 1
f = [1.0,0.0,0.0]
function get_seqs(bot,bodyid,pid,f=[1.0,0.0,0.0],g2=[5.0])

    # Y = RB.build_Y(bot)

    start_sol,parameters0 = RB.get_start_system(bot,RB.DeformMode())
    @show inv.(start_sol.s)
    d0, u0, g0 = parameters0
    u1 = RB.get_cables_restlen(bot)
    @show u0
    @show u1
    # d1 = RB.get_d(bot.st)
    parameters1 = (d=d0, u=u1, g=g0)
    parameters2 = (d=d0, u=u1, g=g2)
    d3 = copy(d0)
    d3[1] = d3[2] = d3[7] = d3[8] = sqrt(0.5)
    parameters3 = (d=d3, u=u1, g=g2)
    parameters_seq = [parameters0,parameters1,parameters2,parameters3]

    F̌ = RB.build_F̌(bot.st,bodyid,pid,f)
    # sys,_ = RB.forward_system(bot.st,RB.DeformMode();F̌)

    # seq = RB.forward_once(bot.st,start_sol,parameters0,parameters1,RB.DeformMode();F)
    # seq = RB.forward_sequence(bot_sys,start_sol,start_parameters,target_parameters)
    seqs = RB.forward_multi_sequence(bot.st,start_sol,parameters_seq,RB.DeformMode();F̌,n=10)
    # @btime RB.forward_multi_sequence($sys,$start_sol,$parameters_seq;n=1)
    # seqs[1]
    # fig = plotstructure(bot,seqs[1])
    # fig = plotstructure(bot,seqs[2])
    # fig = plotstructure(bot,seqs[3])
end
spine2 = spine_true(2,0.105,RotY(0.0)); spine2_video = deepcopy(spine2)
seqs = get_seqs(spine2,bodyid,pid,f)

RB.apply!(spine2,seqs[3])
plot_traj!(spine2;)

# @code_warntype get_seqs(spine2,bodyid,pid,f)
seqs
qs = [seqs[2][begin].q,seqs[3][begin].q,seqs[3][end].q]
RB.set_new_initial!(spine2,qs[1])
RB.get_cables_len(spine2)
h = spine2.st.rigidbodies[2].state.loci_states[pid][3] - spine2.st.rigidbodies[1].state.loci_states[pid][3]
h
nsteps = reduce(vcat,[seqs[i].nstep[2:end] for i = 1:3])
stats_nstep = nsteps |> (x) -> [minimum(x), maximum(x), median(x), mean(x)]
stats_nstep

function set_axis!(ax)
    # hidexdecorations!(ax)
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.zticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.zticklabelsvisible = false
    ax.xlabelvisible = false
    ax.ylabelvisible = false
    ax.zlabelvisible = false
    # hidespines!(ax)
    ax.xspinecolor_1 = :lightgrey
    ax.xspinecolor_3 = :lightgrey
    ax.yspinecolor_1 = :lightgrey
    ax.yspinecolor_3 = :lightgrey
    ax.zspinecolor_1 = :lightgrey
    ax.zspinecolor_3 = :lightgrey

    # ax.xypanelcolor = :lightgrey
    # ax.xzpanelcolor = :lightgrey
    # ax.yzpanelcolor = :lightgrey
    xlims!(ax,-0.10,0.10)
    ylims!(ax,-0.10,0.10)
    zlims!(ax,-0.00,0.25)
end

function plot_force_arrow!(ax,st,bodyid,pid,f,font_pixel)
    rbs = st.rigidbodies
    body = rbs[bodyid]
    rp = body.state.cache.Cp[pid]*body.state.coords.q
    points = [Point3f(rp[1],rp[2],rp[3])]
    directions = [Point3f(f[1],f[2],f[3])]
    size_factor = 0.05
    arrows!(ax, points, directions,
            fxaa = true,
            arrowsize = size_factor.*Vec3f(0.2, 0.2, 0.3),
            linewidth = size_factor*0.05,
            color = :darkred, lengthscale = 0.025)
    text!(ax,["f"],position = points, color = :darkred, textsize = font_pixel,
            align = (:left, :baseline), offset = (-5, 10))
end

function plot2(bot,qs,bodyid,pid,f)
    (;st) = bot
    fig_width = 165.0u"pt"
    fig_height = 90.89u"pt"
    font_size = 10u"pt"
    dpi=300/(1u"inch")
    font_pixel = to_resolution(dpi,font_size)
    resolution = to_resolution.(dpi,(fig_width,fig_height))
    @show resolution
    fig = Figure(;resolution)
    eyepos = Vec3f(0.21623512, 0.5631353, 0.21988256)
    lookat = Vec3f(0.057202004, 0.018000102, 0.12657055)
    azimuth = 1.385530633326996
    elevation = 0.11269908169872413
    ax1 = Axis3(fig[1, 1]; aspect = :data, viewmode = :fit, azimuth, elevation)
    ax2 = Axis3(fig[1, 2]; aspect = :data, viewmode = :fit, azimuth, elevation)
    # ax1 = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    # ax2 = fig[1, 2] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    # ax3 = fig[1, 3] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    # bars1,cables1 = plotstructure!(ax1,st)
    bars2,cables2 = plotstructure!(ax1,st)
    bars3,cables3 = plotstructure!(ax2,st)
    # update_scene!(st,bars1,cables1,qs[1])
    update_scene!(st,bars2,cables2,qs[2])
    plot_force_arrow!(ax1,st,bodyid,pid,f,font_pixel)
    update_scene!(st,bars3,cables3,qs[3])
    plot_force_arrow!(ax2,st,bodyid,pid,f,font_pixel)
    # set_axis!(ax1)
    set_axis!(ax1)
    set_axis!(ax2)
    # cameracontrols(ax.scene).mouse_rotationspeed[] = 0.01
    # update_cam!(ax1.scene, eyepos, lookat)
    # update_cam!(ax2, eyepos, lookat)
    # update_cam!(ax3, eyepos, lookat)
    # sequenceplot(fig,grid_tg,ax,st,bars,cables,update_scene!,seq,dlabels,ulabels,glabels)
    colsize = 500
    rowsize = 480
    for i = 1:2
        colsize!(fig.layout, i, Fixed(colsize))
    end
    rowsize!(fig.layout, 1, Fixed(rowsize))
    colgap!(fig.layout,-250)
    for (ilabel,label) in enumerate(["(a)", "(b)"])
        Label(fig.layout[1, ilabel, TopLeft()], label,
            textsize = font_pixel,
            # font = noto_sans_bold,
            padding = (0, -110, -180, 0),
            halign = :right)
    end
    fig,ax1
end
fig_3D_spine_vis,ax1 = plot2(spine2,qs,bodyid,pid,f)
fig_3D_spine_vis
GLMakie.activate!()
fig_3D_spine_vis
GLMakie.save(texroot*raw"\OneDrive - 中山大学\Papers\ForwardStatics\ieeeconf\images\fig_3D_spine_vis.tiff", fig_3D_spine_vis)

CairoMakie.activate!(type = "pdf")
CairoMakie.save(texroot*raw"\OneDrive - 中山大学\Papers\ForwardStatics\ieeeconf\images\fig_3D_spine_vis.pdf", fig_3D_spine_vis)


function plot_graph(bot,seqs,pid,figname=nothing)
    (;st, traj) = bot

    g = [g[1] for g in seqs[2].g]
    RB.set_new_initial!(bot,seqs[2].q[1])
    append!(traj.q, seqs[2].q[2:end])
    u2 = RB.get_trajectory!(bot,2,pid)

    d = [d[7] for d in seqs[3].d]
    RB.set_new_initial!(bot,seqs[3].q[1])
    append!(traj.q, seqs[3].q[2:end])
    u3 = RB.get_trajectory!(bot,2,pid)

    with_theme(theme_pub;
            resolution = (0.7tw,0.4tw),
            Scatter = (
                markersize = fontsize,
            )
        ) do
        fig = Figure(;)
        ax1, ax2 = [Axis(fig[1,i]) for i = 1:2]
        ax3, ax4 = [Axis(fig[2,i]) for i = 1:2]
        scatter!(ax1,g,u2[1,:],marker=:circle)
        scatter!(ax3,g,u2[3,:],marker=:utriangle)
        # scatter!(ax1,d,u3[1,:])
        dx = scatter!(ax2,d.^2,u3[1,:],marker=:circle,   label=L"x")
        dz = scatter!(ax4,d.^2,u3[3,:],marker=:utriangle,label=L"z")
        # ax1.xreversed = true
        # Legend(fig[1, 3], [dx, dz], ["x", "z"])
        axislegend(ax2,position = :rb)
        axislegend(ax4,position = :rt)
        ylims!(ax1, -.02,.10); ylims!(ax2, -.02,.10)
        ylims!(ax3,  .21,.27); ylims!(ax4,  .21,.27)
        ax1.yticks = -0.04:0.04:0.13; ax2.yticks = -0.04:0.04:0.13
        ax3.yticks =   .22:0.02: .26; ax4.yticks =   .22:0.02: .26
        xlims!(ax1,0,5); xlims!(ax2,0.5,1.0)
        xlims!(ax3,0,5); xlims!(ax4,0.5,1.0)
        ax1.xticks = 0:1:5
        ax3.xticks = 0:1:5;
        ax2.xreversed = true
        ax4.xreversed = true
        # ax1.ylabel = "Coordinate (m)"
        hidex(ax1); hidex(ax2)
        hidey(ax2); hidey(ax4)
        ax3.xlabel = "Force (N)"
        ax4.xlabel = "Area Ratio"
        # colgap!(fig.layout,60)
        # colsize!(fig.layout, 1, Relative(1/2))
        # colsize!(fig.layout, 2, Relative(1/2))
        # rowsize!(fig.layout, 1, Fixed(rowsize))
        for (ilabel,label) in enumerate(["(a)", "(b)"])
            Label(fig.layout[1, ilabel, TopLeft()], label,)
        end
        for (ilabel,label) in enumerate(["(c)", "(d)"])
            Label(fig.layout[2, ilabel, TopLeft()], label,)
        end
        Label(fig[1:2, 0], "Coordinates (m)", rotation = pi/2)
        savefig(fig,figname)
    end
end

GM.activate!(); plot_graph(spine2,seqs,pid)
CM.activate!(); plot_graph(spine2,seqs,pid,"fig_3D_spine_graph")

RB.set_new_initial!(spine2_video,seqs[1][end].q)
seq_sol = seqs[1][end]

start_sol,parameters1 = RB.get_start_system(spine2_video,RB.build_Y(spine2),RB.DeformMode())
parameters2 = deepcopy(parameters1)
parameters2.g[1] = 5.0
F = RB.build_F(spine2.st,bodyid,pid,f)
seq2 = RB.forward_sequence(spine2_video.st,start_sol,parameters1,parameters2,RB.DeformMode();F)
fig = plotstructure(spine2,seq2)

function plot_video(bot,seq_sol,bodyid,pid,f)
    fig = Figure(resolution = (1280, 720))
    azimuth = 1.385530633326996
    elevation = 0.05269908169872413
    # ax1 = Axis3(fig[1, 1]; aspect = :data, viewmode = :fit, azimuth, elevation)
    ax1 = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,cables = plotstructure!(ax1,bot.st)
    eyepos = [0.056068808, 0.5011266, 0.18235058]
    lookat = [0.038956024, 0.0140809715, 0.124568716]
    update_cam!(ax1.scene,eyepos,lookat)
    # start_sol,parameters0 = RB.get_start_system(bot,RB.DeformMode())
    start_sol = (q=seq_sol.q,s=seq_sol.s,λ=seq_sol.λ)
    parameters0 = (d=seq_sol.d,u=seq_sol.u,g=seq_sol.g)
    F = RB.build_F(bot.st,bodyid,pid,f)
    d0,u0,g0 = parameters0
    ls_g = labelslider!(fig, "f:", -5:0.5:5; format = x -> "$(x)N")
    fig[2, 1] = ls_g.layout
    set_startindex(ls_g.slider,0.0)
    # nd = 12
    # dlabels = ["d$i:" for i in 1:nd]
    # dmin, dmax = 0.1, 2.0
    # dranges = Ref([0.0,0.1,0.25,0.5,1.0,1.25,1.5,2.0])
    # formats = Ref(x -> begin @sprintf "%.2fm" x end)
    # ls_d = labelslidergrid!(fig,dlabels,dranges;formats)
    # set_startindex.(ls_d.sliders,d0[1:nd])
    # fig[:,2] = ls_d.layout
    nrad = 2
    rlabels = ["Area $i:" for i in 1:nrad]
    radmin, radmax = 0.5, 1.0
    radranges = Ref(collect(radmin:0.05:radmax))
    formats = Ref(x -> begin @sprintf "%.2f" x end)
    ls_rad = labelslidergrid!(fig,rlabels,radranges;formats)
    set_startindex.(ls_rad.sliders,fill(radmax,2))
    fig[:,2] = ls_rad.layout

    colsize!(fig.layout, 1, Relative(2/3))
    targets_ob = [ls_g.slider.value]
    for sl in ls_rad.sliders
        push!(targets_ob,sl.value)
    end
    onany(targets_ob...) do targets...
        g1 = targets[1]
        rad1 = [rad for rad in targets[2:end]]
        d1 = copy(d0)
        d1[1:2] .= rad1[1]
        d1[7:8] .= rad1[2]
        parameters1 = (d=d1, u=u0, g=[g1])
        result =  RB.forward_once(bot.st,start_sol,parameters0,parameters1,RB.DeformMode();F)
        update_scene!(bot.st,bars,cables,result.q)
        # camera = cameracontrols(ax1.scene)
        # lookat = camera.lookat[]
        # eyepos = camera.eyeposition[]
        # @show eyepos, lookat
    end
    Label(fig[0, :], text = "Forward statics interactive simulation \n of a 2-segment tensegrity spine", textsize = 45)
    fig
end

plot_video(spine2_video,seq_sol,bodyid,pid,f)

# Comparison

includet("../geodesicDR.jl")
function compare(nb)
    pid = 1
    spine_compare = spine_true(nb,0.117,RotY(0.0))
    rs, bs, ss, spine_initial_DR = GDR(spine_compare;β=1e-2,maxiters=10000,ϵ=1e-14,N=50,ξ=1e-14)
    @info "Initial Residual: $(rs[end])"
    f = [5.0,0.0,0.0]
    F = reshape(RB.build_F(spine_initial_DR.st,nb,pid,f),:,1)

    ################### Testing DR #########################
    @info "Testing DR"
    @info "Testing for ϵ=1e-13"
    @info "Precompiling"
    GDR(spine_initial_DR,F;β=1e-2,maxiters=10000,ϵ=1e-13,N=50,ξ=1e-13)
    @info "Benchmarking"
    rs13,_,_,spine_DR13 = @btime $GDR($spine_initial_DR,$F;β=1.5e-2,maxiters=10000,ϵ=1e-13,N=500,ξ=1e-13)
    @info "Step = $(length(rs13)); Residual = $(rs13[end])"
    spine_DR13

    ################### Testing HC #########################
    # @info "Testing HC"
    # start_sol,parameters0 = RB.get_start_system(spine_initial_DR,RB.PrimalMode())
    # u0, g0 = parameters0
    # # @show u0
    # u1 = RB.get_cables_restlen(spine_initial_DR)
    # # @show u1
    # parameters1 = (u=u1, g=g0)
    # parameters2 = (u=u1, g=[1.0])
    # @info "Precompiling"
    # sol1 = RB.forward_once(spine_initial_DR.st,start_sol,parameters0,parameters1,RB.PrimalMode();F)
    # @info "Benchmarking"
    # sol2 = @btime $RB.forward_once($spine_initial_DR.st,$sol1,$parameters1,$parameters2,$RB.PrimalMode();F=$F)
    # @info "Step = $(sol2.nstep); Residual = $(sol2.res)"
end

function compare_loop(irange)
    for i = irange
        @info "Testing $i-segment."
        compare(i)
    end
end

compare_loop(2:5)

@descend_code_warntype compare(2)
plot(rs,axis = (yscale=log10,))
@time GDR(spine_compare;maxiters=900,ϵ=1e-14,N=10,ξ=1e-7)

rs[end]
plotstructure(spine_initial_DR)

RB.set_new_initial!(spine_initial_DR,spine_initial_DR.traj.qs[end])
plotstructure(spine_GDR)


time_GDR.time
rs[end]
plot(rs,axis = (yscale =log10,))
plot(bs,axis = (yscale = log10,))
plot(ss)

# rs,spine_initial_DR = compare(spine_compare,nb,pid)
# rs[end]
#
time_forward = compare(spine_initial_DR,nb,pid,F)
@report_call compare(spine_initial_DR,nb,pid,F)
@descend_code_warntype compare(spine_initial_DR,nb,pid,F)
time_forward.time
time_forward.value.nstep
time_forward.value.res
RB.set_new_initial!(spine_compare,time_forward.value.q)
plotstructure(spine_compare)

(spine_GDR.traj.qs[end] - time_forward.value.q) .|> abs |> maximum

plotstructure(spine_DR)
plot(rs,axis = (yscale = log10,))
@time GDR(spine2.st;maxiters=1000,ϵ=1e-14)
@descend_code_warntype GDR(spine2.st;maxiters=1000,ϵ=1e-14)o
includet("tsff.jl")
RB.get_cables_len(spine2)
errs,Ns = tsff(spine2.st)
plot(errs)
errs[end]
@time tsff(spine2.st)
@descend_code_warntype tsff(spine2.st)
[body.state.ro for rb in spine2.st.rigidbodies]

identical_idx = [1,2,3,4] .== [1,2,3,5]
a = collect(5:8)
a[identical_idx]
.!(identical_idx)
