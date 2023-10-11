#note -- preamble
using OhMyREPL
using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using ElasticArrays
using TypeSortedCollections
using TypedTables
using Rotations
# import GeometricalPredicates as GP
using CoordinateTransformations
using OffsetArrays
import DifferentialEquations as DE
using BenchmarkTools
using RecursiveArrayTools
using CircularArrays
const CA = CircularArray
using Interpolations
# using CubicSplines
# import FLOWMath
using Makie
import GLMakie as GM
import CairoMakie as CM
import WGLMakie as WM
GM.activate!()
using FileIO, MeshIO
using LaTeXStrings
using TexTables
using DataStructures
using Latexify
using PrettyTables
auto_display(false)
using Unitful
using EzXML
using CSV, Tables
using Printf
using StructArrays
using EponymTuples
import GeometryBasics as GB
using Meshing
# using FloatingTableView
import Meshes
using Match
using Cthulhu
using COSMO
using AbbreviatedStackTraces
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true
ENV["JULIA_STACKTRACE_MINIMAL"] = true
using Revise
import TensegrityRobots as TR
cd("examples/ES")
include("../vis.jl")
include("../analysis.jl")
include("../dyn.jl")
include("dyn.jl")
include("def.jl")
include("def3d.jl")
include("../bridge/bridge.jl")

includet("../vis.jl")
includet("../analysis.jl")
includet("../dyn.jl")
includet("dyn.jl")
includet("def.jl")
includet("def3d.jl")
includet("../bridge/bridge.jl")
figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\DynamicTensegrity\ES"
    # figdir::String =raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\dyn"
elseif Sys.isapple()
    figdir::String = raw"."
end
fontsize = 8 |> pt2px
markersize = fontsize
linewidth = 0.5 |> pt2px
cablewidth = 0.75 |> pt2px
barwidth = 1.5 |> pt2px
cw = 469 |> pt2px
tw = 469 |> pt2px
#----- preamble end 

# fontsize = 10 |> pt2px
# cw = lw = tw = 455.24411 |> pt2px

## 1st example, validation and verification
## comparison with Adams
obot = one_bar_one_tri()

prob = TR.SimProblem(obot,(x)->dynfuncs(x;gravity=true))

dt = 1e-3

TR.solve!(prob,TR.Zhong06();dt,tspan=(0.0,4.0),ftol=1e-14)

with_theme(theme_pub;) do 
    plot_traj!(
            obot;
            xlims=(-0.2,0.2),
            ylims=(-0.3,0.1),
            showmesh=false,
            showinfo=false,
            showlabels=false,
            sup! = plot_one_bar_one_tri!,
            dorecord=true,
            showinit=false,
            figname="one_bar_one_tri.mp4"
        )
end
# adams = parse_AdamsResults("./res_1e-5.xml")
ode_sol,ode_M = get_ref_sol(;saveat=obot.traj.t)

function plot_comparison(sol,M,bot,figname=nothing)
    sp = 1; ss = 50

    ω1 = get_angular_velocity!(bot,1)[1,begin:sp:end]
    ω2 = get_angular_velocity!(bot,2)[1,begin:sp:end]

    sol_θ = sol[3:4,begin:ss:end]
    sol_ω = VectorOfArray(
        [
            ode_M(sol[3:4,i]...)\sol[1:2,i]
            for i in 1:ss:length(sol)
        ]
    )

    t = bot.traj.t[begin:sp:end]
    sol_t = sol.t[begin:ss:end]
    rb1rp2 = get_trajectory!(bot,1,2)[begin:sp:end]
    θ1 = atan.(rb1rp2[1,:],-rb1rp2[2,:])#.-sol[3,:]
    rb2rg = get_trajectory!(bot,2,0)[begin:sp:end]
    rb2rg_ = rb2rg .- rb1rp2
    θ2 = atan.(rb2rg_[1,:],-rb2rg_[2,:])#.-sol[4,:]

    with_theme(theme_pub;
            resolution = (0.8cw,0.5cw),
            palette = (
                color = [:black],
                red = [:red]
            ),
            Scatter = (cycle = [:color => :red],),
        ) do
        fig = Figure()
        ax1 = Axis(fig[1,1]; ylabel = L"\theta_1~(\mathrm{Rad})")
        ax2 = Axis(fig[2,1]; ylabel = L"\theta_2~(\mathrm{Rad})")
        ax3 = Axis(fig[3,1]; ylabel = L"\dot{\theta}_1~(\mathrm{Rad/s})")
        ax4 = Axis(fig[4,1]; ylabel = L"\dot{\theta}_2~(\mathrm{Rad/s})", xlabel = tlabel)
        axs = (ax1,ax2,ax3,ax4)
        marker='×'#'○'
        # scales = [-2,-2,-1,-1]
        lines!(ax1,t,θ1;label="Proposed")
        scatter!(ax1,sol_t,sol_θ[1,:];label="Minimal",marker)
        lines!(ax2,t,θ2;label="Proposed")
        scatter!(ax2,sol_t,sol_θ[2,:];label="Minimal",marker)
        lines!(ax3,t,ω1./10;label="Proposed")
        @show length(sol_t), length(sol_ω)
        scatter!(ax3,sol_t,sol_ω[1,:]./10;label="Minimal",marker)
        lines!(ax4,t,ω2./10;label="Proposed")
        scatter!(ax4,sol_t,sol_ω[2,:]./10;label="Minimal",marker)
        # axislegend()

        for ax in axs[1:3]
            hidex(ax)
        end
        for ax in axs
            xlims!(ax,bot.traj.t[begin],bot.traj.t[end])
            ax.ytickformat = "{:.1f}"
        end
        ax4.xlabelpadding = -15
        for (ilabel,label) in enumerate(alphabet[1:4])
            Label(fig.layout[ilabel, 1, TopLeft()], "($label)",
                font = "CMU Serif Bold",
                )
        end
        for (ilabel,label) in zip([3,4],[1,1])
            Label(fig.layout[ilabel, 1, Top()], latexstring("×10^{$label}"),)
        end
        # fig[1, 2] = Legend(fig, ax1, tellwidth=true)
        axislegend(ax1;position=:rb,orientation=:horizontal)
        rowgap!(fig.layout, 0)
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plot_comparison(ode_sol,ode_M,obot)
CM.activate!(); plot_comparison(ode_sol,ode_M,obot,"double_comparison")

## Averaged errors
ode_sol[3,:]

rb1rp2 = get_trajectory!(obot,1,2)
θ1 = atan.(rb1rp2[1,:],-rb1rp2[2,:])
rb2rg = get_trajectory!(obot,2,0)
rb2rg_ = rb2rg .- rb1rp2
θ2 = atan.(rb2rg_[1,:],-rb2rg_[2,:])
err_θ1 = (ode_sol[3,:] .- θ1)./θ1 .|> abs
lines(err_θ1)
mean(err_θ1)

obot.traj.t

## comparison with Alpha
obot_zhong = one_bar_one_tri(); obot_alpha = deepcopy(obot)
prob_zhong = TR.SimProblem(obot_zhong,(x)->dynfuncs(x;gravity=true));
prob_alpha = TR.SimProblem(obot_alpha,(x)->dynfuncs(x;gravity=true));

dt = 1e-3
TR.solve!(prob_zhong,TR.Zhong06();dt,tspan=(0.0,500.0),ftol=1e-10)
TR.solve!(prob_alpha,TR.Alpha(0.7);dt,tspan=(0.0,500.0),ftol=1e-10)

plot_traj!(obot_zhong)

function plot_energy!(bot_zhong,bot_alpha,figname=nothing)
    ms_zhong = TR.mechanical_energy!(bot_zhong;gravity=true)
    ms_alpha = TR.mechanical_energy!(bot_alpha;gravity=true)
    with_theme(theme_pub;
            resolution=(0.6tw,0.3tw),
            figure_padding = (fontsize,fontsize,fontsize,fontsize),
            palette = (color = [:black,:red],),
            cycle = [[:linecolor, :markercolor] => :color,]
        ) do
        fig = Figure()
        scale = -2
        si = 5000
        ax1 = Axis(fig[1,1]; ylabel = L"E~(\mathrm{J})", xlabel = L"t~(\mathrm{s})")
        lines!(ax1,bot_zhong.traj.t[begin:si:end],ms_zhong.E[begin:si:end]./10.0^scale;label="Proposed",linewidth)
        lines!(ax1,bot_alpha.traj.t[begin:si:end],ms_alpha.E[begin:si:end]./10.0^scale;label="Generalized-α",linewidth)
        ax1.xticks = collect(0:100:500)
        t = bot_zhong.traj.t
        xlims!(ax1,t[begin],t[end])
        ylims!(ax1,-6.5,-4.5)
        ax1.xlabelpadding = -15
        # axislegend(ax1,position = :lb)
        axislegend(ax1; position = :lb,)
        # ax1.alignmode = Mixed(;left = -20, right = 20)
        Label(fig.layout[1, 1, Top()], latexstring("×10^{$scale}"))
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plot_energy!(obot_zhong,obot_alpha)
CM.activate!(); plot_energy!(obot_zhong,obot_alpha,"energy")


function plot_comparison_and_energy(sol,M,bot,bot_zhong,bot_alpha,figname=nothing)

    ms_zhong = TR.mechanical_energy!(bot_zhong;gravity=true)
    ms_alpha = TR.mechanical_energy!(bot_alpha;gravity=true)

    sp = 1; ss = 50

    ω1 = get_angular_velocity!(bot,1)[1,begin:sp:end]
    ω2 = get_angular_velocity!(bot,2)[1,begin:sp:end]

    sol_θ = sol[3:4,begin:ss:end]
    sol_ω = VectorOfArray(
        [
            M(sol[3:4,i]...)\sol[1:2,i]
            for i in 1:ss:length(sol)
        ]
    )

    t = bot.traj.t[begin:sp:end]
    sol_t = sol.t[begin:ss:end]
    rb1rp2 = get_trajectory!(bot,1,2)[begin:sp:end]
    θ1 = atan.(rb1rp2[1,:],-rb1rp2[2,:])#.-sol[3,:]
    rb2rg = get_trajectory!(bot,2,0)[begin:sp:end]
    rb2rg_ = rb2rg .- rb1rp2
    θ2 = atan.(rb2rg_[1,:],-rb2rg_[2,:])#.-sol[4,:]

    with_theme(theme_pub;
            resolution = (0.8cw,0.3cw),
            figure_padding = (fontsize,fontsize,0,fontsize),
            # palette = (
            #     color = [:black],
            #     red = [:red]
            # ),            
            palette = (color = [:black,:red], red = [:red]),
            cycle = [[:linecolor, :markercolor] => :color,],
            Scatter = (cycle = [:color => :red], markersize = fontsize),
        ) do
        fig = Figure()
        # ax1 = Axis(fig[1,1]; ylabel = L"\theta_1~(\mathrm{Rad})")
        ax1 = Axis(fig[1,1]; ylabel = L"\theta_2~(\mathrm{Rad})")
        # ax3 = Axis(fig[3,1]; ylabel = L"\dot{\theta}_1~(\mathrm{Rad/s})")
        ax2 = Axis(fig[2,1]; ylabel = L"\dot{\theta}_2~(\mathrm{Rad/s})", xlabel = tlabel)
        axs = (ax1,ax2,)
        marker='×'#'○'
        # scales = [-2,-2,-1,-1]
        # lines!(ax1,t,θ1;label="MSI")
        # scatter!(ax1,sol_t,sol_θ[1,:];label="RK5",marker)
        lines!(ax1,t,θ2;label="MSI")
        scatter!(ax1,sol_t,sol_θ[2,:];label="RK5",marker)
        # lines!(ax3,t,ω1./10;label="MSI")
        # @show length(sol_t), length(sol_ω)
        # scatter!(ax3,sol_t,sol_ω[1,:]./10;label="RK5",marker)
        lines!(ax2,t,ω2./10;label="MSI")
        scatter!(ax2,sol_t,sol_ω[2,:]./10;label="RK5",marker)
        # axislegend()

        ylims!(ax1,-5.5,3.0)
        for ax in axs[1:1]
            hidex(ax)
        end
        for ax in axs
            xlims!(ax,bot.traj.t[begin],bot.traj.t[end])
            ax.ytickformat = "{:.1f}"
        end
        # ax4.xlabelpadding = -15
        for (ilabel,label) in enumerate(alphabet[1:2])
            Label(fig.layout[ilabel, 1, TopLeft()], "($label)",
                font = "CMU Serif Bold",
                )
        end
        # for (ilabel,label) in zip([3,4],[1,1])
        #     Label(fig.layout[ilabel, 1, Top()], latexstring("×10^{$label}"),)
        # end
        # fig[1, 2] = Legend(fig, ax1, tellwidth=true)
        axislegend(ax1;position=:lb,orientation=:horizontal)
        # fig.layout[0, 1,] = Legend(fig,ax1,orientation=:horizontal)
        rowgap!(fig.layout, 0)
        scale = -2
        si = 5000
        ax3 = Axis(fig[:,2]; ylabel = L"E~(\mathrm{J})", xlabel = L"t~(\mathrm{s})")
        lines!(ax3,bot_zhong.traj.t[begin:si:end],ms_zhong.E[begin:si:end]./10.0^scale;label="MSI",linewidth)
        lines!(ax3,bot_alpha.traj.t[begin:si:end],ms_alpha.E[begin:si:end]./10.0^scale;label="Generalized-α",linewidth)
        ax3.xticks = collect(0:100:500)
        t = bot_zhong.traj.t
        xlims!(ax3,t[begin],t[end])
        ylims!(ax3,-6.5,-4.5)
        # ax3.xlabelpadding = -15
        # axislegend(ax3,position = :lb)
        axislegend(ax3; position = :lb,)
        # ax3.alignmode = Mixed(;left = -20, right = 20)
        Label(fig.layout[:, 2, TopLeft()], "($(alphabet[3]))",
            font = "CMU Serif Bold",
        )
        Label(fig.layout[:, 2, Top()], latexstring("×10^{$scale}"))
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plot_comparison_and_energy(ode_sol,ode_M,obot,obot_zhong,obot_alpha,)
CM.activate!(); plot_comparison_and_energy(ode_sol,ode_M,obot,obot_zhong,obot_alpha,"comparison_and_energy")

## simple
gravity = false
c = 0.0
z0 = 0.2
ωz = 5.0
mbar = 0.05
newsim_fixed = simple(;
    c,
    z0,
    ωz,
    mbar
)
GM.activate!(); with_theme(
        theme_pub;
        resolution = (0.5cw,0.25cw),
        figure_padding = (2fontsize,0,0,0),
        Axis3 = (
            azimuth = 5.13952996194027,
            elevation = 0.12269999722606788
        )
    ) do 
    fig = Figure()
    rbs = TR.get_bodies(newsim_fixed)  
    rb = rbs[2]
    gd1 = fig[1,1] = GridLayout(;tellheight = false,)
    ax = Axis3(gd1[1,1];
        aspect=:data,        
        azimuth = 0.0,
        elevation = π/2,
        # tellheight = false,
    )
    Label(
        gd1[1, 1, Top()],
        rich("(a) ", font=:bold),
        justification = :right,
        lineheight = 1.0,
        halign = :center,
        # tellheight = false,
        # alignmode = Mixed(top = 0)
        valign = :bottom,
    )
    # # plot_rigid(
    # #     rbs[2];
    #     fig,
    # #     showpoints = false,
    # #     showmesh = true,
    # # )
    mesh!(ax,build_mesh(rb;update=false);shading = false,)
    xlims!(ax,-0.15,0.15)
    ylims!(ax,-0.15,0.15)
    hidezdecorations!(ax)
    ax.xlabeloffset = 2fontsize

    gd2 = fig[1,2] = GridLayout()
    plot_traj!(newsim_fixed;
        AxisType = Axis3,
        fig = gd2,
        xlims = (-100e-3,100e-3),
        ylims = (-100e-3,100e-3),
        zlims = (-1e-3,400e-3),
        showpoints = false,
        showlabels = false,
        doslide = false,
        showground = false,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(b) ", font=:bold),
                # (@sprintf "t = %.10G (s)" tt)
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            ax.xticklabelsvisible = false
            ax.xlabeloffset = 0
            ax.yticklabelsvisible = false
            ax.ylabeloffset = 0
        end,
    )
    colsize!(fig.layout,1,Relative(0.4))
    rowsize!(gd1,1,Fixed(0.12cw))
    colgap!(fig.layout,0)
    savefig(fig,"simple")
    fig
end

tend = 2.
dt = 1e-4
# Fixed
TR.solve!(
    TR.SimProblem(newsim_fixed,(x)->dynfuncs(x;gravity)),
    TR.Zhong06();
    dt,tspan=(0.0,tend),ftol=1e-14
)

GM.activate!();plot_traj!(newsim_fixed;)
# Freed
newsim_freed = simple(;
    c,
    z0,
    ωz,
    mbar,
    free = true
)
TR.solve!(
    TR.SimProblem(newsim_freed,(x)->dynfuncs(x;gravity)),
    TR.Zhong06();
    dt,tspan=(0.0,tend),ftol=1e-14
)

GM.activate!();plot_traj!(newsim_freed;)

GM.activate!(); with_theme(
        theme_pub;
        resolution = (0.6cw,0.3cw),
        figure_padding = (0,0,0,0),
        Axis3 = (
            azimuth = 5.13952996194027,
            elevation = 0.12269999722606788
        )
    ) do 
    fig = Figure()
    nt = 5
    g1 = fig[1,1] = GridLayout()
    plot_traj!(newsim_fixed;
        AxisType = Axis3,
        fig = g1,
        gridsize = (1,nt),
        attimes = range(0,2,nt),
        xlims = (-100e-3,100e-3),
        ylims = (-100e-3,100e-3),
        zlims = (-1e-3,400e-3),
        showpoints = false,
        showlabels = false,
        doslide = false,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font=:bold),
                rich("t ", font=:math),
                (@sprintf "= %.10G (s)" tt)
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
    g2 = fig[2,1] = GridLayout(;)
    plot_traj!(newsim_freed;
        AxisType = Axis3,
        fig = g2,
        gridsize = (1,nt),
        attimes = range(0,tend,nt),
        xlims = (-100e-3,100e-3),
        ylims = (-100e-3,100e-3),
        zlims = (-1e-3,400e-3),
        showpoints = false,
        showlabels = false,
        doslide = false,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi+nt])) ", font=:bold),
                rich("t ", font=:math),
                (@sprintf "= %.10G (s)" tt)
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
    colgap!(fig.layout,-fontsize)
    rowgap!(fig.layout,0)
    savefig(fig,"simple_snapshots")
    fig
end

GM.activate!(); with_theme(theme_pub;
        resolution = (0.7cw,0.2cw),
        figure_padding = (fontsize/2,fontsize/2,0,0),
        Lines = (
            cycle = Cycle([:linestyle,:color];covary=true),
        )
    ) do 
    fig = Figure()
    (;t) = newsim_freed.traj
    r2p1_fixed = get_trajectory!(newsim_fixed,1,1)
    r2p1_freed = get_trajectory!(newsim_freed,1,1)
    # v2p1_fixed = get_mid_velocity!(newsim_fixed,1,1)
    # v2p1_freed = get_mid_velocity!(newsim_freed,1,1)
    T1_fixed = get_kinetic_energy!(newsim_fixed,1)
    T1_freed = get_kinetic_energy!(newsim_freed,1)
    E_fixed = TR.mechanical_energy!(newsim_fixed;gravity=false)
    E_freed = TR.mechanical_energy!(newsim_freed;gravity=false)
    # TR.kinetic_energy()
    # ax1 = Axis(fig[1,1],xlabel = tlabel, ylabel = L"x~(\mathrm{m})")
    # lines!(ax1,t,r2p1_fixed[1,:],label="Case 1")
    # lines!(ax1,t,r2p1_freed[1,:],label="Case 2")
   
    ax1 = Axis(fig[1,1],xlabel = tlabel, ylabel = L"\mathrm{Energy}~(\mathrm{J})")
    lines!(ax1,t,E_fixed.E,label="E")
    lines!(ax1,t,E_fixed.T,label="T")
    lines!(ax1,t,E_fixed.V,label="V")

    ax2 = Axis(fig[1,2],xlabel = tlabel, ylabel = L"\mathrm{Energy}~(\mathrm{J})")
    lines!(ax2,t,E_freed.E,label="E")
    lines!(ax2,t,E_freed.T,label="T")
    lines!(ax2,t,E_freed.V,label="V")
    
    # ax2 = Axis(fig[1,2],xlabel = tlabel, ylabel = L"x~(\mathrm{m})")
    # lines!(ax2,t,r2p1_fixed[2,:],label="Case 1")
    # lines!(ax2,t,r2p1_freed[2,:],label="Case 2")

    # ax2 = Axis(fig[1,2],xlabel = tlabel, ylabel = L"x~(\mathrm{m})")
    # lines!(ax2,t,T1_fixed,label="Case 1")
    # lines!(ax2,t,T1_freed,label="Case 2") 

    
    # ax3 = Axis(fig[1,3],xlabel = tlabel, ylabel = L"x~(\mathrm{m})")
    # lines!(ax3,t,r2p1_fixed[3,:],label="Case 1")
    # lines!(ax3,t,r2p1_freed[3,:],label="Case 2")

    for ilabel in 1:2
        Label(fig.layout[1, ilabel, TopLeft()], "($(alphabet[ilabel]))",
            # fontsize = fontsize,
            font = "CMU Serif Bold",
            padding = (0, 0, 0, 0),
            # halign = :right
        )
    end
    xlims!(ax1,0,0.5)
    xlims!(ax2,0,0.5)
    # xlims!(ax3,0,0.5)
    Legend(fig[1,3],ax1)
    savefig(fig,"simple_energy")
    fig
end

#-- bridge 3d 
n = 2
newbridge = bridge3d(;n)
GM.activate!();with_theme(theme_pub;
        figure_padding = (2fontsize,2fontsize,2fontsize,fontsize),
        Axis3 = (            
            azimuth = 6.005530633326975,
            elevation = 0.13269908169872408
        )
    ) do    
    fig = Figure(
        resolution = (0.9cw,0.42cw)
    )
    subfig1 = fig[1,1] = GridLayout()
    plot_traj!(newbridge;
        AxisType = Axis3,
        fig = subfig1,
        showinfo = true,
        xlims = (-10,10),
        ylims = (-1,25),
        zlims = (-1e-3,12),
        doslide = false,
        showpoints = false,
        showlabels = false,
        sup! = (ax,tgob,sgi)->begin
            # ax.xlabelvisible=false
            # ax.zlabelvisible=false
            ax.azimuth = 1.5*π - 1e-7
            ax.elevation = 1e-7
            hideydecorations!(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(a) ", font=:bold),
                @sprintf "left view"
            )
        end,
    )
    subfig2 = fig[2,1] = GridLayout()
    plot_traj!(newbridge;
        AxisType = Axis3,
        fig = subfig2,
        showinfo = true,
        xlims = (-10,10),
        ylims = (-1,25),
        zlims = (-1e-3,12),
        doslide = false,
        showpoints = false,
        showlabels = false,
        showground = false,
        sup! = (ax,tgob,sgi)->begin
            # ax.titlevisible=false
            # ax.xlabelvisible=false
            # ax.ylabelvisible=false
            ax.azimuth = 2π
            ax.elevation = π/2
            hidezdecorations!(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(b) ", font=:bold),
                @sprintf "top view"
            )
        end,
    )
    subfig3 = fig[1:2,2] = GridLayout()
    plot_traj!(newbridge;
        AxisType = Axis3,
        fig = subfig3,
        showinfo = true,
        xlims = (-10,10),
        ylims = (-1,25),
        zlims = (-1e-3,12),
        doslide = false,
        showpoints = false,
        showlabels = false,
        sup! = (ax,tgob,sgi)->begin
            # ax.titlevisible=false
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(c) ", font=:bold),
                @sprintf "oblique view"
            )
        end,
    )
    colsize!(fig.layout,1,Relative(0.25))
    # rowsize!(fig.layout,1,Relative(0.5))
    rowgap!(fig.layout,0)
    savefig(fig,"bridge")
    fig
end
μ = TR.inverse_for_restlength(newbridge,newbridge;fmin=1e4,gravity = false,eps_rel=1e-10)

Td = 10.0

F̌ = TR.build_F̌(newbridge.tg,2n+1,2*(2n+1)+1,SVector(0.0,0.0,-1.0))

newbridge_modes = deepcopy(newbridge)

TR.set_restlen!(newbridge_modes.tg,μ)
TR.update!(newbridge_modes.tg; gravity=false)
TR.get_cables_tension(newbridge_modes.tg) |> extrema
TR.GDR!(newbridge_modes;β=1e-4, res=1e-9,gravity= true)
plot_traj!(newbridge_modes)
TR.set_new_initial!(newbridge_modes,newbridge_modes.traj.q[end])
TR.check_static_equilibrium_output_multipliers(newbridge_modes.tg;gravity=true)
TR.undamped_eigen(newbridge_modes.tg;gravity=true)
TR.undamped_eigen!(newbridge_modes;scaling=-0.1,gravity=true)

with_theme(theme_pub;
    figure_padding = (0,0,-fontsize,fontsize),
        Axis3 = (
            azimuth = 5.565530633326985,
            elevation = 0.27269908169872414
        )
    ) do 
    plot_traj!(
        newbridge_modes;
        AxisType=Axis3,
        figsize = (1.0cw,0.18cw),
        gridsize=(1,4),
        doslide=false,
        atsteps=collect(2:5),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font = :bold),
                "Mode $sgi: ",
                rich(
                    (@sprintf "%.3G " tt/2π),
                    font = :math
                ),
                "(Hz)"
            )
        end,
        rowgap=0.0,
        colgap=0.0,
        sup! = (ax,tgob,sgi)->begin
            ax.titlevisible = false
            hidexdecorations!(ax)
            hideydecorations!(ax)
            hidezdecorations!(ax)            
        end,
        xlims = (-10,10),
        ylims = (-1,25),
        zlims = (-1e-3,10),
        showinfo=false,
        showground=false,
        showinit=true,
        showpoints=false,
        showlabels=false,
        slack_linestyle=:solid,
        figname = "bridge_modes"
    )
end

@show newbridge_modes.traj.t./(2π)

plot_traj!(newbridge_modes;showinit=true)
newbridge_modes.traj.t

νs = [
    0.453,
    0.553,
    0.653,
]


tend = 10.0
dstep = 45
dt = 1e-2
newbridge_loads = [
    begin
        bot = deepcopy(newbridge_modes)
        TR.solve!(
            TR.SimProblem(
                bot,
                (x)->dynfuncs(x;actuate=false,gravity=true,(Fˣ!)=(F,t) -> begin
                    F .+= 2e4*(sin(ν*2π*t)+1).*F̌
                end)
            ),
            TR.Zhong06();
            dt,
            # tspan=(0.0,dt*(dstep-1)),
            tspan = (0.0,tend),
            ftol=1e-7,
            exception=false
        )
        bot
    end
    for ν in νs
]

plot_traj!(newbridge_loads)

function plot_bridge_vibes(bots,νs,figname=nothing)
    with_theme(theme_pub;
            resolution = (0.4cw,0.15cw),
            figure_padding = (0,0,0,fontsize/2),
            Lines = (
                cycle = Cycle([:linestyle,:color],covary=true),
                # linewidth = 1 |> pt2px,
            )
        ) do 
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=tlabel, ylabel = L"z~(\mathrm{m})")
        for (bot,ν) in zip(bots,νs)
            p11 = get_trajectory!(bot,5,11)
            lines!(ax,bot.traj.t,p11[3,:],label = latexstring("\\nu=$ν"))
            xlims!(ax,bot.traj.t[begin],bot.traj.t[end])
        end
        fig[1, 2] = Legend(fig, ax,)
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plot_bridge_vibes(newbridge_loads,νs)
CM.activate!(); plot_bridge_vibes(newbridge_loads,νs,"bridge_vibes")


#-- embedding
#without outer and gravity, to obtain rest length
gravity = false

m = 3
α = 2π/m
θ = 1.25α
n = 4
b = 0.14
r = 0.04*sqrt(2)

twotre1 = twotre3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1
) 

prism1 = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
    isprism = true,
)

newembed11 = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
    outer = true,
    # isprism = true,
)

newembed1o = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n,
    outer = true,
)

prism0 = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n = 1,
    isprism = true,
)

newembed0 = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n = 1,
    outer = true,
    # isprism = true,
)

plot_traj!(newembed11)
GM.activate!(); with_theme(theme_pub;
        resolution = (0.6cw,0.2cw),
        figure_padding = (0,0,-fontsize/2,0),
        Axis3 = (
            azimuth = 3.907532586451984,
            elevation = 0.18379626575974045
        )
    ) do
    fig = Figure()
    xlims = (-1.8e-1,1.8e-1)
    ylims = (-1.8e-1,1.8e-1)
    zlims = (-1e-3,2.4e-1)
    doslide = false
    showpoints = false
    showlabels = false
    g0 = fig[1,1] = GridLayout(;tellheight=false,)
    g1 = fig[1,2] = GridLayout(;tellheight=false,)
    # g11 = g1[1,1] = GridLayout()
    # g12 = g1[2,1] = GridLayout()
    g2 = fig[1,3] = GridLayout(;tellheight=false,)
    # g21 = g2[1,1] = GridLayout()
    # g22 = g2[2,1] = GridLayout()
    g3 = fig[1,4] = GridLayout()
    plot_traj!(
        twotre1;
        AxisType = Axis3,
        fig = g0,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(a) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
    plot_traj!(
        prism1;
        AxisType = Axis3,
        fig = g1,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(b) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
        # plot_traj!(
        #     prism0;
        #     AxisType = Axis3,
        #     fig = g12,
        #     xlims,
        #     ylims,
        #     zlims = (-1e-3,1e-1),
        #     doslide,
        #     showpoints,
        #     showlabels,
        #     titleformatfunc = (sgi,tt)-> begin
        #         rich(
        #             rich("(d) ", font=:bold),
        #         )
        #     end,
        #     sup! = (ax,tgob,sgi) -> begin
        #         hidexyz(ax)
        #     end,
        # )
    plot_traj!(newembed11;
        AxisType = Axis3,
        fig = g2,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(c) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
        # plot_traj!(newembed0;
        #     AxisType = Axis3,
        #     fig = g22,
        #     xlims,
        #     ylims,
        #     zlims = (-1e-3,1e-1),
        #     doslide,
        #     showpoints,
        #     showlabels,
        #     titleformatfunc = (sgi,tt)-> begin
        #         rich(
        #             rich("(e) ", font=:bold),
        #         )
        #     end,
        #     sup! = (ax,tgob,sgi) -> begin
        #         hidexyz(ax)
        #     end,
        # )
    plot_traj!(newembed1o;
        AxisType = Axis3,
        fig = g3,
        xlims,
        ylims,
        zlims = (-1e-3,8.5e-1),
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(d) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
    colgap!(fig.layout,0)
    rowsize!(g0,1,Fixed(0.10cw))
    rowsize!(g1,1,Fixed(0.10cw))
    rowsize!(g2,1,Fixed(0.10cw))
    # rowgap!(g1,-fontsize)
    # rowsize!(g1,2,Relative(0.35))
    # rowgap!(g2,-fontsize)
    # rowsize!(g2,2,Relative(0.35))
    savefig(fig,"embedding")
    fig
end
newembed1 = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n,
    outer = false,
)
plot_traj!(newembed1)
μ1 = TR.inverse_for_restlength(
    newembed1,newembed1;gravity,fmin=100.0
)
barplot(μ1)
# (optional) check μ1 
TR.set_restlen!(newembed1.tg,μ1)
TR.update!(newembed1.tg)
f = TR.get_cables_tension(newembed1.tg) 
f |> extrema
TR.check_static_equilibrium_output_multipliers(newembed1.tg;gravity)
TR.undamped_eigen(newembed1.tg;gravity)
TR.undamped_eigen!(newembed1;gravity)
TR.reset!(newembed1)
#-- without outer, see see gravity
TR.GDR!(newembed1;gravity=true)

plot_traj!(newembed1)

#-- folded without outer and gravity, to obtain rest length
newembed0 = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n
)

plot_traj!(newembed0)

μ0 = TR.inverse_for_restlength(
    newembed0,newembed0;gravity,fmin=100.0
)
barplot(μ0)
#-- with outer, using the folded and unfolded 
outer = true
gravity = true
newembed0_outer = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n,outer = true
)
plot_traj!(newembed0_outer)
μ0o = TR.get_cables_restlen(newembed0_outer)
μ0o[begin:3n*m] .= μ0
μ1o = deepcopy(μ0o)
μ1o[begin:3n*m] .= μ1
TR.set_restlen!(newembed0_outer.tg,μ0o)
TR.GDR!(newembed0_outer;gravity)
GM.activate!();plot_traj!(newembed0_outer;)
TR.set_new_initial!(newembed0_outer,newembed0_outer.traj.q[end])
TR.update!(newembed0_outer.tg;gravity)
TR.check_static_equilibrium_output_multipliers(newembed0_outer.tg;gravity)
TR.undamped_eigen(newembed0_outer.tg;gravity)
TR.undamped_eigen!(newembed0_outer;gravity)
μs = [μ0o,μ1o]

barplot(μ0o)
barplot!(μ1o)
pretty_table(
    Dict(
        [
            ("Folded, Cables 1~6", μ0o[1]),
            ("Deployed, Cables 1~6", μ0o[7]),
            ("Folded, Cables 1~6", μ1o[1]),
            ("Deployed, Cables 1~6", μ1o[7]),
            ("Outer Cables 1~6", μ1o[end-2:end]),
        ]
    )
)

Td = 8.0
tend = 10.0

μs = [
    begin
        μ = deepcopy(μ0o)
        μ[begin:3m*j] .= μ1[begin:3m*j]
        μ
    end
    for j = 0:n
]

function make_multi_stages_pres_actor(μs;start = 0.0, stop = 10.0, len = 2, )
    nμ = length(μs[begin])

    function itp(t)
        scaled_itps = extrapolate(
            Interpolations.scale(
                interpolate(
                    reduce(hcat,μs),
                    (NoInterp(),BSpline(Linear()))
                    # (NoInterp(),BSpline(Quadratic(Flat(OnGrid()))))
                ),
                1:nμ, range(;start, stop, length = len,)
            ),
            (Throw(),Flat())
        )
        [scaled_itps(j,t) for j in 1:nμ]
    end

    TR.PrescribedActuator(
        1,
        TR.ManualActuator(1,collect(1:nμ),zeros(nμ),TR.Uncoupled()),
        itp
    )
end
make_multi_stages_pres_actor(μs;start=0.0,stop=Td,len=n+1)

newembed0_outer_deploy = TR.TensegrityRobot(
    deepcopy(newembed0_outer.tg),
    (actuators=[make_multi_stages_pres_actor(μs;start=0.0,stop=Td,len=length(μs))],)
)

TR.solve!(
    TR.SimProblem(
        newembed0_outer_deploy,
        (x)->dynfuncs(x;actuate=true,gravity=true)
    ),
    TR.Zhong06();
    dt=1e-3,
    tspan=(0.0,tend),
    ftol=1e-7,
    exception=true
)

plot_traj!(newembed0_outer_deploy)
GM.activate!();with_theme(theme_pub;
        figure_padding = (0,fontsize,-fontsize/2,fontsize/2),
        # resolution = (0.9cw,0.25cw),
        Axis3 = (
            azimuth = 4.995530633326985,
            elevation = 0.18269908169872415
        )
    ) do
    plot_traj!(
        newembed0_outer_deploy;
        AxisType=Axis3,
        figsize = (0.7cw,0.25cw),
        gridsize = (1,5),
        attimes = range(0,Td,5),
        xlims = (-2e-1,2e-1),
        ylims = (-4e-1,2e-1),
        zlims = (-1e-3,8e-1),
        showinfo = false,
        doslide = false,
        # dorecord = true,
        slack_linestyle = :solid,
        actuate = true,
        showpoints = false,
        showlabels = false,
        sup! = (ax,tgob,sgi) -> begin
            hidex(ax)
            hidey(ax)
            if sgi != 5
                hidez(ax)
            end
        end,
        # rowgap=2fontsize,
        colgap=0,
        figname = "newembed0_outer_deploy"
    )
end
nb = newembed0_outer_deploy.tg.nbodies
ps = [get_trajectory!(newembed0_outer_deploy,nb,j) for j = m+1:m+m]
GM.activate!();with_theme(theme_pub;
        resolution = (0.7cw,0.2cw),
        figure_padding = (0,fontsize,0,0)
    ) do 
    fig = Figure()
    (;t) = newembed0_outer_deploy.traj
    for i = 1:m
        ax = Axis(fig[1,i], 
            xlabel = tlabel, 
            ylabel = latexstring("""$(["x","y","z"][i])~(\\mathrm{m})""")
        )
        lines!(ax,t,ps[2][i,:])
        # lines!(ax,t,ps[2][i,:])
        # lines!(ax,t,ps[3][i,:])
        xlims!(t[begin],t[end])
        Label(
            fig[1,i,TopLeft()],
            "($(alphabet[i]))",
            font = "CMU Serif Bold",
            padding = (0, 0, 0, 0),
            halign = :right,
            valign = :bottom
        )
    end
    savefig(fig,"newembed0_outer_deploy_traj")
    fig
end


## tower3d
## Modal analysis
# tower3dbot = tower3d(;k=500.0,c=10.0,θ=π/6)

## Seismic and deploy
function pres3d!(sysstate,tg;ν=0.5)
    (;t,q̃,q̃̇,q̃̈) = sysstate
    (;bodies,connectivity) = tg
    (;mem2syspres) = connectivity.indexed
    amp = 0.01
    p = ν*2π
    foreach(bodies) do rb
        rbid = rb.prop.id
        if rbid in [1,2,3]
            r0 = zeros(3)
            if rbid == 1
                r0 .= [
                    0.1,
                    0.0,
                    0.0
                ]
                q̃[mem2syspres[rbid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[mem2syspres[rbid]] .=     [amp*p*cos(p*t),0,0]
                q̃̈[mem2syspres[rbid]] .=  [-amp*p^2*sin(p*t),0,0]
            elseif rbid == 2
                r0 .= [
                    -0.05,
                     0.08660254037844387,
                    0.0
                ]
                q̃[mem2syspres[rbid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[mem2syspres[rbid]] .=     [amp*p*cos(p*t),0,0]
                q̃̈[mem2syspres[rbid]] .=  [-amp*p^2*sin(p*t),0,0]
            elseif rbid == 3
                r0 .= [
                    -0.05,
                    -0.08660254037844387,
                    0.0
                ]
                q̃[mem2syspres[rbid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[mem2syspres[rbid]] .=     [amp*p*cos(p*t),0,0]
                q̃̈[mem2syspres[rbid]] .=  [-amp*p^2*sin(p*t),0,0]
            end
        end
    end
end

rest = deleteat!(collect(1:24),13:21)

## inverse statics and modal analysis
tower3dbot0 = tower3d(;k=500.0,k1=1000,c=2.0,α=π/12)
# B,F̃ = TR.build_inverse_statics_for_restlength(tower3dbot0.tg,tower3dbot0.tg;gravity=true)
# rank(B),size(B)
# x0,xn = TR.get_solution_set(B,F̃)
# size(xn)

plot_traj!(tower3dbot0;fontsize=20)
μ0 = TR.inverse_for_restlength(tower3dbot0,tower3dbot0;gravity=true,scale=true,fmin=1.0)
TR.set_restlen!(tower3dbot0.tg,μ0)
TR.update!(tower3dbot0.tg)
TR.get_cables_tension(tower3dbot0.tg) |> extrema
TR.check_static_equilibrium_output_multipliers(tower3dbot0.tg;gravity=true)
TR.undamped_eigen(tower3dbot0.tg;gravity=true)
TR.undamped_eigen!(tower3dbot0;gravity=true)
tower3dbot0.traj.t./(2π)
@show (tower3dbot0.traj.t./(2π))[2]

tower3dbot1 = tower3d(;k=500.0,c=2.0,d = 0.1*√2, r2 = 0.07, α=-π/12)
plot_traj!(tower3dbot1;fontsize=20)
μ1 = TR.inverse_for_restlength(tower3dbot1,tower3dbot1;gravity=true,fmin=1.0)
TR.set_restlen!(tower3dbot1.tg,μ1)
TR.update!(tower3dbot1.tg)
TR.get_cables_tension(tower3dbot1.tg) |> extrema
TR.check_static_equilibrium_output_multipliers(tower3dbot1.tg;gravity=true)
# TR.undamped_eigen(tower3dbot1.tg;gravity=true)
TR.undamped_eigen!(tower3dbot1;gravity=true,scaling=0.2)
tower3dbot1.traj.t./(2π)
@show (tower3dbot1.traj.t./(2π))[2]

twotre1 = twotre3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
    loadmesh = false
) 

spine1 = twotre3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
    loadmesh = false,
    isspine = true,
) 

GM.activate!();with_theme(theme_pub;
        figure_padding = (0,0,-fontsize/2,0),
        resolution = (0.65cw,0.25cw),
        Axis3 = (
            azimuth = 7.20553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    fig = Figure()
    griddecompose_row = fig[1,1] = GridLayout(;tellheight=true)
    # griddecompose_row = griddecompose[1,1] = GridLayout(;tellheight=true)
    griddecompose_twotre = griddecompose_row[1,1] = GridLayout(;tellheight=true)
    griddecompose_spine = griddecompose_row[2,1] = GridLayout(;tellheight=true)
    griddecompose_prism = griddecompose_row[:,2] = GridLayout(;tellheight=false)
    gridfolded = fig[1,2] = GridLayout()
    griddeployed = fig[1,3] = GridLayout()
    xlims=(-0.15,0.15)
    ylims=(-0.15,0.15)
    zlims=(-1e-3,0.5)
    doslide=false
    showpoints=false
    showlabels=false
    showground = false
    plot_traj!(
        twotre1;
        AxisType = Axis3,
        doslide = false,
        xlims,
        ylims,
        zlims=(-1e-3,0.2),
        showpoints,
        showlabels,
        showground,
        fig = griddecompose_twotre,
        sup! = (ax,tgob,sgi)->begin
            hidedecorations!(ax)
            hidespines!(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[1])) ", font=:bold),
            )
        end,
    )
    plot_traj!(
        spine1;
        AxisType = Axis3,
        doslide = false,
        xlims,
        ylims,
        zlims=(-1e-3,0.2),
        showpoints,
        showlabels,
        showground,        
        slack_linestyle = :solid,
        fig = griddecompose_spine,
        sup! = (ax,tgob,sgi)->begin
            hidedecorations!(ax)
            hidespines!(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[2])) ", font=:bold),
            )
        end,
    )
    plot_traj!(
        prism1;
        AxisType = Axis3,
        doslide = false,
        xlims,
        ylims,
        zlims=(-1e-3,0.2),
        showpoints,
        showlabels,
        showground,
        fig = griddecompose_prism,
        sup! = (ax,tgob,sgi)->begin
            hidedecorations!(ax)
            hidespines!(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[3])) ", font=:bold),
            )
        end,
    )
    plot_traj!(
        tower3dbot0;
        AxisType=Axis3,
        fig = gridfolded,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        sup! = (ax,tgob,sgi)->begin
            hidexyz(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[4])) ", font = :bold),
            )
        end,
    )
    plot_traj!(
        tower3dbot1;
        AxisType=Axis3,
        fig = griddeployed,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        sup! = (ax,tgob,sgi)->begin
            hidexyz(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[5])) ", font = :bold),
            )
        end,
    )
    colgap!(griddecompose_row,-4fontsize)
    colsize!(griddecompose_row,1,Fixed(0.12cw))
    colsize!(griddecompose_row,2,Fixed(0.20cw))
    # colsize!(griddecompose_row,2,Fixed(0.20cw))
    rowgap!(griddecompose_row,0)
    colgap!(fig.layout,0)
    colgap!(fig.layout,1,-3fontsize)
    # colsize!(fig.layout,1,Relative(0.25))
    savefig(fig,"tower3d")
    fig
end

WM.activate!();with_theme(theme_pub;
        figure_padding = (0,0,-fontsize,0),
        Axis3 = (
            azimuth = 7.00553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    plot_traj!(
        tower3dbot1;
        AxisType=Axis3,
        figsize = (0.55cw,0.20cw),
        gridsize = (1,3),
        atsteps = 2:4,
        xlims=(-0.2,0.2),
        ylims=(-0.2,0.2),
        zlims=(-1e-3,0.5),
        showmesh=false,
        doslide=false,
        showpoints=false,
        showlabels=false,
        showinfo=false,
        showinit=true,
        dorecord=false,
        sup! = (ax,tgob,sgi)->begin
            hidexyz(ax)
        end,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font = :bold),
                "Mode $sgi",
            )
        end,
        colgap = 0,
        figname="tower3d_modes",
    )
end

@fmt Real = "{:}"
initial_μ = TableCol("Initial", ["$((3i-2,3i-1,3i))" for i in 1:8], ["\\qty{$μ}{\\meter}" for μ in μ0[begin:3:end]])
target_μ = TableCol("Target", ["$((3i-2,3i-1,3i))" for i in 1:8], ["\\qty{$μ}{\\meter}" for μ in μ1[begin:3:end]])
table_μ = hcat(initial_μ,target_μ)
print(TexTables.to_tex(table_μ))

function plot_tensions(bots,figname=nothing)
    fs = [TR.get_cables_tension(bot) for bot in bots]
    # cg = cgrad(:Dark2_6, 6, categorical = true)[[4,3,6]]
    # set_theme!(theme_pub;
    #         palette = (color = cg, ),
    #         Lines = (cycle = [:color],),
    # )
    with_theme(theme_pub;
            # resolution = (0.7tw,0.3tw),
            # font = "Nimbus Rom No9 L",
            # font = "Times New Roman",
            resolution = (2cw,0.7cw),
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1]; 
            xlabel = "Cable No.", 
            ylabel = "Tension (N)"
        )
        x = vec([collect(1:24)';collect(1:24)'])
        grp = repeat([1,2],24)
        display(x)
        height = vec([fs[1]';fs[2]'])
        colors = Makie.wong_colors()[[3,4]]
        barplot!(
            ax1, x, height; 
            dodge = grp,
            color = colors[grp],
            marker = '+', markersize = 1.5fontsize, 
            linewidth, linestyle = :dashdot, 
            # label = "Folded"
        )
        # ax1.yticks = collect(0:15)
        hlines!(ax1,1.0)
        ax1.xticks = collect(1:24)
        # ax1.xlabelpadding = -15
        # axislegend(ax1; position = :lt,)
        # Legend
        labels = ["Folded", "Deployed"]
        elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
        # title = "Groups"

        Legend(fig[0,1], elements, labels,
            orientation=:horizontal,
            tellwidth=false,tellheight=false
        )
        # fig[1,2] = Legend(fig,ax1)
        rowsize!(fig.layout,1,Fixed(0.5cw))
        # rowgap!(fig.layout,0)
        # colsize!(fig.layout,1,Fixed(0.6cw))
        savefig(fig,figname)
    end
end

GM.activate!(); plot_tensions([tower3dbot0,tower3dbot1])

function plot_first_frequencies(bots,figname=nothing)
    n = 5
    frqs = [bot.traj.t[begin+1:end]./(2π) for bot in bots]
    @show frqs[1][1], frqs[2][1]
    unqfrqs = [(round.(frq;sigdigits=4)|> unique)[1:n] for frq in frqs]
    # cg = cgrad(:Dark2_6, 6, categorical = true)[[4,3,6]]
    # set_theme!(theme_pub;
    #         palette = (color = cg, ),
    #         Lines = (cycle = [:color],),
    # )
    with_theme(theme_pub;
            # resolution = (0.7tw,0.3tw),
            # font = "Nimbus Rom No9 L",
            # font = "Times New Roman",
            resolution = (cw,0.5cw),
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1]; 
            xlabel = L"\mathrm{Mode~(Exclude~duplicates)}", 
            ylabel = L"\mathrm{Frequency}~(\mathrm{Hz})"
        )
        x = vec([collect(1:n)';collect(1:n)'])
        grp = repeat([1,2],n)
        display(x)
        height = vec([unqfrqs[1]';unqfrqs[2]'])
        colors = Makie.wong_colors()[[5,6]]
        barplot!(
            ax1, x, height; 
            dodge = grp,
            color = colors[grp],
        )
        ax1.yticks = collect(0:15)
        ax1.xticks = collect(1:n)
        # ax1.xlabelpadding = -15
        # axislegend(ax1; position = :lt,)
        # fig[1,2] = Legend(fig,ax1)
        labels = ["Folded", "Deployed"]
        elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
        # title = "Groups"
        Legend(fig[1,2], elements, labels,
            # orientation=:horizontal,
            tellwidth=false,tellheight=false
        )
        colsize!(fig.layout,1,Fixed(0.6cw))
        savefig(fig,figname)
    end
end

GM.activate!(); plot_first_frequencies([tower3dbot0,tower3dbot1])
CM.activate!(); plot_first_frequencies([tower3dbot0,tower3dbot1],"first_frequencies")

## seismic without deploy
tower3dbot0_nodpl = deepcopy(tower3dbot0)
TR.solve!(TR.SimProblem(tower3dbot0_nodpl,(x)->dynfuncs(x;gravity=true)),
            TR.Zhong06(),
            (
                prescribe! = pres3d!,
                actuate! = nothing
            );
            dt=1e-3,tspan=(0.0,15.0),verbose=false,ftol=1e-10,exception=false)
# plot_tower3d_traj(tower3dbot0_nodpl)
GM.activate!();
with_theme(theme_pub;
        Axis3 = (
            azimuth = 7.00553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    plot_traj!(
        tower3dbot0_nodpl;
        AxisType=Axis3,
        xlims=(-0.2,0.2),
        ylims=(-0.2,0.2),
        zlims=(-1e-3,0.5),
        showmesh=false,
        showpoints=false,
        showlabels=false,
        showinfo=false,
        dorecord=true,
        figname="tower3dbot0_nodpl.mp4",
    )
end

tower3dbot1_nodpl = deepcopy(tower3dbot1)
TR.solve!(TR.SimProblem(tower3dbot1_nodpl,(x)->dynfuncs(x;gravity=true)),
            TR.Zhong06(),
            (
                prescribe! = pres3d!,
                actuate! = nothing
            );
            dt=1e-3,tspan=(0.0,15.0),ftol=1e-14)
# plot_tower3d_traj(tower3dbot1_nodpl)
plot_traj!(tower3dbot1_nodpl;showmesh=false)

with_theme(theme_pub;
        Axis3 = (
            azimuth = 7.00553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    plot_traj!(
        tower3dbot1_nodpl;
        AxisType=Axis3,
        xlims=(-0.2,0.2),
        ylims=(-0.2,0.2),
        zlims=(-1e-3,0.5),
        showmesh=false,
        showpoints=false,
        showlabels=false,
        showinfo=false,
        dorecord=true,
        figname="tower3dbot1_nodpl.mp4",
    )
end

function plot_tower3d_seis_traj(bots,figname=nothing;
        ss = 20,
    )
    cg = cgrad(:seaborn_bright6, 6, categorical = true)
    with_theme(theme_pub;
            # resolution=(0.9tw,0.4tw),
            resolution=(tw,0.7cw),
            figure_padding = (fontsize,fontsize,0,0),
            palette = (color = cg, ),
            Lines = (cycle = [:color],linewidth=2),
        ) do
        fig = Figure(;)
        for i = 1:2
            bot = bots[i]
            t = bot.traj.t[begin:ss:end]
            ax1 = Makie.Axis(fig[i,1];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta x~(\mathrm{m})")
            ax2 = Makie.Axis(fig[i,2];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta y~(\mathrm{m})")
            ax3 = Makie.Axis(fig[i,3];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta z~(\mathrm{m})")
            scalings = [-2,-3,-3]
            rb7r4 = get_trajectory!(bot,7,4,1:ss:length(bot.traj))
            rb8r4 = get_trajectory!(bot,8,4,1:ss:length(bot.traj))
            rb9r3 = get_trajectory!(bot,9,3,1:ss:length(bot.traj))
            # rb9r3 = get_trajectory!(bot,10,3,1:ss:length(bot.traj))
            lines!(ax1,t,(rb7r4[1,:].-rb7r4[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{7,1}")
            lines!(ax1,t,(rb8r4[1,:].-rb8r4[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{8,1}")
            lines!(ax1,t,(rb9r3[1,:].-rb9r3[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{9,1}")
            # lines!(ax1,t,(rb9r3[1,:].-rb9r3[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{10,1}")
            ylims!(ax1,-5,5)
            ax1.yticks = collect(-5:2.5:5)

            lines!(ax2,t,(rb7r4[2,:].-rb7r4[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb8r4[2,:].-rb8r4[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb9r3[2,:].-rb9r3[2,1])./10.0^scalings[2];)
            # lines!(ax2,t,(rb9r3[2,:].-rb9r3[2,1])./10.0^scalings[2];)

            lines!(ax3,t,(rb7r4[3,:].-rb7r4[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb8r4[3,:].-rb8r4[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb9r3[3,:].-rb9r3[3,1])./10.0^scalings[3];)
            # lines!(ax3,t,(rb9r3[3,:].-rb9r3[3,1])./10.0^scalings[3];)
            # @show rb9r3[3,1]
            # ylims!(ax3,1,5)

            function set_ax!(ax)
                # xlims!(ax,0,5)
                # ax.xlabelpadding = -15
                # ax.alignmode =  Mixed(;left = -20, right = 20)
                if i == 1
                    hidex(ax)
                else

                    axislegend(ax1;position=:rb,orientation=:horizontal,)
                end
                xlims!(ax,t[begin],t[end])
                # xlims!(ax,0,10)
            end
            foreach(set_ax!,[ax1,ax2,ax3])
            for (ilabel,label) in enumerate(scalings)
                Label(fig.layout[i, ilabel, Top()], latexstring("×10^{$label}"),
                    # fontsize = fontsize,
                    # # font = "Nimbus Rom No9 L",
                    # padding = (ifelse(ilabel==1,90,55), 0, 0, 0),
                    # halign = :left
                    )
            end
            for (ilabel,label) in enumerate(alphabet[3(i-1)+1:3(i-1)+3])
                Label(fig.layout[i, ilabel, TopLeft()], "($label)",
                    # fontsize = fontsize,
                    font = "CMU Serif Bold",
                    padding = (0, 0, fontsize/2, 0),
                    # halign = :right
                )
            end
        end
        # colgap!(fig.layout,0)
        # rowgap!(fig.layout,0)
        savefig(fig,figname)
    end
end

GM.activate!(); plot_tower3d_seis_traj([tower3dbot0_nodpl,tower3dbot1_nodpl])
CM.activate!(); plot_tower3d_seis_traj([tower3dbot0_nodpl,tower3dbot1_nodpl],"tower3d_seis_traj")

function plot_tower3d_vis_nodpl(bots,figname=nothing)
    fig_width = 1cw
    fig_height = 0.8cw
    with_theme(theme_pub;
        resolution=(fig_width,fig_height),
        figure_padding=(10,0,0,0),
        Axis3 = (
            azimuth = 7.525530633326982,
            elevation = 0.12269908169872408
        )
    ) do
        fig = Figure()
        dt = 1e-3
        ts = [
            2.5,
            3.5
        ]
        nbot = length(bots)
        for i in 1:nbot
            gridi = fig[1,i] = GridLayout()
            boti = bots[i]
            stepj = round(Int,ts[i]/dt)
            TR.goto_step!(boti,1)
            plot_traj!(boti;
                AxisType=Axis3,
                fig=gridi,
                # cablecolor=:slategrey,
                atsteps=[stepj],
                showinit=true,
                tgini=deepcopy(boti.tg),
                showinfo=false,
                doslide=false,
                showlabels=false,
                showpoints=false,
                showground=false,
                xlims = (-0.125,0.125),
                ylims = (-0.125,0.125),
                zlims = (-0.001,0.5),
                titleformatfunc = (sgi,tt) -> "",
                sup! = (ax,tgob,sgi) -> begin
                    rbs = TR.get_bodies(tgob[])
                    for rbid in 7:10
                        # for rbid in [10]
                        if rbid in 7:8
                            pid = 4
                        else
                            pid = 3
                        end
                        meshscatter!(ax,[rbs[rbid].state.rps[pid]], 
                            color = :darkblue, 
                            markersize = 0.008)
                    end                    
                    bjs = Dict(
                        1 => [1],
                        2 => [1],
                        3 => [1],
                        4 => [1,2],
                        5 => [1,2],
                        6 => [1,2],
                        8 => [1]
                    )
                    foreach(bjs) do (key,values)
                        foreach(values) do value
                            meshscatter!(ax,[rbs[key].state.rps[value]], 
                                color = :lightblue, 
                                markersize = 0.008)
                        end
                    end
                    ax.xlabel = "x (m)"
                    ax.ylabel = "y (m)"
                    ax.zlabel = "z (m)"
                    ax.yticks = ax.xticks = [-0.1,0,0.1]
                    ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
                    ax.ylabeloffset = ax.xlabeloffset = 40
                    ax.alignmode =  Mixed(;left = 25, right = -25)               
                end
            )
        end
        titles = ["Folded configuration","Deployed configuration"]
        for (ilabel,label) in enumerate(alphabet[1:nbot])
            Label(fig.layout[1, ilabel, Bottom()],
                rich(
                    rich("($label) ", font = "CMU Serif Bold"),
                    "$(titles[ilabel])",
                ),
                padding = (70, 0, 0, 80),
                halign = :left,
                valign = :top)
        end
        # for (ilabel,label) in enumerate(alphabet[1:nbot])
        #     Label(fig.layout[1, ilabel, Bottom()],
        #         latexstring("t= $(ts[ilabel])(s)"),
        #         padding = (
        #             ifelse(
        #                 ilabel==1,
        #                 9.5fontsize, 
        #                 10.5fontsize, 
        #             ),
        #             0, 0, 83
        #         ),
        #         halign = :left,
        #         valign = :top)
        # end
        colgap!(fig.layout,0)
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plot_tower3d_vis_nodpl([tower3dbot0_nodpl,tower3dbot1_nodpl])
GM.activate!(); plot_tower3d_vis_nodpl([tower3dbot0_nodpl,tower3dbot1_nodpl], "tower3d_vis_nodpl")

## deploy
Td = [10.0]
function do_noseis(;Td,tend=15.0)
    bot   = TR.TensegrityRobot(deepcopy(tower3dbot0.tg),(actuators=[make_pres_actor(μ0,μ1,0.0,Td)],))
    TR.solve!(TR.SimProblem(bot,(x)->dynfuncs(x;actuate=true,gravity=true)),
                TR.Zhong06();
                dt=1e-3,tspan=(0.0,tend),ftol=1e-14)
    bot
end
function do_seis(;ν,Td,tend=15.0)
    bot   = TR.TensegrityRobot(deepcopy(tower3dbot0.tg),(actuators=[make_pres_actor(μ0,μ1,0.0,Td)],))
    TR.solve!(TR.SimProblem(bot,(x)->dynfuncs(x;actuate=true,gravity=true)),
                TR.Zhong06(),
                (
                    prescribe! = (sysstate,tg)->pres3d!(sysstate,tg;ν),
                    actuate! = nothing
                );
                dt=1e-3,tspan=(0.0,tend),ftol=1e-14)
    bot
end
tower3d_noseis_10 = do_noseis(;Td=Td[1])
tower3d_seis_10 = do_seis(;ν=1.0,Td=Td[1])
tower3d_seis_lo_10 = do_seis(;ν=1.5,Td=Td[1])
tower3d_seis_hi_10 = do_seis(;ν=0.5,Td=Td[1])

plot_traj!(tower3d_seis_10;actuate=true)

GM.activate!();with_theme(theme_pub;
        figure_padding = (0,0,-fontsize,0),
        Axis3 = (
            azimuth = 7.80553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    plot_traj!(
        tower3d_seis_10;
        AxisType=Axis3,
        figsize = (0.9cw,0.2cw),
        gridsize = (1,5),
        attimes = range(1,10,5),
        xlims=(-0.3,0.2),
        ylims=(-0.2,0.2),
        zlims=(-1e-3,0.5),
        # showmesh=false,
        showpoints=false,
        showlabels=false,
        showinfo=false,
        doslide=false,
        # dorecord=true,
        actuate=true,
        figname="tower3d_seis_10",
        sup! = (ax,tgob,sgi)->begin
            hidexyz(ax)
        end,
    )
end

plot_traj!(tower3d_seis_10;actuate=true)
plot_traj!(tower3d_noseis_15;actuate=true)
plot_traj!(tower3d_seis_hi;actuate=true)
me_series = TR.mechanical_energy!(tower3d_seis;actuate=true,gravity=true)
lines(tower3d_seis.traj.t,me_series.E)

function plot_tower3d_dpl_traj(bots,figname=nothing;
        ss=10,
    )
    cg = cgrad(:seaborn_bright6, 6, categorical = true)
    with_theme(theme_pub;
            # resolution = (cw,0.8cw),
            resolution = (0.8cw,0.20cw),
            # palette = (color = cg, ),
            figure_padding = (0,0,0,0),
            Lines = (cycle = Cycle([:color,]), linewidth = 4),
        ) do
        nrow = length(bots)
        fig = Figure(;)
        ats = [10.5,14.0,18.0]
        Tds = [10.0,15.0,20.0]
        axs = [
            [
                Makie.Axis(fig[i,1];xlabel=L"t~(\mathrm{s})",ylabel=L"x~(\mathrm{m})"),
                Makie.Axis(fig[i,2];xlabel=L"t~(\mathrm{s})",ylabel=L"z~(\mathrm{m})"),
                # Makie.Axis(fig[i,3];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta z~(\mathrm{m})")
            ]
            for i = 1:nrow
        ]
        scalings = [-1,-1]

        for i = 1:nrow
            bot0, bot1, botlo, bothi = bots[i]
            at = ats[i]
            Td = Tds[i]
            ax1,ax2 = axs[i]
            t0 = bot0.traj.t[begin:ss:end]
            step_collapse = findfirst((x)->x>=at,bot1.traj.t)
            step_Td = findfirst((x)->x>=Td,bot1.traj.t)
            step_Tdp1 = findfirst((x)->x>=Td+1,bot1.traj.t)
            tsteps = 1:ss:step_collapse
            t1 = bot1.traj.t[tsteps]
            # bot0
            bot0rb9r3 = get_trajectory!(bot0,9,3,1:ss:length(bot0.traj))
            # bot1
            bot1rb9r3 = get_trajectory!(bot1,9,3,tsteps)
            # botlo
            botlorb9r3 = get_trajectory!(botlo,9,3,1:ss:length(botlo.traj))
            # bothi
            bothirb9r3 = get_trajectory!(bothi,9,3,1:ss:length(bothi.traj))

            # x
            lines!(ax1,t0,bot0rb9r3[1,:]./10.0^scalings[1];)
            lines!(ax1,t0,botlorb9r3[1,:]./10.0^scalings[1];)
            lines!(ax1,t1,bot1rb9r3[1,:]./10.0^scalings[1];)
            lines!(ax1,t0,bothirb9r3[1,:]./10.0^scalings[1];)
            ylims!(ax1,-1,2.0)
            # ax1.yticks = collect(-5:2.5:5)
            xlims!(ax1,t0[begin],t0[end])

            # # y
            # lines!(ax2,t0,bot0rb9r3[2,:]./10.0^scalings[2];label=L"\mathrm{Static~ground}")
            # lines!(ax2,t0,botlorb9r3[2,:]./10.0^scalings[2];label=L"\mathrm{Seismic~ground},~\nu = 0.2~(\mathrm{Hz})")
            # lines!(ax2,t0,bothirb9r3[2,:]./10.0^scalings[2];label=L"\mathrm{Seismic~ground},~\nu = 1.0~(\mathrm{Hz})")
            # lines!(ax2,t1,bot1rb9r3[2,:]./10.0^scalings[2];label=L"\mathrm{Seismic~ground},~\nu = 0.5~(\mathrm{Hz})")
            # # ylims!(ax2,-0.5,7.5)
            # xlims!(ax2,t0[begin],t0[end])
            # ax2.yticks = collect(-1.0:1.0:7.5)

            # z
            lines!(ax2,t0,bot0rb9r3[3,:]./10.0^scalings[2];label=L"\mathrm{Static~ground}")
            lines!(ax2,t0,botlorb9r3[3,:]./10.0^scalings[2];label=L"\mathrm{Seismic~ground},~\nu = 0.2~(\mathrm{Hz})")
            lines!(ax2,t1,bot1rb9r3[3,:]./10.0^scalings[2];label=L"\mathrm{Seismic~ground},~\nu = 0.5~(\mathrm{Hz})")
            lines!(ax2,t0,bothirb9r3[3,:]./10.0^scalings[2];label=L"\mathrm{Seismic~ground},~\nu = 1.0~(\mathrm{Hz})")
            # ylims!(ax2,3.5,5.5)
            xlims!(ax2,t0[begin],t0[end])

            # rb9r3_z0 = 0.49497474683058335
            # lines!(ax3,t0,(bot0rb9r3[3,:].-rb9r3_z0)./10.0^scalings[3];)
            # lines!(ax3,t0,(botlorb9r3[3,:].-rb9r3_z0)./10.0^scalings[3];)
            # lines!(ax3,t1,(bot1rb9r3[3,:].-rb9r3_z0)./10.0^scalings[3];)
            # lines!(ax3,t0,(bothirb9r3[3,:].-rb9r3_z0)./10.0^scalings[3];)
            # if i == 1
            #     ylims!(ax3,-1.5,1.5)
            # else
            #     ylims!(ax3,-.8,1.0)
            # end
            # xlims!(ax3,Td,Td+2)

            # if i == 2
            fig[:,3] = Legend(fig,ax2;orientation=:vertical,)
            # end labelsize=pt2px(6.5)
            function set_ax!(ax)
                # xlims!(ax,0,5)
                # ax.alignmode =  Mixed(;left = -20, right = 20)
                vlines!(ax,[Td];linewidth,color=:slategrey)
                # ax.xticks = collect(0:1:15)
                ax.xminorticks = collect(0:1:25)
                ax.xminorgridvisible = true
                if i != nrow
                    hidex(ax)
                end
            end
            foreach(set_ax!,axs[i])
            leftpads = [90,80,80]
            for (ilabel,label) in enumerate(scalings)
                Label(fig.layout[i, ilabel, Top()], latexstring("×10^{$label}"),)
            end
            for (ilabel,label) in enumerate(alphabet[[2i-1,2i]])
                Label(fig.layout[i, ilabel, TopLeft()], "($label)",
                    fontsize = fontsize,
                    font = "CMU Serif Bold",
                    padding = (0, 0, fontsize/2, 0),
                    # halign = :right
                )
            end
        end
        # colgap!(fig.layout,0)
        # rowgap!(fig.layout,0)
        # colsize!(fig.layout,3,Relative(0.15))
        savefig(fig,figname)
    end
end

GM.activate!(); plot_tower3d_dpl_traj(
    [
        [tower3d_noseis_10,tower3d_seis_10,tower3d_seis_lo_10,tower3d_seis_hi_10],
    ];
    ss=20
)

CM.activate!(); plot_tower3d_dpl_traj(
    [
        [tower3d_noseis_10,tower3d_seis_10,tower3d_seis_lo_10,tower3d_seis_hi_10],
    ],
    "tower3d_dpl_traj";
    ss=20
)


function plot_tower3d_vis(bot0,bot1,figname=nothing)
    fig_width = 0.95tw
    fig_height = 0.8cw
    with_theme(theme_pub;
            resolution=(fig_width,fig_height),
            figure_padding=(0,0,0,0),            
            Axis3 = (
                azimuth = 7.525530633326982,
                elevation = 0.12269908169872408
            )
        ) do
        fig = Figure()
        dt = 1e-3
        ts = [1,5.5,10]

        for i in 1:3        
            gridi = fig[1,i] = GridLayout()
            stepj = round(Int,ts[i]/dt)
            TR.goto_step!(bot0,stepj;actuate=true)
            TR.goto_step!(bot1,stepj;actuate=true)
            plot_traj!(bot1;
                AxisType=Axis3,
                fig=gridi,
                # cablecolor=:slategrey,
                atsteps=[stepj],
                showinit=true,
                tgini=deepcopy(bot0.tg),
                showinfo=false,
                doslide=false,
                showlabels=false,
                showpoints=false,
                showground=false,
                xlims = (-0.125,0.125),
                ylims = (-0.125,0.125),
                zlims = (-0.001,0.5),
                titleformatfunc = (sgi,tt) -> "",
                sup! = (ax,tgob,sgi) -> begin
                    rbs = TR.get_bodies(tgob[])
                    # for rbid in 7:10
                    for rbid in [10]
                        if rbid in 7:8
                            pid = 4
                        else
                            pid = 3
                        end
                        meshscatter!(ax,[rbs[rbid].state.rps[pid]], 
                            color = :darkblue, 
                            markersize = 0.008)
                    end
                    ax.xlabel = "x (m)"
                    ax.ylabel = "y (m)"
                    ax.zlabel = "z (m)"
                    ax.yticks = ax.xticks = [-0.1,0,0.1]
                    ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
                    ax.ylabeloffset = ax.xlabeloffset = 40
                    ax.alignmode =  Mixed(;left = 25, right = -25)               
                end
            )

        end

        for (ilabel,label) in enumerate(alphabet[1:3])
            Label(fig.layout[1, ilabel, Bottom()], "($label)",
                fontsize = fontsize,
                font = "CMU Serif Bold",
                padding = (-110, 0, 0, 70),
                halign = :center,
                valign = :top)
        end
        for (ilabel,label) in enumerate(alphabet[1:3])
            Label(fig.layout[1, ilabel, Bottom()], latexstring("t=$(ts[ilabel])~(\\mathrm{s})"),
                fontsize = fontsize,
                padding = (110, 0, 0, 80),
                halign = :center,
                valign = :top)
        end
        colgap!(fig.layout,0)
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plot_tower3d_vis(tower3d_noseis_10,tower3d_seis_10,)
GM.activate!(); plot_tower3d_vis(tower3d_noseis_10,tower3d_seis_10,"tower3d_vis")

function plot_actuate(bot)
    (;t) = bot.traj
    (;actuators) = bot.hub
    act1 = actuators[1]
    itp = act1.pres
    fig_width = columnwidth
    fig_height = 0.4columnwidth
    fig = Figure(;resolution=(fig_width,fig_height),figure_padding= (0,20,0,20))
    ax1 = Axis(fig[1,1];xlabel=L"t~(\mathrm{s})",ylabel=L"\tau")
    μ1 = VectorOfArray(itp.(t))[1,:]
    μ1 .-= μ1[begin]
    τ = μ1./μ1[end]
    lines!(ax1,t,τ;linewidth,color=:blue)
    xlims!(ax1,t[begin],t[end])
    ylims!(ax1,0,1.5)
    ax1.xlabelpadding = -15
    colsize!(fig.layout,1,Fixed(0.6columnwidth))
    fig
end

fig_actuate = plot_actuate(tower3d_noseis)
GM.activate!()
CM.activate!()
CM.save(texroot*raw"\OneDrive - 中山大学\Papers\DynamicTensegrity\CS\images\actuate.pdf", fig_actuate)

fig_tower3d_traj = plot_tower3d_traj(tower3d_noseis)
fig_tower3d_traj = plot_tower3d_traj(tower3d_seis)

GM.activate!()
CM.activate!()
CM.save(texroot*raw"\OneDrive - 中山大学\Papers\DynamicTensegrity\CS\images\tower3d_traj.pdf", fig_tower3d_traj)

tower3dbot0 = tower3d(;k=500.0,c=2.0)

tower3d_seis   = TR.TensegrityRobot(deepcopy(tower3dbot0.tg),(actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
TR.solve!(TR.SimProblem(tower3d_seis,(x)->dynfuncs(x;actuate=true,gravity=true)),
            TR.Zhong06(),(sysstate,tg)->pres3d!(sysstate,tg;T=0.5);
            dt=1e-3,tspan=(0.0,15.0),ftol=1e-14)

tower3dact_i = TR.TensegrityRobot(
                    tower3d(;k=500.0,c=2.0,ijkl=1).tg,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
tower3dact_j = TR.TensegrityRobot(
                    tower3d(;k=500.0,c=2.0,ijkl=2).tg,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
tower3dact_k = TR.TensegrityRobot(
                    tower3d(;k=500.0,c=2.0,ijkl=3).tg,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
tower3dact_l = TR.TensegrityRobot(
                    tower3d(;k=500.0,c=2.0,ijkl=4).tg,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))

dt = 1e-3
TR.solve!(TR.SimProblem(tower3dact_i,(x)->dynfuncs(x;actuate=true,gravity=true)),
            TR.Zhong06(),(sysstate,tg)->pres3d!(sysstate,tg;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)
TR.solve!(TR.SimProblem(tower3dact_j,(x)->dynfuncs(x;actuate=true,gravity=true)),
            TR.Zhong06(),(sysstate,tg)->pres3d!(sysstate,tg;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)
TR.solve!(TR.SimProblem(tower3dact_k,(x)->dynfuncs(x;actuate=true,gravity=true)),
            TR.Zhong06(),(sysstate,tg)->pres3d!(sysstate,tg;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)
TR.solve!(TR.SimProblem(tower3dact_l,(x)->dynfuncs(x;actuate=true,gravity=true)),
            TR.Zhong06(),(sysstate,tg)->pres3d!(sysstate,tg;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)

function plot_compare_ijkl(boti,botj,botk,botl)
    irb9rg = get_trajectory!(boti,9,0)
    jrb9rg = get_trajectory!(botj,9,0)
    krb9rg = get_trajectory!(botk,9,0)
    lrb9rg = get_trajectory!(botl,9,0)
    fig_width = 0.7columnwidth
    fig_height = 0.4columnwidth
    cg = [:black,:red,:blue]
    set_theme!(theme_pub;
            palette = (color = cg, ),
            Lines = (cycle = [:color],),
    )
    fig = Figure(;resolution=(fig_width,fig_height))
    ax1 = Axis(fig[1,1];xlabel=L"t~(\mathrm{s})",ylabel=L"ϵ~(\mathrm{m})")
    # lines!(ax1,boti.traj.t,irb9rg[1,:])
    scaling = -13
    lines!(ax1,botj.traj.t,(jrb9rg[2,:].-irb9rg[2,:])./10.0^scaling,label=L"\mathbf{q}_{rrvw}")
    lines!(ax1,botk.traj.t,(krb9rg[2,:].-irb9rg[2,:])./10.0^scaling,label=L"\mathbf{q}_{rrrw}")
    lines!(ax1,botl.traj.t,(lrb9rg[2,:].-irb9rg[2,:])./10.0^scaling,label=L"\mathbf{q}_{rrrr}")
    # xlims!(ax1,0,5)
    # ylims!(ax1,-7,7)
    ax1.xlabelpadding = -15
    Label(fig.layout[1, 1, Top()], latexstring("×10^{$scaling}"),
        fontsize = fontsize,
        font = "CMU Serif Bold",
        padding = (0, 0, 0, 0),
        halign = :left
    )
    fig[1,2] = Legend(fig,ax1)
    fig
end
fig_compare_ijkl = plot_compare_ijkl(tower3dact_i,tower3dact_j,tower3dact_k,tower3dact_l)
GM.activate!()

me_series = TR.mechanical_energy!(tower3dact;actuate=true,gravity=true)
lines(tower3dact.traj.t,me_series.E)
lines(tower3dact.traj.t,me_series.T)
lines(tower3dact.traj.t,me_series.V)

newspine = newspine3d(2,)

plot_traj!(newspine;showground=false)


## lander
newlander = lander()
plot_traj!(newlander;)

k,μ=TR.optimize_for_stiffness_and_restlen(newlander;fmin=1.0,gravity=false)
k,μ

newlander_opt = lander(;k)
TR.set_restlen!(newlander_opt.tg,μ)
TR.update!(newlander_opt.tg)
TR.get_cables_tension(newlander_opt.tg) |> extrema
TR.check_static_equilibrium_output_multipliers(newlander_opt.tg;gravity=false)
TR.undamped_eigen(newlander_opt.tg;gravity=false)
TR.undamped_eigen!(newlander_opt;gravity=false)
plot_traj!(newlander_opt;)

#- scissor_tower
gravity = true
newst = scissor_tower(θ = π/5)
newst_opt = deepcopy(newst)

plot_traj!(newst;showmesh=false,showground=false)
with_theme(theme_pub;) do 
    plot_traj!(
        newst;
        showmesh=false,
        showinfo=false,
        showlabels=true,
        showpoints=true,
        sup! = plot_one_bar_one_tri!,
        # dorecord=true,
        showinit=false,
        # figname="one_bar_one_tri.mp4"
    )
end

μ = TR.inverse_for_restlength(newst,newst;gravity,fmin=0.0)

@show newst_opt.tg.state.system

k,μ=TR.optimize_for_stiffness_and_restlen(newst;fmin = 10.0,gravity)
k .= 1000.0
TR.set_restlen!(newst_opt.tg,μ)
TR.update!(newst_opt.tg;gravity)

f = TR.get_cables_tension(newst_opt.tg) 

TR.check_static_equilibrium_output_multipliers(newst_opt.tg;gravity=true)
TR.undamped_eigen(newst_opt.tg;gravity=true)
TR.undamped_eigen!(newst_opt;gravity=true)
plot_traj!(newst_opt;showmesh=false,)

newst_folded = scissor_tower(θ = π/5;k)
plot_traj!(newst_folded;showmesh=false,showlabels=true,showpoints=true,showground=false)


k,μ=TR.optimize_for_stiffness_and_restlen(newst_folded;gravity=true)

μ_folded = TR.inverse_for_restlength(newst_folded,newst_folded;gravity=false,fmin=0.0)
TR.set_restlen!(newst_folded.tg,μ_folded)
TR.update!(newst_folded.tg)
TR.get_cables_tension(newst_folded.tg) |> extrema

#-- three 3-prism
prism3 = prism_modules(;n=1,p=3)
bot = prism3

with_theme(theme_pub;
    Poly = (
        transparency=true,
    )
    ) do
    plot_traj!(prism3;
        showlabels=false,
        show_cable_labels=true,
        show_node_labels=false,
        showground=false,
    )
end
k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)
f = TR.get_cables_tension(bot)
q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
# Ň_ = TR.nullspace(Ǎ)
# Ň = modified_gram_schmidt(Ň_)
Ň = TR.make_intrinsic_nullspace(bot.tg,q)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

rank(Ň)
Ǎ*Ň |> norm
# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ

S,D = TR.static_kinematic_determine(ℬᵀ)
ns = size(S,2)
nk = size(D,2)

_,multipliers=  TR.check_static_equilibrium_output_multipliers(bot.tg)
multipliers
ω²,δq̌ = TR.undamped_eigen(bot.tg)
sqrt.(ω²[7:end])./(2π)

# bar
A = ((5e-3)^2-(4e-3)^2)*π
ρ = 1800.0
b = 1.01
m = A*b*ρ

# cable 
A = (0.32e-3)^2*π
ρ = 1435
E = 131e9
# cable diagonal
l = 0.6489
# cable horizontal
l = 0.5464

k = E*A/l
17/k

m = A*l*ρ
