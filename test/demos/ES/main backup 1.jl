#-- preamble
using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using ElasticArrays
using TypeSortedCollections
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
GM.activate!()
using LaTeXStrings
using TexTables
using Latexify
auto_display(true)
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
import Rible as RB
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
figdir::String = raw"C:\Users\luo22\OneDrive\Papers\DynamicTensegrity\ES"
# figdir::String =raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\dyn"

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


## first example, validation and verification
## comparison with Adams
obot = one_bar_one_tri()

prob = RB.DynamicsProblem(obot,(x)->dynfuncs(x;gravity=true))

dt = 1e-3

RB.solve!(prob,RB.Zhong06();dt,tspan=(0.0,4.0),ftol=1e-14)

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
    rb1rp2 = RB.get_trajectory!(bot,1,2)[begin:sp:end]
    θ1 = atan.(rb1rp2[1,:],-rb1rp2[2,:])#.-sol[3,:]
    rb2rg = RB.get_trajectory!(bot,2,0)[begin:sp:end]
    rb2rg_ = rb2rg .- rb1rp2
    θ2 = atan.(rb2rg_[1,:],-rb2rg_[2,:])#.-sol[4,:]

    with_theme(theme_pub;
            size = (0.8cw,0.5cw),
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

rb1rp2 = RB.get_trajectory!(obot,1,2)
θ1 = atan.(rb1rp2[1,:],-rb1rp2[2,:])
rb2rg = RB.get_trajectory!(obot,2,0)
rb2rg_ = rb2rg .- rb1rp2
θ2 = atan.(rb2rg_[1,:],-rb2rg_[2,:])
err_θ1 = (ode_sol[3,:] .- θ1)./θ1 .|> abs
lines(err_θ1)
mean(err_θ1)

obot.traj.t

## comparison with Alpha
obot_zhong = one_bar_one_tri(); obot_alpha = deepcopy(obot)
prob_zhong = RB.DynamicsProblem(obot_zhong,(x)->dynfuncs(x;gravity=true));
prob_alpha = RB.DynamicsProblem(obot_alpha,(x)->dynfuncs(x;gravity=true));

dt = 1e-3
RB.solve!(prob_zhong,RB.Zhong06();dt,tspan=(0.0,500.0),ftol=1e-14)
RB.solve!(prob_alpha,RB.Alpha(0.7);dt,tspan=(0.0,500.0),ftol=1e-14)

plot_traj!(obot_zhong)

function plot_energy!(bot_zhong,bot_alpha,figname=nothing)
    ms_zhong = RB.mechanical_energy!(bot_zhong;gravity=true)
    ms_alpha = RB.mechanical_energy!(bot_alpha;gravity=true)
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

    ms_zhong = RB.mechanical_energy!(bot_zhong;gravity=true)
    ms_alpha = RB.mechanical_energy!(bot_alpha;gravity=true)

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
    rb1rp2 = RB.get_trajectory!(bot,1,2)[begin:sp:end]
    θ1 = atan.(rb1rp2[1,:],-rb1rp2[2,:])#.-sol[3,:]
    rb2rg = RB.get_trajectory!(bot,2,0)[begin:sp:end]
    rb2rg_ = rb2rg .- rb1rp2
    θ2 = atan.(rb2rg_[1,:],-rb2rg_[2,:])#.-sol[4,:]

    with_theme(theme_pub;
            size = (0.8cw,0.3cw),
            figure_padding = (fontsize,fontsize,0,0.5fontsize),
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
        # lines!(ax1,t,θ1;label="Proposed")
        # scatter!(ax1,sol_t,sol_θ[1,:];label="Minimal",marker)
        lines!(ax1,t,θ2;label="Proposed")
        scatter!(ax1,sol_t,sol_θ[2,:];label="Minimal",marker)
        # lines!(ax3,t,ω1./10;label="Proposed")
        # @show length(sol_t), length(sol_ω)
        # scatter!(ax3,sol_t,sol_ω[1,:]./10;label="Minimal",marker)
        lines!(ax2,t,ω2./10;label="Proposed")
        scatter!(ax2,sol_t,sol_ω[2,:]./10;label="Minimal",marker)
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
        lines!(ax3,bot_zhong.traj.t[begin:si:end],ms_zhong.E[begin:si:end]./10.0^scale;label="Proposed",linewidth)
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
## 2nd example, 2D tower
## Natural frequency comparison
tower2dbot = tower2d(;k=100.0,ratio=0.9,ratio1=0.9,slack=false)
RB.undamped_eigen!(tower2dbot;scaling=0.035)
tower2d_frequency = tower2dbot.traj.t[begin+1:end]./2pi

adams_frequency = parse_Adams_frequency("./tower2d_frequency_3.xml")

function plot_frequency_comparison(ref_frequency,bot_frequency,figname=nothing)
    frequency_rel_err = abs.(tower2d_frequency.-adams_frequency)./tower2d_frequency
    @show frequency_rel_err, frequency_rel_err |> maximum
    # fig_width = cw
    # fig_height = 0.4cw
    with_theme(theme_pub;
            # resolution=(0.6tw,0.3tw),
            resolution=(cw,0.4cw),
            palette = (color = [:black,:red],),
            cycle = [[:linecolor, :markercolor] => :color,],
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1]; ylabel = L"\mathrm{Frequency}~(\mathrm{Hz})", xlabel = L"\mathrm{Mode}")
        scatterlines!(bot_frequency;label="Proposed",marker=:+,linewidth,linestyle=:dashdot)
        scatterlines!(ref_frequency;label="Adams",marker='×',linewidth=0)
        ylims!(ax1,0,30)
        ax1.yticks = collect(0:5:30)
        # ax1.xlabelpadding = -15
        # axislegend(ax1;position=:lt)
        fig[1,2] = Legend(fig,ax1)
        colsize!(fig.layout,1,Fixed(0.6cw))
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plot_frequency_comparison(adams_frequency,tower2d_frequency)
CM.activate!(); plot_frequency_comparison(adams_frequency,tower2d_frequency,"frequency")

## Mode shape
function plot_tower2d_mode_shape(bot,figname=nothing)
    (;traj,st) = bot
    (;q,t) = traj
    shapes = q
    frequencies = t[begin+1:end]./2pi
    # @sprintf "Mode %s \n %.4f (Hz)"  1  tower2d_frequency[1]
    with_theme(theme_pub;
            # resolution=(tw,0.3tw),
            resolution=(0.8tw,0.5cw),
            Axis = (
                titlefont = "CMU Serif",
            )
        ) do
        fig = Figure(;)
        axs = [
            begin
                ax = Axis(
                    fig[1,i];
                    # visible = true,
                    aspect = DataAspect(),
                    limits = (nothing,(0,0.4)),
                    titlegap = -40,
                    title = @sprintf "Mode %s \n %.4f (Hz)"  i  frequencies[i]
                )
                hidedecorations!(ax)
                hidespines!(ax)
                ax
            end
            for i = 1:5
        ]
        RB.update_bodies!(st,shapes[1])
        plot_tower2d!(axs[1],st;isref=true)
        plot_tower2d!(axs[2],st;isref=true)
        plot_tower2d!(axs[3],st;isref=true)
        plot_tower2d!(axs[4],st;isref=true)
        plot_tower2d!(axs[5],st;isref=true)

        RB.update_bodies!(st,shapes[2])
        plot_tower2d!(axs[1],st)
        RB.update_bodies!(st,shapes[3])
        plot_tower2d!(axs[2],st)
        RB.update_bodies!(st,shapes[4])
        plot_tower2d!(axs[3],st)
        RB.update_bodies!(st,shapes[5])
        plot_tower2d!(axs[4],st)
        RB.update_bodies!(st,shapes[6])
        plot_tower2d!(axs[5],st)
        colgap!(fig.layout,0)
        savefig(fig,figname)
    end
end
GM.activate!(); plot_tower2d_mode_shape(tower2dbot)
CM.activate!(); plot_tower2d_mode_shape(tower2dbot,"mode_shapes")

plot_traj!(tower2dbot)

## Natural frequency with varying rest lengths
function vary_restlengths(ratio_range)
    VectorOfArray([
        begin
            bot = tower2d(;k=100.0,ratio,slack=false)
            RB.undamped_eigen!(bot)
            _frequency = bot.traj.t[begin+1:end]./2pi
        end
        for ratio in ratio_range
    ])
end
ratio_range = collect(0.845:0.001:0.999)
frequencies = vary_restlengths(ratio_range)
adams_frequencies = VectorOfArray(
    [
        parse_Adams_frequency("./tower2d_frequency_$i.xml")
        for i = 1:6
    ]
)

adams_ratio = [
    parse_Adams_static("./tower2d_frequency_$i.xml").VARIABLE_ratio
    for i = 1:6
]

function plot_frequency_varying_restlengths(figname=nothing)
    cg = cgrad(:Dark2_6, 6, categorical = true)[[1,2,3,4,6]]
    mk = ['⨉','○','◇','✳','◻']
    with_theme(theme_pub;
            # resolution=(0.7tw,0.4tw),
            resolution=(cw,0.5cw),
            palette = (color = cg, marker = mk),
            Lines = (cycle = [:color],),
            Scatter = (cycle = Cycle([:color, :marker, :strokecolor=>:color],covary=true),),
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1];
                ylabel = L"\mathrm{Frequency}~(\mathrm{Hz})",
                xlabel = L"\alpha~\mathrm{and}~\beta",
                )
        nmode = size(frequencies,1)
        for i = 1:nmode
            lines!(ax1,ratio_range,frequencies[i,:];linewidth,label="Proposed")
            scatter!(ax1,adams_ratio,adams_frequencies[i,:];markersize,strokewidth=0,label="Adams (Mode $i)")
        end
        ylims!(ax1,0,30)
        ax1.yticks = collect(0:5:30)
        xlims!(ax1,0.84,1.0)
        ax1.xlabelpadding = -15
        legend_elems = [
            [
                LineElement(;color = cg[i],linewidth),
                MarkerElement(;color = cg[i],marker=mk[i],markersize,strokewidth=0,strokecolor=cg[i])
            ]
            for i in nmode:-1:1
        ]
        legend_names = [
            "Mode $i" for i in length(legend_elems):-1:1
        ]

        fig[1, 2] = Legend(fig, legend_elems, legend_names)
        # fig[1, 2] = Legend(fig, ax1, nbanks = 2)
        # colsize!(fig.layout,1,Fixed(0.6columnwidth))
        savefig(fig,figname)
    end
end

GM.activate!(); plot_frequency_varying_restlengths()
CM.activate!(); plot_frequency_varying_restlengths("frequency_varying_restlengths")

## seismic simulation

function pres2d!(sysstate,st;aux=true)
    (;t,q̃,q̃̇,q̃̈) = sysstate
    (;bodies,connectivity) = st
    (;bodyid2sys_pres_coords) = connectivity.indexed
    amp = 0.01
    p = 6π
    foreach(bodies) do body
        bodyid = body.prop.id
        if bodyid == 1
            q̃[bodyid2sys_pres_coords[bodyid]] .= [     amp*sin(p*t),0]
            q̃̇[bodyid2sys_pres_coords[bodyid]] .= [   p*amp*cos(p*t),0]
            q̃̈[bodyid2sys_pres_coords[bodyid]] .= [-p^2*amp*sin(p*t),0]
        end
        if bodyid == 2
            q̃[bodyid2sys_pres_coords[bodyid]] .= [ 0.1+amp*sin(p*t),0]
            q̃̇[bodyid2sys_pres_coords[bodyid]] .= [   p*amp*cos(p*t),0]
            q̃̈[bodyid2sys_pres_coords[bodyid]] .= [-p^2*amp*sin(p*t),0]
        end
        if aux
            if bodyid == 8
                q̃[bodyid2sys_pres_coords[bodyid][3:4]] .= [-0.1+amp*sin(p*t),0]
                q̃̇[bodyid2sys_pres_coords[bodyid][3:4]] .= [   p*amp*cos(p*t),0]
                q̃̈[bodyid2sys_pres_coords[bodyid][3:4]] .= [-p^2*amp*sin(p*t),0]
                q̃[bodyid2sys_pres_coords[bodyid][5:6]] .= [ 0.2+amp*sin(p*t),0]
                q̃̇[bodyid2sys_pres_coords[bodyid][5:6]] .= [   p*amp*cos(p*t),0]
                q̃̈[bodyid2sys_pres_coords[bodyid][5:6]] .= [-p^2*amp*sin(p*t),0]
            end
        end
    end
end

ratio1 = 0.998
dt = 1e-3
## undamped noslack
tower2dbot_undamped_noslack = tower2d(;k=100.0,c=0.0,ratio=0.95,ratio1,slack=false)
RB.solve!(RB.DynamicsProblem(tower2dbot_undamped_noslack,(x)->dynfuncs(x;gravity=true)),
          RB.Zhong06(),
          (
            prescribe! = (sysstate,st)->pres2d!(sysstate,st;aux=true),
            actuate! = nothing
          );
          dt,tspan=(0.0,2.0),ftol=1e-14,verbose=true,exception=false)

adams_undamped_noslack = parse_Adams_dynamic("tower2d_undamped_noslack.xml")

plot_traj!(tower2dbot_undamped_noslack)


## damped noslack
tower2dbot_damped_noslack = tower2d(;k=100.0,c=1.0,ratio=0.95,ratio1,slack=false)
RB.solve!(RB.DynamicsProblem(tower2dbot_damped_noslack,(x)->dynfuncs(x;gravity=true)),
          RB.Zhong06(),
          (
            prescribe! = (sysstate,st)->pres2d!(sysstate,st;aux=true),
            actuate! = nothing
          );
          dt,tspan=(0.0,2.0),ftol=1e-14)

adams_damped_noslack = parse_Adams_dynamic("tower2d_damped_noslack.xml")

## undamped slack
tower2dbot_undamped_slack = tower2d(;k=100.0,c=0.0,ratio=0.95,ratio1,slack=true)
f = RB.get_cables_tension(tower2dbot_undamped_slack)
l = RB.get_cables_len(tower2dbot_undamped_slack)
μ = RB.get_cables_restlen(tower2dbot_undamped_slack)
k = RB.get_cables_stiffness(tower2dbot_undamped_slack)
x1,y1 = approx_slack(;k=k[1],μ=μ[1],l0=l[1],filename="xy1")
x2,y2 = approx_slack(;k=k[3],μ=μ[3],l0=l[3],filename="xy2")
x3,y3 = approx_slack(;k=k[5],μ=μ[5],l0=l[5],filename="xy3")
x4,y4 = approx_slack(;k=k[11],μ=μ[11],l0=l[11],filename="xy4")

x0,y0 = approx_slack(;k=100.0,μ=0.0,l0=0.0)
function plot_data_points(x,y,figname=nothing)
    with_theme(theme_pub;
        size = (cw,0.8cw),
        Scatter = (
            marker=:xcross,
        )
    ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1],
            xlabel = L"\Delta l_j~(\mathrm{m})",
            ylabel = L"f_j~(\mathrm{N})",
        )
        ax2 = Axis(fig[2,1],
            # xlabel = L"l_j-\mu_j",
            # ylabel = L"f_j",
            tellwidth = false,
            tellheight = false,
            width = 0.225cw,
            height = 0.16cw,
            xlabelsize = 5 |> pt2px,
            ylabelsize = 5 |> pt2px,
            xticklabelsize = 5 |> pt2px,
            yticklabelsize = 5 |> pt2px,

        )
        lines!(ax1, x, y; label = "ideal curve")
        scatter!(ax1, x, y; markersize = fontsize*0.8, color = :red, label = "data points")
        ylims!(ax1,-0.5,10.5)
        lines!(ax2, x, y, )
        scatter!(ax2, x, y; markersize = fontsize*0.6, color = :red,)
        xlims!(ax2, -0.002, 0.002)
        ylims!(ax2,-0.02,0.20)
        fig[1,2] = Legend(fig,ax1;)
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plot_data_points(x0,y0)
CM.activate!(); plot_data_points(x0,y0,"f_curve")

RB.solve!(RB.DynamicsProblem(tower2dbot_undamped_slack,(x)->dynfuncs(x;gravity=true)),
          RB.Zhong06(),
          (
            prescribe! = (sysstate,st)->pres2d!(sysstate,st;aux=true),
            actuate! = nothing
          );
          dt,tspan=(0.0,3.0),ftol=1e-14,exception=false)

adams_undamped_slack = parse_Adams_dynamic("tower2d_undamped_slack.xml")
plot_traj!(tower2dbot_undamped_slack;showmesh)

## undamped slack with fix aux
tower2dbot_undamped_slack_noaux = tower2d(;k=100.0,c=0.0,ratio=0.95,ratio1,slack=true)
RB.solve!(RB.DynamicsProblem(tower2dbot_undamped_slack_noaux,(x)->dynfuncs(x;gravity=true)),
          RB.Zhong06(),
          (
            prescribe! = (sysstate,st)->pres2d!(sysstate,st;aux=false),
            actuate! = nothing
          );
          dt,tspan=(0.0,3.0),ftol=1e-14)
adams_undamped_slack_aux = parse_Adams_dynamic("tower2d_undamped_slack_noaux.xml")

function plot_tower2d_point_traj(ba,figname=nothing;at=1.23)
    cg = cgrad(:seaborn_bright6, 6, categorical = true)
    mk = ['○','□','×','+']
    with_theme(theme_pub;
            resolution=(tw,0.5tw),
            palette = (color = cg, marker = mk),
            Lines = (cycle = [:color],),
            Scatter = (cycle = Cycle([:color, :marker, :strokecolor=>:color],covary=true),),
    ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1]; xlabel = L"t~(\mathrm{s})", ylabel = L"x~(\mathrm{m})")
        linewidth = 0.75 |> pt2px
        markersizes = [8,8,10,10] .|> pt2px
        begins = [1,21,31,41]
        cases = ["UN","DN","US","US-AUX"]
        si = 50
        for (bot,adams,markersize,be,case) in zip(ba.is...,markersizes,begins,cases)
            rg = RB.get_trajectory!(bot,3,0)
            if case == "US"
                step_collapse = findfirst((x)->x>=at,bot.traj.t)
                lines!(ax1,bot.traj.t[1:step_collapse],rg[1,1:step_collapse];linewidth,label="Proposed($case)")
            else
                lines!(ax1,bot.traj.t,rg[1,:];linewidth,label="Proposed($case)")
            end
            if case in ["UN","DN"]
                strokewidth = 0
            else
                strokewidth = 0
            end
            scatter!(ax1,adams.time[be:si:end],
                    adams.PART_3_cm_x[be:si:end];
                    markersize,
                    strokewidth,
                    label="Adams($case)"
            )
        end
        vlines!(ax1,[at],color=:slategrey)
        text!(ax1,latexstring("t=$at~(\\mathrm{Collapse})"), fontsize = fontsize, position = (at+0.02, 0.015), align = (:left, :center))
        xlims!(ax1,0.,2.)
        ylims!(ax1,-0.005,0.092)
        ax1.yticks = collect(0.00:0.01:0.09)
        # ax1.xlabelpadding = -15
        axislegend(ax1,position=:lb,
                    orientation=:horizontal,
                    nbanks = 2,
                    labelsize=round(Int,0.8fontsize))
        savefig(fig,figname)
    end
end

GM.activate!(); plot_tower2d_point_traj(
    zip(
        [tower2dbot_undamped_noslack,tower2dbot_damped_noslack,tower2dbot_undamped_slack,tower2dbot_undamped_slack_noaux],
        [adams_undamped_noslack,adams_damped_noslack,adams_undamped_slack,adams_undamped_slack_aux]
    )
)

CM.activate!(); plot_tower2d_point_traj(
    zip(
        [tower2dbot_undamped_noslack,tower2dbot_damped_noslack,tower2dbot_undamped_slack,tower2dbot_undamped_slack_noaux],
        [adams_undamped_noslack,adams_damped_noslack,adams_undamped_slack,adams_undamped_slack_aux]
    ),
    "tower2d_point_traj"
)

function plot_slack_cables(bots,figname=nothing;at=1.23)
    cg = cgrad(:mk_12, 12, categorical = true)
    with_theme(theme_pub;
            resolution=(0.7tw,0.5tw),
            # figure_padding=(0,50,0,0),
            palette = (color = cg, ),
            Lines = (cycle = [:color],),
            Scatter = (
                markersize = fontsize/3,
            )
        ) do
        fig = Figure(;)
        linewidth = 1 |> pt2px
        for i = 1:2
            ax = Axis(fig[i,1]; xlabel = tlabel, ylabel = L"\mathrm{Cable}~j")
            vlines!(ax,[at],color=:slategrey)
            if i == 1
                ax.xlabelvisible = false
                ax.xticklabelsvisible = false
                ax.yticks = collect(1:12)
                ylims!(ax,0.5,12.5)
                text!(ax,latexstring("t=$at~(\\mathrm{Collapse})"), fontsize = fontsize,
                            position = (at+0.055, 4.5),
                            align = (:left, :center)
                )
            else
                ax.yticks = collect(1:4)
                ylims!(ax,0.5,4.5)
                at = 2.0
            end
            bot = bots[i]
            (;cables) = bot.st.tensiles
            ncables = length(cables)
            f = [get_tension!(bot,j) for j in collect(1:12)]
            t_mids = get_time_mids(bot)
            step_at = findfirst((x)->x>=at, t_mids)
            noslack = [findall((x)->x>0,fj[1:step_at]) for (j,fj) in enumerate(f)]
            slack   = [findall((x)->x<=0,fj[1:step_at]) for (j,fj) in enumerate(f)]
            for (i,(ns,s)) in enumerate(zip(noslack,slack))
                scatter!(ax,t_mids[1:step_at][s],fill(i,length(s));linewidth)
            end
            xlims!(ax,0,2)
            # ax.xlabelpadding = -15
            # ax.alignmode = Mixed(;left = -20, right = 20)
            # ax.xticks = [0.275]
        end
        for (ilabel,label) in enumerate(alphabet[1:2])
            Label(fig.layout[ilabel, 1, TopLeft()], "($label)",
                fontsize = fontsize,
                font = "CMU Serif Bold",
                padding = (10, 0, -10, 0),
                halign = :left)
        end
        rowsize!(fig.layout, 2, Relative(1/4))
        rowgap!(fig.layout, 0)
        savefig(fig,figname)
    end
end
GM.activate!(); plot_slack_cables(
    [
        tower2dbot_undamped_slack,
        tower2dbot_undamped_slack_noaux
    ]
)
CM.activate!(); plot_slack_cables(
    [
        tower2dbot_undamped_slack,
        tower2dbot_undamped_slack_noaux
    ],
    "slack_cable_tensions"
)

## At time
function plot_tower2d_at_time(bots,figname=nothing;at=1.23)
    with_theme(theme_pub;
            resolution=(0.8tw,0.25tw),
        ) do
        fig = Figure(;)

        for i = 1:4
            ax = Axis(
                fig[1,i];
                # visible = true,
                aspect = DataAspect(),
                limits = ((-0.12,0.22),(0,0.36)),
            )
            hidedecorations!(ax)
            hidespines!(ax)
            RB.goto_step!(bots[i],1)
            plot_tower2d!(ax,bots[i].st;isref=true,markit=true)
            istep = findfirst((x)->x>=at, bots[i].traj.t)
            @show istep
            RB.goto_step!(bots[i],istep)
            # @show RB.get_cables_tension(bots[i].st)
            plot_tower2d!(ax,bots[i].st,markit=true)
        end
        for (ilabel,label) in enumerate(alphabet[1:4])
            Label(fig.layout[1, ilabel, TopLeft()], "($label)",
                fontsize = fontsize,
                font = "CMU Serif Bold",
                padding = (0, -50, -50, 0),
                halign = :right)
        end
        colgap!(fig.layout,0)
        rowgap!(fig.layout,0)
        savefig(fig,figname)
    end
end

GM.activate!(); plot_tower2d_at_time(
    [
        tower2dbot_undamped_noslack,
        tower2dbot_damped_noslack,
        tower2dbot_undamped_slack,
        tower2dbot_undamped_slack_noaux
    ]
)
CM.activate!(); plot_tower2d_at_time(
    [
        tower2dbot_undamped_noslack,
        tower2dbot_damped_noslack,
        tower2dbot_undamped_slack,
        tower2dbot_undamped_slack_noaux
    ],
    "tower2d_at_time"
)

ms_tower2d = RB.mechanical_energy!(tower2dbot_undamped_noslack;gravity=true)

## simple
gravity = false
newsim = simple(;
    c=0.0,
    z0 = 0.2,
    ωz = 5.0,
    mbar = 0.05,
    free = true
)

plot_traj!(newsim)

RB.solve!(
    RB.DynamicsProblem(newsim,(x)->dynfuncs(x;gravity)),
    RB.Zhong06();
    dt=1e-3,tspan=(0.0,15.0),ftol=1e-14
)
# plot_tower3d_traj(newsim)
plot_traj!(newsim;)

me = RB.mechanical_energy!(newsim;gravity)
me.E |> lines



## tower3d
## Modal analysis
# tower3dbot = tower3d(;k=500.0,c=10.0,θ=π/6)

## Seismic and deploy
function pres3d!(sysstate,st;ν=0.5)
    (;t,q̃,q̃̇,q̃̈) = sysstate
    (;bodies,connectivity) = st
    (;bodyid2sys_pres_coords) = connectivity.indexed
    amp = 0.01
    p = ν*2π
    foreach(bodies) do body
        bodyid = body.prop.id
        if bodyid in [1,2,3]
            r0 = zeros(3)
            if bodyid == 1
                r0 .= [
                    0.1,
                    0.0,
                    0.0
                ]
                q̃[bodyid2sys_pres_coords[bodyid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[bodyid2sys_pres_coords[bodyid]] .=     [amp*p*cos(p*t),0,0]
                q̃̈[bodyid2sys_pres_coords[bodyid]] .=  [-amp*p^2*sin(p*t),0,0]
            elseif bodyid == 2
                r0 .= [
                    -0.05,
                     0.08660254037844387,
                    0.0
                ]
                q̃[bodyid2sys_pres_coords[bodyid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[bodyid2sys_pres_coords[bodyid]] .=     [amp*p*cos(p*t),0,0]
                q̃̈[bodyid2sys_pres_coords[bodyid]] .=  [-amp*p^2*sin(p*t),0,0]
            elseif bodyid == 3
                r0 .= [
                    -0.05,
                    -0.08660254037844387,
                    0.0
                ]
                q̃[bodyid2sys_pres_coords[bodyid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[bodyid2sys_pres_coords[bodyid]] .=     [amp*p*cos(p*t),0,0]
                q̃̈[bodyid2sys_pres_coords[bodyid]] .=  [-amp*p^2*sin(p*t),0,0]
            end
        end
    end
end

rest = deleteat!(collect(1:24),13:21)

## inverse statics and modal analysis
tower3dbot0 = tower3d(;k=500.0,k1=1000,c=2.0,α=π/12)
# B,F̃ = RB.build_inverse_statics_for_restlength(tower3dbot0.st,tower3dbot0.st;gravity=true)
# rank(B),size(B)
# x0,xn = RB.get_solution_set(B,F̃)
# size(xn)

plot_traj!(tower3dbot0;fontsize=20)
μ0 = RB.inverse_for_restlength(tower3dbot0,tower3dbot0;gravity=true,scale=true,fmin=1.0,verbose=true)
RB.set_restlen!(tower3dbot0.st,μ0)
RB.update!(tower3dbot0.st)
RB.get_cables_tension(tower3dbot0.st) |> extrema
RB.check_static_equilibrium_output_multipliers(tower3dbot0.st;gravity=true)
RB.undamped_eigen(tower3dbot0.st;gravity=true)
RB.undamped_eigen!(tower3dbot0;gravity=true)
tower3dbot0.traj.t./(2π)
@show (tower3dbot0.traj.t./(2π))[2]

tower3dbot1 = tower3d(;k=500.0,c=2.0,d = 0.1*√2, r2 = 0.07, α=-0.0)
plot_traj!(tower3dbot1;fontsize=20)
μ1 = RB.inverse_for_restlength(tower3dbot1,tower3dbot1;gravity=true,fmin=1.0)
RB.set_restlen!(tower3dbot1.st,μ1)
RB.update!(tower3dbot1.st)
RB.get_cables_tension(tower3dbot1.st) |> extrema
RB.check_static_equilibrium_output_multipliers(tower3dbot1.st;gravity=true)
# RB.undamped_eigen(tower3dbot1.st;gravity=true)
RB.undamped_eigen!(tower3dbot1;gravity=true,scaling=0.2)
tower3dbot1.traj.t./(2π)
@show (tower3dbot1.traj.t./(2π))[2]

GM.activate!();
with_theme(theme_pub;
        Axis3 = (
            azimuth = 7.00553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    plot_traj!(
        tower3dbot1;
        AxisType=Axis3,
        xlims=(-0.2,0.2),
        ylims=(-0.2,0.2),
        zlims=(-1e-3,0.5),
        showmesh=false,
        showpoints=false,
        showlabels=false,
        showinfo=false,
        showinit=true,
        dorecord=false,
        # figname="tower3d_seis_hi_15.mp4",
    )
end

@fmt Real = "{:}"
initial_μ = TableCol("Initial", ["$((3i-2,3i-1,3i))" for i in 1:8], ["\\qty{$μ}{\\meter}" for μ in μ0[begin:3:end]])
target_μ = TableCol("Target", ["$((3i-2,3i-1,3i))" for i in 1:8], ["\\qty{$μ}{\\meter}" for μ in μ1[begin:3:end]])
table_μ = hcat(initial_μ,target_μ)
print(TexTables.to_tex(table_μ))

fig_decompose_tower3d = plot_decompose_tower3d(tower3dbot0,tower3dbot1)
GM.activate!()
CM.activate!()
CM.save(texroot*raw"\OneDrive - 中山大学\Papers\DynamicTensegrity\CS\images\decompose_tower3d.pdf", fig_decompose_tower3d)

fig_compose_tower3d = plot_compose_tower3d(tower3dbot0,tower3dbot1)
GM.activate!()
CM.activate!()
CM.save(texroot*raw"\OneDrive - 中山大学\Papers\DynamicTensegrity\CS\images\compose_tower3d.pdf", fig_compose_tower3d)

function plot_restlengths(bots,figname=nothing)
    friction_coefficients = [RB.get_cables_restlen(bot) for bot in bots]
    # cg = cgrad(:Dark2_6, 6, categorical = true)[[4,3,6]]
    # set_theme!(theme_pub;
    #         palette = (color = cg, ),
    #         Lines = (cycle = [:color],),
    # )
    with_theme(theme_pub;
            # size = (0.7tw,0.3tw),
            # font = "Nimbus Rom No9 L",
            # font = "Times New Roman",
            size = (2cw,0.7cw),
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1]; 
            xlabel = "Cable No.", 
            ylabel = "Rest length (m)"
        )
        x = vec([collect(1:24)';collect(1:24)'])
        grp = repeat([1,2],24)
        display(x)
        height = vec([friction_coefficients[1]';friction_coefficients[2]'])
        colors = Makie.wong_colors()
        barplot!(
            ax1, x, height; 
            dodge = grp,
            color = colors[grp],
            marker = '+', markersize = 1.5fontsize, 
            linewidth, linestyle = :dashdot, 
            # label = "Folded"
        )
        # ax1.yticks = collect(0:15)
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

GM.activate!(); plot_restlengths([tower3dbot0,tower3dbot1])

function plot_tensions(bots,figname=nothing)
    fs = [RB.get_cables_tension(bot) for bot in bots]
    # cg = cgrad(:Dark2_6, 6, categorical = true)[[4,3,6]]
    # set_theme!(theme_pub;
    #         palette = (color = cg, ),
    #         Lines = (cycle = [:color],),
    # )
    with_theme(theme_pub;
            # size = (0.7tw,0.3tw),
            # font = "Nimbus Rom No9 L",
            # font = "Times New Roman",
            size = (2cw,0.7cw),
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
            # size = (0.7tw,0.3tw),
            # font = "Nimbus Rom No9 L",
            # font = "Times New Roman",
            size = (cw,0.5cw),
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
RB.solve!(RB.DynamicsProblem(tower3dbot0_nodpl,(x)->dynfuncs(x;gravity=true)),
            RB.Zhong06(),
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
RB.solve!(RB.DynamicsProblem(tower3dbot1_nodpl,(x)->dynfuncs(x;gravity=true)),
            RB.Zhong06(),
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
            rb7r4 = RB.get_trajectory!(bot,7,4,1:ss:length(bot.traj))
            rb8r4 = RB.get_trajectory!(bot,8,4,1:ss:length(bot.traj))
            rb9r3 = RB.get_trajectory!(bot,9,3,1:ss:length(bot.traj))
            # rb10r3 = RB.get_trajectory!(bot,10,3,1:ss:length(bot.traj))
            lines!(ax1,t,(rb7r4[1,:].-rb7r4[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{7,1}")
            lines!(ax1,t,(rb8r4[1,:].-rb8r4[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{8,1}")
            lines!(ax1,t,(rb9r3[1,:].-rb9r3[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{9,1}")
            # lines!(ax1,t,(rb10r3[1,:].-rb10r3[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{10,1}")
            ylims!(ax1,-5,5)
            ax1.yticks = collect(-5:2.5:5)

            lines!(ax2,t,(rb7r4[2,:].-rb7r4[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb8r4[2,:].-rb8r4[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb9r3[2,:].-rb9r3[2,1])./10.0^scalings[2];)
            # lines!(ax2,t,(rb10r3[2,:].-rb10r3[2,1])./10.0^scalings[2];)

            lines!(ax3,t,(rb7r4[3,:].-rb7r4[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb8r4[3,:].-rb8r4[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb9r3[3,:].-rb9r3[3,1])./10.0^scalings[3];)
            # lines!(ax3,t,(rb10r3[3,:].-rb10r3[3,1])./10.0^scalings[3];)
            # @show rb10r3[3,1]
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
            RB.goto_step!(boti,1)
            plot_traj!(boti;
                AxisType=Axis3,
                fig=gridi,
                # cablecolor=:slategrey,
                atsteps=[stepj],
                showinit=true,
                tgini=deepcopy(boti.st),
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
                    rbs = RB.get_bodies(tgob[])
                    for bodyid in 7:10
                        # for bodyid in [10]
                        if bodyid in 7:8
                            pid = 4
                        else
                            pid = 3
                        end
                        meshscatter!(ax,[rbs[bodyid].state.loci_states[pid]], 
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
                            meshscatter!(ax,[rbs[key].state.loci_states[value]], 
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
Td = [5.0,10.0,15.0,20.0]
function do_noseis(;Td,tend=25.0)
    bot   = RB.Robot(deepcopy(tower3dbot0.st),(actuators=[make_pres_actor(μ0,μ1,0.0,Td)],))
    RB.solve!(RB.DynamicsProblem(bot,(x)->dynfuncs(x;actuate=true,gravity=true)),
                RB.Zhong06();
                dt=1e-3,tspan=(0.0,tend),ftol=1e-14)
    bot
end
function do_seis(;ν,Td,tend=25.0)
    bot   = RB.Robot(deepcopy(tower3dbot0.st),(actuators=[make_pres_actor(μ0,μ1,0.0,Td)],))
    RB.solve!(RB.DynamicsProblem(bot,(x)->dynfuncs(x;actuate=true,gravity=true)),
                RB.Zhong06(),
                (
                    prescribe! = (sysstate,st)->pres3d!(sysstate,st;ν),
                    actuate! = nothing
                );
                dt=1e-3,tspan=(0.0,tend),ftol=1e-14)
    bot
end
tower3d_noseis_5 = do_noseis(;Td=Td[1])
tower3d_seis_5 = do_seis(;ν=0.5,Td=Td[1])
tower3d_seis_lo_5 = do_seis(;ν=0.2,Td=Td[1])
tower3d_seis_hi_5 = do_seis(;ν=1.0,Td=Td[1],tend=7.0)
tower3d_noseis_10 = do_noseis(;Td=Td[2])
tower3d_seis_10 = do_seis(;ν=0.5,Td=Td[2])
tower3d_seis_lo_10 = do_seis(;ν=0.2,Td=Td[2])
tower3d_seis_hi_10 = do_seis(;ν=1.0,Td=Td[2])
tower3d_noseis_15 = do_noseis(;Td=Td[3])
tower3d_seis_15 = do_seis(;ν=0.5,Td=Td[3])
tower3d_seis_lo_15 = do_seis(;ν=0.2,Td=Td[3])
tower3d_seis_hi_15 = do_seis(;ν=1.0,Td=Td[3])
tower3d_noseis_20 = do_noseis(;Td=Td[4])
tower3d_seis_20 = do_seis(;ν=0.5,Td=Td[4],tend=21.7)
tower3d_seis_lo_20 = do_seis(;ν=0.2,Td=Td[4])
tower3d_seis_hi_20 = do_seis(;ν=1.0,Td=Td[4])

GM.activate!();
with_theme(theme_pub;
        Axis3 = (
            azimuth = 7.00553063332698,
            elevation = 0.22269908169872407
        )
    ) do
    plot_traj!(
        tower3d_seis_15;
        AxisType=Axis3,
        xlims=(-0.2,0.2),
        ylims=(-0.2,0.2),
        zlims=(-1e-3,0.5),
        showmesh=false,
        showpoints=false,
        showlabels=false,
        showinfo=false,
        dorecord=true,
        actuate=true,
        figname="tower3d_seis_15.mp4",
    )
end

plot_traj!(tower3d_seis_hi_15;actuate=true)
plot_traj!(tower3d_noseis_15;actuate=true)
plot_traj!(tower3d_seis_hi;actuate=true)
me_series = RB.mechanical_energy!(tower3d_seis;actuate=true,gravity=true)
lines(tower3d_seis.traj.t,me_series.E)

rb10r3 = RB.get_trajectory!(tower3d_seis_hi_20,10,3,1:10:length(tower3d_seis_hi_20.traj))

rb10r3
function plot_tower3d_dpl_traj(bots,figname=nothing;
        ss=10,
    )
    cg = cgrad(:seaborn_bright6, 6, categorical = true)
    with_theme(theme_pub;
            # size = (tw,0.8tw),
            size = (tw,1.25cw),
            palette = (color = cg, ),
            Lines = (cycle = [:color], linewidth = 2),
        ) do
        fig = Figure(;)
        ats = [10.5,14.0,18.0]
        Tds = [10.0,15.0,20.0]
        axs = [
            [
                Makie.Axis(fig[i,1];xlabel=L"t~(\mathrm{s})",ylabel=L"x~(\mathrm{m})"),
                Makie.Axis(fig[i,2];xlabel=L"t~(\mathrm{s})",ylabel=L"z~(\mathrm{m})"),
                # Makie.Axis(fig[i,3];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta z~(\mathrm{m})")
            ]
            for i = 1:3
        ]
        scalings = [-1,-1]

        for i = 1:3
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
            bot0rb10r3 = RB.get_trajectory!(bot0,10,3,1:ss:length(bot0.traj))
            # bot1
            bot1rb10r3 = RB.get_trajectory!(bot1,10,3,tsteps)
            # botlo
            botlorb10r3 = RB.get_trajectory!(botlo,10,3,1:ss:length(botlo.traj))
            # bothi
            bothirb10r3 = RB.get_trajectory!(bothi,10,3,1:ss:length(bothi.traj))

            # x
            lines!(ax1,t0,bot0rb10r3[1,:]./10.0^scalings[1];linewidth)
            lines!(ax1,t0,botlorb10r3[1,:]./10.0^scalings[1];linewidth)
            lines!(ax1,t1,bot1rb10r3[1,:]./10.0^scalings[1];linewidth)
            lines!(ax1,t0,bothirb10r3[1,:]./10.0^scalings[1];linewidth)
            ylims!(ax1,0,2.0)
            # ax1.yticks = collect(-5:2.5:5)
            xlims!(ax1,t0[begin],t0[end])

            # # y
            # lines!(ax2,t0,bot0rb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Static~ground}")
            # lines!(ax2,t0,botlorb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.2~(\mathrm{Hz})")
            # lines!(ax2,t0,bothirb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 1.0~(\mathrm{Hz})")
            # lines!(ax2,t1,bot1rb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.5~(\mathrm{Hz})")
            # # ylims!(ax2,-0.5,7.5)
            # xlims!(ax2,t0[begin],t0[end])
            # ax2.yticks = collect(-1.0:1.0:7.5)

            # z
            lines!(ax2,t0,bot0rb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Static~ground}")
            lines!(ax2,t0,botlorb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.2~(\mathrm{Hz})")
            lines!(ax2,t1,bot1rb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.5~(\mathrm{Hz})")
            lines!(ax2,t0,bothirb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 1.0~(\mathrm{Hz})")
            ylims!(ax2,3.5,5.5)
            xlims!(ax2,t0[begin],t0[end])

            # rb10r3_z0 = 0.49497474683058335
            # lines!(ax3,t0,(bot0rb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            # lines!(ax3,t0,(botlorb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            # lines!(ax3,t1,(bot1rb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            # lines!(ax3,t0,(bothirb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            # if i == 1
            #     ylims!(ax3,-1.5,1.5)
            # else
            #     ylims!(ax3,-.8,1.0)
            # end
            # xlims!(ax3,Td,Td+2)

            # if i == 2
            fig[4,1:2] = Legend(fig,ax2;orientation=:horizontal,)
            # end labelsize=pt2px(6.5)
            function set_ax!(ax)
                # xlims!(ax,0,5)
                # ax.alignmode =  Mixed(;left = -20, right = 20)
                vlines!(ax,[Td];linewidth,color=:slategrey)
                # ax.xticks = collect(0:1:15)
                ax.xminorticks = collect(0:1:25)
                ax.xminorgridvisible = true
                if i !== 4
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
        # [tower3d_noseis_5,tower3d_seis_5,tower3d_seis_lo_5,tower3d_seis_hi_5],
        [tower3d_noseis_10,tower3d_seis_10,tower3d_seis_lo_10,tower3d_seis_hi_10],
        [tower3d_noseis_15,tower3d_seis_15,tower3d_seis_lo_15,tower3d_seis_hi_15],
        [tower3d_noseis_20,tower3d_seis_20,tower3d_seis_lo_20,tower3d_seis_hi_20],
    ];
    ss=20
)

CM.activate!(); plot_tower3d_dpl_traj(
    [
        [tower3d_noseis_5,tower3d_seis_5,tower3d_seis_lo_5,tower3d_seis_hi_5],
        [tower3d_noseis_10,tower3d_seis_10,tower3d_seis_lo_10,tower3d_seis_hi_10],
        [tower3d_noseis_15,tower3d_seis_15,tower3d_seis_lo_15,tower3d_seis_hi_15],
        [tower3d_noseis_20,tower3d_seis_20,tower3d_seis_lo_20,tower3d_seis_hi_20],
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
            RB.goto_step!(bot0,stepj;actuate=true)
            RB.goto_step!(bot1,stepj;actuate=true)
            plot_traj!(bot1;
                AxisType=Axis3,
                fig=gridi,
                # cablecolor=:slategrey,
                atsteps=[stepj],
                showinit=true,
                tgini=deepcopy(bot0.st),
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
                    rbs = RB.get_bodies(tgob[])
                    # for bodyid in 7:10
                    for bodyid in [10]
                        if bodyid in 7:8
                            pid = 4
                        else
                            pid = 3
                        end
                        meshscatter!(ax,[rbs[bodyid].state.loci_states[pid]], 
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

tower3d_seis   = RB.Robot(deepcopy(tower3dbot0.st),(actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
RB.solve!(RB.DynamicsProblem(tower3d_seis,(x)->dynfuncs(x;actuate=true,gravity=true)),
            RB.Zhong06(),(sysstate,st)->pres3d!(sysstate,st;T=0.5);
            dt=1e-3,tspan=(0.0,15.0),ftol=1e-14)

tower3dact_i = RB.Robot(
                    tower3d(;k=500.0,c=2.0,ijkl=1).st,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
tower3dact_j = RB.Robot(
                    tower3d(;k=500.0,c=2.0,ijkl=2).st,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
tower3dact_k = RB.Robot(
                    tower3d(;k=500.0,c=2.0,ijkl=3).st,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))
tower3dact_l = RB.Robot(
                    tower3d(;k=500.0,c=2.0,ijkl=4).st,
                    (actuators=[make_pres_actor(μ0,μ1,0.0,10.0)],))

dt = 1e-3
RB.solve!(RB.DynamicsProblem(tower3dact_i,(x)->dynfuncs(x;actuate=true,gravity=true)),
            RB.Zhong06(),(sysstate,st)->pres3d!(sysstate,st;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)
RB.solve!(RB.DynamicsProblem(tower3dact_j,(x)->dynfuncs(x;actuate=true,gravity=true)),
            RB.Zhong06(),(sysstate,st)->pres3d!(sysstate,st;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)
RB.solve!(RB.DynamicsProblem(tower3dact_k,(x)->dynfuncs(x;actuate=true,gravity=true)),
            RB.Zhong06(),(sysstate,st)->pres3d!(sysstate,st;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)
RB.solve!(RB.DynamicsProblem(tower3dact_l,(x)->dynfuncs(x;actuate=true,gravity=true)),
            RB.Zhong06(),(sysstate,st)->pres3d!(sysstate,st;T=1.0);dt=dt,tspan=(0.0,10.0),ftol=1e-14)

function plot_compare_ijkl(boti,botj,botk,botl)
    irb9rg = RB.get_trajectory!(boti,9,0)
    jrb9rg = RB.get_trajectory!(botj,9,0)
    krb9rg = RB.get_trajectory!(botk,9,0)
    lrb9rg = RB.get_trajectory!(botl,9,0)
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

me_series = RB.mechanical_energy!(tower3dact;actuate=true,gravity=true)
lines(tower3dact.traj.t,me_series.E)
lines(tower3dact.traj.t,me_series.T)
lines(tower3dact.traj.t,me_series.V)

newspine = newspine3d(2,)

plot_traj!(newspine;showground=false)

#--sldkr 
gravity = false
n = 2
newbridge_tag = bridge3d(;n)
plot_traj!(newbridge_tag)

μ = RB.inverse_for_restlength(newbridge_tag,newbridge_tag;fmin=10.0,gravity)
k,μ = RB.optimize_for_stiffness_and_restlen(newbridge_tag;fmin=100.0)
μ
function pres3d!(sysstate,st;Td=10.0,q_ini,q_tag)
    (;t,q̃,q̃̇,q̃̈) = sysstate
    (;connectivity) = st
    (;sys_pres_coords_idx) = connectivity.indexed
    q̃_ini = @view q_ini[sys_pres_coords_idx]
    q̃_tag = @view q_tag[sys_pres_coords_idx]
    Δq̃ = q̃_tag .- q̃_ini
    τ = t/Td
    q̃ .= q̃_ini .+ τ.*Δq̃
    q̃̇ .= Δq̃./Td
    q̃̈ .= 0.0
end
Td = 10.0

F̌ = RB.build_F̌(newbridge_tag.st,2n+1,2*(2n+1)+1,SVector(0.0,0.0,-1.0))

newbridge_tag_dyn = deepcopy(newbridge_tag)

RB.set_restlen!(newbridge_tag_dyn.st,μ)
RB.update!(newbridge_tag_dyn.st)
RB.get_cables_tension(newbridge_tag_dyn.st) |> extrema
RB.GDR!(newbridge_tag_dyn;gravity)
plot_traj!(newbridge_tag_dyn)
RB.set_new_initial!(newbridge_tag_dyn,newbridge_tag_dyn.traj.q[end])
RB.check_static_equilibrium_output_multipliers(newbridge_tag_dyn.st;gravity)
RB.undamped_eigen(newbridge_tag_dyn.st;gravity)
RB.undamped_eigen!(newbridge_tag_dyn;gravity)

plot_traj!(newbridge_tag_dyn)
newbridge_tag_dyn.traj.t

tend = 10.0
dstep = 45
dt = 1e-4
RB.solve!(
    RB.DynamicsProblem(
        newbridge_tag_dyn,
        (x)->dynfuncs(x;actuate=false,gravity,(Fˣ!)=(F,t) -> begin
            F .+= F̌*sin(t)*1e2
        end)
    ),
    RB.Zhong06();
    dt=1e-3,
    # tspan=(0.0,dt*(dstep-1)),
    tspan = (0.0,tend),
    ftol=1e-13,
    exception=false
)
plot_traj!(newbridge_tag_dyn)

bridge_deploy = RB.Robot(
    deepcopy(newbridge_tag_deploy.st),
    (actuators = [make_pres_actor(μ_tag,μ_tag,0.0,Td)],)
)
    
RB.check_static_equilibrium_output_multipliers(bridge_deploy.st;gravity=true)
q_tag = RB.get_coords(newbridge_tag.st)

RB.solve!(
    RB.DynamicsProblem(
        bridge_deploy,(x)->dynfuncs(x;actuate=false,gravity=true)),
        RB.Zhong06(),
        # (
        #     prescribe! = (sysstate,st)->pres3d!(sysstate,st;Td,q_ini,q_tag=q_ini),
        #     actuate! = nothing
        # );
        dt=1e-4,tspan=(0.0,2.0),ftol=1e-7
)
plot_traj!(bridge_deploy)
## lander
newlander = lander()
plot_traj!(newlander;)

k,μ=RB.optimize_for_stiffness_and_restlen(newlander;fmin=1.0,gravity=false)
k,μ

newlander_opt = lander(;k)
RB.set_restlen!(newlander_opt.st,μ)
RB.update!(newlander_opt.st)
RB.get_cables_tension(newlander_opt.st) |> extrema
RB.check_static_equilibrium_output_multipliers(newlander_opt.st;gravity=false)
RB.undamped_eigen(newlander_opt.st;gravity=false)
RB.undamped_eigen!(newlander_opt;gravity=false)
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

μ = RB.inverse_for_restlength(newst,newst;gravity,fmin=0.0)

@show newst_opt.st.state.system

k,μ=RB.optimize_for_stiffness_and_restlen(newst;fmin = 10.0,gravity)
k .= 1000.0
RB.set_restlen!(newst_opt.st,μ)
RB.update!(newst_opt.st;gravity)

f = RB.get_cables_tension(newst_opt.st) 

RB.check_static_equilibrium_output_multipliers(newst_opt.st;gravity=true)
RB.undamped_eigen(newst_opt.st;gravity=true)
RB.undamped_eigen!(newst_opt;gravity=true)
plot_traj!(newst_opt;showmesh=false,)

newst_folded = scissor_tower(θ = π/5;k)
plot_traj!(newst_folded;showmesh=false,showlabels=true,showpoints=true,showground=false)


k,μ=RB.optimize_for_stiffness_and_restlen(newst_folded;gravity=true)

μ_folded = RB.inverse_for_restlength(newst_folded,newst_folded;gravity=false,fmin=0.0)
RB.set_restlen!(newst_folded.st,μ_folded)
RB.update!(newst_folded.st)
RB.get_cables_tension(newst_folded.st) |> extrema


#-- without outer and gravity, to obtain rest length
gravity = false
m = 3
α = 2π/m
θ = 1.25α
n = 4
b = 1.4
r = 0.4*sqrt(2)
newembed1 = embed3d(;
    r1= 0.3*sqrt(2),
    r,b,m,α,θ,n
)

q = RB.get_coords(newembed1.st)
Φ = RB.make_cstr_function(newembed1.st)
Φ(q) |> norm

μ1 = RB.inverse_for_restlength(
    newembed1,newembed1;gravity,fmin=100.0
)

# (optional) check μ1 
RB.set_restlen!(newembed1.st,μ1)
RB.update!(newembed1.st)
f = RB.get_cables_tension(newembed1.st) 
f |> extrema
RB.check_static_equilibrium_output_multipliers(newembed1.st;gravity)
RB.undamped_eigen(newembed1.st;gravity)
RB.undamped_eigen!(newembed1;gravity)
RB.reset!(newembed1)
#-- without outer, see see gravity
RB.GDR!(newembed1;gravity=true)

plot_traj!(newembed1)

#-- folded without outer and gravity, to obtain rest length
newembed0 = embed3d(;
    r1= 0.6*sqrt(2),
    r,b,m,α,θ,n
)

plot_traj!(newembed0)

μ0 = RB.inverse_for_restlength(
    newembed0,newembed0;gravity,fmin=100.0
)

#-- with outer, using the folded and unfolded 
outer = true
newembed0_outer = embed3d(;
    r1= 0.6*sqrt(2),
    r,b,m,α,θ,n,outer
)
plot_traj!(newembed0_outer)
μ0o = RB.get_cables_restlen(newembed0_outer)
μ0o[begin:end-n*m] .= μ0
μ1o = deepcopy(μ0o)
μ1o[begin:end-n*m] .= μ1
RB.set_restlen!(newembed0_outer.st,μ0o)
RB.GDR!(newembed0_outer)
plot_traj!(newembed0_outer;)
RB.set_new_initial!(newembed0_outer,newembed0_outer.traj.q[end])
RB.update!(newembed0_outer.st)
RB.check_static_equilibrium_output_multipliers(newembed0_outer.st;gravity)
RB.undamped_eigen(newembed0_outer.st;gravity)
RB.undamped_eigen!(newembed0_outer;gravity)
friction_coefficients = [μ0o,μ1o]

Td = 15.0
tend = 20.0

friction_coefficients = [
    begin
        μ = deepcopy(μ0o)
        μ[begin:3m*j] .= μ1[begin:3m*j]
        μ
    end
    for j = 0:n
]

function make_multi_stages_pres_actor(friction_coefficients;start = 0.0, stop = 10.0, len = 2, )
    nμ = length(friction_coefficients[begin])

    function itp(t)
        scaled_itps = extrapolate(
            Interpolations.scale(
                interpolate(
                    reduce(hcat,friction_coefficients),
                    (NoInterp(),BSpline(Linear()))
                    # (NoInterp(),BSpline(Quadratic(Flat(OnGrid()))))
                ),
                1:nμ, range(;start, stop, length = len,)
            ),
            (Throw(),Flat())
        )
        [scaled_itps(j,t) for j in 1:nμ]
    end

    RB.PrescribedActuator(
        1,
        RB.ManualActuator(1,collect(1:nμ),zeros(nμ),RB.Uncoupled()),
        itp
    )
end
make_multi_stages_pres_actor(friction_coefficients;start=0.0,stop=Td,len=n+1)

newembed1_outer_deploy = RB.Robot(
    deepcopy(newembed1_outer.st),
    (actuators=[make_multi_stages_pres_actor(friction_coefficients;start=0.0,stop=Td,len=length(friction_coefficients))],)
)

RB.solve!(
    RB.DynamicsProblem(
        newembed1_outer_deploy,
        (x)->dynfuncs(x;actuate=true,gravity=true)
    ),
    RB.Zhong06();
    dt=5e-3,
    tspan=(0.0,tend),
    ftol=1e-6,
    exception=false
)

with_theme(theme_pub;
        figure_padding = (2fontsize,fontsize,0,fontsize),
        Axis3 = (
            azimuth = 3.995530633326985,
            elevation = 0.18269908169872415
        )
    ) do
    plot_traj!(
        newembed1_outer_deploy;
        figsize = (820,1280),
        AxisType=Axis3,
        xlims = (-2,2),
        ylims = (-2,2),
        zlims = (-1e-3,7),
        showinfo = false,
        doslide = false,
        dorecord = true,
        actuate = true,
        showpoints = false,
        showlabels = false,
        figname = "newembed1_outer_deploy.mp4"
    )
end