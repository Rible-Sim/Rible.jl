using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using TypeSortedCollections
using Rotations
# import GeometricalPredicates as GP
using CoordinateTransformations
using OffsetArrays
import DifferentialEquations as DE
using BenchmarkTools
using RecursiveArrayTools
using Interpolations
# using CubicSplines
# import FLOWMath
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using TexTables
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
using Revise
import TensegrityRobots as TR
cd("examples/test_prescribed")
include("../vis.jl")
include("../analysis.jl")
include("../dyn.jl")
include("dyn.jl")
include("def.jl")
include("def3d.jl")

includet("../vis.jl")
includet("../analysis.jl")
includet("../dyn.jl")
includet("dyn.jl")
includet("def.jl")
includet("def3d.jl")
# figdir::String = raw"C:\Users\luo22\OneDrive\Papers\DynamicTensegrity\CS"
figdir::String =raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\dyn"

fontsize = 8 |> pt2px
markersize = fontsize
linewidth = 0.5 |> pt2px
cablewidth = 0.75 |> pt2px
barwidth = 1.5 |> pt2px
cw = 252 |> pt2px
tw = 522 |> pt2px

# fontsize = 10 |> pt2px
# cw = lw = tw = 455.24411 |> pt2px

# set_theme!(;font = "Nimbus Rom No9 L", fontsize)

## 1st example, validation and verification
## comparison with Adams
obot = one_bar_one_tri()

prob = TR.SimProblem(obot,(x)->dynfuncs(x;gravity=true))

dt = 1e-3

TR.solve!(prob,TR.Zhong06();dt,tspan=(0.0,4.0),ftol=1e-14)

plot_traj!(obot;showmesh=false,showinfo=false,showground=false)
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

    fig_width = cw
    fig_height = cw
    with_theme(theme_pub;
            resolution = (0.8tw,0.5tw),
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
TR.solve!(prob_zhong,TR.Zhong06();dt,tspan=(0.0,500.0),ftol=1e-14)
TR.solve!(prob_alpha,TR.Alpha(0.7);dt,tspan=(0.0,500.0),ftol=1e-14)

plot_traj!(obot_zhong)

function plot_energy!(bot_zhong,bot_alpha,figname=nothing)
    ms_zhong = TR.mechanical_energy!(bot_zhong;gravity=true)
    ms_alpha = TR.mechanical_energy!(bot_alpha;gravity=true)
    fig_width = 0.9columnwidth
    fig_height = 0.35columnwidth
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


## 2nd example, 2D tower
## Natural frequency comparison
tower2dbot = tower2d(;k=100.0,ratio=0.9,ratio1=0.9,slack=false)
TR.undamped_eigen!(tower2dbot;scaling=0.035)
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
    (;traj,tg) = bot
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
                    # visible = false,
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
        TR.update_rigids!(tg,shapes[1])
        plot_tower2d!(axs[1],tg;isref=true)
        plot_tower2d!(axs[2],tg;isref=true)
        plot_tower2d!(axs[3],tg;isref=true)
        plot_tower2d!(axs[4],tg;isref=true)
        plot_tower2d!(axs[5],tg;isref=true)

        TR.update_rigids!(tg,shapes[2])
        plot_tower2d!(axs[1],tg)
        TR.update_rigids!(tg,shapes[3])
        plot_tower2d!(axs[2],tg)
        TR.update_rigids!(tg,shapes[4])
        plot_tower2d!(axs[3],tg)
        TR.update_rigids!(tg,shapes[5])
        plot_tower2d!(axs[4],tg)
        TR.update_rigids!(tg,shapes[6])
        plot_tower2d!(axs[5],tg)
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
            TR.undamped_eigen!(bot)
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

function pres2d!(sysstate,tg;aux=true)
    (;t,q̃,q̃̇,q̃̈) = sysstate
    (;bodies,connectivity) = tg
    (;mem2syspres) = connectivity.indexed
    amp = 0.01
    p = 6π
    foreach(bodies) do rb
        rbid = rb.prop.id
        if rbid == 1
            q̃[mem2syspres[rbid]] .= [     amp*sin(p*t),0]
            q̃̇[mem2syspres[rbid]] .= [   p*amp*cos(p*t),0]
            q̃̈[mem2syspres[rbid]] .= [-p^2*amp*sin(p*t),0]
        end
        if rbid == 2
            q̃[mem2syspres[rbid]] .= [ 0.1+amp*sin(p*t),0]
            q̃̇[mem2syspres[rbid]] .= [   p*amp*cos(p*t),0]
            q̃̈[mem2syspres[rbid]] .= [-p^2*amp*sin(p*t),0]
        end
        if aux
            if rbid == 8
                q̃[mem2syspres[rbid][3:4]] .= [-0.1+amp*sin(p*t),0]
                q̃̇[mem2syspres[rbid][3:4]] .= [   p*amp*cos(p*t),0]
                q̃̈[mem2syspres[rbid][3:4]] .= [-p^2*amp*sin(p*t),0]
                q̃[mem2syspres[rbid][5:6]] .= [ 0.2+amp*sin(p*t),0]
                q̃̇[mem2syspres[rbid][5:6]] .= [   p*amp*cos(p*t),0]
                q̃̈[mem2syspres[rbid][5:6]] .= [-p^2*amp*sin(p*t),0]
            end
        end
    end
end


ratio1 = 0.998
dt = 1e-3
## undamped noslack
tower2dbot_undamped_noslack = tower2d(;k=100.0,c=0.0,ratio=0.95,ratio1,slack=false)
TR.solve!(TR.SimProblem(tower2dbot_undamped_noslack,(x)->dynfuncs(x;gravity=true)),
          TR.Zhong06(),
          (
            prescribe! = (sysstate,tg)->pres2d!(sysstate,tg;aux=true),
            actuate! = nothing
          );
          dt,tspan=(0.0,2.0),ftol=1e-14,verbose=true,exception=false)

adams_undamped_noslack = parse_Adams_dynamic("tower2d_undamped_noslack.xml")

plot_traj!(tower2dbot_undamped_noslack)


## damped noslack
tower2dbot_damped_noslack = tower2d(;k=100.0,c=1.0,ratio=0.95,ratio1,slack=false)
TR.solve!(TR.SimProblem(tower2dbot_damped_noslack,(x)->dynfuncs(x;gravity=true)),
          TR.Zhong06(),
          (
            prescribe! = (sysstate,tg)->pres2d!(sysstate,tg;aux=true),
            actuate! = nothing
          );
          dt,tspan=(0.0,2.0),ftol=1e-14)

adams_damped_noslack = parse_Adams_dynamic("tower2d_damped_noslack.xml")

## undamped slack
tower2dbot_undamped_slack = tower2d(;k=100.0,c=0.0,ratio=0.95,ratio1,slack=true)
f = TR.get_cables_tension(tower2dbot_undamped_slack)
l = TR.get_cables_len(tower2dbot_undamped_slack)
μ = TR.get_cables_restlen(tower2dbot_undamped_slack)
k = TR.get_cables_stiffness(tower2dbot_undamped_slack)
x1,y1 = approx_slack(;k=k[1],μ=μ[1],l0=l[1],filename="xy1")
x2,y2 = approx_slack(;k=k[3],μ=μ[3],l0=l[3],filename="xy2")
x3,y3 = approx_slack(;k=k[5],μ=μ[5],l0=l[5],filename="xy3")
x4,y4 = approx_slack(;k=k[11],μ=μ[11],l0=l[11],filename="xy4")

x0,y0 = approx_slack(;k=100.0,μ=0.0,l0=0.0)
function plot_data_points(x,y,figname=nothing)
    with_theme(theme_pub;
        resolution = (cw,0.8cw),
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

dy1 = map(xs1) do x
    CubicSplines.gradient(spline1, x, 1)
end
lines(xs1,dy1)

spline1 = FLOWMath.Akima(x1, y1)
ys1 = spline1(xs1)
lines(xs1, ys1)
scatter!(x1,y1)
dydx1 = FLOWMath.gradient(spline1, xs1)
lines(xs1, dydx1)

scatter(x1,y1)
TR.solve!(TR.SimProblem(tower2dbot_undamped_slack,(x)->dynfuncs(x;gravity=true)),
          TR.Zhong06(),
          (
            prescribe! = (sysstate,tg)->pres2d!(sysstate,tg;aux=true),
            actuate! = nothing
          );
          dt,tspan=(0.0,3.0),ftol=1e-14)

adams_undamped_slack = parse_Adams_dynamic("tower2d_undamped_slack.xml")
plot_traj!(tower2dbot_undamped_slack)

## undamped slack with fix aux
tower2dbot_undamped_slack_noaux = tower2d(;k=100.0,c=0.0,ratio=0.95,ratio1,slack=true)
TR.solve!(TR.SimProblem(tower2dbot_undamped_slack_noaux,(x)->dynfuncs(x;gravity=true)),
          TR.Zhong06(),
          (
            prescribe! = (sysstate,tg)->pres2d!(sysstate,tg;aux=false),
            actuate! = nothing
          );
          dt,tspan=(0.0,3.0),ftol=1e-14)
adams_undamped_slack_aux = parse_Adams_dynamic("tower2d_undamped_slack_noaux.xml")

function plot_tower2d_point_traj(ba,figname=nothing;at=1.23)
    fig_width = 0.9testwidth
    fig_height = 0.75columnwidth
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
            rg = get_trajectory!(bot,3,0)
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
    fig_width = 0.95columnwidth
    fig_height = 0.7columnwidth
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
            (;cables) = bot.tg.tensiles
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
    fig_width = 1.5columnwidth
    fig_height = 0.45columnwidth
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
            TR.goto_step!(bots[i],1)
            plot_tower2d!(ax,bots[i].tg;isref=true,markit=true)
            istep = findfirst((x)->x>=at, bots[i].traj.t)
            @show istep
            TR.goto_step!(bots[i],istep)
            # @show TR.get_cables_tension(bots[i].tg)
            plot_tower2d!(ax,bots[i].tg,markit=true)
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

ms_tower2d = TR.mechanical_energy!(tower2dbot_undamped_noslack;gravity=true)


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
tower3dbot0 = tower3d(;k=500.0,c=2.0)
TR.get_cables_stiffness(tower3dbot0)
GM.activate!(); plot_traj!(tower3dbot0;showmesh=false,fontsize=20)

ℓ0 = TR.get_cables_len(tower3dbot0.tg)
ℓ0[begin:3:end]
μ0 = zero(ℓ0)
μ0[13:15] .= 0.03
μ0[16:18] .= 0.10
μ0[19:21] .= 0.03

B,F̃ = TR.build_inverse_statics_for_restlength(tower3dbot0.tg,tower3dbot0.tg;gravity=true)
F̃ .-= B*μ0
μx,μb = TR.get_solution_set(B[:,rest],F̃)
μ0[rest] .= μx
@show μ0[begin:3:end]
(ℓ0 - μ0)[begin:3:end]
TR.set_restlen!(tower3dbot0.tg,μ0)
TR.update!(tower3dbot0.tg)
TR.check_static_equilibrium_output_multipliers(tower3dbot0.tg;gravity=true)
TR.undamped_eigen(tower3dbot0.tg;gravity=true)
TR.undamped_eigen!(tower3dbot0;gravity=true)
tower3dbot0.traj.t./(2π)


tower3dbot1 = tower3d(;k=500.0,c=2.0,d = 0.1*√2, r2 = 0.07)
TR.get_cables_stiffness(tower3dbot1)[begin:3:end]
plot_traj!(tower3dbot1;fontsize=20)
μ1 = zero(TR.get_cables_restlen(tower3dbot1.tg))
μ1[13:21] .= μ0[13:21]
B,F̃ = TR.build_inverse_statics_for_restlength(tower3dbot1.tg,tower3dbot1.tg;gravity=true)
F̃ .-= B*μ1
μx,μb = TR.get_solution_set(B[:,rest],F̃)
μ1[rest] .= μx
ℓ1 = TR.get_cables_len(tower3dbot1.tg)
(ℓ1 - μ1)[begin:3:end]
TR.set_restlen!(tower3dbot1.tg,μ1)
TR.update!(tower3dbot1.tg)
TR.check_static_equilibrium_output_multipliers(tower3dbot1.tg;gravity=true)
TR.undamped_eigen(tower3dbot1.tg;gravity=true)
TR.undamped_eigen!(tower3dbot1;gravity=true)
tower3dbot1.traj.t./(2π)

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



function plot_rest_lengths(bots)
    μs = [TR.get_cables_restlen(bot.tg)[begin:3:end] for bot in bots]
    fig_width = 0.8columnwidth
    fig_height = 0.5columnwidth
    fig = Figure(;resolution=(fig_width,fig_height))
    ax1 = Makie.Axis(fig[1,1]; xlabel = L"\mathrm{Cable~Group}", ylabel = L"\mathrm{μ}~(\mathrm{m})")
    scatterlines!(ax1, μs[1]; marker = '+', markersize, linewidth, color = :blue, linestyle = :dashdot, label = "initial")
    scatterlines!(ax1, μs[2]; marker = '×', markersize, linewidth, color = :red,  linestyle = :dashdot, label = "target")
    axislegend(ax1,position=:rt)
    ylims!(ax1,0,0.2)
    ax1.xticks = collect(1:8)
    ax1.xlabelpadding = -15
    fig
end

fig_rest_lengths = plot_rest_lengths([tower3dbot0,tower3dbot1])
GM.activate!()
CM.activate!()
CM.save(texroot*raw"\OneDrive - 中山大学\Papers\DynamicTensegrity\CS\images\rest_lengths.pdf", fig_rest_lengths)

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
            resolution = (cw,0.5cw),
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1]; xlabel = L"\mathrm{Mode~(Exclude~duplicates)}", ylabel = L"\mathrm{Frequency}~(\mathrm{Hz})")
        scatterlines!(ax1, unqfrqs[1]; marker = '+', markersize = 1.5fontsize, linewidth, linestyle = :dashdot, color = :blue, label = "initial")
        scatterlines!(ax1, unqfrqs[2]; marker = '×', markersize = 1.5fontsize, linewidth, linestyle = :dashdot, color = :red, label = "target")
        ax1.yticks = collect(0:15)
        ax1.xticks = collect(1:n)
        # ax1.xlabelpadding = -15
        # axislegend(ax1; position = :lt,)
        fig[1,2] = Legend(fig,ax1)
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
            dt=1e-3,tspan=(0.0,15.0),verbose=false,ftol=1e-14)
# plot_tower3d_traj(tower3dbot0_nodpl)
plot_traj!(tower3dbot0_nodpl)

tower3dbot1_nodpl = deepcopy(tower3dbot1)
TR.solve!(TR.SimProblem(tower3dbot1_nodpl,(x)->dynfuncs(x;gravity=true)),
            TR.Zhong06(),
            (
                prescribe! = pres3d!,
                actuate! = nothing
            );
            dt=1e-3,tspan=(0.0,15.0),ftol=1e-14)
# plot_tower3d_traj(tower3dbot1_nodpl)
plot_traj!(tower3dbot1_nodpl)

function plot_tower3d_seis_traj(bots,figname=nothing)
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
            (;t) = bot.traj
            ax1 = Makie.Axis(fig[i,1];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta x~(\mathrm{m})")
            ax2 = Makie.Axis(fig[i,2];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta y~(\mathrm{m})")
            ax3 = Makie.Axis(fig[i,3];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta z~(\mathrm{m})")
            scalings = [-2,-3,-3]
            rb7r4 = get_trajectory!(bot,7,4)
            rb8r4 = get_trajectory!(bot,8,4)
            rb9r3 = get_trajectory!(bot,9,3)
            rb10r3 = get_trajectory!(bot,10,3)
            lines!(ax1,t,(rb7r4[1,:].-rb7r4[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{7,1}")
            lines!(ax1,t,(rb8r4[1,:].-rb8r4[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{8,1}")
            lines!(ax1,t,(rb9r3[1,:].-rb9r3[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{9,1}")
            lines!(ax1,t,(rb10r3[1,:].-rb10r3[1,1])./10.0^scalings[1];label=L"\mathbf{r}_{10,1}")
            ylims!(ax1,-5,5)
            ax1.yticks = collect(-5:2.5:5)

            lines!(ax2,t,(rb7r4[2,:].-rb7r4[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb8r4[2,:].-rb8r4[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb9r3[2,:].-rb9r3[2,1])./10.0^scalings[2];)
            lines!(ax2,t,(rb10r3[2,:].-rb10r3[2,1])./10.0^scalings[2];)

            lines!(ax3,t,(rb7r4[3,:].-rb7r4[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb8r4[3,:].-rb8r4[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb9r3[3,:].-rb9r3[3,1])./10.0^scalings[3];)
            lines!(ax3,t,(rb10r3[3,:].-rb10r3[3,1])./10.0^scalings[3];)
            @show rb10r3[3,1]
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
    fig = Figure(;resolution=(fig_width,fig_height),figure_padding=(10,0,0,0))
    dt = 1e-3
    ts = [
        [2.5],
        [3.5]
    ]
    nbot = length(bots)
    for i in 1:nbot
        ax = Axis3(fig[1,i],aspect=:data,)
        boti = bots[i]
        TR.goto_step!(boti,1;)
        plot_tower3d!(ax,boti.tg,cablecolor=:slategrey,barcolor=:slategrey,tetcolor=:slategrey)
        for t in ts[i]
            stepj = round(Int,t/dt)
            TR.goto_step!(boti,stepj;)
            # @show TR.get_cables_tension(boti)
            plot_tower3d!(ax,boti.tg;markit = true)
        end
    end

    for (ilabel,label) in enumerate(alphabet[1:nbot])
        Label(fig.layout[1, ilabel, Bottom()],
            "($label)",
            fontsize = fontsize,
            font = "CMU Serif Bold",
            padding = (-110, 0, 0, 80),
            halign = :center,
            valign = :top)
    end
    titles = ["initial","target"]
    for (ilabel,label) in enumerate(alphabet[1:nbot])
        Label(fig.layout[1, ilabel, Bottom()],
            latexstring("\\mathrm{$(titles[ilabel])}"),
            fontsize = fontsize,
            padding = (50, 0, 0, 90),
            halign = :center,
            valign = :top)
    end
    colgap!(fig.layout,0)
    savefig(fig,figname)
    fig
end

GM.activate!(); plot_tower3d_vis_nodpl([tower3dbot0_nodpl,tower3dbot1_nodpl])
CM.activate!(); plot_tower3d_vis_nodpl([tower3dbot0_nodpl,tower3dbot1_nodpl], "tower3d_vis_nodpl")

## deploy
Td = [5.0,10.0,15.0,20.0]
function do_noseis(;Td,tend=25.0)
    bot   = TR.TensegrityRobot(deepcopy(tower3dbot0.tg),(actuators=[make_pres_actor(μ0,μ1,0.0,Td)],))
    TR.solve!(TR.SimProblem(bot,(x)->dynfuncs(x;actuate=true,gravity=true)),
                TR.Zhong06();
                dt=1e-3,tspan=(0.0,tend),ftol=1e-14)
    bot
end
function do_seis(;ν,Td,tend=25.0)
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
tower3d_noseis_5 = do_noseis(;Td=Td[1])
tower3d_seis_5 = do_seis(;ν=0.5,Td=Td[1])
tower3d_seis_lo_5 = do_seis(;ν=0.2,Td=Td[1])
tower3d_seis_hi_5 = do_seis(;ν=1.0,Td=Td[1])
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

plot_traj!(tower3d_seis_20;actuate=true)
plot_traj!(tower3d_seis_lo;actuate=true)
plot_traj!(tower3d_seis_hi;actuate=true)
me_series = TR.mechanical_energy!(tower3d_seis;actuate=true,gravity=true)
lines(tower3d_seis.traj.t,me_series.E)
function plot_tower3d_dpl_traj(bots,figname=nothing)
    cg = cgrad(:seaborn_bright6, 6, categorical = true)
    with_theme(theme_pub;
            # resolution = (tw,0.8tw),
            resolution = (tw,1.25cw),
            palette = (color = cg, ),
            Lines = (cycle = [:color], linewidth = 2),
        ) do
        fig = Figure(;)
        ats = [16.5,10.5,14.0,19.5]
        Tds = [5.0,10.0,15.0,20.0]
        axs = [
            [
                Makie.Axis(fig[i,1];xlabel=L"t~(\mathrm{s})",ylabel=L"x~(\mathrm{m})"),
                Makie.Axis(fig[i,2];xlabel=L"t~(\mathrm{s})",ylabel=L"z~(\mathrm{m})"),
                Makie.Axis(fig[i,3];xlabel=L"t~(\mathrm{s})",ylabel=L"\Delta z~(\mathrm{m})")
            ]
            for i = 1:4
        ]
        scalings = [-1,-1,-2]

        for i = 1:4
            bot0, bot1, botlo, bothi = bots[i]
            at = ats[i]
            Td = Tds[i]
            ax1,ax2,ax3 = axs[i]
            t0 = bot0.traj.t
            step_collapse = findfirst((x)->x>=at,bot1.traj.t)
            step_Td = findfirst((x)->x>=Td,bot1.traj.t)
            step_Tdp1 = findfirst((x)->x>=Td+1,bot1.traj.t)
            tsteps = 1:step_collapse
            t1 = bot1.traj.t[tsteps]
            # bot0
            bot0rb10r3 = get_trajectory!(bot0,10,3)
            # bot1
            bot1rb10r3 = get_trajectory!(bot1,10,3)
            # botlo
            botlorb10r3 = get_trajectory!(botlo,10,3)
            # bothi
            bothirb10r3 = get_trajectory!(bothi,10,3)

            # x
            lines!(ax1,t0,bot0rb10r3[1,:]./10.0^scalings[1];linewidth)
            lines!(ax1,t0,botlorb10r3[1,:]./10.0^scalings[1];linewidth)
            lines!(ax1,t1,bot1rb10r3[1,tsteps]./10.0^scalings[1];linewidth)
            lines!(ax1,t0,bothirb10r3[1,:]./10.0^scalings[1];linewidth)
            ylims!(ax1,0,2.0)
            # ax1.yticks = collect(-5:2.5:5)
            xlims!(ax1,t0[begin],t0[end])

            # # y
            # lines!(ax2,t0,bot0rb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Static~ground}")
            # lines!(ax2,t0,botlorb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.2~(\mathrm{Hz})")
            # lines!(ax2,t0,bothirb10r3[2,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 1.0~(\mathrm{Hz})")
            # lines!(ax2,t1,bot1rb10r3[2,tsteps]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.5~(\mathrm{Hz})")
            # # ylims!(ax2,-0.5,7.5)
            # xlims!(ax2,t0[begin],t0[end])
            # ax2.yticks = collect(-1.0:1.0:7.5)

            # z
            lines!(ax2,t0,bot0rb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Static~ground}")
            lines!(ax2,t0,botlorb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.2~(\mathrm{Hz})")
            lines!(ax2,t1,bot1rb10r3[3,tsteps]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 0.5~(\mathrm{Hz})")
            lines!(ax2,t0,bothirb10r3[3,:]./10.0^scalings[2];linewidth,label=L"\mathrm{Seismic~ground},~\nu = 1.0~(\mathrm{Hz})")
            ylims!(ax2,3.5,5.5)
            xlims!(ax2,t0[begin],t0[end])

            rb10r3_z0 = 0.49497474683058335
            lines!(ax3,t0,(bot0rb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            lines!(ax3,t0,(botlorb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            lines!(ax3,t1,(bot1rb10r3[3,tsteps].-rb10r3_z0)./10.0^scalings[3];linewidth)
            lines!(ax3,t0,(bothirb10r3[3,:].-rb10r3_z0)./10.0^scalings[3];linewidth)
            if i == 1
                ylims!(ax3,-1.5,1.5)
            else
                ylims!(ax3,-.8,1.0)
            end
            xlims!(ax3,Td,Td+2)

            # if i == 2
            fig[5,1:3] = Legend(fig,ax2;orientation=:horizontal,)
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
        colsize!(fig.layout,3,Relative(0.15))
        savefig(fig,figname)
    end
end

GM.activate!(); plot_tower3d_dpl_traj(
    [
        [tower3d_noseis_5,tower3d_seis_5,tower3d_seis_lo_5,tower3d_seis_hi_5],
        [tower3d_noseis_10,tower3d_seis_10,tower3d_seis_lo_10,tower3d_seis_hi_10],
        [tower3d_noseis_15,tower3d_seis_15,tower3d_seis_lo_15,tower3d_seis_hi_15],
        [tower3d_noseis_20,tower3d_seis_20,tower3d_seis_lo_20,tower3d_seis_hi_20],
    ]
)


CM.activate!(); plot_tower3d_dpl_traj(
    [
        [tower3d_noseis_5,tower3d_seis_5,tower3d_seis_lo_5,tower3d_seis_hi_5],
        [tower3d_noseis_10,tower3d_seis_10,tower3d_seis_lo_10,tower3d_seis_hi_10],
        [tower3d_noseis_15,tower3d_seis_15,tower3d_seis_lo_15,tower3d_seis_hi_15],
        [tower3d_noseis_20,tower3d_seis_20,tower3d_seis_lo_20,tower3d_seis_hi_20],
    ],
    "tower3d_dpl_traj"
)


function plot_tower3d_vis(bot0,bot1,figname=nothing)
    fig_width = 0.95tw
    fig_height = 0.8cw
    fig = Figure(;resolution=(fig_width,fig_height),figure_padding=(0,0,0,0))
    dt = 1e-3
    ts = [1,5.5,10]

    for i in 1:3
        ax = Axis3(fig[1,i],aspect=:data,)
        stepi = round(Int,ts[i]/dt)
        TR.goto_step!(bot0,stepi;actuate=true)
        plot_tower3d!(ax,bot0.tg,cablecolor=:slategrey,barcolor=:slategrey,tetcolor=:slategrey)
        TR.goto_step!(bot1,stepi;actuate=true)
        plot_tower3d!(ax,bot1.tg;markit = true)
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
GM.activate!(); plot_tower3d_vis(tower3d_noseis_10,tower3d_seis_10,)
CM.activate!(); plot_tower3d_vis(tower3d_noseis_10,tower3d_seis_10,"tower3d_vis")

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
