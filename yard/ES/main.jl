#note -- preamble
using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/tensegrity.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
import Rible as RB

figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\DynamicTensegrity\ES" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/DynamicTensegrity/ES" #src
end #src

include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

## include("../bridge/bridge.jl"); includet("../bridge/bridge.jl")

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
RB.solve!(prob_zhong,RB.Zhong06();dt,tspan=(0.0,500.0),ftol=1e-10)
RB.solve!(prob_alpha,RB.Alpha(0.7);dt,tspan=(0.0,500.0),ftol=1e-10)

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

#---- simple
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
        size = (0.5cw,0.25cw),
        figure_padding = (2fontsize,0,0,0),
        Axis3 = (
            azimuth = 5.13952996194027,
            elevation = 0.12269999722606788
        )
    ) do 
    fig = Figure()
    rbs = RB.get_bodies(newsim_fixed)  
    body = rbs[2]
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
    mesh!(ax,build_mesh(body;update=false);shading = false,)
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
RB.solve!(
    RB.DynamicsProblem(newsim_fixed,(x)->dynfuncs(x;gravity)),
    RB.Zhong06();
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
RB.solve!(
    RB.DynamicsProblem(newsim_freed,(x)->dynfuncs(x;gravity)),
    RB.Zhong06();
    dt,tspan=(0.0,tend),ftol=1e-14
)

GM.activate!();plot_traj!(newsim_freed;)

GM.activate!(); with_theme(
        theme_pub;
        size = (0.6cw,0.3cw),
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
        size = (0.7cw,0.2cw),
        figure_padding = (fontsize/2,fontsize/2,0,0),
        Lines = (
            cycle = Cycle([:linestyle,:color];covary=true),
        )
    ) do 
    fig = Figure()
    (;t) = newsim_freed.traj
    r2p1_fixed = RB.get_trajectory!(newsim_fixed,1,1)
    r2p1_freed = RB.get_trajectory!(newsim_freed,1,1)
    # v2p1_fixed = get_mid_velocity!(newsim_fixed,1,1)
    # v2p1_freed = get_mid_velocity!(newsim_freed,1,1)
    T1_fixed = get_kinetic_energy!(newsim_fixed,1)
    T1_freed = get_kinetic_energy!(newsim_freed,1)
    E_fixed = RB.mechanical_energy!(newsim_fixed;gravity=false)
    E_freed = RB.mechanical_energy!(newsim_freed;gravity=false)
    # RB.kinetic_energy()
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
        size = (0.9cw,0.42cw)
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
μ = RB.inverse_for_restlength(newbridge,newbridge;fmin=1e4,gravity = false,eps_rel=1e-10)

Td = 10.0

F̌ = RB.build_F̌(newbridge.st,2n+1,2*(2n+1)+1,SVector(0.0,0.0,-1.0))

newbridge_modes = deepcopy(newbridge)

RB.set_restlen!(newbridge_modes.st,μ)
RB.update!(newbridge_modes.st; gravity=false)
RB.get_cables_tension(newbridge_modes.st) |> extrema
RB.GDR!(newbridge_modes;β=1e-4, res=1e-9,gravity= true)
plot_traj!(newbridge_modes)
RB.set_new_initial!(newbridge_modes,newbridge_modes.traj.q[end])
RB.check_static_equilibrium_output_multipliers(newbridge_modes.st;gravity=true)
RB.undamped_eigen(newbridge_modes.st;gravity=true)
RB.undamped_eigen!(newbridge_modes;scaling=-0.1,gravity=true)

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
        RB.solve!(
            RB.DynamicsProblem(
                bot,
                (x)->dynfuncs(x;actuate=false,gravity=true,(Fˣ!)=(F,t) -> begin
                    F .+= 2e4*(sin(ν*2π*t)+1).*F̌
                end)
            ),
            RB.Zhong06();
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
            size = (0.4cw,0.15cw),
            figure_padding = (0,0,0,fontsize/2),
            Lines = (
                cycle = Cycle([:linestyle,:color],covary=true),
                # linewidth = 1 |> pt2px,
            )
        ) do 
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=tlabel, ylabel = L"z~(\mathrm{m})")
        for (bot,ν) in zip(bots,νs)
            p11 = RB.get_trajectory!(bot,5,11)
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
μ0 = RB.inverse_for_restlength(tower3dbot0,tower3dbot0;gravity=true,scale=true,fmin=1.0)
RB.set_restlen!(tower3dbot0.st,μ0)
RB.update!(tower3dbot0.st)
RB.get_cables_tension(tower3dbot0.st) |> extrema
RB.check_static_equilibrium_output_multipliers(tower3dbot0.st;gravity=true)
RB.undamped_eigen(tower3dbot0.st;gravity=true)
RB.undamped_eigen!(tower3dbot0;gravity=true)
tower3dbot0.traj.t./(2π)
@show (tower3dbot0.traj.t./(2π))[2]

tower3dbot1 = tower3d(;k=500.0,c=2.0,d = 0.1*√2, r2 = 0.07, α=-π/12)
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
        size = (0.65cw,0.25cw),
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
            xlabel = "DistanceSpringDamper No.", 
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
            # rb9r3 = RB.get_trajectory!(bot,10,3,1:ss:length(bot.traj))
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
Td = [10.0]
function do_noseis(;Td,tend=15.0)
    bot   = RB.Robot(deepcopy(tower3dbot0.st),(actuators=[make_pres_actor(μ0,μ1,0.0,Td)],))
    RB.solve!(RB.DynamicsProblem(bot,(x)->dynfuncs(x;actuate=true,gravity=true)),
                RB.Zhong06();
                dt=1e-3,tspan=(0.0,tend),ftol=1e-14)
    bot
end
function do_seis(;ν,Td,tend=15.0)
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
me_series = RB.mechanical_energy!(tower3d_seis;actuate=true,gravity=true)
lines(tower3d_seis.traj.t,me_series.E)

function plot_tower3d_dpl_traj(bots,figname=nothing;
        ss=10,
    )
    cg = cgrad(:seaborn_bright6, 6, categorical = true)
    with_theme(theme_pub;
            # size = (cw,0.8cw),
            size = (0.8cw,0.20cw),
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
            bot0rb9r3 = RB.get_trajectory!(bot0,9,3,1:ss:length(bot0.traj))
            # bot1
            bot1rb9r3 = RB.get_trajectory!(bot1,9,3,tsteps)
            # botlo
            botlorb9r3 = RB.get_trajectory!(botlo,9,3,1:ss:length(botlo.traj))
            # bothi
            bothirb9r3 = RB.get_trajectory!(bothi,9,3,1:ss:length(bothi.traj))

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

#-- three 3-prism
prism3 = prism_modules(;n=1,p=1)
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
_,λ=  RB.check_static_equilibrium_output_multipliers(bot.st)
λ
k = RB.get_cables_stiffness(bot.st)
l = RB.get_cables_len(bot.st)
f = RB.get_cables_tension(bot)
q = RB.get_coords(bot.st)
q̌ = RB.get_free_coords(bot.st)
Ǎ = RB.make_cstr_jacobian(bot.st)(q)
# Ň_ = RB.nullspace(Ǎ)
# Ň = RB.modified_gram_schmidt(Ň_)
Ň = RB.make_intrinsic_nullspace(bot.st,q)
Q̃ = RB.build_Q̃(bot.st)
L̂ = RB.build_L̂(bot.st)
M̌ = RB.build_M̌(bot.st)
Ǩ = RB.build_Ǩ(bot.st,λ)
ℳ = transpose(Ň)*M̌*Ň
𝒦 = transpose(Ň)*Ǩ*Ň
# @show ℳ, 𝒦
# ω²,_ = RB.undamped_eigen(bot.st;)
ω²,Ξ = eigen(Symmetric(𝒦),Symmetric(ℳ))
ω²[1:6] .= 0
ω = sqrt.(ω²)
frq = ω./(2π)
M = transpose(Ξ)*ℳ*Ξ
M[findall((x)->abs(x)<1e-14,M)] .= 0.0
f_input = range(1.0,30.0;step=0.1)
function Hd(f)
    [
        ξ[3]*ξ[3]/(-(2π*f)^2+ω[i]^2)
        for (i,ξ) in enumerate(eachcol(Ξ))
    ] |> sum
end

function Hv(f)
    im*(2π*f)*Hd(f) |> abs
end

function Ha(f)
    -(2π*f)*Hd(f)
end
    
fig = Figure()
ax1 = Axis(fig[1,1];
    xlabel="Mode",
    ylabel="Frequency (Hz)",
)
barplot!(ax1,frq[1:18];width=0.1)
ax2 = Axis(fig[1,2];
    xscale=log10,
    xlabel="Frequency (Hz)",
    ylabel="Response"
)
lines!(ax2,f_input,Hd.(f_input))
ax3 = Axis(fig[2,1];
    xscale=log10,
    xlabel="Frequency (Hz)",
    ylabel="Response"
)
lines!(ax3,f_input,Hv.(f_input))
ax4 = Axis(fig[2,2];
    xscale=log10,
    xlabel="Frequency (Hz)",
    ylabel="Response"
)
lines!(ax4,f_input,Ha.(f_input))
Makie.DataInspector(fig)
fig

rank(Ň)
Ǎ*Ň |> norm
# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ

S,D = RB.static_kinematic_determine(ℬᵀ)
ns = size(S,2)
nk = size(D,2)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
Ǩa = RB.cstr_forces_jacobian(bot.st,λ)
𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)

Ǩm = RB.build_material_stiffness_matrix!(bot.st,q,k)
Ǩg = RB.build_geometric_stiffness_matrix!(bot.st,q,f)
vec𝒦ps = [
    begin
        si = S[:,i]
        # s = S\f
        # @show s
        λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
        # @show f,λ
        Ǩai = - RB.cstr_forces_jacobian(bot.st,λi)

        Ǩgi = RB.build_geometric_stiffness_matrix!(bot.st,q,si)

        𝒦pi = transpose(Ň)*(Ǩgi.+Ǩai)*Ň |> Symmetric 
        # vec𝒦pi = SymmetricPacked(𝒦pi).tri
        vec𝒦pi = vec(𝒦pi)
    end
    for i = 1:ns
]

mat𝒦ps = reduce(hcat,vec𝒦ps)

𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric
vec𝒦m = vec(𝒦m)
vecI = vec(Matrix(1.0I,size(𝒦m)))

𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 
𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
𝒦 = 𝒦m.+ 𝒦p |> Symmetric

vals_𝒦m,vecs_𝒦m = eigen(𝒦m)
sort(vals_𝒦m)
vm = vecs_𝒦m[:,1:nk]

vals_𝒦p,vecs_𝒦p = eigen(𝒦p)
sort(vals_𝒦p)

vals_𝒦,vecs_𝒦 = eigen(𝒦)
sort(vals_𝒦)

v = vecs_𝒦[:,1:6]
v'*𝒦*v


Ňv = Ň*nullspace(v')

r𝒦m = transpose(Ňv)*(Ǩm)*Ňv |> Symmetric 
# vecr𝒦m = SymmetricPacked(r𝒦m).tri
rd = nullspace(r𝒦m)
vecr𝒦m = vec(r𝒦m)

# vecI = SymmetricPacked(Matrix(1.0I,size(r𝒦m))).tri
vecI = vec(Matrix(1.0I,size(r𝒦m)))
r𝒦m |> issymmetric

r𝒦g = transpose(Ňv)*(Ǩg)*Ňv |> Symmetric 
r𝒦a = transpose(Ňv)*(Ǩa)*Ňv |> Symmetric 
r𝒦p = r𝒦g .+ r𝒦a
vals_r𝒦p,vecs_r𝒦p = eigen(r𝒦p)
@myshow sort(vals_𝒦p)

vals_rd𝒦pd,vecs_rd𝒦pd = eigen(rd'*r𝒦p*rd)
@myshow sort(vals_rd𝒦pd)

vecr𝒦ps = [
    begin
        si = S[:,i]
        # s = S\f
        # @show s
        λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
        # @show f,λ
        Ǩai = - RB.cstr_forces_jacobian(bot.st,λi)

        Ǩgi = RB.build_geometric_stiffness_matrix!(bot.st,q,si)

        r𝒦pi = transpose(Ňv)*(Ǩgi.+Ǩai)*Ňv |> Symmetric 
        # vecr𝒦pi = SymmetricPacked(r𝒦pi).tri
        vecr𝒦pi = vec(r𝒦pi)
    end
    for i = 1:ns
]

matr𝒦ps = reduce(hcat,vecr𝒦ps)

ᾱ = [1.0]
A = hcat(
    -Matrix(1.0I,ns,ns),
    ᾱ,
    zero(ᾱ)
)
b = [0.0]
nx = ns+2
result_max = RB.optimize_maximum_stiffness(matr𝒦ps,vecr𝒦m,vecI,A,b,nx)
σ_max = result_max.x[end-1]
ρ_max = result_max.x[end]

𝒦_max = 𝒦m + σ_max*reshape(mat𝒦ps*ᾱ,size(𝒦m))
vals_𝒦_max, vecs_𝒦_max = eigen(𝒦_max)

vals, vecs = eigen(𝒦_max - ρ_max*I)
@myshow vals

result_zero = RB.optimize_zero_stiffness(matr𝒦ps,vecr𝒦m,vecI,
    hcat(
        -Matrix(1.0I,ns,ns),
        ᾱ,
    ),
    [0.0],
    ns+1,
    result_max.x[1:end-1]
)
σ_zero = result_zero.x[end]

𝒦_zero = 𝒦m + σ_zero*reshape(mat𝒦ps*ᾱ,size(𝒦m))
vals_𝒦_zero, vecs_𝒦_zero = eigen(𝒦_zero)
ρ_zero = vals_𝒦_zero[1]
maxminmodes = hcat(
    vecs_𝒦_max[:,1],
    vecs_𝒦_zero[:,1:3],
)

 
σs = LinRange(0,40000,400)
rρs =  [
    begin
        r𝒦 = r𝒦m + σ*reshape(matr𝒦ps*ᾱ,size(r𝒦m))
        vals_r𝒦, vecs_r𝒦 = eigen(r𝒦)
        vals_r𝒦
    end
    for σ in σs
] |> VectorOfArray

size(rρs,1)
ρs =  [
    begin
        𝒦 = 𝒦m + σ*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦, vecs_𝒦 = eigen(𝒦)
        vals_𝒦[begin+1]
    end
    for σ in σs
]

with_theme(theme_pub;
        size = (0.3tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{(\mathrm{1})}"
    )
    lines!(ax,σs,rρs[1,:],)
    # xlims!(ax,0,5500)
    # ylims!(ax,-400,600)
    # for i = axes(rρs,1)
    #     lines!(ax,σs,rρs[i,:],)
    # end
    scatter!(
        ax,
        [σ_max,σ_zero],
        [ρ_max,ρ_zero]
    )
    text!([σ_max], [ρ_max], 
        text = [L"\rho_{(1),\mathrm{max}}"],
        align = (:center,:bottom),
        offset = (0, fontsize/4)
    )
    text!([σ_zero], [ρ_zero], 
        text = [L"\sigma_{\mathrm{max}}"],
        align = (:right,:center),
        offset = (-fontsize/2, 0)
    )
    # text!(x, y, text = string.(aligns), align = aligns)
    # savefig(fig,"superball_curve")
    fig
end

RB.undamped_eigen!(bot;scaling=0.1)
with_theme(theme_pub;
        fontsize = 6 |> pt2px,
        Axis3 = (
            azimuth = 7.045530633326983,
            elevation = 0.7926990816987238
        )
    ) do 
    plot_traj!(
        bot,
        AxisType=Axis3,
        gridsize=(2,3),
        atsteps = collect(2:7).+6,
        showinit = true,
        doslide=false,
        showinfo=false,
        showground=false,
        showpoints=false,
        showlabels=false,
        showcables=true,
        meshcolor=:black,
        xlims = (-0.7,0.7),
        ylims = (-0.7,0.7),
        # xlims = (-1.5,0.7),
        # ylims = (-1.2,1.2),
        zlims = (-1e-4,0.5),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font=:bold),
                (@sprintf "f = %.4G (Hz)" tt/(2π))
            )
        end,
    )
end
sqrt.(ω²[7:end])./(2π)
ω²,δq̌ =
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
