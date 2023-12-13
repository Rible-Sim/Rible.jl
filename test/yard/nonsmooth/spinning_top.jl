#-- deps
using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../test/yard/nonsmooth/deps.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
figdir::String = joinpath(pathof(RB),"../../tmp") #jl
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME" #src
end #src
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl
tw::Float64 = 455.24411 #src
scalefactor = 2 #nb
include(joinpath(pathof(RB),"../../examples/robots/spinningtop.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/spinningtop.jl"))

# ## First scenario
# Contact Surfaces
planes = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([0,0,1.0],[0,0,0.0]),
    ]
)

# Initial conditions
origin_position = [0,0,0.5]
R = RotX(π/24)
origin_velocity = [1.0,0.0,0.0]
Ω = [0.0,0.0,200.0]
# parameters
μ = 0.95
e = 0.5
# time 
tspan = (0.0,1.8)
h = 1e-4

topq = make_top(origin_position,R,origin_velocity,Ω;μ,e,loadmesh=true)

RB.solve!(
    RB.DynamicsProblem(
        topq,
        planes,
        RB.RestitutionFrictionCombined(
            RB.NewtonRestitution(),
            RB.CoulombFriction(),
        )
    ),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    tspan,
    dt=h,
    ftol=1e-12,
    maxiters=50,exception=true,verbose_contact=true
)

me = RB.mechanical_energy!(topq)
lines(me.E)
v5 = RB.get_velocity!(topq,1,5) 
lines(v5[3,:])
GM.activate!(;scalefactor)
plot_traj!(
    topq;
    showinfo=false,
    # rigidcolor=:white,
    showwire=true,
    showarrows=false,
)

topn = make_top(origin_position,R,origin_velocity,Ω,RB.NCF.NC;μ,e,loadmesh=true)
RB.solve!(
    RB.DynamicsProblem(
        topn,
        planes,
        RB.RestitutionFrictionCombined(
            RB.NewtonRestitution(),
            RB.CoulombFriction(),
        )
    ),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    tspan,
    dt=h,
    ftol=1e-14,
    maxiters=50,exception=false,verbose_contact=false
)

v5 = RB.get_velocity!(topn,1,5)
lines(v5[3,:])
plot_traj!(
    topn;
    showinfo=false,
    # rigidcolor=:white,
    showwire=false,
    showarrows=false,
    # figsize=(0.6tw,0.6tw)
)

GM.activate!(;scalefactor);with_theme(theme_pub;
        size = (1.0tw,0.3tw),
        Axis3 = (
            azimuth = 5.1155306333269825,
            elevation = 0.1426990816987241,
        ),
        Poly = (
            cycle = [:patchcolor=>:color],
            transparency = true,
        )
    ) do
    bot = topn
    (;t) = bot.traj
    rp5 = RB.get_trajectory!(bot,1,5)
    vp5 = RB.get_velocity!(bot,1,5)
    me = RB.mechanical_energy!(bot)
    steps = 1:1000:15000
    nstep = length(steps)
    alphas = fill(0.1,nstep)
    alphas[1:3] = [1,0.4,0.2]
    alphas[end] = 1
    cg = cgrad(:winter, nstep, categorical = true)
    fig = Figure()
    gd1 = fig[1,1] = GridLayout()
    gd2 = fig[1,2] = GridLayout()
    plot_traj!(
        bot;
        AxisType=Axis3,
        fig = gd1,
        doslide = false,
        showinfo=false,
        # rigidcolor=:white,
        showwire=false,
        showpoints=false,
        showlabels=false,
        showarrows=false,
        showcables=false,
        showmesh=false,
        showtitle=false,
        xlims = (-0.1,1.0),
        ylims = (-0.1,0.2),
        zlims = (-1e-6,0.6),
        sup! = (ax,tgob,sgi) -> begin
            RB.hidey(ax)
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                (;r,b,g) = cg[istep]
                RB.viz!(ax,Observable(suptg);
                    showlabels=false,
                    showarrows=false,
                    showpoints=false,
                    meshcolor=Makie.RGBA(r,b,g,alphas[istep])
                )
            end
            lines!(ax,rp5)
        end
    )
    ax1 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = L"\acute{v}_{n}~(\mathrm{m/s})",
    )
    lines!(ax1,t,vp5[3,:])
    ax2 = Axis(gd2[1:2,2],
        xlabel = tlabel,
        ylabel = "Energy (J)",
    )
    lines!(ax2,t,me.E, label="E")
    lines!(ax2,t,me.T, label="T")
    lines!(ax2,t,me.V.+6, label="V")
    axislegend(ax2,position=:rt)
    ax3 = Axis(gd2[2,1],
        xlabel = tlabel,
        ylabel = L"\acute{v}_{tb}~(\mathrm{m/s})",
    )
    lines!(ax3,t,[norm(vp5[2:3,i]) for i = 1:length(t)])
    RB.hidex(ax1)
    xlims!(ax1,0,1.8)
    xlims!(ax2,0,1.8)
    xlims!(ax3,0,1.8)
    Label(
        fig[1,1,TopLeft()],
        rich("($(alphabet[1]))",font=:bold),
        justification = :right,
        padding = (0,-5fontsize,0,0)
        
    )
    Label(
        gd2[1,1,TopLeft()],
        rich("($(alphabet[2]))",font=:bold)
    )
    Label(
        gd2[2,1,TopLeft()],
        rich("($(alphabet[3]))",font=:bold)
    )
    Label(
        gd2[1,2,TopLeft()],
        rich("($(alphabet[4]))",font=:bold)
    )
    colsize!(fig.layout,2,Fixed(0.50tw))
    colgap!(fig.layout,1.5fontsize)
    rowgap!(gd2,0)
    colgap!(gd2,1.5fontsize)
    # savefig(fig,"spinningtop_drop")
    fig
end

topq_longtime = deepcopy(topq)
RB.set_new_initial!(topq_longtime,topq.traj.q[end],topq.traj.q̇[end])
RB.solve!(
    RB.DynamicsProblem(topq_longtime,top_contact_dynfuncs),
    RB.ZhongQCCP();
    tspan = (0.0,50.0),
    dt=2e-3,
    ftol=1e-14,
    maxiters=100,exception=false,verbose=false
)
me = RB.mechanical_energy!(topq_longtime)
lines(me.E)
rp5 = RB.get_trajectory!(topq_longtime,1,5)
lines(topq_longtime.traj.t,-rp5[3,:].-(-rp5[3,1]))

topn_longtimes = [
    begin
        h = 2*0.01897941
        θ = π/18
        origin_position = h.*[0,sin(θ),cos(θ)]
        R = RotX(θ)
        origin_velocity = [0.0,0.0,0.0]
        Ω = [0.0,0.0,200.0]
        bot = make_top(origin_position,R,origin_velocity,Ω,RB.NCF.NC;μ,e,loadmesh=true)
        RB.solve!(
            RB.DynamicsProblem(bot,(x)->top_contact_dynfuncs(x;checkpersist=check)),
            RB.ZhongCCP();
            tspan = (0.0,500.0),
            dt=2e-3,
            ftol=1e-14,
            maxiters=100,exception=false,verbose=false
        )
    end
    for check in [true,false]
]

plot_traj!(
    topn_longtimes[1];
    showinfo=false,
    # rigidcolor=:white,
    showwire=false,
    showarrows=false,
    # figsize=(0.6tw,0.6tw)
)

plotsave_contact_persistent(topn_longtimes[1])
true_me = RB.mechanical_energy!(topn_longtimes[1])
lines(true_me.E)
false_me = RB.mechanical_energy!(topn_longtimes[2])
lines!(false_me.E)
true_rp5 = RB.get_trajectory!(topn_longtimes[1],1,5)
lines(true_rp5)
false_rp5 = RB.get_trajectory!(topn_longtimes[2],1,5)
lines!(false_rp5)
plot_traj!(topn_longtimes[1])
# vp5 = RB.get_velocity!(topn,1,5)
with_theme(theme_pub;
        size = (1tw,0.25tw),
        figure_padding = (0,fontsize,0,fontsize/2)
    ) do
    bot = topn_longtimes[1]
    # bot = topq_longtime
    (;t) = bot.traj
    fig = Figure()
    ax1 = Axis(fig[1,1];xlabel=tlabel, ylabel="Rel. Err.")
    ax2 = Axis(fig[1,2];xlabel=tlabel, ylabel="Abs. Err. (m)")
    ax3 = Axis(fig[1,3];xlabel="x (m)", ylabel="y (m)", aspect=DataAspect())
    skipstep = 2000
    startstep = time2step(1.5,t)
    @myshow length(t[startstep:skipstep:end])
    mo=5
    scaling = 10.0^(-mo)
    Label(fig[1,1,Top()],latexstring("\\times 10^{-$(mo)}"))
    lines!(ax1,
        t[startstep:skipstep:end],
        (true_me.E[startstep:skipstep:end].-true_me.E[startstep])./true_me.E[startstep]./scaling
    )
    lines!(ax1,
        t[startstep:skipstep:end],
        (false_me.E[startstep:skipstep:end].-false_me.E[startstep])./false_me.E[startstep]./scaling
    )
    # ylims!(ax1,-1e-4,1e-4)
    xlims!(ax1,extrema(t)...)
    @myshow true_rp5[3,startstep]
    mo=5
    scaling = 10.0^(-mo)
    Label(fig[1,2,Top()],latexstring("\\times 10^{-$(mo)}"))
    lines!(ax2,
        t[startstep:skipstep:end],
        ((-true_rp5[3,startstep:skipstep:end]).-(-true_rp5[3,startstep]))./scaling
    )
    lines!(ax2,
        t[startstep:skipstep:end],
        ((-false_rp5[3,startstep:skipstep:end]).-(-false_rp5[3,startstep]))./scaling
    )
    xlims!(ax2,extrema(t)...)

    scatterlines!(ax3,
        true_rp5[1:2,end:end];markersize=fontsize/2
    )
    lines!(ax3,
        true_rp5[1:2,begin:skipstep:end];label="Classified CCP"
    )
    lines!(ax3,
        false_rp5[1:2,begin:skipstep:end];label="Unclassified CCP"
    )
    ax3.xticks = [-0.001,0,0.001]
    Legend(fig[1,4],ax3)
    Label(fig[1,1,TopLeft()], rich("($(alphabet[1]))",font=:bold))
    Label(fig[1,2,TopLeft()], rich("($(alphabet[2]))",font=:bold))
    Label(fig[1,3,TopLeft()], rich("($(alphabet[3]))",font=:bold))
    savefig(fig,"spinningtop_longtime")
    fig    
end

# sliding

R = RotX(π/18)
origin_position = [0,0,0.037]
Ω = [0,0,50.0]
dts = [1e-3,1e-2]
checks = [true,false]
tops_e0 = [
    begin   
        topbot = make_top(origin_position,R,origin_velocity,Ω,RB.NCF.NC; μ = 0.01, e = 0.0,loadmesh=true)
        RB.solve!(
            RB.DynamicsProblem(topbot,(x)->top_contact_dynfuncs(x;checkpersist=check)),
            RB.ZhongCCP();
            tspan=(0.0,2.0),
            dt,ftol=1e-14,maxiters=50,exception=false,#verbose_contact=false
        )
    end
    for dt in dts, check in checks
]

plot_traj!(tops_e0[1])
vp5 = RB.get_velocity!(tops_e0[1],1,5)

GM.activate!();with_theme(theme_pub;
        size = (1.0tw,0.45tw),
        figure_padding = (0,1fontsize,0,0),
        Axis3 = (
            azimuth = 4.695530633326983,
            elevation = 0.07269908169872409,
            protrusions = 0.0,
        ),
        Poly = (
            cycle = [:patchcolor=>:color],
            transparency = true,
        )
    ) do
    bot = tops_e0[1]
    stepstart = 1
    t = bot.traj.t[stepstart:end]
    rp5 = RB.get_trajectory!(bot,1,5)
    rp1 = RB.get_trajectory!(bot,1,1)
    vp5 = RB.get_velocity!(bot,1,5)[stepstart:end]
    me = RB.mechanical_energy!(bot)[stepstart:end]
    steps = 1:100:1800
    nstep = length(steps)
    alphas = fill(0.1,nstep)
    alphas[1:3] = [1,0.4,0.2]
    alphas[end] = 1
    cg = cgrad(:winter, nstep, categorical = true)
    fig = Figure()
    gd1 = fig[1,1:2] = GridLayout()
    gd2 = fig[2,1] = GridLayout()
    gd3 = fig[2,2] = GridLayout(;tellheight=false)
    plot_traj!(
        bot;
        AxisType=Axis3,
        fig = gd1,
        doslide = false,
        showinfo=true,
        # rigidcolor=:white,
        showwire=false,
        showpoints=false,
        showlabels=false,
        showarrows=false,
        showcables=false,
        showmesh=false,
        showtitle=false,
        xlims = (-0.1,1.65),
        ylims = (-0.1,0.2),
        zlims = (-1e-6,0.1),
        sup! = (ax,tgob,sgi) -> begin
            hidey(ax)
            ax.xlabeloffset = 0.0
            ax.zlabeloffset = 2fontsize
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                (;r,g,b) = cg[istep]
                viz!(ax,Observable(suptg);
                    showlabels=false,
                    showarrows=false,
                    showpoints=false,
                    show_nodes=true,
                    meshcolor=Makie.RGBA(r,g,b,alphas[istep])
                )
            end
            lines!(ax,rp5)
        end
    )
    ax2 = Axis(gd1[2,1],
        xlabel = L"x~(\mathrm{m})",
        ylabel = L"y~(\mathrm{m})",
        tellwidth = false,
    )
    # ax2.xlabelpadding = -fontsize
    lines!(ax2,rp5[1:2,:],label="tip")
    lines!(ax2,rp1[1:2,:],label="edge")
    ylims!(ax2,-0.05,0.05)
    xlims!(ax2,-0.05,2.0)
    axislegend(ax2;position=:rc)
    as3 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = L"\alpha-\pi~(\mathrm{Rad})",
    )
    mo_α=8
    scaling = 10.0^(-mo_α)
    Label(gd2[1,1,Top()],latexstring("\\times 10^{-$(mo_α)}"))
    contacts_traj_voa = VectorOfArray(bot.contacts_traj)[:,stepstart:end]
    c1s = contacts_traj_voa[end,:]
    idx_per = findall((x)->RB.doespersist(x;Λtol=0),c1s) #∩ idx_sli
    # @show idx_per
    α_per = map(c1s[idx_per]) do c
        RB.get_contact_angle(c;Λtol=0)
    end
    lines!(as3,t[idx_per],abs.(α_per.-π)./scaling;)
    ax4 = Axis(gd3[1,1:2],
        xlabel = tlabel,
        ylabel = "Energy (J)",
    )
    lines!(ax4,t,me.E, label="E")
    lines!(ax4,t,me.T, label="T")
    lines!(ax4,t,me.V, label="V")
    Legend(gd3[1,3],ax4,orientation=:vertical,tellwidth=false,tellheight=false)
    # axislegend(ax4,position=:rt)
    xlims!(as3,0,2.0)
    xlims!(ax4,0,2.0)
    Label(
        gd1[1,1,TopLeft()],
        rich("($(alphabet[1]))",font=:bold),
        justification = :right,
        padding = (5fontsize,0,0,0),
    )
    Label(
        gd1[2,1,TopLeft()],
        rich("($(alphabet[2]))",font=:bold),
        justification = :right,
        padding = (fontsize,0,0,0)
    )
    Label(
        gd2[1,1,TopLeft()],
        rich("($(alphabet[3]))",font=:bold)
    )
    Label(
        gd3[1,1,TopLeft()],
        rich("($(alphabet[4]))",font=:bold)
    )
    # colsize!(fig.layout,1,Fixed(0.40tw))
    colgap!(fig.layout,1,3fontsize)
    rowgap!(gd1,0)
    rowgap!(fig.layout,0)
    rowsize!(fig.layout,2,Fixed(0.11tw))
    # rowsize!(gd1,2,Fixed(0.05tw))
    # rowgap!(gd1,0)
    # rowsize!(gd3,2,0.1tw)
    savefig(fig,"spinningtop_sliding")
    fig
end

with_theme(theme_pub;
        size = (0.95tw, 0.18tw),
        figure_padding = (0,fontsize/2,0,0)
    ) do
    fig = Figure()
    ax1 = Axis(fig[1,1],
            xlabel = tlabel,
            ylabel = "Abs. Err. (m)",
        )
    ax2 = Axis(fig[1,2],
        xlabel = tlabel,
        ylabel = "Abs. Err. (m)",
    )
    mo_rp5=6
    scaling = 10.0^(-mo_rp5)
    Label(fig[1,1,Top()],latexstring("\\times 10^{-$(mo_rp5)}"))
    for bot in tops_e0[1,:]
        t = bot.traj.t
        rp5 = RB.get_trajectory!(bot,1,5)
        lines!(ax1,t,((-rp5[3,begin:end]).-(-rp5[3,begin]))./scaling)
    end
    mo_rp5=4
    scaling = 10.0^(-mo_rp5)
    Label(fig[1,2,Top()],latexstring("\\times 10^{-$(mo_rp5)}"))
    labels = [
        "Classified CCP",
        "Unclassified CCP"
    ]
    for (label,bot) in zip(labels,tops_e0[2,:])
        t = bot.traj.t
        rp5 = RB.get_trajectory!(bot,1,5)
        lines!(ax2,t,((-rp5[3,begin:end]).-(-rp5[3,begin]))./scaling;label)
    end
    # hidex(ax1)
    Legend(fig[1,3],ax2)
    Label(
        fig[1,1,TopLeft()],
        rich("($(alphabet[1]))",font=:bold)
    )
    Label(
        fig[1,2,TopLeft()],
        rich("($(alphabet[2]))",font=:bold)
    )
    xlims!(ax1,0,2)
    xlims!(ax2,0,2)
    savefig(fig,"spinningtop_sliding_cccp")
    fig
end

plotsave_contact_persistent(tops_e0[1],cid=5,tol=0)
#dt

dts = [1e-2,5e-3,3e-3,2e-3,1e-3,5e-4,3e-4,2e-4,1e-4,1e-5]
tops_dt = [
 begin
    top = deepcopy(tops_e0[1])
     RB.solve!(RB.DynamicsProblem(top,top_contact_dynfuncs),
             RB.ZhongCCP();
             tspan=(0.0,0.1),dt,ftol=1e-14,maxiters=50,exception=true)
 end
 for dt in dts
]

me = RB.mechanical_energy!(bot)
lines(me.E)
fig = Figure()
ax = Axis3(fig[1,1])
for top in tops_dt
    lines!(ax,RB.get_trajectory!(top,1,1),label="$(top.traj.t[2])")
end
Legend(fig[1,2],ax)
fig
_,traj_err_avg = get_err_avg(tops_dt;bid=1,pid=1,di=1,field=:traj)
_,vel_err_avg = get_err_avg(tops_dt;bid=1,pid=1,di=1,field=:vel)
GM.activate!();with_theme(theme_pub;
        size = (0.7tw,0.2tw),
        figure_padding = (0,fontsize/2,0,0)
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1]; ylabel = "Traj. Err.")
    ax2 = Axis(fig[1,2]; ylabel = "Vel. Err.")
    plot_convergence_order!(ax1,dts[begin:end-1],traj_err_avg;show_orders=true)
    plot_convergence_order!(ax2,dts[begin:end-1],vel_err_avg;show_orders=true)
    # ax1.xticks = [1e-4,1e-3,1e-2]
    xlims!(ax1,5e-5,2e-2)
    xlims!(ax2,5e-5,2e-2)
    Legend(fig[1,3],ax2)
    Label(fig[1,1,TopLeft()],"($(alphabet[1]))",font=:bold)
    Label(fig[1,2,TopLeft()],"($(alphabet[2]))",font=:bold)
    # savefig(fig,"spinningtop_order")
    fig
end





