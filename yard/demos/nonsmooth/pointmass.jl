#-- deps
include(joinpath(@__DIR__,"demos/nonsmooth/deps.jl"))
import Rible as RB
include(joinpath(@__DIR__,"vis.jl"))
includet(joinpath(@__DIR__,"vis.jl")) #jl
#-- deps end

#-- point Mass
include(joinpath(@__DIR__,"../examples/robots/pointmass.jl"))
includet(joinpath(@__DIR__,"../examples/robots/pointmass.jl")) #jl

# pm_contact_dynfuncs
function pm_contact_dynfuncs(bot;θ=0.0)
    a = tan(θ)
    n = [-a,0,1]
    ## horizontal plane
    inclined_plane = RB.Plane(n,zeros(3))
    RB.contact_dynfuncs(bot;flatplane = inclined_plane)
end

# time
tspan = (0.0,0.455)
tspan = (0.0,1.5)
h = 1e-3

restitution_coefficients = [0.5]
v0s = [1.0]

# pointmass
pm = new_pointmass(;
    e = restitution_coefficients[1],
    μ=0.1,
    origin_velocity = [v0s[1],0,0]
);

# simulation
prob = RB.SimProblem(pm,pm_contact_dynfuncs)
RB.solve!(
    prob,
    RB.ZhongCCP();
    tspan,dt=h,
    ftol=1e-14,
    maxiters=50,
    exception=false
);

# visualize
GM.activate!()
with_theme(RB.theme_pub;
            resolution = (0.8tw,0.2tw),
            figure_padding = (0,fontsize,0,0),
            Axis3 = (
                azimuth = 4.575530633326986,
                elevation = 0.13269908169872408
            ),
            color = :red,
            palette = (
                color = [:red, :blue],
            ),
            Scatter = (
                color = :red,
                cycle = [:color ],
            )
    ) do
    bot = pm
    fig = Figure()
    gd = fig[1,1] = GridLayout()
    gd1 = fig[1,2] = GridLayout()
    colsize!(fig.layout,2,Fixed(0.5tw))
    steps = 1:60:1000
    cg = cgrad(:winter, length(steps), categorical = true)
    plot_traj!(bot;
        doslide=false,
        AxisType = Axis3,
        fig = gd,
        ## gridsize=(1,4),
        xlims=(-1e-3,1.0),
        ylims=(-0.4,0.4),
        zlims=(-1e-3,1.0),
        showinfo =false,
        showmesh=false,
        showwire=false,
        showlabels=false,
        showcables=false,
        showpoints=false,
        showtitle=false,
        sup! = (ax,tgob,sgi) -> begin
            RB.hidey(ax)
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                RB.viz!(ax,Observable(suptg);
                    showlabels=false,
                    showmesh=false,
                    showcables=false,
                    showwire=false,
                    showpoints=true,
                    pointcolor=cg[istep]
                )
            end
        end
        ## figname="pointmass_e05_v00",
        ## figname="pointmass_e05_v10"
    )
    Label(gd[1,1,TopLeft()], 
        rich("($(alphabet[1]))",font=:bold)
    )
    (;t) = bot.traj
    ax1 = Axis(gd1[1,1], xlabel = tlabel, ylabel = L"\dot{z}~(\mathrm{m/s})")
    vp1 = RB.get_velocity!(bot,1,1)
    lines!(ax1,t,vp1[3,:])
    xlims!(ax1,t[begin],t[end])
    ## ylims!(ax1,-6,6)
    Label(gd1[1,1,TopLeft()], 
        rich("($(alphabet[2]))",font=:bold)
    )
    hlines!(ax1,[0],color=:gray)
    ax2 = Axis(gd1[1,2], xlabel = tlabel, ylabel = L"\dot{x}~(\mathrm{m/s})")
    ## me = RB.mechanical_energy!(bot)
    ## lines!(ax2,t,me.E)
    lines!(ax2,t,vp1[1,:])
    xlims!(ax2,t[begin],t[end])
    ## ylims!(ax2,-1,11)
    ## ax2.yticks = [0,0.3,0.7,1.0]
    ## ax2.xticks = collect(1:xtickmaxs[botid])
    Label(gd1[1,2,TopLeft()], 
        rich("($(alphabet[3]))",font=:bold)
    )
    savefig(fig,"pointmass_bouncing")
    fig
end

# inclined_plane
θ = 15 |> deg2rad
inclined_plane = RB.Plane([-tan(θ),0,1],zeros(3))
origin_position = [0.0,0,-1e-7]
origin_velocity = [2.0cos(θ),0,2.0sin(θ)]
## origin_velocity = [0,-2.0,0]

# analytical
g = 9.81
μ=0.3
vo = norm(origin_velocity)
a = -μ*g*cos(θ)-g*sin(θ)
## μ*g*cos(θ)-g*sin(θ)
tf = -vo/a
## d(t) -> vo*t+1/2*a*t^2

# simulation
tspan = (0.0,0.6)
pm = new_pointmass(;e=0.0, μ, origin_position, origin_velocity)

prob = RB.SimProblem(pm,(x)->pm_contact_dynfuncs(x;θ))
RB.solve!(prob,RB.ZhongCCP();tspan,dt=1e-3,ftol=1e-14,maxiters=50,exception=false);

# post processing
rp1 = RB.get_trajectory!(pm,1,1)
ṙp1 = RB.get_velocity!(pm,1,1)
dp1 = rp1.u .|> norm
vl1 = [u ⋅ normalize(origin_velocity) for u in ṙp1]
## overshoot!
scatterlines(vl1)

# visualize
GM.activate!(); with_theme(RB.theme_pub;
        ## fontsize = 6 |> pt2px,
        resolution = (1tw,0.2tw),
        figure_padding = (fontsize,fontsize,0,0),
        Axis3 = (
            azimuth = 4.575530633326984,
            elevation = 0.16269908169872405,
            ## zlabelvisible = false,
            ## yticklabelsvisible = false,
            zlabeloffset = 2.5fontsize,
        )
    ) do
    steps = 1:50:600
    cg = cgrad(:winter, length(steps), categorical = true)
    bot = pm
    (;t) = bot.traj
    fig = Figure()
    gd1 = fig[1,1] = GridLayout()
    gd2 = fig[1,2] = GridLayout()
    Label(gd1[1,1,TopLeft()], 
        rich("($(alphabet[1]))",font=:bold)
    )
    colsize!(fig.layout,2,Fixed(0.55tw))
    colgap!(fig.layout,fontsize/2)
    plot_traj!(bot;
        doslide=false,
        AxisType=Axis3,
        fig = gd1,
        xlims=(-8e-3,0.5),
        ylims=(-0.1,0.1),
        zlims=(-8e-3,0.2),
        showinfo =false,
        showmesh=false,
        showwire=false,
        showlabels=false,
        showcables=false,
        showpoints=false,
        showtitle=false,
        ground=inclined_plane,
        sup! = (ax,tgob,sgi) -> begin
            RB.hidey(ax)
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                RB.viz!(ax,Observable(suptg);
                    showlabels=false,
                    showmesh=false,
                    showcables=false,
                    showwire=false,
                    showpoints=true,
                    pointcolor=cg[istep]
                )
            end
        end,
        ## figname="pointmass_sliding"
    )
    ax1 = Axis(gd2[1,1], xlabel = tlabel, ylabel = "disp. (m)")
    Label(gd2[1,1,TopLeft()], 
        rich("($(alphabet[2]))",font=:bold)
    )
    ax2 = Axis(gd2[1,2], xlabel = tlabel, ylabel = "disp. (m)")
    Label(gd2[1,2,TopLeft()], 
        rich("($(alphabet[3]))",font=:bold)
    )
    rp1 = RB.get_trajectory!(bot,1,1)
    ṙp1 = RB.get_velocity!(bot,1,1)
    dp1 = rp1.u .|> norm
    vl1 = ṙp1.u .|> norm
    @myshow findmax(dp1)
    stopstep = RB.time2step(tf,t)
    lines!(ax1,t,dp1)
    xlims!(ax1,0,0.6)
    scatterlines!(ax2,t[1:stopstep],(t) -> vo*t+1/2*a*t^2,color=:red,label="Analytic")
    scatterlines!(ax2,t[stopstep:end],(t) -> vo*tf+1/2*a*tf^2,color=:red)
    scatter!(ax2,t,dp1,color=:blue,marker=:diamond,label="NMSI")
    ## lines!(ax2,t,vl1)
    axislegend(ax2;position=:rt,orientation=:horizontal,)
    xlims!(ax2,0.369,0.389)
    ylims!(ax2,3.7161e-1,3.7164e-1)
    vlines!(ax1,t[stopstep])
    vlines!(ax2,t[stopstep])
    text!(ax1,latexstring("t_s=0.372(s)"), 
        position = (tf+0.02, 0.19),
        fontsize = 6 |> pt2px,
        align = (:left, :center)
    )
    text!(ax2,latexstring("t_s=0.372(s)"), 
        position = (tf+0.02, 1.0),
        fontsize = 6 |> pt2px,
        align = (:left, :center)
    )
    @show extrema(dp1), extrema(vl1)
    savefig(fig,"pointmass_sliding")
    fig
end
