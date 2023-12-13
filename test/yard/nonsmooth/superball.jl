
#-------- SUPERBall ---------

GM.activate!(); with_theme(theme_pub;
        Scatter = (markersize = fontsize,)	
    ) do
    plot_traj!(ballbot;
        AxisType = Axis3,
        figsize = (0.6tw,0.5tw),
        zlims = [1.0,3.0].-2.0,
        doslide = false,
        showinfo = false,
        # showlabels = false,
        showground = false,
        sup! = (ax,_,_) -> begin
            ax.zlabeloffset = 2fontsize
            ax.azimuth = 3.315530633326984
            ax.elevation = 0.2326990816987242
        end,
        # figname = "superball"
    )
end

GM.activate!(); plot_traj!(ballbot;
    xlims = [-1,10],
    ylims = [-1,3],
    zlims = [-1e-3,2.4],
)

#-- testing
l = 1.7/2
d = l/2
ballbot = superball(
    0.0;
    origin_velocity = SVector(2.0,1.0,0),
    ω = SVector(0.0,0.0,1.0),
    μ = 0.05,
    e = 0.0,
    l,d,
    z0 = l^2/(sqrt(5)*d) - 1e-3,
    visible = true,
    loadmesh = false,
)

function ball_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0])
    )
end

# testing
tspan = (0.0,5.0)
h = 1e-2
prob = RB.DynamicsProblem(ballbot,ball_dynfuncs)
@time RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=100,exception=false)
#-- end testing
GM.activate!(); plotsave_contactpoints(ballbot)

plot_traj!(ballbot;)

me = RB.mechanical_energy!(ballbot)
me.E |> lines

step_start = time2step(1.6,ballbot.traj.t)
step_stop = time2step(2.5,ballbot.traj.t)
r2p1 = RB.get_trajectory!(ballbot,2,1)
r1p2 = RB.get_trajectory!(ballbot,1,2)
r6p2 = RB.get_trajectory!(ballbot,6,2)
lines(r2p1)
lines(r1p2)
lines(r6p2)


dts = [1e-2,3e-3,1e-3,3e-4,1e-4,1e-5]
superballs_dt = [
    begin
        ballbot_dt = deepcopy(ballbot)
        # RB.set_new_initial!(ballbot_dt,ballbot.traj.q[step_start],ballbot.traj.q̇[step_start])
        prob = RB.DynamicsProblem(ballbot_dt,ball_dynfuncs)
        RB.solve!(prob,
            RB.ZhongCCP();
            tspan=(0.0,0.1),dt,ftol=1e-14,
            maxiters=500,exception=false
        )
    end
    for dt in dts
]


_,err_avg = get_err_avg(superballs_dt;bid=2,pid=1,di=1)

GM.activate!(); with_theme(theme_pub;
        figure_padding = (0,0.5fontsize,0,0),
        size = (1.0tw,0.45tw),
        Axis3 = (
            azimuth = 4.825530633326982,
            elevation = 0.6726990816987243
        )
    ) do
    fig = Figure()
    gd2 = fig[1,2] = GridLayout()
    gd3 = fig[2,2] = GridLayout()
    gd1 = fig[:,1] = GridLayout()
    steps = 1:100:501    
    nstep = length(steps)
    alphas = fill(0.15,nstep)
    alphas[1:3] = [1,0.2,0.2]
    alphas[end] = 1
    cg = cgrad(:winter, length(steps), categorical = true, rev = true)
    r1p2 = RB.get_trajectory!(ballbot,1,2)
    r6p2 = RB.get_trajectory!(ballbot,6,2)
    r2p1 = RB.get_trajectory!(ballbot,2,1)
    me = RB.mechanical_energy!(ballbot,)
    plot_traj!(ballbot;
        AxisType = Axis3,
        fig = gd1,
        xlims = [-1,6],
        ylims = [-1,3],
        zlims = [-1e-3,2.4],
        doslide = false,
        showinfo = true,
        showpoints = false,
        showlabels = false,
        showtitle = false,
        showcables = false,
        showmesh = false,
        showwire = false,
        sup! = (ax,_,_) -> begin
            for (i,step) in enumerate(steps)
                RB.goto_step!(ballbot,step)
                tgvis = deepcopy(ballbot.st)
                (;r,g,b) = cg[i]
                db = Makie.parse(Makie.RGBA,"deepskyblue")
                viz!(ax,tgvis;
                showcables=true,
                cablecolor=Makie.RGBAf(db.r,db.g,db.b,Makie.N0f8(alphas[i])),
                meshcolor = Makie.RGBAf(r,g,b,alphas[i]))
            end
            lines!(ax,r2p1)
        end
        # figname = "ballbot"
    )
    ax2 = Axis(gd2[1,1],xlabel=tlabel,ylabel = "Energy (J)")
    lines!(ax2,me.E,label="E")
    lines!(ax2,me.T,label="T")
    lines!(ax2,me.V,label="V")
    Legend(gd2[1,2],ax2,)
    ax3 = Axis(gd3[1,1],ylabel="Traj. Err.")
    plot_convergence_order!(ax3,dts[begin:end-1],err_avg;show_orders=true)
    Legend(gd3[1,2],ax3)
    Label(
        gd1[1,1,TopLeft()],"($(alphabet[1]))",font=:bold
    )
    Label(
        gd2[1,1,TopLeft()],"($(alphabet[2]))",font=:bold
    )
    Label(
        gd3[1,1,TopLeft()],"($(alphabet[3]))",font=:bold
    )
    colsize!(fig.layout,1,0.55tw)
    colgap!(fig.layout,2fontsize)
    # savefig(fig,"ballbot_sliding")
    fig
end

CM.activate!(); plotsave_contactpoints(ballbot,"ballbot_contactpoints")

function plotsave_velocity_restitution(bots,showlegend,
        figname=nothing;
        cps,
        ymids = 0.5ones(length(cps)),
        size = (tw,0.5th),
    )
    with_theme(theme_pub;
            resolution,
            figure_padding = (fontsize,fontsize,0,0),
            fontsize = 6.5 |> pt2px
        ) do
        fig = Figure()
        grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
        for (botid,bot) in enumerate(bots)
            ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})",)
            bodyid,pid = divrem(cps[botid]-1,2) .+ 1
            ṙpc_mid = get_mid_velocity!(bot,bodyid,pid)
            t_mids = get_time_mids(bot)
            lines!(ax1,t_mids,ṙpc_mid[3,:], label="v́ₙ")
            lines!(ax1,t_mids,ṙpc_mid[1,:], label="v́ₜ" )
            xlims!(ax1,0,10.0)
            # ylims!(ax1,-6,ymax)
            Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
            Label(grids[botid,1][1,1,Top()], "Contact No.$(cps[botid])", halign = :center)

            ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e_{\mathrm{eff}}")
            vz = RB.get_velocity!(bot,bodyid,pid)[3,:]

            contacts_traj_voa = VectorOfArray(bot.contacts_traj)
            c1_traj = contacts_traj_voa[cps[botid],:]
            # steps = 1:length(c1_traj)
            idx_imp = findall(isimpact, c1_traj)
            @show idx_imp |> length
            e_eff = -vz[idx_imp]./vz[idx_imp.-1]
            scatter!(ax2,e_eff)
            @show mean(e_eff)
            ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
            # ax2.yticks = [0,0.3,0.7,1.0]
            ax2.xticks = collect(1:length(e_eff))
            xlims!(ax2,0.5,length(e_eff)+0.5)
            Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")

            Label(grids[botid,2][1,1,Top()], "Contact No.$(cps[botid])", halign = :center)

            if botid == 3 && showlegend
                axislegend(ax1;position=:cb,orientation=:horizontal)
            end
        end
        rowgap!(fig.layout,0)
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plotsave_velocity_restitution(
    [ballbot,ballbot,ballbot,ballbot],true;
    cps = [2,3,6,12],
)
CM.activate!(); plotsave_velocity_restitution(
    [ballbot,ballbot,ballbot,ballbot],true,"ballbot_restitution";
    cps = [2,3,6,12],
)

# rolling
ballbot = superball(
    0.0;
    origin_velocity = SVector(7.0,2.0,-7.0),
    ω = SVector(0.0,0.0,0.0),
    μ = 0.9,
    e = 0.8,
    l,d,
    z0 = l^2/(sqrt(5)*d) + 2.0,
    visible = true,
    loadmesh = false,
)

# test rolling
tspan = (0.0,5.0)
h = 1e-2
prob = RB.DynamicsProblem(ballbot,ball_dynfuncs)
@time RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=200,exception=false)

GM.activate!(); plotsave_contactpoints(ballbot)

plot_traj!(ballbot;auto=true)

GM.activate!(); with_theme(theme_pub;
        figure_padding = (0,0.5fontsize,0,0),
        size = (1.0tw,0.40tw),
        Axis3 = (
            azimuth = 4.7855306333269805,
            elevation = 0.03269908169872391
        )
    ) do
    fig = Figure()
    bot = ballbot
    (;t) = bot.traj
    gd1 = fig[1,1] = GridLayout()
    gd23 = fig[2,1] = GridLayout()
    gd2 = gd23[1,1] = GridLayout()
    gd3 = gd23[1,2] = GridLayout()
    imptimes = [0.25,0.29,0.30,0.34]
    impstep = time2step(imptimes[1],bot.traj.t)
    steps = vcat(1,impstep,collect(impstep:50:length(t)))
    nstep = length(steps)
    alphas = fill(0.15,nstep)
    alphas[1:3] = [1,0.2,0.2]
    alphas[end] = 1
    cg = cgrad(:winter, nstep, categorical = true, rev = true)
    step_start = time2step(0.1,bot.traj.t)
    step_stop = time2step(0.35,bot.traj.t)
    v2p1 = RB.get_velocity!(bot,2,1,step_start:step_stop)
    v1p2 = RB.get_velocity!(bot,1,2,step_start:step_stop)
    v6p2 = RB.get_velocity!(bot,6,2,step_start:step_stop)
    r2p1 = RB.get_trajectory!(bot,2,1)
    me = RB.mechanical_energy!(bot,)
    plot_traj!(bot;
        AxisType = Axis3,
        fig = gd1,
        xlims = [-1,20],
        ylims = [-1,8],
        zlims = [-1e-3,3.0],
        doslide = false,
        showinfo = true,
        showpoints = false,
        showlabels = false,
        showmesh = false,
        showwire = false,
        showtitle = false,
        showcables  = false,
        sup! = (ax,_,_) -> begin
            for (i,step) in enumerate(steps)
                RB.goto_step!(bot,step)
                tgvis = deepcopy(bot.structure)
                (;r,g,b) = cg[i]
                db = Makie.parse(Makie.RGBA,"deepskyblue")
                viz!(ax,tgvis;
                showcables=true,
                cablecolor=Makie.RGBAf(db.r,db.g,db.b,Makie.N0f8(alphas[i])),
                meshcolor = Makie.RGBAf(r,g,b,alphas[i]))
            end
            hidey(ax)
            lines!(ax,r2p1)
        end
        # figname = "bot"
    )
    ax31 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = L"\dot{z}~(\mathrm{m/s})",
        limits = (t[step_start]+0.06,t[step_stop],-11.6,11.6)
    )
    vlines!(ax31,imptimes[1:2],linestyle=:dash)
    lines!(ax31,t[step_start:step_stop],v2p1[3,:])
    ax32 = Axis(gd2[1,2],
        xlabel = tlabel,
        ylabel = L"\dot{z}~(\mathrm{m/s})",
        limits = (t[step_start]+0.06,t[step_stop],-11.6,11.6)
    )
    vlines!(ax32,imptimes[1:2],linestyle=:dash)
    lines!(ax32,t[step_start:step_stop],v1p2[3,:])
    ax33 = Axis(gd2[1,3],
        xlabel = tlabel,
        ylabel = L"\dot{z}~(\mathrm{m/s})",
        limits = (t[step_start]+0.06,t[step_stop],-11.6,11.6)
    )
    vlines!(ax33,imptimes[[1,3,4]],linestyle=:dash)
    lines!(ax33,t[step_start:step_stop],v6p2[3,:])
    ax2 = Axis(gd3[1,1],
        xlabel = tlabel,
        ylabel = "Energy (J)",
    )
    lines!(ax2,t,me.E,label="E")
    lines!(ax2,t,me.T,label="T")
    lines!(ax2,t,me.V,label="V")
    xlims!(ax2,extrema(t)...)
    
    # (gd2[1,2],ax3)
    Legend(gd3[1,2],
        ax2;
        # position=:rt,
        # orientation=:horizontal,
        tellheight=false
    )

    Label(
        gd1[1,1,TopLeft()],"($(alphabet[1]))",font=:bold
    )
    Label(
        gd2[1,1,TopLeft()],"($(alphabet[2]))",font=:bold
    )
    Label(
        gd2[1,2,TopLeft()],"($(alphabet[3]))",font=:bold
    )
    Label(
        gd2[1,3,TopLeft()],"($(alphabet[4]))",font=:bold
    )
    Label(
        gd3[1,1,TopLeft()],"($(alphabet[5]))",font=:bold
    )
    colsize!(gd23,1,0.55tw)
    rowsize!(fig.layout,1,0.17tw)
    savefig(fig,"ballbot_rolling")
    fig
end


c6= VectorOfArray(ballbot.contacts_traj)[6,:]
c6_1 = c6[4163]
c6_1.state.Λ ⋅ c6_1.state.v
c6_1.state.v
check_Coulomb(1,c6_1)

GM.activate!(); plotsave_friction_direction(
        [ballbot,ballbot,ballbot],L"\mathrm{DistanceSpringDamper~No.}",
        [3,6,12]; 
        size = (tw,0.5tw),
        mo_α = 8,
        mo = 8,
        vtol=1e-7,
)

CM.activate!(); plotsave_friction_direction(
        [ballbot,ballbot,ballbot],L"\mathrm{DistanceSpringDamper~No.}",
        [3,6,12], 
        "ballbot_friction_direction"; 
        size = (tw,0.5tw),
        mo_α = 8,
        mo = 8,
        vtol=1e-7,
)