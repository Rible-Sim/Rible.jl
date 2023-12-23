
#----------- unibot ---------------
unibot = uni(0.0;
            μ = 0.9,
            e = 0.0,
            z0 = 0.2
)

plot_traj!(unibot;
    # zlims = (0,1),
    showground = true
)

function uni_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0])
    )
end

tspan = (0.0,1.5)
h = 1e-3

prob = RB.DynamicsProblem(unibot,uni_dynfuncs)

RB.solve!(prob,RB.Zhong06();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false)

RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=true)

plot_traj!(unibot)

GM.activate!()
with_theme(theme_pub;
        markersize = 0.4fontsize,
        fontsize = 0.8fontsize,
        figure_padding = (fontsize,0,1.5fontsize,fontsize),
    ) do 
    plot_traj!(unibot;
        figsize = (0.9tw,0.7tw),
        xlims = (-0.01,0.15),
        ylims = (-0.15,0.05),
        zlims = (-0.01,0.3),
        gridsize = (2,3),
        attimes = [0,0.179,0.255,0.412,0.694,1.177],
        AxisType=Axis3,
        doslide =false,
        showinfo = false,
        showlabels = false,
        showpoints = false,
        showground = true,
        sup! = (ax,_) -> begin
            # ax.zlabeloffset = 2fontsize
            # ax.xlabeloffset = 1.5fontsizew
            # ax.azimuth = 4.335530633326986
            # ax.elevation = 0.27269908169872403
        end,
        # figname = "unibot",
    )
end

GM.activate!(); plotsave_energy(unibot)
CM.activate!(); plotsave_energy(unibot,"unibot_energy")

unibots = [
    begin
        RB.solve!(
            RB.DynamicsProblem(
                uni(0.0; μ = 0.9, e, z0 = 0.2),
                uni_dynfuncs
            ),
            RB.ZhongCCP();
            tspan,dt=h,ftol=1e-14,maxiters=50,exception=false
        )
    end
    for e = [1.0,0.7,0.3,0.0]
]

plot_traj!(unibots[1])
GM.activate!(); plotsave_energy(unibots[3])
function plotsave_velocity_restitution(bots,showlegend,
        figname=nothing;
        ymids = [1.0,0.7,0.3,0.0]
    )
    with_theme(theme_pub;
            size = (tw,0.4th),
            figure_padding = (fontsize,fontsize,0,0),

        ) do
        fig = Figure()
        grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
        for (botid,bot) in enumerate(bots)
            ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})")
            ṙp2_mid = get_mid_velocity!(bot,1,2)
            t_mids = get_mid_times(bot)
            lines!(ax1,t_mids,ṙp2_mid[3,:], label="v́ₙ")
            lines!(ax1,t_mids,ṙp2_mid[1,:], label="v́ₜ" )
            if botid == 4 && showlegend
                axislegend(ax1;position=:cb,orientation=:horizontal)
            end
            xlims!(ax1,0,1.5)
            # ylims!(ax1,-6,ymax)
            Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
            ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e_{\mathrm{eff}}")
            vz = RB.get_velocity!(bot,1,2)[3,:]

            contacts_traj_voa = VectorOfArray(bot.contacts_traj)
            c1_traj = contacts_traj_voa[2,:]
            # steps = 1:length(c1_traj)
            idx_imp = findall(isimpact, c1_traj)
            @show idx_imp |> length
            e_eff = -vz[idx_imp]./vz[idx_imp.-1]
            scatter!(ax2,e_eff)
            # @show e_eff
            @show mean(e_eff)
            ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
            # ax2.yticks = [0,0.3,0.7,1.0]
            # ax2.xticks = collect(1:length(e_eff))
            Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
        end
        rowgap!(fig.layout,-fontsize)
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plotsave_velocity_restitution(unibots,true)
CM.activate!(); plotsave_velocity_restitution(unibots,true,"unibot_restitution")

unibot_e5 = RB.solve!(
    RB.DynamicsProblem(
        uni(0.0; μ = 0.01, e=0.5, z0 = 0.2),
        uni_dynfuncs
    ),
    RB.ZhongCCP();
    tspan,dt=0.5h,ftol=1e-14,maxiters=50,exception=false,verbose=true,
)
GM.activate!(); plot_traj!(unibot_e5)

me = RB.mechanical_energy!(unibot_e5)
me.E |> lines

GM.activate!(); plotsave_friction_direction(
    [unibot_e5],L"e",[0.5]; 
    size = (0.9tw,0.28tw),
    cid=2
    )

CM.activate!(); plotsave_friction_direction(
        [unibot_e5],L"e=",[0.5], 
        "unibot_friction_direction"; 
        size = (0.9tw,0.28tw)
)

function plotsave_point_traj_vel(bots,figname=nothing)
    with_theme(theme_pub;
        resolution=(1.0tw,0.4tw),
    ) do
        fig = Figure()
        nbot = length(bots)
        for (i,bot) in enumerate(bots)
            ax1 = Axis(
                fig[i,1],
                xlabel=L"x~(\mathrm{m})",
                ylabel=L"y~(\mathrm{m})",
                aspect=DataAspect(),
            )
            ax2 = Axis(
                fig[i,2],
                xlabel=tlabel,
                ylabel=L"\theta",
            )
            Label(fig[i,1,TopLeft()], "($(alphabet[2i-1]))")
            Label(fig[i,2,TopLeft()], "($(alphabet[2i]))")
            ax3 = Axis(fig[i,2],
                        xlabel=tlabel,
                        ylabel=L"v~(\mathrm{m/s})",
                        yticklabelcolor = :red,
                        yaxisposition = :right
                )
            hidespines!(ax3)
            hidexdecorations!(ax3)

            (;t) = bot.traj
            rp2 = RB.get_trajectory!(bot,1,2)
            rpx = rp2[1,:]
            rpy = rp2[2,:]
            lines!(ax1,rpx,rpy)
            # xlims!(ax1,0,13.0)
            # ylims!(ax1,-0.3,0.45)

            ṙp2 = get_mid_velocity!(bot,1,2)
            ṙpx = ṙp2[1,:]
            ṙpy = ṙp2[2,:]
            t_mids = get_mid_times(bot)
            θ_mids = atan.(ṙpy,ṙpx)
            ṙpxy = map(ṙp2.u) do ṙ
                norm(ṙ[1:2])
            end
            lines!(ax2,t_mids,θ_mids)
            lines!(ax3,t_mids,ṙpxy,color=:red)
            xlims!(ax2,extrema(bot.traj.t)...)
            xlims!(ax3,extrema(bot.traj.t)...)
            ylims!(ax2,-π,π)

            if i !== nbot
                hidex(ax1)
                hidex(ax2)
            end

        end
        colsize!(fig.layout,1,Relative(0.35))
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plotsave_point_traj_vel([unibot_e5])
CM.activate!(); plotsave_point_traj_vel([unibot_e5],"unibot_contact_point_traj_vel")

#dt
dts = [1e-2,3e-3,1e-3,3e-4,1e-4,1e-5]
# dts = [1e-2]

unibot_z0 = [
    RB.solve!(
        RB.DynamicsProblem(
            uni(0.0; μ = 0.01, e=0.0, z0 = 0.2-0.13468-1e-5, ωz = 50.0, mbar = 1.0),
            uni_dynfuncs
        ),
        RB.ZhongCCP();
        tspan=(0.0,1.0),dt,ftol=1e-14,maxiters=50,exception=true
    )
    for dt in dts
]

plot_traj!(unibot_z0[6];
    # zlims = (0,1),
    showground = true
)

GM.activate!(); plotsave_error(unibot_z0,dts;size = (0.4tw,0.3tw),)
CM.activate!(); plotsave_error(unibot_z0,dts,"unibot_err";size = (0.4tw,0.3tw),)

GM.activate!(); plotsave_velocity_restitution([unibot_z0],true; ymids = [0.0])
plot_traj!(unibot_z0)

GM.activate!(); plotsave_point_traj_vel([unibot_z0])

me = RB.mechanical_energy!(unibot_z0)
me.E |> lines

GM.activate!(); plot_traj!(unibot_z0)
#-- uni bot end ---