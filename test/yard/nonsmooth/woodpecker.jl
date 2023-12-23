using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../test/yard/nonsmooth/deps.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME" #src
end #src
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

tw = 455.8843 #pt |> pt2px
scalefactor = 4
#--  slider crank 
include(joinpath(pathof(RB),"../../examples/robots/woodpecker.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/woodpecker.jl"))#jl

# Contact Surfaces
halfspaces = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([0,-1.0,0],[0,-0.0025,0]),
        RB.HalfSpace([0, 1.0,0],[0, 0.0025,0])
    ]
)

# ## Limit cylces
wt = woodpecker() 
wt = woodpecker(;case=3,frictional_coefﬁcient = 0.35) #src
dt = 1e-4
tspan = (0.0,2.0)
prob = RB.DynamicsProblem(
    wt,
    halfspaces,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    )
)

# Zhong06
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-10,maxiters=50,verbose=false,
    exception=false,progress=false,
    max_restart=3
)

b1vg = RB.get_mid_velocity!(wt,1,0) #src
lines(b1vg[3,:]) #src

GM.activate!(;scalefactor);with_theme(theme_pub;
        size = (1tw,0.24tw),
        figure_padding = (0,fontsize,0,0)
    ) do 

    b1g = RB.get_trajectory!(wt,1,0)
    b1ω = RB.get_mid_angular_velocity!(wt,1)
    b1θ = RB.get_mid_orientation!(wt,1)
    b2ω = RB.get_mid_angular_velocity!(wt,2)
    b2θ = RB.get_mid_orientation!(wt,2)
    fig = Figure()
    ax1 = Axis(fig[1,1],xlabel = tlabel, ylabel = L"z~(\mathrm{m})")
    ax2 = Axis(fig[1,2],xlabel = tlabel, ylabel = L"z~(\mathrm{m})")
    ax3 = Axis(fig[1,3],xlabel = L"\theta_S~(\mathrm{rad})", ylabel = L"\omega_S~(\mathrm{rad}/s)",)
    ax4 = Axis(fig[1,4],xlabel = L"\theta_M~(\mathrm{rad})", ylabel = L"\omega_M~(\mathrm{rad}/s)",)

    period = 0.13054830287206265
    @show f = 1/period
    stiction_times = collect(0.05:period:2.0)
    stiction_steps = RB.time2step.(stiction_times,Ref(wt.traj.t))
    falling_distances = b1g[3,stiction_steps]
    @show length(falling_distances) - 1
    @show falling_distances[begin+1:end] .- falling_distances[begin:end-1] |> mean
    hlines!(ax1,falling_distances,color=:slategrey,linestyle=:dash)
    lines!(ax1,wt.traj.t,b1g[3,:])
    xlims!(ax1,0,2.0)
    hlines!(ax2,falling_distances,color=:slategrey,linestyle=:dash)
    lines!(ax2,wt.traj.t,b1g[3,:])
    xlims!(ax2,0,0.18)
    ylims!(ax2,-0.023,0)

    lines!(ax3,b2θ[3,:],b2ω[1,:])
    xlims!(ax3,-0.6,0.2)
    ylims!(ax3,-10,20)
    lines!(ax4,b1θ[3,:],b1ω[1,:],)
    ylims!(ax4,-25,35)

    Label(fig[1,1,TopLeft()],"($(alphabet[1]))",font=:bold)
    Label(fig[1,2,TopLeft()],"($(alphabet[2]))",font=:bold)
    Label(fig[1,3,TopLeft()],"($(alphabet[3]))",font=:bold)
    Label(fig[1,4,TopLeft()],"($(alphabet[4]))",font=:bold)
    colsize!(fig.layout,2,Relative(0.2))
    colgap!(fig.layout,0.5fontsize)
    savefig(fig,"limit_cycles";backend=CM)
    fig
end

GM.activate!(;scalefactor=1);plot_traj!(wt;showmesh=false,showground=false) #src

# ## Convergence analysis
dts = vcat(
    [10^(-x) for x in range(3,4;length=7)],
    1e-7
)
tspan = (0.0,0.02)
function get_stats(tspan,dts)
    [
        begin
            @timed RB.solve!(
                RB.DynamicsProblem(
                    woodpecker(),
                    halfspaces,
                    RB.RestitutionFrictionCombined(
                        RB.NewtonRestitution(),
                        RB.CoulombFriction(),
                    )
                ),
                RB.DynamicsSolver(
                    solver,
                    RB.InnerLayerContactSolver(
                        RB.InteriorPointMethod()
                    )
                );
                dt,tspan,ftol=ifelse(dt==1e-7,1e-14,1e-12),maxiters=50,verbose=false,
                exception=true,progress=true,
                max_restart=3
            ).prob.bot
        end
        for dt in dts, solver in (RB.Zhong06(), RB.Moreau(0.5))
    ];
end
stats_wt_dts = get_stats(tspan, dts)
wt_dts = map((x)->x.value,stats_wt_dts)

fig = Figure() #src
ax = Axis(fig[1,1]) #src
for bot in wt_dts[begin+1:end,1] #src
    b1r2 = RB.get_mid_velocity!(bot,2,2) #src
    lines!(ax,RB.get_mid_times(bot),b1r2[3,:]) #src
end #src
fig #src

GM.activate!(;scalefactor); with_theme(theme_pub;
        size = (1tw,0.2tw)
    ) do
    fig = Figure()
    ax1 = Axis(fig[1,1],ylabel = "Abs. Err. (rad)")
    ax2 = Axis(fig[1,2],ylabel = "Abs. Err. (rad/s)")
    ax3 = Axis(fig[1,3],xlabel = "Abs. Err. (rad)", ylabel = "Walltime (s)")
    _,err_avg_nmsi_traj = RB.get_err_avg(vcat(wt_dts[begin:end-1,1],wt_dts[end,2]);bid=2,pid=2,di=2,field=:traj)
    plot_convergence_order!(ax1,dts[begin:end-1],err_avg_nmsi_traj;orders=[2],label="NMSI")
    _,err_avg_moreau_traj = RB.get_err_avg(wt_dts[:,2];bid=2,pid=2,di=2,field=:traj)
    plot_convergence_order!(ax1,dts[begin:end-1],err_avg_moreau_traj;orders=[1],marker=:utriangle,color=:blue,label="Moreau")
    _,err_avg_nmsi_vel = RB.get_err_avg(vcat(wt_dts[begin:end-1,1],wt_dts[end,2]);bid=2,pid=2,di=3,field=:midvel)
    plot_convergence_order!(ax2,dts[begin:end-1],err_avg_nmsi_vel;orders=[2],label="NMSI")
    _,err_avg_moreau_vel = RB.get_err_avg(wt_dts[:,2];bid=2,pid=2,di=3,field=:midvel)
    plot_convergence_order!(ax2,dts[begin:end-1],err_avg_moreau_vel;orders=[1],marker=:utriangle,color=:blue,label="Moreau")
    Legend(fig[1,4],ax2,)
    moreau_time = map((x)->x.time-x.gctime,stats_wt_dts[begin:end-1,2])
    nmsi_time = map((x)->x.time-x.gctime,stats_wt_dts[begin:end-1,1])
    scatterlines!(ax3,err_avg_moreau_traj,moreau_time;marker=:utriangle,color=:blue,)
    scatterlines!(ax3,err_avg_nmsi_traj,nmsi_time;marker=:rect,color=:red)
    ax3.xscale = Makie.log10
    ax3.xminorticksvisible = true 
    ax3.xminorgridvisible = true 
    ax3.xminorticks = IntervalsBetween(8)
    ax3.yscale = Makie.log10
    ax3.yminorticksvisible = true 
    ax3.yminorgridvisible = true 
    ax3.yminorticks = IntervalsBetween(4)
    Label(fig[1,1,TopLeft()],"($(alphabet[1]))",font=:bold)
    Label(fig[1,2,TopLeft()],"($(alphabet[2]))",font=:bold)
    Label(fig[1,3,TopLeft()],"($(alphabet[3]))",font=:bold)
    savefig(fig,"woodpecker_convergence";backend=CM)
    fig
end