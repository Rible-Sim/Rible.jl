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
scalefactor = 4

include(joinpath(pathof(RB),"../../examples/bodies/hammer.jl"))
includet(joinpath(pathof(RB),"../../examples/bodies/hammer.jl"))
include(joinpath(pathof(RB),"../../examples/bodies/cable_ancf.jl"))
includet(joinpath(pathof(RB),"../../examples/bodies/cable_ancf.jl"))
include(joinpath(pathof(RB),"../../examples/robots/meteor_hammer.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/meteor_hammer.jl"))

# set up the inclined plane
R = RotXY(deg2rad(0),deg2rad(-0))
inclined_plane = RB.Plane(R*[0,0,1.0],[-1.039230,-0.0,0.781])
p1 =  [-1.039230,0,0.780]
p2 = [-1.173952,-0.11,0.8646278]
p3 = [-0.901447,-0.112500,0.762216]
n = (p2-p1) × (p3 - p1 )
inclined_plane = RB.Plane(n,p3.+[0,0,0.001])
inclined_plane_env = RB.StaticContactSurfaces(
    [
        inclined_plane
    ]
)

# use dynamic relax to find a static equilibrium, fixing two ends of the cable.
meteor_hammer_DR = make_meteor_hammer(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
RB.GDR!(meteor_hammer_DR;β=2e-4,maxiters=2e5,verbose=false);

# set the initial conditions
meteor_hammer = make_meteor_hammer(;
    ri  = [ 0.0, 0.0, 1.5],
    rix = [ 0.0, 0.0,-1.0],
    rj  = [-0.6*√3, 0.0, 1.2],
    rjx = [ 0.0,-1.0, 0.0],
    e = 0.0,
    μ=0.1,L=1.2,nx=5,R=RotY(-π/2),
    ω=-200n,
)
meteor_hammer_DR.traj.q[end][end-8:end] .= meteor_hammer.traj.q[end][end-8:end]
meteor_hammer_DR.traj.q̇[end][end-8:end] .= meteor_hammer.traj.q̇[end][end-8:end]
RB.set_new_initial!(meteor_hammer,meteor_hammer_DR.traj.q[end],meteor_hammer_DR.traj.q̇[end])


with_theme(theme_pub;
        Poly= (
            transparency = true,
        ) 
    ) do 
    plot_traj!(
        meteor_hammer;
        ground=inclined_plane,
        xlims=(-1.2,1.0),
        ylims=(-0.5,0.5),
        zlims=(-0.4,1.8),
    )
end

# sliding and avg err 
tspan = (0.0,1.5)
h = 2e-4

RB.solve!(
    RB.DynamicsProblem(
        meteor_hammer,
        inclined_plane_env,
        RB.RestitutionFrictionCombined(
            RB.NewtonRestitution(),
            RB.CoulombFriction(),
        )
    ),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        ),
    );
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

plotsave_contactpoints(meteor_hammer,) #src

rp2 = RB.get_trajectory!(meteor_hammer,6,10)#src
lines(rp2[2,:])#src

me = RB.mechanical_energy!(meteor_hammer)#src
lines((me.E.-me.E[begin])./me.E[begin])#src

# Convergence analysis
# dts = [1e-1,3e-2,1e-2,3e-3,1e-3,1e-4] #src
dts = [1e-2,5e-3,2e-3,0.99e-3,5e-4,2e-4,2e-5]
meteor_hammers_dt = [
    begin
        meteor_hammer_dt = deepcopy(meteor_hammer)
        RB.solve!(
            RB.DynamicsProblem(
                meteor_hammer_dt,
                inclined_plane_env,
                RB.RestitutionFrictionCombined(
                    RB.NewtonRestitution(),
                    RB.CoulombFriction(),
                )
            ),
            RB.DynamicsSolver(
                RB.Zhong06(),
                RB.InnerLayerContactSolver(
                    RB.InteriorPointMethod()
                ),
            );
            tspan=(0.0,0.12),dt,ftol=1e-14,maxiters=50,exception=false
        ).prob.bot
    end
    for dt in dts
]

_, err_avg = RB.get_err_avg(meteor_hammers_dt;bid=6,pid=10,di=1)
GM.activate!(;scalefactor);with_theme(theme_pub;
        size = (0.9tw,0.45tw),
        figure_padding = (fontsize,fontsize,0,0),
        Axis3=(
            azimuth = 7.595530633326987,
            elevation = 0.14269908169872403
        )
    ) do
    bot = meteor_hammer
    (;t) = bot.traj
    bots = meteor_hammers_dt
    rp1 = RB.get_trajectory!(bot,6,1)
    rp12 = RB.get_trajectory!(bot,6,12)
    rp10 = RB.get_trajectory!(bot,6,10)
    fig = Figure()
    gd2 = fig[1,2] = GridLayout()
    gd3 = fig[2,2] = GridLayout()
    gd1 = fig[:,1] = GridLayout()
    steps = 1:1000:length(bot.traj)
    laststep = last(steps)
    cg = cgrad(:winter, length(steps), categorical = true, rev = true)
    nstep = length(steps)
    alphas = fill(0.2,nstep)
    alphas[1:3] = [1, 0.2, 0.2]
    alphas[end] = 1
    plot_traj!(
        bot;
        AxisType=Axis3,
        fig = gd1,
        showtitle=false,
        showpoints=false,
        showlabels=false,
        showarrows=false,
        showcables=false,
        showwire=false,
        showmesh=false,
        xlims=(-1.2,1.0),
        ylims=(-0.3,0.8),
        zlims=(-0.0,1.8),
        doslide=false,
        showinfo=false,
        ground=inclined_plane,
        sup! = (ax,_,_) -> begin
            for (i,step) in enumerate(steps)
                RB.goto_step!(bot,step)
                tgvis = deepcopy(bot.structure)
                (;r,g,b) = cg[i]
                RB.viz!(ax,tgvis;
                    meshcolor=Makie.RGBAf(r,g,b,alphas[i])
                )
            end
            lines!(ax,rp1[begin:laststep],color=:red)
            lines!(ax,rp12[begin:laststep],color=:blue)
            lines!(ax,rp10[begin:laststep],color=:purple)
            handlemesh = load(joinpath(RB.assetpath(),"把柄.STL")) |> RB.make_patch(;
                scale = 1/400,
                trans = [ 0.0, 0.0, 1.55],
                rot = RotY(-π/2)
            )
            RB.hidey(ax)
            mesh!(ax,handlemesh)
        end
    )
    ax2 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = "Energy (J)",
        # aspect = DataAspect()
    )
    lines!(ax2,t,me.E,label="E")
    lines!(ax2,t,me.T,label="T")
    lines!(ax2,t,me.V,label="V")
    # axislegend(ax2)
    Legend(gd2[1,2],ax2)
    ax3 = Axis(gd3[1,1])
    plot_convergence_order!(ax3,dts[begin:end-1],err_avg;show_orders=true)
    Legend(gd3[1,2],ax3)
    colsize!(fig.layout,1,0.45tw)
    rowgap!(fig.layout,0)
    Label(
        gd1[1,1,TopLeft()],"($(alphabet[1]))",font = :bold,
        padding = (0,0,0,0),
        justification = :right,
        lineheight = 1.0,
        halign = :center,
        valign = :bottom,        
    )
    Label(
        gd2[1,1,TopLeft()],"($(alphabet[2]))",font = :bold
    )
    Label(
        gd3[1,1,TopLeft()],"($(alphabet[3]))",font = :bold
    )
    savefig(fig,"meteor_sliding")
    fig
end

# First Scenario
meteor_hammer_DR = make_meteor_hammer(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
RB.GDR!(meteor_hammer_DR;β=2e-4,maxiters=2e5,verbose=false);
meteor_hammer = make_meteor_hammer(;
    ri  = [ 0.0, 0.0, 1.5],
    rix = [ 0.0, 0.0,-1.0],
    rj  = [-0.6*√3, 0.0, 1.2],
    rjx = [ 0.0,-1.0, 0.0],
    e = 0.5,
    μ = 0.9,
    L=1.2,nx=5,R=RotY(-π/2)
)
RB.set_new_initial!(meteor_hammer,meteor_hammer_DR.traj.q[end],meteor_hammer_DR.traj.q̇[end])

plane_normal = RotZY(deg2rad(-7.5),deg2rad(-90))*[0,0,1.0]
inclined_plane = RB.Plane(plane_normal,[0.05,0,0])

inclined_plane_env = RB.StaticContactSurfaces(
    [
        inclined_plane
    ]
)

GM.activate!();with_theme(theme_pub;
        Poly= (
            transparency = true,
        ) 
    ) do 
    plot_traj!(meteor_hammer;
    ground=inclined_plane
    )
end

h = 2e-4
tspan = (0.0,1.5)
RB.solve!(
    RB.DynamicsProblem(
        meteor_hammer,
        inclined_plane_env,
        RB.RestitutionFrictionCombined(
            RB.NewtonRestitution(),
            RB.CoulombFriction(),
        )
    ),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        ),
    );
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

me = RB.mechanical_energy!(meteor_hammer)
lines(meteor_hammer.traj.t,(me.E.-me.E[begin])./me.E[begin])
rp2 = RB.get_trajectory!(meteor_hammer,6,2)
vp2 = RB.get_velocity!(meteor_hammer,6,2)
pnvp2 = [plane_normal'*v for v in vp2 ]
lines(pnvp2)
plotsave_contact_persistent(meteor_hammer,tol=0) #src

with_theme(theme_pub;
        size = (0.95tw,0.5tw),
        figure_padding = (0,fontsize,0,0),
        Axis3=(
            azimuth = 8.125530633326981,
            elevation = 0.18269908169872404
        ),
        Poly=(
            transparency=true,
        )
    ) do
    bot = meteor_hammer
    (;t) = bot.traj
    fig = Figure()
    gd2 = fig[1,2] = GridLayout()
    gd3 = fig[2,2] = GridLayout()
    gd1 = fig[:,1] = GridLayout()
    steps = 1:500:4001
    alphas = [1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, 0.2, 1]
    cg = cgrad(:winter, length(steps), categorical = true, rev = true)
    plot_traj!(
        bot;
        AxisType=Axis3,
        fig = gd1,
        showtitle=false,
        showpoints=false,
        showlabels=false,
        showarrows=false,
        showcables=false,
        showmesh=false,
        showwire=false,
        xlims=(-1.2,0.1),
        ylims=(-0.5,0.5),
        zlims=(-0.1,1.8),
        doslide=false,
        showinfo=false,
        ground=inclined_plane,
        sup! = (ax,_,_) -> begin
            for (i,step) in enumerate(steps)
                RB.goto_step!(bot,step)
                tgvis = deepcopy(bot.structure)
                (;r,g,b) = cg[i]
                RB.viz!(ax,tgvis;meshcolor=Makie.RGBA(r,g,b,alphas[i]))
            end
            lines!(ax,rp2[:,begin:4001],color=:blue)
            handlemesh = load(joinpath(RB.assetpath(),"把柄.STL")) |> RB.make_patch(;
                scale = 1/400,
                trans = [ 0.0, 0.0, 1.55],
                rot = RotY(-π/2)
            )
            mesh!(ax,handlemesh)
        end
    )
    ax2 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = L"\acute{v}_n~(\mathrm{m/s})",
    )
    impact_time = 2995*h
    step_before_impact = RB.time2step(impact_time,t)
    @myshow step_before_impact
    vminus, vplus = pnvp2[step_before_impact:step_before_impact+1]
    @myshow vminus, vplus, vplus/(-vminus)
    vlines!(ax2,[impact_time];linestyle=:dash)
    lines!(ax2,t,pnvp2,color=:red)
    ylims!(ax2,-7,7)
    xlims!(ax2,extrema(t)...)
    poly!(ax2,Point2f.([[0.45,-7],[0.45,7],[0.75,7],[0.75,-7]]),
        color=Makie.RGBAf(0,0,0,0),
        strokecolor=:black,
        strokewidth=1,
        )
    ax22 = Axis(gd2[1,2],
        xlabel = tlabel,
        ylabel = L"\acute{v}_n~(\mathrm{m/s})",
    )
    vlines!(ax22,[impact_time];linestyle=:dash)
    lines!(ax22,t,pnvp2,color=:red)
    ylims!(ax22,-7,7)
    xlims!(ax22,0.45,0.75)
    ax3 = Axis(gd3[1,1],
        xlabel = tlabel,
        ylabel = "Energy (J)"
    )
    vlines!(ax3,[impact_time];linestyle=:dash)
    lines!(ax3,t,me.E,label="E")
    lines!(ax3,t,me.T,label="T")
    lines!(ax3,t,me.V,label="V")
    Legend(gd3[1,2],ax3,orientation=:horizontal,tellheight=false)
    xlims!(ax3,extrema(t)...)
    Label(
        gd1[1,1,TopLeft()],"($(alphabet[1]))",font = :bold,
    )
    Label(
        gd2[1,1,TopLeft()],"($(alphabet[2]))",font = :bold
    )
    Label(
        gd2[1,2,TopLeft()],"($(alphabet[3]))",font = :bold
    )
    Label(
        gd3[1,1,TopLeft()],"($(alphabet[4]))",font = :bold
    )
    colsize!(fig.layout,1,0.35tw)
    colsize!(gd3,1,0.3tw)
    rowgap!(fig.layout,0)
    colgap!(fig.layout,0)
    savefig(fig,"meteor_swing")
    fig
end
