#-- preamble
using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using Rotations
# using Parameters
import GeometryBasics as GB
using OffsetArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using BlockDiagonals
using RecursiveArrayTools
using Interpolations
using EponymTuples
using CircularArrays
using TypeSortedCollections
using DataStructures
using PrettyTables
using Printf
using CoordinateTransformations
using Meshing
using ForwardDiff
using BenchmarkTools
using IterTools
using Unitful
using Match
using FileIO
using Cthulhu
using JET
using TypedTables
using Revise
using AbbreviatedStackTraces
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true
ENV["JULIA_STACKTRACE_MINIMAL"] = true
import Rible as RB
import Meshes
cd(@__DIR__)
include("def.jl"); includet("def.jl")
# includet("plotting.jl")
include("../../vis.jl"); includet("../../vis.jl")
include("../../dyn.jl"); includet("../../dyn.jl")
figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME"
elseif Sys.isapple()
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME"
end
#-- preamble end

#-- point Mass

function new_pointmass(;
        origin_position = [0.0,0.0,1.0],
        origin_velocity = zero(origin_position),
        m = 1.0,
        μ = 0.3,
        e = 0.9
    )
    movable = false
    constrained = true
    Ia = SMatrix{3,3}(Matrix(m*I,3,3))
    mass_locus  = SVector{3}([ 0.0, 0.0, 0.0])
    r̄p1 = SVector{3}([ 0.0, 0.0, 0.0])
    loci = [r̄p1]
    axes = [SVector{3}([ 1.0, 0.0, 0.0])]
    friction_coefficients = [μ]
    restitution_coefficients = [e]
    prop = RB.RigidBodyProperty(
        1,movable,m,Ia,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients,
        ;constrained=constrained
    )
    ω = zero(origin_position)
    R = RotX(0.0)
    loci = Ref(origin_position) .+ Ref(R).*loci
    nmcs = RB.NCF.NC3D1P(loci[1],)
    ci = Int[]
    cstr_idx = Int[]
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω,ci,cstr_idx)
    rb1 = RB.RigidBody(prop,state)

    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)
    cnt = RB.Connectivity(numberedpoints,indexedcoords,tensioned)
    st = RB.Structure(rbs,tensiles,cnt,)
    bot = RB.Robot(st)
end

function pm_contact_dynfuncs(bot;θ=0.0)
    a = tan(θ)
    n = [-a,0,1]
    inclined_plane = RB.Plane(n,zeros(3))
    contact_dynfuncs(bot;flatplane = inclined_plane)
end

tspan = (0.0,0.455)
tspan = (0.0,1.5)
h = 1e-3

# horizontal plane
restitution_coefficients = [0.5]
v0s = [1.0]
pms_hp = [
    begin
        pm = new_pointmass(;
                e,μ=0.1,origin_velocity = [v0,0,0]
            )
        RB.solve!(
            RB.SimProblem(pm,pm_contact_dynfuncs),
            RB.ZhongCCP();
            tspan,dt=h,
            ftol=1e-14,
            maxiters=50,
            exception=false
        )
    end
    for v0 in v0s for e in restitution_coefficients
]

GM.activate!();with_theme(theme_pub;
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
    bot = pms_hp[1]
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
        # gridsize=(1,4),
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
            hidey(ax)
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                viz!(ax,Observable(suptg);
                    showlabels=false,
                    showmesh=false,
                    showcables=false,
                    showwire=false,
                    showpoints=true,
                    pointcolor=cg[istep]
                )
            end
        end
        # figname="pointmass_e05_v00",
        # figname="pointmass_e05_v10"
    )
    Label(gd[1,1,TopLeft()], 
        rich("($(alphabet[1]))",font=:bold)
    )
    (;t) = bot.traj
    ax1 = Axis(gd1[1,1], xlabel = tlabel, ylabel = L"\dot{z}~(\mathrm{m/s})")
    vp1 = RB.get_velocity!(bot,1,1)
    lines!(ax1,t,vp1[3,:])
    xlims!(ax1,t[begin],t[end])
    # ylims!(ax1,-6,6)
    Label(gd1[1,1,TopLeft()], 
        rich("($(alphabet[2]))",font=:bold)
    )
    hlines!(ax1,[0],color=:gray)
    ax2 = Axis(gd1[1,2], xlabel = tlabel, ylabel = L"\dot{x}~(\mathrm{m/s})")
    # me = RB.mechanical_energy!(bot)
    # lines!(ax2,t,me.E)
    lines!(ax2,t,vp1[1,:])
    xlims!(ax2,t[begin],t[end])
    # ylims!(ax2,-1,11)
    # ax2.yticks = [0,0.3,0.7,1.0]
    # ax2.xticks = collect(1:xtickmaxs[botid])
    Label(gd1[1,2,TopLeft()], 
        rich("($(alphabet[3]))",font=:bold)
    )
    savefig(fig,"pointmass_bouncing")
    fig
end

θ = 15 |> deg2rad
inclined_plane = RB.Plane([-tan(θ),0,1],zeros(3))
origin_position = [0.0,0,-1e-7]
origin_velocity = [2.0cos(θ),0,2.0sin(θ)]
# origin_velocity = [0,-2.0,0]

# analytical
g = 9.81
μ=0.3
vo = norm(origin_velocity)
a = -μ*g*cos(θ)-g*sin(θ)
# μ*g*cos(θ)-g*sin(θ)
tf = -vo/a
# d(t) -> vo*t+1/2*a*t^2
tspan = (0.0,0.6)
pm = new_pointmass(;e=0.0, μ, origin_position, origin_velocity)

prob = RB.SimProblem(pm,(x)->pm_contact_dynfuncs(x;θ))
RB.solve!(prob,RB.ZhongCCP();tspan,dt=1e-3,ftol=1e-14,maxiters=50,exception=false)

rp1 = RB.get_trajectory!(pm,1,1)
ṙp1 = RB.get_velocity!(pm,1,1)
dp1 = rp1.u .|> norm
vl1 = [u ⋅ normalize(origin_velocity) for u in ṙp1]
# overshoot!
scatterlines(vl1)
GM.activate!(); with_theme(theme_pub;
        # fontsize = 6 |> pt2px,
        resolution = (1tw,0.2tw),
        figure_padding = (fontsize,fontsize,0,0),
        Axis3 = (
            azimuth = 4.575530633326984,
            elevation = 0.16269908169872405,
            # zlabelvisible = false,
            # yticklabelsvisible = false,
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
            hidey(ax)
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                viz!(ax,Observable(suptg);
                    showlabels=false,
                    showmesh=false,
                    showcables=false,
                    showwire=false,
                    showpoints=true,
                    pointcolor=cg[istep]
                )
            end
        end,
        # figname="pointmass_sliding"
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
    stopstep = time2step(tf,t)
    lines!(ax1,t,dp1)
    xlims!(ax1,0,0.6)
    scatterlines!(ax2,t[1:stopstep],(t) -> vo*t+1/2*a*t^2,color=:red,label="Analytic")
    scatterlines!(ax2,t[stopstep:end],(t) -> vo*tf+1/2*a*tf^2,color=:red)
    scatter!(ax2,t,dp1,color=:blue,marker=:diamond,label="NMSI")
    # lines!(ax2,t,vl1)
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

#-- point mass end

#--  Spinning top

function make_top(origin_position = [0.0,0.0,0.0],
        R = one(RotMatrix{3}),
        origin_velocity = [0.0,0.0,0.0],
        Ω = [0.0,0.0,5.0],
        cT = RB.QCF.QC;
        μ = 0.5,
        e = 0.9,
        color = :slategray,
        constrained=false,
        loadmesh=false,
    )
    ω = R*Ω
    movable = true
    if constrained
        pres_idx = [1,2,3]
    else
        pres_idx = Int[]
    end

    m =  0.58387070
    mass_locus = @SVector zeros(3)
    # Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
    Ī = SMatrix{3,3}(
        Diagonal(SA[
            0.00022129,
            0.00022129,
            0.00030207
            ])
        )
    
    # h = 0.02292582
    radius = 0.044/√2
    h = 2*0.01897941
    loci = [radius.*[1,1,0] for i = 1:4]
    push!(loci,[0,0,-h])
    axes = [SVector(1.0,0,0) for i = 1:5]
    friction_coefficients = [μ for i = 1:5]
    restitution_coefficients = [e for i = 1:5]
    if loadmesh
        topmesh = load(
            RB.assetpath("Toupise2.STL")
        ) |> make_patch(;
            scale=1/1000,
            color,
        )
    else
        pts = Point3.(loci)
        fcs = GB.TriangleFace.([
            [5,1,2],
            [5,4,3],
            [5,3,1],
            [5,2,4],
            [1,4,2],
            [4,1,3],
            [3,2,1],
            [2,3,4]
        ])
        nls = GB.normals(pts,fcs)
        topmesh = GB.Mesh(GB.meta(pts,normals=nls),fcs)
    end
    prop = RB.RigidBodyProperty(
        1,movable,m,Ī,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients;
        constrained
    )
    ri = origin_position+R*loci[5]
    @myshow ri
    if cT == RB.QCF.QC
        nmcs = RB.QCF.QC(m,Ī)
    else
        nmcs = RB.NCF.NC1P3V(ri,origin_position,R,origin_velocity,ω)
    end
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω,pres_idx)
    rb1 = RB.RigidBody(prop,state,topmesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)
    cnt = RB.Connectivity(numberedpoints,indexedcoords,tensioned)
    st = RB.Structure(rbs,tensiles,cnt,)
    bot = RB.Robot(st)
end

function top_contact_dynfuncs(bot;checkpersist=true,)
    contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0]),
        checkpersist,
    )
end

origin_position = [0,0,0.5]
R = RotX(0.0)
origin_velocity = [1.0,0.0,0.0]
Ω = [0.0,0.0,200.0]
# R = rand(RotMatrix3)
# origin_velocity = rand(3)
# Ω = rand(3)

# tspan = (0.0,1.0)
μ = 0.95
e = 0.5
tspan = (0.0,1.8)
h = 1e-4

topq = make_top(origin_position,R,origin_velocity,Ω;μ,e,loadmesh=true)
#note subsequent iteration slow convergence 
#note initial guess can not improve it?
RB.solve!(
    RB.SimProblem(topq,top_contact_dynfuncs),
    RB.ZhongQCCP();
    tspan,
    dt=h,
    ftol=1e-10,
    maxiters=50,exception=false,verbose_contact=true
)
me = RB.mechanical_energy!(topq)
lines(me.E)
rp5 = RB.get_trajectory!(topq,1,5)

plot_traj!(
    topq;
    showinfo=false,
    # rigidcolor=:white,
    showwire=true,
    showarrows=false,
)

topn = make_top(origin_position,R,origin_velocity,Ω,RB.NCF.NC;μ,e,loadmesh=true)
RB.solve!(
    RB.SimProblem(topn,top_contact_dynfuncs),
    RB.ZhongCCP();
    tspan,
    dt=h,
    ftol=1e-14,
    maxiters=50,exception=false,verbose_contact=false
)

plot_traj!(
    topn;
    showinfo=false,
    # rigidcolor=:white,
    showwire=false,
    showarrows=false,
    # figsize=(0.6tw,0.6tw)
)

GM.activate!();with_theme(theme_pub;
        resolution = (1.0tw,0.3tw),
        figure_padding = (0,fontsize,0,0),
        Axis3 = (
            azimuth = 5.1155306333269825,
            elevation = 0.1426990816987241,
            # protrusions = 0.0,
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
    colsize!(fig.layout,2,Fixed(0.50tw))
    colgap!(fig.layout,1.5fontsize)
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
            hidey(ax)
            for (istep,step) in enumerate(steps)
                suptg = deepcopy(bot.structure)
                suptg.state.system.q .= bot.traj.q[step]
                RB.update!(suptg)
                (;r,b,g) = cg[istep]
                viz!(ax,Observable(suptg);
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
    hidex(ax1)
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
    savefig(fig,"spinningtop_drop")
    fig
end

topq_longtime = deepcopy(topq)
RB.set_new_initial!(topq_longtime,topq.traj.q[end],topq.traj.q̇[end])
RB.solve!(
    RB.SimProblem(topq_longtime,top_contact_dynfuncs),
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
            RB.SimProblem(bot,(x)->top_contact_dynfuncs(x;checkpersist=check)),
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
        resolution = (1tw,0.25tw),
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
            RB.SimProblem(topbot,(x)->top_contact_dynfuncs(x;checkpersist=check)),
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
        resolution = (1.0tw,0.45tw),
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
        resolution = (0.95tw, 0.18tw),
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
     RB.solve!(RB.SimProblem(top,top_contact_dynfuncs),
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
        resolution = (0.7tw,0.2tw),
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

#-- spinning top

#--  painleve's Bar
function make_bar(;
        μ = 0.1,
        e = 0.0
    )
    movable = true
    constrained = false
    m = 0.402
    l = 1.0
    Ixx = 1/12*m*l^2
    Ī = SMatrix{3,3}(Diagonal([Ixx,0.0,0.0]))
    mass_locus = SVector{3}(0.0,0.0,0.0)
    loci = SVector{3}.([
        [ -l/2, 0.0, 0.0],
        [  l/2, 0.0, 0.0]
    ])
    θ = π/4
    zoffset = 1e-6
    xoffset = 1e-6
    ri = [     0.0-xoffset,0.0,sin(θ)*l-zoffset]
    rj = [cos(θ)*l-xoffset,0.0,     0.0-zoffset]
    origin_position = (ri + rj)./2
    origin_velocity = zero(origin_position)
    ω = zero(origin_position)
    R = Matrix(RotY(θ))

    # b = l/2*ω[2]*sin(θ)-9.81
    # p⁺ = 1 + 3*cos(θ)^2-μ*cos(θ)*sin(θ)
    # @show b,p⁺

    prop = RB.RigidBodyProperty(1,movable,m,Ī,mass_locus,loci;constrained)

    nmcs = RB.NCF.NC3D2P(ri,rj,origin_position,R,origin_velocity,ω)
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω)

    p1 = Meshes.Point(loci[1])
    p2 = Meshes.Point(loci[2])
    s = Meshes.Segment(p1,p2)
    cyl_bar = Meshes.Cylinder(Meshes.length(s)/40,s)
    cylsurf_bar = Meshes.boundary(cyl_bar)
    cyl_bar_simple = Meshes.triangulate(cylsurf_bar)
    cyl_bar_mesh = cyl_bar_simple |> simple2mesh
    rb1 = RB.RigidBody(prop,state,cyl_bar_mesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    connections = RB.connect(rbs,zeros(Int,0,0))
    cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
    st = RB.Structure(rbs,tensiles,cnt,)
    bot = RB.Robot(st)
end

bar = make_bar(;)

plot_traj!(bar;)

function bar_contact_dynfuncs(bot)
end


tspan = (0.0,0.4)
h = 1e-3


prob = RB.SimProblem(bar,bar_contact_dynfuncs)

RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=50,exception=false)

contacts_traj_voa = VectorOfArray(bar.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,1:end])
 check_Coulomb(i,c)
end

plot_traj!(bar;showinfo=false)

me = RB.mechanical_energy!(bar)
lines(me.E)

[c[1].state.active for c in bar.contacts_traj] |> lines
[c[1].state.persistent for c in bar.contacts_traj] |> lines

[c[2].state.active for c in bar.contacts_traj] |> lines
[c[2].state.persistent for c in bar.contacts_traj] |> lines
#-- painleve's bar

#-- meteor hammer

cablemesh = sample(ancs,e,1000)

mesh(cablemesh,transparency=false)

function make_hammer(id,r̄ijkl,origin_position,R,ri,rj=nothing,rk=nothing,rl=nothing;
        μ,e,
        ω = zero(origin_position),
        movable = true,
        constrained = false,
        pres_idx = Int[],
        cstr_idx = collect(1:6),
    )
    # free_idx = collect(1:6)
    m = 4.58794901
    Īg = SMatrix{3,3}(
        Matrix(Diagonal([
            0.00522676,
            0.00512073,
            0.00512073,
            ])
        )
    )
    mass_locus = SVector(-0.20,0,0)
    # loci = deepcopy(r̄ijkl)
    loci = vcat(
        [
            SVector(0.22*cos(i*π/2),0.0,0.22*sin(i*π/2)) + mass_locus
            for i = -2:1
        ],
        [
            RotY( π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + mass_locus
            for i = 1:6
        ],
        [
            RotY(-π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + mass_locus
            for i = [1,2,4,5]
        ]
    )
    axes = [
        SVector(1.0,0,0) 
        for _ in eachindex(loci)
    ]
    friction_coefficients = fill(μ,length(loci))
    restitution_coefficients = fill(e,length(loci))
    @show m,diag(Īg),mass_locus,length(loci),friction_coefficients,restitution_coefficients
    prop = RB.RigidBodyProperty(
                id,movable,m,Īg,
                mass_locus,loci,axes,
                friction_coefficients,restitution_coefficients;
                constrained=constrained
                )
    origin_velocity = zero(origin_position)
    if rj isa Nothing
        nmcs = RB.NCF.NC1P3V(ri,origin_position,R,origin_velocity,ω)
    elseif rk isa Nothing
        nmcs = RB.NCF.NC2P2V(ri,rj,origin_position,R,origin_velocity,ω)
    elseif rl isa Nothing
        nmcs = RB.NCF.NC3P1V(ri,rj,rk,origin_position,R,origin_velocity,ω)
    else
        nmcs = RB.NCF.NC4P(ri,rj,rk,rl,origin_position,R,origin_velocity,ω)
    end
    # cf = RB.NCF.CoordinateFunctions(nmcs,q0,pres_idx,free_idx)
    # @show typeof(nmcs)
       # trimesh = Meshes.discretize(box,) |> simple2mesh
    meteormesh = load(joinpath(RB.assetpath("流星锤.STL"))) |> make_patch(;
        scale = 1/400,
        trans = [-0.20,0,0],
    )
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω,pres_idx,cstr_idx)
    body = RB.RigidBody(prop,state,meteormesh)
end

function cable_ancf(pres_idx, 𝐞, L = 1.0) 
    radius = 2.0e-3
    # ancs = RB.ANCF.ANC3DRURU(8.96e3;E=110e6,L,radius)
    # mat_cable = filter(
    #     row->row.name == "Nylon 66",
    #     material_properties
    # )[1]
    # ρ = ustrip(Unitful.kg/Unitful.m^3,mat_cable.density)
    # E = ustrip(Unitful.Pa,mat_cable.modulus_elas)
    ρ = 1.03e3
    E = 0.2e9
    @myshow ρ, E
    ancs = RB.ANCF.ANC3DRURU(ρ;E,L,radius)
    mass = RB.ANCF.build_mass(ancs)
    @show mass
    T = typeof(L)
    mass_locus = SVector(L/2,T(0),T(0))
    r̄p1 = SVector(T(0),T(0),T(0))
    # r̄p2 = SVector(L,T(0),T(0))
    loci = [
        r̄p1,
    ]
    prop = RB.FlexibleBodyProperty(
        1,
        :cable,
        mass,
        mass_locus,
        # length(loci),
        loci
    )
    # cache = RB.get_CoordinatesCache(prop,ancs,𝐞)
    state = RB.FlexibleBodyState(prop,ancs,𝐞;pres_idx)
    fb = RB.FlexibleBody(prop,state)
end

function make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.5*√3, 0.0, 0.5-1e-2],
        rjx = [ 0.0,-1.0, 0.0],
        doDR=false,
        μ = 0.5,
        e = 0.9,
        L = 1.0,
        nx = 2,
        R = RotY(deg2rad(-60)),
        ω = zero(ri)
    )
    if doDR
        fb_pres_idx = [7,8,9]
        rb_pres_idx = [1,2,3]
    else
        fb_pres_idx = Int[]
        rb_pres_idx = Int[]
    end
    𝐞 = [ri;rix;rj;rjx]
    fb1 = cable_ancf(fb_pres_idx, 𝐞, L)
    # @show fb1.prop.mass
    subfbs,subsm = RB.subdivide(fb1,nx)
    r̄ijkl = SVector{3,Float64}.(
        0.1 .*[
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
        ]
    )
    rigidbody = make_hammer(
        nx+1,
        r̄ijkl,
        rj,
        R,
        rj;
        μ,e,
        ω,
        pres_idx=rb_pres_idx,
        constrained=ifelse(!isempty(rb_pres_idx),true,false)
    )
    fbs = TypeSortedCollection(vcat(subfbs,rigidbody))
    # fbs = TypeSortedCollection([fb1,])
    numberedpoints = RB.number(fbs)
    # indexedcoords = RB.index(fbs,)
    sm = zeros(Int,size(subsm,1)+3,size(subsm,2)+1)
    sm[1:size(subsm,1),1:size(subsm,2)] = subsm
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)  ] = 7:9
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)+1] = 1:3
    # display(sm)
    indexedcoords = RB.index(fbs,sm)
    ss = Int[]
    tensiles = (cables = ss,)
    connected = RB.connect(fbs,)
    tensioned = @eponymtuple(connected,)
    cst1 = RB.FixedIndicesConstraint(1,[1,2,3],ri)
    jointed = RB.join((cst1,),indexedcoords)
    cnt = RB.Connectivity(numberedpoints,indexedcoords,tensioned,jointed)
    st = RB.Structure(fbs,tensiles,cnt,)
    bot = RB.Robot(st)
end

function flexcable_contact_dynfuncs(bot,ground_plane)
    contact_dynfuncs(bot;
        flatplane = ground_plane
    )
end

# don't bother
R = RotXY(deg2rad(0),deg2rad(-0))
inclined_plane = RB.Plane(R*[0,0,1.0],[-1.039230,-0.0,0.781])
p1 =  [-1.039230,0,0.780]
p2 = [-1.173952,-0.11,0.8646278]
p3 = [-0.901447,-0.112500,0.762216]
n = (p2-p1) × (p3 - p1 )

inclined_plane = RB.Plane(n,p3.+[0,0,0.001])


flexcable_DR = make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
RB.GDR!(flexcable_DR;β=2e-4,maxiters=2e5,verbose=false)
flexcable = make_flexcable(;
    ri  = [ 0.0, 0.0, 1.5],
    rix = [ 0.0, 0.0,-1.0],
    rj  = [-0.6*√3, 0.0, 1.2],
    rjx = [ 0.0,-1.0, 0.0],
    e = 0.0,
    μ=0.1,L=1.2,nx=5,R=RotY(-π/2),
    ω=-200n,
)

# flexcable_DR = make_flexcable(;R = RotY(-π/2),L=1.5,nx=5,doDR=true)
# RB.GDR!(flexcable_DR;β=1e-3,maxiters=1e5)
# flexcable = make_flexcable(;R = RotY(-π/2),e=0.0,μ=0.01,L=1.5,nx=5)

flexcable_DR.traj.q[end][end-8:end] .= flexcable.traj.q[end][end-8:end]
flexcable_DR.traj.q̇[end][end-8:end] .= flexcable.traj.q̇[end][end-8:end]

RB.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])


with_theme(theme_pub;
        Poly= (
            transparency = true,
        ) 
    ) do 
    plot_traj!(
        flexcable;
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
    RB.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
    RB.ZhongCCP();
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

plotsave_contactpoints(flexcable,)

rp2 = RB.get_trajectory!(flexcable,6,10)
lines(rp2[2,:])

me = RB.mechanical_energy!(flexcable)
lines((me.E.-me.E[begin])./me.E[begin])

#dt
# dts = [1e-1,3e-2,1e-2,3e-3,1e-3,1e-4]
dts = [1e-2,5e-3,2e-3,0.99e-3,5e-4,2e-4,2e-5]
flexcables_dt = [
    begin
        flexcable_dt = deepcopy(flexcable)
        RB.solve!(
            RB.SimProblem(flexcable_dt,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
                RB.ZhongCCP();
                tspan=(0.0,0.12),dt,ftol=1e-14,maxiters=50,exception=false)
    end
    for dt in dts
]

_, err_avg = get_err_avg(flexcables_dt;bid=6,pid=10,di=1)

with_theme(theme_pub;
        resolution = (0.9tw,0.45tw),
        figure_padding = (fontsize,fontsize,0,0),
        Axis3=(
            azimuth = 7.595530633326987,
            elevation = 0.14269908169872403
        )
    ) do
    bot = flexcable
    (;t) = bot.traj
    bots = flexcables_dt
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
        # dorecord=true,
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
                viz!(ax,tgvis;
                    meshcolor=Makie.RGBAf(r,g,b,alphas[i])
                )
            end
            lines!(ax,rp1[begin:laststep],color=:red)
            lines!(ax,rp12[begin:laststep],color=:blue)
            lines!(ax,rp10[begin:laststep],color=:purple)
            # lines!(ax,rp2[:,506:end],color=:blue)
            handlemesh = load(joinpath(RB.assetpath(),"把柄.STL")) |> make_patch(;
                scale = 1/400,
                trans = [ 0.0, 0.0, 1.55],
                rot = RotY(-π/2)
            )
            hidey(ax)
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
    # savefig(fig,"meteor_sliding")
    fig
end

# swing
flexcable_DR = make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
RB.GDR!(flexcable_DR;β=2e-4,maxiters=2e5,verbose=true)
flexcable = make_flexcable(;
    ri  = [ 0.0, 0.0, 1.5],
    rix = [ 0.0, 0.0,-1.0],
    rj  = [-0.6*√3, 0.0, 1.2],
    rjx = [ 0.0,-1.0, 0.0],
    e = 0.5,
    μ = 0.9,
    L=1.2,nx=5,R=RotY(-π/2)
)
# flexcable_DR.traj.q[end][end-8:end] .= flexcable.traj.q[end][end-8:end]

RB.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])

plane_normal = RotZY(deg2rad(-7.5),deg2rad(-90))*[0,0,1.0]
inclined_plane = RB.Plane(plane_normal,[0.05,0,0])

plane_normal'*rp+inclined_plane.d

GM.activate!();with_theme(theme_pub;
        Poly= (
            transparency = true,
        ) 
    ) do 
    plot_traj!(flexcable;
    ground=inclined_plane
    )
end

h = 2e-4
tspan = (0.0,1.5)
RB.solve!(
    RB.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
    RB.ZhongCCP();
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

me = RB.mechanical_energy!(flexcable)
lines(flexcable.traj.t,(me.E.-me.E[begin])./me.E[begin])
rp2 = RB.get_trajectory!(flexcable,6,2)
vp2 = RB.get_velocity!(flexcable,6,2)
pnvp2 = [plane_normal'*v for v in vp2 ]
lines(pnvp2)

plotsave_contact_persistent(flexcable,tol=0)

with_theme(theme_pub;
        resolution = (0.95tw,0.5tw),
        figure_padding = (0,fontsize,0,0),
        Axis3=(
            azimuth = 8.125530633326981,
            elevation = 0.18269908169872404
        ),
        Poly=(
            transparency=true,
        )
    ) do
    bot = flexcable
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
        # dorecord=true,
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
                viz!(ax,tgvis;meshcolor=Makie.RGBA(r,g,b,alphas[i]))
            end
            lines!(ax,rp2[:,begin:4001],color=:blue)
            handlemesh = load(joinpath(RB.assetpath(),"把柄.STL")) |> make_patch(;
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
        # aspect = DataAspect()
    )
    impact_time = 2995*h
    step_before_impact = time2step(impact_time,t)
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
        # aspect = DataAspect()
    )
    vlines!(ax22,[impact_time];linestyle=:dash)
    lines!(ax22,t,pnvp2,color=:red)
    ylims!(ax22,-7,7)
    xlims!(ax22,0.45,0.75)
    # axislegend(ax2)
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
    # axislegend(ax3,position=:lc)
    colsize!(fig.layout,1,0.35tw)
    colsize!(gd3,1,0.3tw)
    rowgap!(fig.layout,0)
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
    savefig(fig,"meteor_swing")
    fig
end

#-- meteor end

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

prob = RB.SimProblem(unibot,uni_dynfuncs)

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
            RB.SimProblem(
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
            resolution = (tw,0.4th),
            figure_padding = (fontsize,fontsize,0,0),

        ) do
        fig = Figure()
        grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
        for (botid,bot) in enumerate(bots)
            ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})")
            ṙp2_mid = get_mid_velocity!(bot,1,2)
            t_mids = get_time_mids(bot)
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
    RB.SimProblem(
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
    resolution = (0.9tw,0.28tw),
    cid=2
    )

CM.activate!(); plotsave_friction_direction(
        [unibot_e5],L"e=",[0.5], 
        "unibot_friction_direction"; 
        resolution = (0.9tw,0.28tw)
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
            t_mids = get_time_mids(bot)
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
        RB.SimProblem(
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

GM.activate!(); plotsave_error(unibot_z0,dts;resolution = (0.4tw,0.3tw),)
CM.activate!(); plotsave_error(unibot_z0,dts,"unibot_err";resolution = (0.4tw,0.3tw),)

GM.activate!(); plotsave_velocity_restitution([unibot_z0],true; ymids = [0.0])
plot_traj!(unibot_z0)

GM.activate!(); plotsave_point_traj_vel([unibot_z0])

me = RB.mechanical_energy!(unibot_z0)
me.E |> lines

GM.activate!(); plot_traj!(unibot_z0)
#-- uni bot end ---

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
    constrained = false,
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
prob = RB.SimProblem(ballbot,ball_dynfuncs)
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
        prob = RB.SimProblem(ballbot_dt,ball_dynfuncs)
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
        resolution = (1.0tw,0.45tw),
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
        resolution = (tw,0.5th),
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
    constrained = false,
    loadmesh = false,
)

# test rolling
tspan = (0.0,5.0)
h = 1e-2
prob = RB.SimProblem(ballbot,ball_dynfuncs)
@time RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=200,exception=false)

GM.activate!(); plotsave_contactpoints(ballbot)

plot_traj!(ballbot;auto=true)

GM.activate!(); with_theme(theme_pub;
        figure_padding = (0,0.5fontsize,0,0),
        resolution = (1.0tw,0.40tw),
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
        [ballbot,ballbot,ballbot],L"\mathrm{Cable~No.}",
        [3,6,12]; 
        resolution = (tw,0.5tw),
        mo_α = 8,
        mo = 8,
        vtol=1e-7,
)

CM.activate!(); plotsave_friction_direction(
        [ballbot,ballbot,ballbot],L"\mathrm{Cable~No.}",
        [3,6,12], 
        "ballbot_friction_direction"; 
        resolution = (tw,0.5tw),
        mo_α = 8,
        mo = 8,
        vtol=1e-7,
)


#----------- quadruped ---------------
quadbot = quad(10.0)
quadbot.st.connectivity.numbered

plot_traj!(quadbot;
    zlims = (0,1),
    showground = false
)

function quad_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0])
    )
end

tspan = (0.0,1.91)
tspan = (0.0,1.0)
h = 1e-3

prob = RB.SimProblem(quadbot,quad_dynfuncs)

RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true)

RB.solve!(prob,RB.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(quadbot;
    zlims = (0,1),
    showground = false
)

me = RB.mechanical_energy!(quadbot)

me.E |> lines
contacts_traj_voa = VectorOfArray(quadbot.contacts_traj)

#-- cube
function make_cube(origin_position = [0.0,0.0,0.5],
        R = one(RotMatrix{3}),
        origin_velocity = [2.0,0.0,0.0],
        ω = [0.0,0.0,5.0];
        μ = 0.9,
        e = 0.0
    )
    m = 1.0
    mass_locus = zeros(3)
    inertia = SMatrix{3,3}(Matrix(1.0I,3,3))
    h = 5.0
    loci = [h.*[x,y,z] for z = [-1,1] for y = [-1,1] for x = [-1,1]]
    axes = [SVector(1.0,0,0) for _ = 1:length(loci)]
    friction_coefficients = [μ for _ = 1:length(loci)]
    restitution_coefficients = [e for _ = 1:length(loci)]
    pts = Point3.(loci)
    fcs = GB.QuadFace.([
        [1,3,4,2],
        [3,7,8,4],
        [7,5,6,8],
        [5,1,2,6],
        [2,4,8,6],
        [1,5,7,3]
    ])
    nls = GB.normals(pts,fcs)
    cube_mesh = GB.Mesh(GB.meta(pts,normals=nls),fcs)
    prop = RB.RigidBodyProperty(
        1,true,m,inertia,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients,
        )
    ri = copy(origin_position)
    nmcs = RB.NCF.NC1P3V(ri,origin_position,R,origin_velocity,ω)
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω)
    rb1 = RB.RigidBody(prop,state,cube_mesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    hub = nothing
    #
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected)
    cnt = RB.Connectivity(numberedpoints,indexedcoords,tensioned)
    st = RB.Structure(rbs,tensiles,cnt,)
    bot = RB.Robot(st,nothing)
end

cubebot = make_cube()
plot_traj!(cubebot)

function cube_contact_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0])
    )
end

ros = [[0,0,z] for z = 10.0:-0.1:0.0]
qs = [
    begin
        cube = make_cube(origin_position,RotXY(π/4,π/4))
        q = RB.get_coords(cube.st)
    end
    for origin_position in ros
]
# p = RB.generalized_α(1.0)
Ro = [
    0.368869  -0.824063   0.429949;
    0.769524   0.011314  -0.638518;
    0.521314   0.566385   0.63831
]
origin_velocity = [0.0,0.0,0.0]
ωo = [3.0,4.0,5.0]
# cube = make_cube()
# cube = make_cube(ros[1])
cube = make_cube(ros[2],Ro,origin_velocity,0.1.*ωo)
tspan = (0.0,1.91)
tspan = (0.0,5.0)
dt = 1e-3

prob = RB.SimProblem(cube,cube_contact_dynfuncs)
RB.solve!(prob,
    RB.ZhongCCP();tspan,dt,ftol=1e-10,maxiters=50,exception=false
)

plot_traj!(cube;
    rigidcolor=:white,
    xlims = [-10,10],
    ylims = [-10,10],
    zlims = [-1e-3,10],
    showinfo=false,
    showwire=true
)

ME = RB.mechanical_energy!(cube)
lines(cube.traj.t,ME.E)

ts,qs,vs,as,ṽ̇s = RB.NSSFC.nssfc(nq,nλ,q0,v0,p,h,contact_dynfuncs(cube),tspan;tol=1e-6,imax=10)

ts,cs,qs,vs,ps,λs,friction_coefficients = RB.NSSFC.nhsolve(nq,nλ,nμ,q0,v0,contact_dynfuncs(cube),tspan;
            dt=h,ftol=1e-6,imax=50,exception=true)

ts,cs,qs,vs,ps,λs,friction_coefficients = RB.NSSFC.ipsolve(nq,nλ,nμ,q0,v0,contact_dynfuncs(cube),tspan;
            dt=h,ftol=1e-10,imax=50,exception=false)

ts,qs,vs,as,ṽ̇s = RB.NSSFC.nssfc(n,b,q26026,v26026,p,h,contact_dynfuncs(cube),tspan;tol=1e-10,imax=100)

cuberef = deepcopy(cube)
prob = RB.SimProblem(cuberef,dynfuncs,tspan)
RB.solve!(prob,RB.Zhong06();dt=h,ftol=1e-14)

vis(cuberef,contact_dynamics)

cube = make_cube(ros[1],Ro,origin_velocity,ωo)
# cube = make_cube(ros[1])
# cube = make_cube()

RB.prepare_traj!(cube.traj;tspan,dt=h,restart=true)

cube.traj.t[2:end] = ts[2:end]
cube.traj.q[2:end] = qs[2:end]
cube.traj.q̇[2:end] = vs[2:end]

vis(cube,contact_dynamics)

vis(cube,contact_dynamics;do_record=true)

plot(cube.traj.ts,VectorOfArray(cube.traj.qs)[7,:],color="yellow")
plot!(cuberef.traj.ts,VectorOfArray(cuberef.traj.qs)[7,:],color="blue")

ME = RB.mechanical_energy!(cube;gravity=true,)

ME.E
scatter(ME.E)
kes,epes,gpes,restitution_coefficients,es_err = analyse_energy(cube;gravity=true,elasticity=false)
es_off = OffsetArray(restitution_coefficients,0:length(restitution_coefficients)-1)
kes_off = OffsetArray(kes,0:length(restitution_coefficients)-1)
gpes_off = OffsetArray(gpes,0:length(restitution_coefficients)-1)
steprange = 0:length(restitution_coefficients)-1
steprange = 558:560
steprange = 5574:5600
steprange = 55740:56000
es1 = es_off[steprange[1]]
plot(steprange,(es_off[steprange].-es1)./es1)
plot(steprange,kes_off[steprange])
plot(steprange,gpes_off[steprange])
smoothsteps = findall((x)->x==0,cs)
plotsmoothsteps = [x for x in smoothsteps if x∈steprange]
plot!(plotsmoothsteps,es_off[plotsmoothsteps],color=:blue)

-0.0016049537616469219-(-0.0016048102389008441)

plot(epes[6400:6500])
ylims!(78,80)

plot(ts,kes[begin+1:end])
plot(ts,gpes[begin+1:end])

ros = [q[1:3] for q in qs]
ṙos = [q̇[1:3] for q̇ in vs]
r̈os = [ṽ̇[1:3] for ṽ̇ in ṽ̇s]
aos = [a[1:3] for a in as]

plot(ts,norm.(ros))

plot(ts,norm.(ṙos))

plot!(ts,norm.(ros))

plot(ts,norm.(r̈os))


plot(ts,norm.(aos))

plot!(ts,gpes[begin+1:end])

norm(ṙos[1393])/norm(ṙos[1392])

y_split = [[0.09510319720143898, -0.008672570657933804, 0.09470694080223288], [0.23134216351905035, 0.1781735286559687, 0.10655619736781423], [0.25492807263830397, 0.18779000020080971, -0.14810414657700796]]
[yi[1] - norm(yi[2:3]) for yi in y_split]
[norm(yi[2:3]) for yi in y_split]
Λ_split = [[6.077211920457062e-5, 5.541879908804313e-6, -6.051890646476349e-5], [1.6851978930926635e-17, -1.4431551427742453e-17, -7.662097843905011e-18], [1.0182657945568905e-17, -7.963039923870281e-18, 5.867035717865434e-18]]
[Λi[1] - norm(Λi[2:3]) for Λi in Λ_split]
[Λi[1] for Λi in Λ_split]
[norm(Λi[2:3]) for Λi in Λ_split]
