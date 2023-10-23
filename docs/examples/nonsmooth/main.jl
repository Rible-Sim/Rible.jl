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
using Test
using IterTools
using Unitful
using Match
using FileIO
using Cthulhu
using JET
using TypedTables
using Revise
using AbbreviatedStackTraces
import TensegrityRobots as TR
import Meshes
cd(@__DIR__)
include("def.jl"); includet("def.jl")
# includet("plotting.jl")
include("../analysis.jl"); includet("../analysis.jl")
include("../vis.jl"); includet("../vis.jl")
include("../dyn.jl"); includet("../dyn.jl")
figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME"
elseif Sys.isapple()
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME"
end
#-- preamble end

#-- point Mass

function new_pointmass(;
        ro = [0.0,0.0,1.0],
        ṙo = zero(ro),
        m = 1.0,
        μ = 0.3,
        e = 0.9
    )
    movable = false
    constrained = true
    Ia = SMatrix{3,3}(Matrix(m*I,3,3))
    r̄g  = SVector{3}([ 0.0, 0.0, 0.0])
    r̄p1 = SVector{3}([ 0.0, 0.0, 0.0])
    r̄ps = [r̄p1]
    ās = [SVector{3}([ 1.0, 0.0, 0.0])]
    μs = [μ]
    es = [e]
    prop = TR.RigidBodyProperty(
        1,movable,m,Ia,
        r̄g,r̄ps,ās,
        μs,es,
        ;constrained=constrained
    )
    R = RotX(0.0)
    ω = zero(ro)
    rps = Ref(ro) .+ Ref(R).*r̄ps
    lncs, _ = TR.NCF.NC3D1P(rps[1],)
    ci = Int[]
    Φi = Int[]
    state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,ci,Φi)
    rb1 = TR.RigidBody(prop,state)

    rbs = TypeSortedCollection((rb1,))
    numberedpoints = TR.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = TR.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    connected = TR.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)
    cnt = TR.Connectivity(numberedpoints,indexedcoords,tensioned)
    contacts = [TR.Contact(i,μ,e) for i = [1]]
    tg = TR.TensegrityStructure(rbs,tensiles,cnt,contacts)
    bot = TR.TensegrityRobot(tg)
end

function pm_contact_dynfuncs(bot;θ=0.0)
    a = tan(θ)
    n = [-a,0,1]
    inclined_plane = TR.Plane(n,zeros(3))
    contact_dynfuncs(bot;flatplane = inclined_plane)
end

tspan = (0.0,0.455)
tspan = (0.0,1.5)
h = 1e-3

# horizontal plane
es = [0.5]
v0s = [1.0]
pms_hp = [
    begin
        pm = new_pointmass(;
                e,μ=0.1,ṙo = [v0,0,0]
            )
        TR.solve!(
            TR.SimProblem(pm,pm_contact_dynfuncs),
            TR.ZhongCCP();
            tspan,dt=h,
            ftol=1e-14,
            maxiters=50,
            exception=false
        )
    end
    for v0 in v0s for e in es
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
                suptg = deepcopy(bot.tg)
                suptg.state.system.q .= bot.traj.q[step]
                TR.update!(suptg)
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
    vp1 = get_velocity!(bot,1,1)
    lines!(ax1,t,vp1[3,:])
    xlims!(ax1,t[begin],t[end])
    # ylims!(ax1,-6,6)
    Label(gd1[1,1,TopLeft()], 
        rich("($(alphabet[2]))",font=:bold)
    )
    hlines!(ax1,[0],color=:gray)
    ax2 = Axis(gd1[1,2], xlabel = tlabel, ylabel = L"\dot{x}~(\mathrm{m/s})")
    # me = TR.mechanical_energy!(bot)
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

function plotsave_pointmass_xz_energy(bots,figname=nothing)
    with_theme(theme_pub;
            resolution = (0.9tw,0.45tw),
        ) do
        fig = Figure()
        grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
        for (botid,bot) in enumerate(bots)
            (;t) = bot.traj
            ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = "Coords. (m)")
            rp1 = get_trajectory!(bot,1,1)
            lines!(ax1,t,rp1[3,:],label=L"z")
            lines!(ax1,t,rp1[1,:],label=L"x")
            if botid == 2
                axislegend(ax1,position=:ct,orientation=:horizontal)
            end
            xlims!(ax1,t[begin],t[end])
            Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
            ax2 = Axis(grids[botid,2][1,1], xlabel = tlabel, ylabel = L"\mathrm{Energy}~(\mathrm{J})")
            me = TR.mechanical_energy!(bot)
            lines!(ax2,t,me.E)
            xlims!(ax2,t[begin],t[end])
            ylims!(ax2,-1,11)
            if botid !== 4
                hidex(ax1); hidex(ax2)
            end
            Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
        end
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plotsave_pointmass_xz_energy(reshape(pms_hp,1,1)[:,2])
CM.activate!(); plotsave_pointmass_xz_energy(reshape(pms_hp,3,2)[:,2],"pointmass_xz_energy")

function plotsave_pointmass_velocity_restitution(bots,figname=nothing)
    with_theme(theme_pub;
            resolution = (0.9tw,0.6tw),
        ) do
        fig = Figure()
        grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
        ymids = es
        xtickmaxs = [4,5,1]
        for (botid,bot) in enumerate(bots)
            (;t) = bot.traj
            ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})")
            ṙp1 = get_velocity!(bot,1,1)
            lines!(ax1,t,ṙp1[3,:])
            xlims!(ax1,t[begin],t[end])
            ylims!(ax1,-6,6)
            Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
            ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e")
            vz = ṙp1[3,:]
            pos_pks_raw, _ = findmaxima(vz)
            pos_pks, _ = peakproms(pos_pks_raw,vz,minprom=0.2)
            neg_pks_raw, _ = findminima(vz)
            neg_pks, _ = peakproms(neg_pks_raw,vz,minprom=0.2)
            if (botid == length(bots))
                pos_pks = neg_pks.+1
            elseif length(neg_pks) > length(pos_pks)
                pop!(neg_pks)
            end
            ratios = -vz[pos_pks]./vz[neg_pks]
            scatter!(ax2,ratios)
            ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
            # ax2.yticks = [0,0.3,0.7,1.0]
            ax2.xticks = collect(1:xtickmaxs[botid])
            Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
        end
        rowgap!(fig.layout,0)
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,1])
CM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,1],"pointmass_v00_velocity_restitution")
GM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,2])
CM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,2],"pointmass_v10_velocity_restitution")

θ = 15 |> deg2rad
inclined_plane = TR.Plane([-tan(θ),0,1],zeros(3))
ro = [0.0,0,-1e-7]
ṙo = [2.0cos(θ),0,2.0sin(θ)]
# ṙo = [0,-2.0,0]

# analytical
g = 9.81
μ=0.3
vo = norm(ṙo)
a = -μ*g*cos(θ)-g*sin(θ)
# μ*g*cos(θ)-g*sin(θ)
tf = -vo/a
# d(t) -> vo*t+1/2*a*t^2
tspan = (0.0,0.6)
pm = new_pointmass(;e=0.0, μ, ro, ṙo)

prob = TR.SimProblem(pm,(x)->pm_contact_dynfuncs(x;θ))
TR.solve!(prob,TR.ZhongCCP();tspan,dt=1e-3,ftol=1e-14,maxiters=50,exception=false)

rp1 = get_trajectory!(pm,1,1)
ṙp1 = get_velocity!(pm,1,1)
dp1 = rp1.u .|> norm
vl1 = [u ⋅ normalize(ṙo) for u in ṙp1]
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
                suptg = deepcopy(bot.tg)
                suptg.state.system.q .= bot.traj.q[step]
                TR.update!(suptg)
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
    rp1 = get_trajectory!(bot,1,1)
    ṙp1 = get_velocity!(bot,1,1)
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
GM.activate!(); plotsave_dis_vel(pm)
CM.activate!(); plotsave_dis_vel(pm,"pointmass_dis_vel")

fig = Figure()
ax = Axis3(fig[1,1])
#dt
dts = [1e-1,3e-2,1e-2,3e-3,1e-3,1e-4]
pms = [
    begin
        pm = new_pointmass(;e=0.0, μ=0.3, ro, ṙo)
        TR.solve!(
            TR.SimProblem(pm,(x)->pm_contact_dynfuncs(x;θ)),
            TR.ZhongCCP();tspan,dt,ftol=1e-14,maxiters=50,exception=false
        )
    end
    for dt in dts
]

GM.activate!(); plotsave_error(pms,dts,pid=1,di=1)

GM.activate!(); plotsave_contact_persistent(pm)


me = TR.mechanical_energy!(pm); me.E[500:end] |> lines

me = TR.mechanical_energy!(pm); @show me.E[end]

contacts_traj_voa = VectorOfArray(pm.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,:])
check_Coulomb(i,c)
end

# contact_friction
function plotsave_pointmass_energy_friction(bot,figname=nothing)
    with_theme(theme_pub;
        resolution = (0.9tw,0.3tw)
    ) do
        (;t) = bot.traj
        fig = Figure()
        ax1 = Axis(fig[1,1], xlabel = tlabel, ylabel = "Energy (J)")
        Label(fig[1,1,TopLeft()], "($(alphabet[1]))")
        ax2 = Axis(fig[1,2], xlabel = tlabel, ylabel = L"\lambda~(\mathrm{N})")
        Label(fig[1,2,TopLeft()], "($(alphabet[2]))")
        me = TR.mechanical_energy!(bot)
        (;T,V,E) = me
        lines!(ax1,t,me.T, label="Kinetic")
        lines!(ax1,t,me.V, label="Potential")
        lines!(ax1,t,me.E, label="Total")
        xlims!(ax1,extrema(t)...)
        axislegend(ax1)

        c1_traj = VectorOfArray(bot.contacts_traj)[1,5:end]
        Λn = [c.state.Λ[1]./c.μ for c in c1_traj]
        Λt = [c.state.Λ[2:3] |> norm for c in c1_traj]
        Λt_mid = (Λt[2:2:end]+Λt[1:2:end-1])./2
        th = t[5:end]
        lines!(ax2,th,Λn,label=L"\lambda_{n}")
        lines!(ax2,(th[2:2:end]+th[1:2:end-1])./2,Λt_mid,label=L"\lambda_{t}")
        # lines!(ax2, th, Λt ,label=L"\lambda_{t}")
        @show Λn[end], Λt_mid[begin], Λt_mid[end]
        axislegend(ax2,position=:rc)
        xlims!(ax2,extrema(t)...)
        ylims!(ax2,-5,10)
        savefig(fig,figname)
        fig
    end
end
GM.activate!(); plotsave_pointmass_energy_friction(pm)
CM.activate!(); plotsave_pointmass_energy_friction(pm,"pointmass_energy_friction")
#-- point mass end

#--  Spinning top

function make_top(ro = [0.0,0.0,0.0],
        R = one(RotMatrix{3}),
        ṙo = [0.0,0.0,0.0],
        Ω = [0.0,0.0,5.0],
        cT = TR.QBF.QC;
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
    r̄g = @SVector zeros(3)
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
    r̄ps = [radius.*[1,1,0] for i = 1:4]
    push!(r̄ps,[0,0,-h])
    ās = [SVector(1.0,0,0) for i = 1:5]
    μs = [μ for i = 1:5]
    es = [e for i = 1:5]
    if loadmesh
        topmesh = load(
            TR.assetpath("Toupise2.STL")
        ) |> make_patch(;
            scale=1/1000,
            color,
        )
    else
        pts = Point3.(r̄ps)
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
    prop = TR.RigidBodyProperty(
        1,movable,m,Ī,
        r̄g,r̄ps,ās,
        μs,es;
        constrained
    )
    ri = ro+R*r̄ps[5]
    @myshow ri
    if cT == TR.QBF.QC
        nmcs = TR.QBF.QC(m,Ī)
    else
        nmcs, _ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
    end
    state = TR.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,pres_idx)
    rb1 = TR.RigidBody(prop,state,topmesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = TR.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = TR.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    connected = TR.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)
    cnt = TR.Connectivity(numberedpoints,indexedcoords,tensioned)
    tg = TR.TensegrityStructure(rbs,tensiles,cnt,)
    bot = TR.TensegrityRobot(tg)
end

function top_contact_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = TR.Plane([0,0,1.0],[0,0,0.0])
    )
end

ro = [0,0,0.5]
R = RotX(0.0)
ṙo = [1.0,0.0,0.0]
Ω = [0.0,0.0,200.0]
# R = rand(RotMatrix3)
# ṙo = rand(3)
# Ω = rand(3)

# tspan = (0.0,1.0)
μ = 0.95
e = 0.5
tspan = (0.0,1.8)
h = 1e-4

topq = make_top(ro,R,ṙo,Ω;μ,e,loadmesh=true)
#note subsequent iteration slow convergence 
#note initial guess can not improve it?
TR.solve!(
    TR.SimProblem(topq,top_contact_dynfuncs),
    TR.ZhongQCCP();
    tspan,
    dt=h,
    ftol=1e-10,
    maxiters=50,exception=false,verbose_contact=true
)
me = TR.mechanical_energy!(topq)
lines(me.E)
rp5 = get_trajectory!(topq,1,5)

plot_traj!(
    topq;
    showinfo=false,
    # rigidcolor=:white,
    showwire=true,
    showarrows=false,
)

topn = make_top(ro,R,ṙo,Ω,TR.NCF.LNC;μ,e,loadmesh=true)
TR.solve!(
    TR.SimProblem(topn,top_contact_dynfuncs),
    TR.ZhongCCP();
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
    rp5 = get_trajectory!(bot,1,5)
    vp5 = get_velocity!(bot,1,5)
    me = TR.mechanical_energy!(bot)
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
                suptg = deepcopy(bot.tg)
                suptg.state.system.q .= bot.traj.q[step]
                TR.update!(suptg)
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
TR.set_new_initial!(topq_longtime,topq.traj.q[end],topq.traj.q̇[end])
TR.solve!(
    TR.SimProblem(topq_longtime,top_contact_dynfuncs),
    TR.ZhongQCCP();
    tspan = (0.0,50.0),
    dt=2e-3,
    ftol=1e-14,
    maxiters=100,exception=false,verbose=false
)
me = TR.mechanical_energy!(topq_longtime)
lines(me.E)
rp5 = get_trajectory!(topq_longtime,1,5)
lines(topq_longtime.traj.t,-rp5[3,:].-(-rp5[3,1]))

topn_longtime = deepcopy(topn)
TR.set_new_initial!(topn_longtime,topn.traj.q[end],topn.traj.q̇[end])
TR.solve!(
    TR.SimProblem(topn_longtime,top_contact_dynfuncs),
    TR.ZhongCCP();
    tspan = (0.0,500.0),
    dt=2e-3,
    ftol=1e-14,
    maxiters=100,exception=false,verbose=false
)
plotsave_contact_persistent(topn_longtime)
me = TR.mechanical_energy!(topn_longtime)
lines(me.E)
rp5 = get_trajectory!(topn_longtime,1,5)
lines(rp5)
plot_traj!(topn_longtime)
# vp5 = get_velocity!(topn,1,5)
with_theme(theme_pub;
        resolution = (0.7tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize/2)
    ) do
    bot = topn_longtime
    # bot = topq_longtime
    (;t) = bot.traj
    fig = Figure()
    ax1 = Axis(fig[1,1];xlabel=tlabel, ylabel="Rel. Err.")
    ax2 = Axis(fig[1,2];xlabel=tlabel, ylabel="Abs. Err. (m)")
    skipstep = 500
    startstep = time2step(12.1,t)
    @myshow length(t[startstep:skipstep:end])
    lines!(ax1,
        t[startstep:skipstep:end],
        (me.E[startstep:skipstep:end].-me.E[startstep])./me.E[startstep]
    )
    ylims!(ax1,-1e-4,1e-4)
    xlims!(ax1,extrema(t)...)
    @myshow rp5[3,startstep]
    lines!(ax2,
        t[startstep:skipstep:end],
        (-rp5[3,startstep:skipstep:end]).-(-rp5[3,startstep])
    )
    xlims!(ax2,extrema(t)...)
    Label(fig[1,1,TopLeft()], rich("($(alphabet[1]))",font=:bold))
    Label(fig[1,2,TopLeft()], rich("($(alphabet[2]))",font=:bold))
    savefig(fig,"spinningtop_longtime")
    fig    
end


c1_topq = get_trajectory!(topq,1,5)
c1_topn = get_trajectory!(topn,1,5)
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax,c1_topq[3,:])
lines!(ax,c1_topn[3,:])
fig

# no contact
prob = TR.SimProblem(top,(x)->dynfuncs(x;gravity=true))

TR.solve!(prob,TR.Zhong06Q();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true,verbose=false)

plot_traj!(top;showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw))

me = TR.mechanical_energy!(top)
me.E |> lines

contacts_traj_voa = VectorOfArray(top.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,end-10:end])
 check_Coulomb(i,c)
end


# contacts_traj = TR.solve!(prob,TR.AlphaCCP(0.95);tspan,dt=h,ftol=1e-8,maxiters=50,exception=true)

with_theme(theme_try;
     Axis3 = (
         azimuth = 4.2555306333269835,
         elevation = 0.2326990816987238
     )
 ) do
 plot_traj!(top;
     AxisType=Axis3,
     doslide=false,
     gridsize=(2,3),
     attimes=[0,2,4,6,8,10],
     xlims=(-0.5,1.5),
     zlims=(-1e-3,3),
     rigidcolor=:white,
     showwire=true,
     showlabels=false,
     fontsize=8 |> pt2px,
     figsize=(0.9tw,0.7tw),
     savefig=true,
     figname="spinningtop_traj"
 )
end

GM.activate!(); plotsave_contact_persistent(top)
CM.activate!(); plotsave_contact_persistent(top,"spinningtop_contact_persistent")

# friction_direction

R = RotX(π/18)
μs = [
    0.01,
]
ro = [0,0,0.037]
tops_e0 = [
    begin
        Ω = [0,0,50.0]
        top = make_top(ro,R,ṙo,Ω,TR.NCF.LNC; μ, e = 0.0,loadmesh=true)
        TR.solve!(
            TR.SimProblem(top,top_contact_dynfuncs),
            TR.ZhongCCP();
            tspan=(0.0,2.0),
            dt=1e-3,ftol=1e-14,maxiters=50,exception=false,#verbose_contact=false
        )
    end
    for μ in μs
]

plot_traj!(tops_e0[1])
vp5 = get_velocity!(tops_e0[1],1,5)

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
    rp5 = get_trajectory!(bot,1,5)
    rp1 = get_trajectory!(bot,1,1)
    vp5 = get_velocity!(bot,1,5)[stepstart:end]
    me = TR.mechanical_energy!(bot)[stepstart:end]
    steps = 1:100:1800
    nstep = length(steps)
    alphas = fill(0.1,nstep)
    alphas[1:3] = [1,0.4,0.2]
    alphas[end] = 1
    cg = cgrad(:winter, nstep, categorical = true)
    fig = Figure()
    gd1 = fig[1,1:3] = GridLayout()
    gd2 = fig[2,1:2] = GridLayout()
    gd3 = fig[2,3] = GridLayout(;tellheight=false)
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
                suptg = deepcopy(bot.tg)
                suptg.state.system.q .= bot.traj.q[step]
                TR.update!(suptg)
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
    ax3 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = "Abs. Err. (m)",
    )
    mo_rp5=13
    scaling = 10.0^(-mo_rp5)
    Label(gd2[1,1,Top()],latexstring("\\times 10^{-$(mo_rp5)}"))
    @myshow -rp5[3,stepstart]
    lines!(ax3,t,((-rp5[3,stepstart:end]).-(-rp5[3,stepstart]))./scaling)
    # hidex(ax3)
    ax4 = Axis(gd2[1,2],
        xlabel = tlabel,
        ylabel = L"\alpha-\pi~(\mathrm{Rad})",
    )
    mo_α=8
    scaling = 10.0^(-mo_α)
    Label(gd2[1,2,Top()],latexstring("\\times 10^{-$(mo_α)}"))
    contacts_traj_voa = VectorOfArray(bot.contacts_traj)[:,stepstart:end]
    c1s = contacts_traj_voa[end,:]
    idx_per = findall((x)->doespersist(x;Λtol=0),c1s) #∩ idx_sli
    # @show idx_per
    α_per = map(c1s[idx_per]) do c
        get_contact_angle(c;Λtol=0)
    end
    lines!(ax4,t[idx_per],abs.(α_per.-π)./scaling;)
    ax5 = Axis(gd3[1,1],
        xlabel = tlabel,
        ylabel = "Energy (J)",
    )
    lines!(ax5,t,me.E, label="E")
    lines!(ax5,t,me.T, label="T")
    lines!(ax5,t,me.V, label="V")
    Legend(gd3[1,2],ax5,orientation=:vertical,tellheight=false)
    # axislegend(ax5,position=:rt)
    xlims!(ax3,0,2.0)
    xlims!(ax4,0,2.0)
    xlims!(ax5,0,2.0)
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
        gd2[1,2,TopLeft()],
        rich("($(alphabet[4]))",font=:bold)
    )
    Label(
        gd3[1,1,TopLeft()],
        rich("($(alphabet[5]))",font=:bold)
    )
    # colsize!(fig.layout,1,Fixed(0.40tw))
    colgap!(fig.layout,1,3fontsize)
    rowgap!(gd1,0)
    rowgap!(fig.layout,0)
    rowsize!(fig.layout,2,Fixed(0.13tw))
    # rowsize!(gd1,2,Fixed(0.05tw))
    # rowgap!(gd1,0)
    # rowsize!(gd3,2,0.1tw)
    savefig(fig,"spinningtop_sliding")
    fig
end

function plotsave_friction_direction_error(bot,figname=nothing;
        resolution = (0.9tw,0.4tw),
        mo = 8,
        mo_α = mo,
        vtol=1e-5,
        Λtol=1e-5,
    )
    with_theme(theme_pub;
        resolution
    ) do
    fig = Figure()
    ax1 = Axis(fig[1,1],
                xlabel=tlabel,
                ylabel=L"\delta\alpha~(\mathrm{Rad})",
                # title=latexstring("$x=$(xs[1])")
        )
    ax2 = Axis(fig[1,2],
                xlabel=tlabel,
                ylabel=L"\delta\theta~(\mathrm{Rad})",
                # title=latexstring("$x=$(xs[1])")
        )
    # ax3 = Axis(fig[1,2],
    #         xlabel=L"\delta\theta~(\mathrm{Rad})",
    #         title=latexstring("$x=$(xs[1])"),
    #         xtickformat = values -> [latexstring("10^{$(value)}") for value in values]
            # xscale = Makie.log10,
            # xticks = [-1e-7, -1e-8, -1e-9, 0, 1e-9, 1e-8, 1e-7],
            # limits = ((1e-14, 1e-7), nothing),
    # )
    # Label(fig[1,1,TopLeft()],"($(alphabet[2botid-1]))")
    # Label(fig[1,2,TopLeft()],"($(alphabet[2botid]))")
    tstart = 309
    t = bot.traj.t[tstart:end]
    contacts_traj_voa = VectorOfArray(bot.contacts_traj)[:,tstart:end]
    c1s = contacts_traj_voa[1,:]
    # idx_sli = findall(c1s) do c
    #     issliding(c;vtol)
    # end
    idx_imp = findall((x)->isimpact(x;Λtol),c1s) #∩ idx_sli
    δα_imp = map(c1s[idx_imp]) do c
        get_contact_angle(c;Λtol)
    end
    idx_per = findall((x)->doespersist(x;Λtol=0),c1s) #∩ idx_sli
    # @show idx_per[begin]
    δα_per = map(c1s[idx_per]) do c
        get_contact_angle(c;Λtol)
    end
    rp5 = get_trajectory!(bot,1,5,tstart:length(bot.traj))
    scaling = 10.0^(-mo_α)
    # Label(fig[1,1,Top()],latexstring("\\times 10^{-$(mo_α)}"))
    markersize = fontsize
    lines!(ax1,δα_per.-π;)
    lines!(ax2,t[idx_per],rp5[3,idx_per].-rp5[3,begin])
    # xlims!(ax2,0.65,2)
    # scatter!(ax2,idx_per,ones(length(idx_per)))
    # scatter!(ax2,idx_imp,ones(length(idx_imp)))
    # rp5[3,idx_per])
    # # scatter!(ax1,t[idx_imp],δα_imp./scaling;
    # #             marker=:xcross, markersize)
    # ylims!(ax1,-1.0,1.0)

    # @myshow δα_per

    # θ_imp = map(c1s[idx_imp]) do c
    #         get_friction_direction(c)
    #     end
    # θ_per = map(c1s[idx_per]) do c
    #         get_friction_direction(c)
    #     end
    scaling = 10.0^(-mo)
    # Label(fig[1,2,Top()],latexstring("\\times 10^{-$(mo)}"))
    # scatter!(ax2,t[idx_per],(θ_per.-π/4)./scaling, label="Persistent";
    #             marker=:diamond, markersize)
    # # scatter!(ax2,t[idx_imp],(θ_imp.-π/4)./scaling, label="Impact";
    # #             marker=:xcross, markersize)
    # xlims!(ax2,extrema(t)...)
    # ylims!(ax2,-1.0,1.0)
    # hist!(ax3,collect(log10.(abs.(skipmissing(δα_per))));normalization  = :pdf)
    # hidex(ax3)
    # Legend(fig[1,:],ax3;
    #     orientation=:horizontal
    # )
    rowgap!(fig.layout,fontsize/2)
    savefig(fig,figname)
    fig
    end
end

GM.activate!(); plotsave_friction_direction_error(tops_e0[1];vtol=1e-7)

plotsave_contact_persistent(tops_e0[1],tol=0)
#dt

dts = [1e-2,5e-3,3e-3,2e-3,1e-3,5e-4,3e-4,2e-4,1e-4,1e-5]
tops_dt = [
 begin
    top = deepcopy(tops_e0[1])
     TR.solve!(TR.SimProblem(top,top_contact_dynfuncs),
             TR.ZhongCCP();
             tspan=(0.0,0.1),dt,ftol=1e-14,maxiters=50,exception=true)
 end
 for dt in dts
]

me = TR.mechanical_energy!(bot)
lines(me.E)
fig = Figure()
ax = Axis3(fig[1,1])
for top in tops_dt
    lines!(ax,get_trajectory!(top,1,1),label="$(top.traj.t[2])")
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
    
# contact point trajectory

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
                title=latexstring("\\mu = $(μs[i])")
            )
            ax2 = Axis(
                fig[i,2],
                xlabel=tlabel,
                ylabel=L"\theta",
                title=latexstring("\\mu = $(μs[i])")
            )
            Label(fig[i,1,TopLeft()], "($(alphabet[2i-1]))")
            Label(fig[i,2,TopLeft()], "($(alphabet[2i]))")
            ax3 = Axis(fig[i,2],
                xlabel=tlabel,
                ylabel=L"v~(\mathrm{m/s})",
                yticklabelcolor = :red,
                yaxisposition = :right
            )
            ax4 = Axis(fig[i,3],
                xlabel=tlabel,
                ylabel="Energy (J)",
            )
            hidespines!(ax3)
            hidexdecorations!(ax3)

            (;t) = bot.traj
            rp5 = get_trajectory!(bot,1,5)
            rpx = rp5[1,:]
            rpy = rp5[2,:]
            lines!(ax1,rpx,rpy)
            # xlims!(ax1,0,13.0)
            # ylims!(ax1,-0.3,0.45)

            ṙp5 = get_mid_velocity!(bot,1,5)
            ṙpx = ṙp5[1,:]
            ṙpy = ṙp5[2,:]
            t_mids = get_time_mids(bot)
            θ_mids = atan.(ṙpy,ṙpx)
            ṙpxy = map(ṙp5.u) do ṙ
                norm(ṙ[1:2])
            end
            lines!(ax2,t_mids,θ_mids)
            lines!(ax3,t_mids,ṙpxy,color=:red)
            # xlims!(ax2,extrema(bot.traj.t)...)
            # xlims!(ax3,extrema(bot.traj.t)...)
            # ylims!(ax2,-π,π)

            me = TR.mechanical_energy!(bot)
            lines!(ax4,t[309:end],me.E[309:end])
            lines!(ax4,t[309:end],me.T[309:end])
            lines!(ax4,t[309:end],me.V[309:end].+6)
            if i !== nbot
                hidex(ax1)
                hidex(ax2)
                hidex(ax4)
            end
            # xlims!(ax4,0.308,t[end])

            # vlines!(ax1,t[554],)
            # vlines!(ax2,t[554],)
        end
        savefig(fig,figname)
        fig
    end
end

GM.activate!(); plotsave_point_traj_vel(tops_e0)
CM.activate!(); plotsave_point_traj_vel(tops_e0,"contact_point_traj_vel")

ωs = [5.0,10.0,20.0]
tops_ω = [
 begin
     top = make_top([0,0,cos(π/24)*0.5-1e-4],R,[0,0,0.0],[0,0,ω]; μ=0.9, e = 0.0)
     TR.solve!(TR.SimProblem(top,top_contact_dynfuncs),
             TR.ZhongCCP();
             tspan=(0.0,500.0),dt=1e-2,ftol=1e-14,maxiters=50,exception=false)
 end
 for ω in ωs
]

plot_traj!(tops_ω[1];showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw),
 xlims=(-1.0,2.0),
 ylims=(-1.0,1.0),
 zlims=(-1e-3,1.0),
)


r1v5 = get_mid_velocity!(tops_ω[1],1,5)
r1v5[3,:] |> scatter

r1p5 = get_trajectory!(tops_ω[1],1,5)
r1p5[3,:] |> scatter

# energy conserving
function plotsave_energy_conserving(bots,figname=nothing)
 with_theme(theme_pub;
         resolution = (0.7tw,0.5tw)
     ) do
     fig = Figure()
     # mos = [4,5,6]
     mos = 12 .*[1,1,1]
     for (botid,bot) in enumerate(bots)
         ax1 = Axis(
             fig[botid,1],
             xlabel=tlabel,
             ylabel="Rel. Err.",
             title=latexstring("\\omega_z=$(ωs[botid])")
         )
         Label(fig[botid,1,TopLeft()],"($(alphabet[botid]))")
         if botid !== 3
             # hidex(ax1)
         end
         (;t) = bot.traj
         mo = mos[botid]
         scaling = 10.0^(-mo)
         Label(fig[botid,1,Top()],latexstring("\\times 10^{-$mo}"))
         me = TR.mechanical_energy!(bot)
         lines!(ax1,t,(me.E.-me.E[begin])./me.E[begin]./scaling)
         xlims!(ax1,extrema(t)...)
         ylims!(ax1,-1,1)
     end
     savefig(fig,figname)
     fig
 end
end
GM.activate!(); plotsave_energy_conserving(tops_ω)
CM.activate!(); plotsave_energy_conserving(tops_ω,"energy_conserving")

top_ω5 = make_top([0,0,cos(π/24)*0.5-1e-6],RotX(π/24),[0,0,0.0],[0,0,5.0]; μ=0.9, e = 0.0)
TR.solve!(TR.SimProblem(top_ω5,top_contact_dynfuncs),
     TR.ZhongCCP();
     tspan=(0.0,50.0),dt=1e-2,ftol=1e-14,maxiters=50,exception=false)

findall(isactive,VectorOfArray(top_ω5.contacts_traj)[1,:])

plot_traj!(top_ω5;showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw),
 xlims=(-1.0,2.0),
 ylims=(-1.0,1.0),
 zlims=(-1e-3,1.0),
)

# DAE
top_fix = make_top([0,0,cos(π/24)*0.5],RotX(π/24),[0,0,0.0],[0,0,5.0]; μ=0.9, e = 0.0, constrained=true)
TR.solve!(TR.SimProblem(top_fix,top_contact_dynfuncs),
     TR.Zhong06();
     tspan=(0.0,50.0),dt=1e-2,ftol=1e-14,maxiters=50,exception=false)

plot_traj!(top_fix;showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw),
 xlims=(-1.0,2.0),
 ylims=(-1.0,1.0),
 zlims=(-1e-3,1.0),
)

function plotsave_compare_traj(top1,top2,figname=nothing)
 ro_top1 = get_trajectory!(top1,1,0)
 ro_top2 = get_trajectory!(top2,1,0)
 ylabels = ["x","y","z"]
 with_theme(theme_pub;
         resolution = (0.8tw,0.5tw),
     ) do
     fig = Figure()
     for i = 1:3
         ax = Axis(fig[i,1]; xlabel = tlabel, ylabel = latexstring("$(ylabels[i])~(\\mathrm{m})"))
         Label(fig[i,1,TopLeft()],"($(alphabet[i]))")
         lines!(ax,top1.traj.t, ro_top1[i,:],label="NMSI")
         scatter!(ax,top2.traj.t[begin:25:end], ro_top2[i,begin:25:end],label="MSI",
                     marker='⨉', strokewidth=0, color = :red, strokecolor=:red)
         xlims!(ax,extrema(top2.traj.t)...)
         if i !== 3
             hidex(ax)
         else
             axislegend(ax;position=:rb,orientation=:horizontal)
         end
     end
     savefig(fig,figname)
     fig
 end
end
GM.activate!(); plotsave_compare_traj(top_fix,top_ω5)
CM.activate!(); plotsave_compare_traj(top_fix,top_ω5,"compare_traj")




plot_traj!(tops_dt_v[end];showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw),
 xlims=(-1.0,2.0),
 ylims=(-1.0,1.0),
 zlims=(-1e-3,1.0),
)

r1p5 = get_trajectory!(tops_dt[end],1,5)
lines(r1p5[1,:])
# friction_direction
# Dare you try
GM.activate!(); plotsave_friction_direction(tops_e0)

# contact_friction
function plotsave_contact_friction(bots,figname=nothing)
 with_theme(theme_pub;
     resolution = (0.9tw,0.5tw)
 ) do
     fig = Figure()
     nbots = length(bots)

     for (botid,bot) in enumerate(bots)
         ax = Axis(
             fig[botid,1],
             xlabel=tlabel,
             ylabel=L"\Lambda~(\mathrm{N}\cdot\mathrm{s})",
         )
         if botid !== nbots
             hidex(ax)
         end
         Label(fig[botid,1,TopLeft()],"($(alphabet[botid]))")
         (;t) = bot.traj
         c1s = VectorOfArray(bot.contacts_traj)[1,:]
         idx_per = findall(doespersist, c1s)
         Λt1 = [c.state.Λ[2] for c in c1s]
         Λt2 = [c.state.Λ[3] for c in c1s]
         lines!(ax,t[idx_per],Λt1[idx_per],label=L"\Lambda_{t_1}")
         lines!(ax,t[idx_per],Λt2[idx_per],label=L"\Lambda_{t_2}")
         # pos_normal_idx = findall((x)->x[1]>1e-4,c1_states.Λ)
         # active_pos_idx = intersect(idx_per,pos_normal_idx)
         # t_mids = 0.5 .* (t[active_pos_idx] .+ t[active_pos_idx .+ 1])
         # Λt1_mids = 0.5 .* (Λt1[active_pos_idx] .+ Λt1[active_pos_idx .+ 1])
         # Λt2_mids = 0.5 .* (Λt2[active_pos_idx] .+ Λt2[active_pos_idx .+ 1])
         # lines!(ax,t_mids,Λt1_mids,label=L"\Lambda_{t_1}")
         # lines!(ax,t_mids,Λt2_mids,label=L"\Lambda_{t_2}")
         # xlims!(ax,t[556],t[end])
     end
     # Legend(fig[:,3],axs[1])
     savefig(fig,figname)
     fig
 end
end
GM.activate!(); plotsave_contact_friction(tops_dt[8:8])
CM.activate!(); plotsave_contact_friction(tops_e0,"contact_friction")

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
 r̄g = SVector{3}(0.0,0.0,0.0)
 r̄ps = SVector{3}.([
     [ -l/2, 0.0, 0.0],
     [  l/2, 0.0, 0.0]
 ])
 θ = π/4
 zoffset = 1e-6
 xoffset = 1e-6
 ri = [     0.0-xoffset,0.0,sin(θ)*l-zoffset]
 rj = [cos(θ)*l-xoffset,0.0,     0.0-zoffset]
 ro = (ri + rj)./2
 ṙo = zero(ro)
 ω = zero(ro)
 R = Matrix(RotY(θ))

 # b = l/2*ω[2]*sin(θ)-9.81
 # p⁺ = 1 + 3*cos(θ)^2-μ*cos(θ)*sin(θ)
 # @show b,p⁺

 prop = TR.RigidBodyProperty(1,movable,m,Ī,r̄g,r̄ps;constrained)

 lncs, _ = TR.NCF.NC3D2P(ri,rj,ro,R,ṙo,ω)
 state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω)

 p1 = Meshes.Point(r̄ps[1])
 p2 = Meshes.Point(r̄ps[2])
 s = Meshes.Segment(p1,p2)
 cyl_bar = Meshes.Cylinder(Meshes.length(s)/40,s)
 cylsurf_bar = Meshes.boundary(cyl_bar)
 cyl_bar_simple = Meshes.triangulate(cylsurf_bar)
 cyl_bar_mesh = cyl_bar_simple |> simple2mesh
 rb1 = TR.RigidBody(prop,state,cyl_bar_mesh)
 rbs = TypeSortedCollection((rb1,))
 numberedpoints = TR.number(rbs)
 matrix_sharing = zeros(Int,0,0)
 indexedcoords = TR.index(rbs,matrix_sharing)
 ss = Int[]
 tensiles = (cables = ss,)
 connections = TR.connect(rbs,zeros(Int,0,0))
 contacts = [TR.Contact(i,μ,e) for i = 1:2]
 cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
 tg = TR.TensegrityStructure(rbs,tensiles,cnt,contacts)
 bot = TR.TensegrityRobot(tg)
end

bar = make_bar(;)

plot_traj!(bar;)

function bar_contact_dynfuncs(bot)
 (;tg) = bot
 function F!(F,q,q̇,t)
     TR.clear_forces!(tg)
     TR.update_rigids!(tg,q,q̇)
     TR.update_cables_apply_forces!(tg)
     TR.apply_gravity!(tg)
     F .= TR.generate_forces!(tg)
 end
 function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
     ∂F∂q̌ .= 0
     ∂F∂q̌̇ .= 0
     TR.clear_forces!(tg)
     TR.update_rigids!(tg,q,q̇)
     TR.update_cables_apply_forces!(tg)
     TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
     TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
 end

 rbs = TR.get_bodies(tg)
 function prepare_contacts!(contacts, q)
     TR.update_rigids!(tg,q)
     rb = rbs[1]
     dirs = [1,3]
     for pid = 1:2
         gap = rb.state.rps[pid][dirs[pid]]
         TR.activate!(contacts[pid],gap)
     end
     active_contacts = filter(contacts) do c
         c.state.active
     end
     na = length(active_contacts)
     D = Matrix{eltype(q)}(undef,3na,length(q))
     inv_μ_vec = ones(eltype(q),3na)
     for (i,ac) in enumerate(active_contacts)
         (;id,state) = ac
         n = float.(1:3 .== dirs[id])
         state.frame = TR.SpatialFrame(n)
         (;n,t1,t2) = state.frame
         rbid = 1
         C = rbs[rbid].state.cache.Cps[id]
         CT = C*TR.build_T(tg,rbid)
         Dn = Matrix(transpose(n)*CT)
         Dt1 = Matrix(transpose(t1)*CT)
         Dt2 = Matrix(transpose(t2)*CT)
         D[3(i-1)+1,:] = Dn
         D[3(i-1)+2,:] = Dt1
         D[3(i-1)+3,:] = Dt2
         inv_μ_vec[3(i-1)+1] = 1/ac.μ
     end
     es = [ac.e for ac in active_contacts]
     gaps = [ac.state.gap for ac in active_contacts]
     H = Diagonal(inv_μ_vec)
     active_contacts, na, gaps, D, H, es
 end
 F!,Jac_F!,prepare_contacts!
end


tspan = (0.0,0.4)
h = 1e-3


prob = TR.SimProblem(bar,bar_contact_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=50,exception=false)

contacts_traj_voa = VectorOfArray(bar.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,1:end])
 check_Coulomb(i,c)
end

plot_traj!(bar;showinfo=false)

me = TR.mechanical_energy!(bar)
lines(me.E)

[c[1].state.active for c in bar.contacts_traj] |> lines
[c[1].state.persistent for c in bar.contacts_traj] |> lines

[c[2].state.active for c in bar.contacts_traj] |> lines
[c[2].state.persistent for c in bar.contacts_traj] |> lines
#-- painleve's bar

#-- meteor hammer
L = 1.0
radius = 0.01
ancs = TR.ANCF.ANC3DRURU(8.96e3;E=110e6,L=1.0,radius)
ancs
S = TR.ANCF.make_S(ancs)(0.5)
Sₓ = TR.ANCF.make_Sₓ(ancs)(0.5)
Sₓₓ = TR.ANCF.make_Sₓₓ(ancs)(0.5)
ri = rand(3)
riₓ = rand(3)
rj = rand(3)
rjₓ = rand(3)
e = [ri;riₓ;rj;rjₓ]
r = TR.ANCF.make_r(ancs,e)(0.5)
rₓ = TR.ANCF.make_rₓ(ancs,e)(0.5)
rₓₓ = TR.ANCF.make_rₓₓ(ancs,e)(0.5)
κ = TR.ANCF.make_κ(ancs,e)(0.5)
_κ = TR.ANCF.make_κ(ancs,e)
function ∂κ∂eᵀ_forwarddiff!(out,ancs,x,e)
    function κ(e)
        TR.ANCF.make_κ(ancs,e)(x)
    end
    ForwardDiff.gradient!(out,κ,e)
end
ne = length(e)
eT = eltype(e)
out = zeros(eT,ne)
∂κ∂eᵀ_forwarddiff!(out,ancs,0.5L,e)
TR.ANCF.make_∂κ∂eᵀ(ancs,e)(0.5) .- ∂κ∂eᵀ_forwarddiff!(out,ancs,0.5,e) |> norm 


@btime ∂κ∂eᵀ_forwarddiff!($out,$ancs,0.5,$e)

@btime TR.ANCF.make_S($ancs)(0.5)
@btime TR.ANCF.make_Sₓ($ancs)(0.5)
@btime TR.ANCF.make_Sₓₓ($ancs)(0.5)
@btime TR.ANCF.make_r($ancs,$e)(0.5)
@btime TR.ANCF.make_rₓ($ancs,$e)(0.5)
@btime TR.ANCF.make_rₓₓ($ancs,$e)(0.5)
@btime TR.ANCF.make_κ($ancs,$e)(0.5)
@btime TR.ANCF.make_∂κ∂eᵀ($ancs,$e)(0.5)
Q = TR.ANCF.make_Q(ancs)

@btime TR.ANCF.build_M($ancs)

M = TR.ANCF.build_M(ancs)
G = TR.ANCF.build_G(ancs) 

# ne = length(e)
# eT = eltype(e)
# out = zeros(eT,ne)
# ∂Q∂e_forwarddiff!(out,ancs,0.5L,e)
# @btime ForwardDiff.jacobian!($out,$Q,$e)
# ∂Q∂e |> issymmetric

cablemesh = sample(ancs,e,1000)

mesh(cablemesh,transparency=false)

function make_hammer(id,r̄ijkl,ro,R,ri,rj=nothing,rk=nothing,rl=nothing;
        μ,e,
        ω = zero(ro),
        movable = true,
        constrained = false,
        pres_idx = Int[],
        Φ_mask = collect(1:6),
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
    r̄g = SVector(-0.20,0,0)
    # r̄ps = deepcopy(r̄ijkl)
    r̄ps = vcat(
        [
            SVector(0.22*cos(i*π/2),0.0,0.22*sin(i*π/2)) + r̄g
            for i = -2:1
        ],
        [
            RotY( π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + r̄g
            for i = 1:6
        ],
        [
            RotY(-π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + r̄g
            for i = [1,2,4,5]
        ]
    )
    ās = [
        SVector(1.0,0,0) 
        for _ in eachindex(r̄ps)
    ]
    μs = fill(μ,length(r̄ps))
    es = fill(e,length(r̄ps))
    @show m,diag(Īg),r̄g,length(r̄ps),μs,es
    prop = TR.RigidBodyProperty(
                id,movable,m,Īg,
                r̄g,r̄ps,ās,
                μs,es;
                constrained=constrained
                )
    ṙo = zero(ro)
    if rj isa Nothing
        lncs,_ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
    elseif rk isa Nothing
        lncs,_ = TR.NCF.NC2P2V(ri,rj,ro,R,ṙo,ω)
    elseif rl isa Nothing
        lncs,_ = TR.NCF.NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
    else
        lncs,_ = TR.NCF.NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
    end
    # cf = TR.NCF.CoordinateFunctions(lncs,q0,pres_idx,free_idx)
    # @show typeof(lncs)
       # trimesh = Meshes.discretize(box,) |> simple2mesh
    meteormesh = load(joinpath(TR.assetpath("流星锤.STL"))) |> make_patch(;
        scale = 1/400,
        trans = [-0.20,0,0],
    )
    state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,pres_idx,Φ_mask)
    rb = TR.RigidBody(prop,state,meteormesh)
end

function cable_ancf(pres_idx, 𝐞, L = 1.0) 
    radius = 2.0e-3
    # ancs = TR.ANCF.ANC3DRURU(8.96e3;E=110e6,L,radius)
    # mat_cable = filter(
    #     row->row.name == "Nylon 66",
    #     material_properties
    # )[1]
    # ρ = ustrip(Unitful.kg/Unitful.m^3,mat_cable.density)
    # E = ustrip(Unitful.Pa,mat_cable.modulus_elas)
    ρ = 1.03e3
    E = 0.2e9
    @myshow ρ, E
    ancs = TR.ANCF.ANC3DRURU(ρ;E,L,radius)
    mass = TR.ANCF.build_mass(ancs)
    @show mass
    T = typeof(L)
    r̄g = SVector(L/2,T(0),T(0))
    r̄p1 = SVector(T(0),T(0),T(0))
    # r̄p2 = SVector(L,T(0),T(0))
    r̄ps = [
        r̄p1,
    ]
    prop = TR.FlexibleBodyProperty(
        1,
        :cable,
        mass,
        r̄g,
        # length(r̄ps),
        r̄ps
    )
    # cache = TR.get_CoordinatesCache(prop,ancs,𝐞)
    state = TR.FlexibleBodyState(prop,ancs,𝐞;pres_idx)
    fb = TR.FlexibleBody(prop,state)
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
    subfbs,subsm = TR.subdivide(fb1,nx)
    r̄ijkl = SVector{3,Float64}.(
        0.1 .*[
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
        ]
    )
    rb = make_hammer(
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
    fbs = TypeSortedCollection(vcat(subfbs,rb))
    # fbs = TypeSortedCollection([fb1,])
    numberedpoints = TR.number(fbs)
    # indexedcoords = TR.index(fbs,)
    sm = zeros(Int,size(subsm,1)+3,size(subsm,2)+1)
    sm[1:size(subsm,1),1:size(subsm,2)] = subsm
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)  ] = 7:9
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)+1] = 1:3
    # display(sm)
    indexedcoords = TR.index(fbs,sm)
    ss = Int[]
    tensiles = (cables = ss,)
    connected = TR.connect(fbs,)
    tensioned = @eponymtuple(connected,)
    cst1 = TR.FixedIndicesConstraint(1,[1,2,3],ri)
    jointed = TR.join((cst1,),indexedcoords)
    cnt = TR.Connectivity(numberedpoints,indexedcoords,tensioned,jointed)
    tg = TR.TensegrityStructure(fbs,tensiles,cnt,)
    bot = TR.TensegrityRobot(tg)
end

function flexcable_contact_dynfuncs(bot,ground_plane)
    contact_dynfuncs(bot;
        flatplane = ground_plane
    )
end

# don't bother
R = RotXY(deg2rad(0),deg2rad(-0))
inclined_plane = TR.Plane(R*[0,0,1.0],[-1.039230,-0.0,0.781])
p1 =  [-1.039230,0,0.780]
p2 = [-1.173952,-0.11,0.8646278]
p3 = [-0.901447,-0.112500,0.762216]
n = (p2-p1) × (p3 - p1 )

inclined_plane = TR.Plane(n,p3.+[0,0,0.001])


flexcable_DR = make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
TR.GDR!(flexcable_DR;β=2e-4,maxiters=2e5,verbose=false)
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
# TR.GDR!(flexcable_DR;β=1e-3,maxiters=1e5)
# flexcable = make_flexcable(;R = RotY(-π/2),e=0.0,μ=0.01,L=1.5,nx=5)

flexcable_DR.traj.q[end][end-8:end] .= flexcable.traj.q[end][end-8:end]
flexcable_DR.traj.q̇[end][end-8:end] .= flexcable.traj.q̇[end][end-8:end]

TR.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])


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

TR.solve!(
    TR.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
    TR.ZhongCCP();
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

plotsave_contactpoints(flexcable,)

rp2 = get_trajectory!(flexcable,6,10)
lines(rp2[2,:])

me = TR.mechanical_energy!(flexcable)
lines((me.E.-me.E[begin])./me.E[begin])

#dt
# dts = [1e-1,3e-2,1e-2,3e-3,1e-3,1e-4]
dts = [1e-2,5e-3,2e-3,0.99e-3,5e-4,2e-4,2e-5]
flexcables_dt = [
    begin
        flexcable_dt = deepcopy(flexcable)
        TR.solve!(
            TR.SimProblem(flexcable_dt,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
                TR.ZhongCCP();
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
    rp1 = get_trajectory!(bot,6,1)
    rp12 = get_trajectory!(bot,6,12)
    rp10 = get_trajectory!(bot,6,10)
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
                TR.goto_step!(bot,step)
                tgvis = deepcopy(bot.tg)
                (;r,g,b) = cg[i]
                viz!(ax,tgvis;
                    meshcolor=Makie.RGBAf(r,g,b,alphas[i])
                )
            end
            lines!(ax,rp1[begin:laststep],color=:red)
            lines!(ax,rp12[begin:laststep],color=:blue)
            lines!(ax,rp10[begin:laststep],color=:purple)
            # lines!(ax,rp2[:,506:end],color=:blue)
            handlemesh = load(joinpath(TR.assetpath(),"把柄.STL")) |> make_patch(;
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
TR.GDR!(flexcable_DR;β=2e-4,maxiters=2e5,verbose=true)
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

TR.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])

plane_normal = RotZY(deg2rad(-7.5),deg2rad(-90))*[0,0,1.0]
inclined_plane = TR.Plane(plane_normal,[0.05,0,0])

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
TR.solve!(
    TR.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
    TR.ZhongCCP();
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

me = TR.mechanical_energy!(flexcable)
lines(flexcable.traj.t,(me.E.-me.E[begin])./me.E[begin])
rp2 = get_trajectory!(flexcable,6,2)
vp2 = get_velocity!(flexcable,6,2)
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
                TR.goto_step!(bot,step)
                tgvis = deepcopy(bot.tg)
                (;r,g,b) = cg[i]
                viz!(ax,tgvis;meshcolor=Makie.RGBA(r,g,b,alphas[i]))
            end
            lines!(ax,rp2[:,begin:4001],color=:blue)
            handlemesh = load(joinpath(TR.assetpath(),"把柄.STL")) |> make_patch(;
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
        flatplane = TR.Plane([0,0,1.0],[0,0,0.0])
    )
end

tspan = (0.0,1.5)
h = 1e-3

prob = TR.SimProblem(unibot,uni_dynfuncs)

TR.solve!(prob,TR.Zhong06();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=true)

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
        TR.solve!(
            TR.SimProblem(
                uni(0.0; μ = 0.9, e, z0 = 0.2),
                uni_dynfuncs
            ),
            TR.ZhongCCP();
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
            vz = get_velocity!(bot,1,2)[3,:]

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

unibot_e5 = TR.solve!(
    TR.SimProblem(
        uni(0.0; μ = 0.01, e=0.5, z0 = 0.2),
        uni_dynfuncs
    ),
    TR.ZhongCCP();
    tspan,dt=0.5h,ftol=1e-14,maxiters=50,exception=false,verbose=true,
)
GM.activate!(); plot_traj!(unibot_e5)

me = TR.mechanical_energy!(unibot_e5)
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
            rp2 = get_trajectory!(bot,1,2)
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
    TR.solve!(
        TR.SimProblem(
            uni(0.0; μ = 0.01, e=0.0, z0 = 0.2-0.13468-1e-5, ωz = 50.0, mbar = 1.0),
            uni_dynfuncs
        ),
        TR.ZhongCCP();
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

me = TR.mechanical_energy!(unibot_z0)
me.E |> lines

GM.activate!(); plot_traj!(unibot_z0)
#-- uni bot end ---

#-------- SUPERBall ---------


plot_traj!(ballbot)

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
    ṙo = SVector(2.0,1.0,0),
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
        flatplane = TR.Plane([0,0,1.0],[0,0,0.0])
    )
end

# testing
tspan = (0.0,5.0)
h = 1e-2
prob = TR.SimProblem(ballbot,ball_dynfuncs)
@time TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=100,exception=false)
#-- end testing
GM.activate!(); plotsave_contactpoints(ballbot)

plot_traj!(ballbot;)

me = TR.mechanical_energy!(ballbot)
me.E |> lines

step_start = time2step(1.6,ballbot.traj.t)
step_stop = time2step(2.5,ballbot.traj.t)
r2p1 = get_trajectory!(ballbot,2,1)
r1p2 = get_trajectory!(ballbot,1,2)
r6p2 = get_trajectory!(ballbot,6,2)
lines(r2p1)
lines(r1p2)
lines(r6p2)


dts = [1e-2,3e-3,1e-3,3e-4,1e-4,1e-5]
superballs_dt = [
    begin
        ballbot_dt = deepcopy(ballbot)
        # TR.set_new_initial!(ballbot_dt,ballbot.traj.q[step_start],ballbot.traj.q̇[step_start])
        prob = TR.SimProblem(ballbot_dt,ball_dynfuncs)
        TR.solve!(prob,
            TR.ZhongCCP();
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
    r1p2 = get_trajectory!(ballbot,1,2)
    r6p2 = get_trajectory!(ballbot,6,2)
    r2p1 = get_trajectory!(ballbot,2,1)
    me = TR.mechanical_energy!(ballbot,)
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
                TR.goto_step!(ballbot,step)
                tgvis = deepcopy(ballbot.tg)
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
            rbid,pid = divrem(cps[botid]-1,2) .+ 1
            ṙpc_mid = get_mid_velocity!(bot,rbid,pid)
            t_mids = get_time_mids(bot)
            lines!(ax1,t_mids,ṙpc_mid[3,:], label="v́ₙ")
            lines!(ax1,t_mids,ṙpc_mid[1,:], label="v́ₜ" )
            xlims!(ax1,0,10.0)
            # ylims!(ax1,-6,ymax)
            Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
            Label(grids[botid,1][1,1,Top()], "Contact No.$(cps[botid])", halign = :center)

            ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e_{\mathrm{eff}}")
            vz = get_velocity!(bot,rbid,pid)[3,:]

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
    ṙo = SVector(7.0,2.0,-7.0),
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
prob = TR.SimProblem(ballbot,ball_dynfuncs)
@time TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=200,exception=false)

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
    v2p1 = get_velocity!(bot,2,1,step_start:step_stop)
    v1p2 = get_velocity!(bot,1,2,step_start:step_stop)
    v6p2 = get_velocity!(bot,6,2,step_start:step_stop)
    r2p1 = get_trajectory!(bot,2,1)
    me = TR.mechanical_energy!(bot,)
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
                TR.goto_step!(bot,step)
                tgvis = deepcopy(bot.tg)
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
quadbot.tg.connectivity.numbered

plot_traj!(quadbot;
    zlims = (0,1),
    showground = false
)

function quad_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = TR.Plane([0,0,1.0],[0,0,0.0])
    )
end

tspan = (0.0,1.91)
tspan = (0.0,1.0)
h = 1e-3

prob = TR.SimProblem(quadbot,quad_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true)

TR.solve!(prob,TR.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(quadbot;
    zlims = (0,1),
    showground = false
)

me = TR.mechanical_energy!(quadbot)

me.E |> lines
contacts_traj_voa = VectorOfArray(quadbot.contacts_traj)

#-- cube
function make_cube(ro = [0.0,0.0,0.5],
        R = one(RotMatrix{3}),
        ṙo = [2.0,0.0,0.0],
        ω = [0.0,0.0,5.0];
        μ = 0.9,
        e = 0.0
    )
    m = 1.0
    r̄g = zeros(3)
    inertia = SMatrix{3,3}(Matrix(1.0I,3,3))
    h = 5.0
    r̄ps = [h.*[x,y,z] for z = [-1,1] for y = [-1,1] for x = [-1,1]]
    ās = [SVector(1.0,0,0) for _ = 1:length(r̄ps)]
    μs = [μ for _ = 1:length(r̄ps)]
    es = [e for _ = 1:length(r̄ps)]
    pts = Point3.(r̄ps)
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
    prop = TR.RigidBodyProperty(
        1,true,m,inertia,
        r̄g,r̄ps,ās,
        μs,es,
        )
    ri = copy(ro)
    lncs,q,q̇ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
    state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω)
    rb1 = TR.RigidBody(prop,state,cube_mesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = TR.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = TR.index(rbs,matrix_sharing)
    ss = Int[]
    tensiles = (cables = ss,)
    hub = nothing
    #
    connected = TR.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected)
    cnt = TR.Connectivity(numberedpoints,indexedcoords,tensioned)
    contacts = [TR.Contact(i,μ,e) for i = 1:8]
    tg = TR.TensegrityStructure(rbs,tensiles,cnt,contacts)
    bot = TR.TensegrityRobot(tg,nothing)
end

cubebot = make_cube()
plot_traj!(cubebot)

function cube_contact_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = TR.Plane([0,0,1.0],[0,0,0.0])
    )
end

ros = [[0,0,z] for z = 10.0:-0.1:0.0]
qs = [
    begin
        cube = make_cube(ro,RotXY(π/4,π/4))
        q = TR.get_q(cube.tg)
    end
    for ro in ros
]
# p = TR.generalized_α(1.0)
Ro = [
    0.368869  -0.824063   0.429949;
    0.769524   0.011314  -0.638518;
    0.521314   0.566385   0.63831
]
ṙo = [0.0,0.0,0.0]
ωo = [3.0,4.0,5.0]
# cube = make_cube()
# cube = make_cube(ros[1])
cube = make_cube(ros[2],Ro,ṙo,0.1.*ωo)
tspan = (0.0,1.91)
tspan = (0.0,5.0)
dt = 1e-3

prob = TR.SimProblem(cube,cube_contact_dynfuncs)
TR.solve!(prob,
    TR.ZhongCCP();tspan,dt,ftol=1e-10,maxiters=50,exception=false
)

plot_traj!(cube;
    rigidcolor=:white,
    xlims = [-10,10],
    ylims = [-10,10],
    zlims = [-1e-3,10],
    showinfo=false,
    showwire=true
)

ME = TR.mechanical_energy!(cube)
lines(cube.traj.t,ME.E)

ts,qs,vs,as,ṽ̇s = TR.NSSFC.nssfc(nq,nλ,q0,v0,p,h,contact_dynfuncs(cube),tspan;tol=1e-6,imax=10)

ts,cs,qs,vs,ps,λs,μs = TR.NSSFC.nhsolve(nq,nλ,nμ,q0,v0,contact_dynfuncs(cube),tspan;
            dt=h,ftol=1e-6,imax=50,exception=true)

ts,cs,qs,vs,ps,λs,μs = TR.NSSFC.ipsolve(nq,nλ,nμ,q0,v0,contact_dynfuncs(cube),tspan;
            dt=h,ftol=1e-10,imax=50,exception=false)

ts,qs,vs,as,ṽ̇s = TR.NSSFC.nssfc(n,b,q26026,v26026,p,h,contact_dynfuncs(cube),tspan;tol=1e-10,imax=100)

cuberef = deepcopy(cube)
prob = TR.SimProblem(cuberef,dynfuncs,tspan)
TR.solve!(prob,TR.Zhong06();dt=h,ftol=1e-14)

vis(cuberef,contact_dynamics)

cube = make_cube(ros[1],Ro,ṙo,ωo)
# cube = make_cube(ros[1])
# cube = make_cube()

TR.prepare_traj!(cube.traj;tspan,dt=h,restart=true)

cube.traj.t[2:end] = ts[2:end]
cube.traj.q[2:end] = qs[2:end]
cube.traj.q̇[2:end] = vs[2:end]

vis(cube,contact_dynamics)

vis(cube,contact_dynamics;do_record=true)

plot(cube.traj.ts,VectorOfArray(cube.traj.qs)[7,:],color="yellow")
plot!(cuberef.traj.ts,VectorOfArray(cuberef.traj.qs)[7,:],color="blue")

ME = TR.mechanical_energy!(cube;gravity=true,)

ME.E
scatter(ME.E)
kes,epes,gpes,es,es_err = analyse_energy(cube;gravity=true,elasticity=false)
es_off = OffsetArray(es,0:length(es)-1)
kes_off = OffsetArray(kes,0:length(es)-1)
gpes_off = OffsetArray(gpes,0:length(es)-1)
steprange = 0:length(es)-1
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
