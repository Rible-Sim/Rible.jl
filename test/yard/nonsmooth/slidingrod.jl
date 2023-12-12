using LinearAlgebra
using StaticArrays
using RecursiveArrayTools
using BenchmarkTools
using GeometryBasics
using GLMakie
import Meshes
import MeshViz
using RecursiveArrayTools
using EponymTuples
using TypeSortedCollections
using Rotations
using Revise
import Rible as RB
cd("examples/nonsmooth")
# includet("plotting.jl")
includet("../analysis.jl")
includet("../vis.jl")


function make_rod()
    contactable = true
    visible = true
    m = 0.402
    mass_locus = @SVector zeros(3)
    ṙo = zero(mass_locus)
    ωo = zero(mass_locus)
    l = 1.0
    rad = 0.002
    Ixx = 1/2*m*rad^2
    Iyy = Izz = 1/12*m*(3rad^2+l^2)
    Ī = SMatrix{3,3}(Diagonal([Ixx,Iyy,Izz]))
    r = 0.0
    loci = SVector{3}.([
        [-l/2, 0.0, 0.0],
        [ l/2, 0.0, 0.0],
        # [-l/2, 0.0,  -r],
        # [ l/2, 0.0,  -r]
    ])
    θ = π/4
     leftcontact_point = [0.0,0.0,l/2/sin(θ)]
    rightcontact_point = [l/2/sin(θ),0.0,0.0]
    offset = r*[sin(θ),0.0,cos(θ)]
    r1 =  leftcontact_point + offset
    r2 = rightcontact_point + offset
    ro =(r1 + r2)/2
    ri = ro
    # @show ri,rj,ro
    Ro = Matrix(RotY(θ))
    prop = RB.RigidBodyProperty(1,contactable,m,Ī,mass_locus,loci;visible)

    nmcs = RB.NCF.NC1P3V(ri,ro,Ro)
    state = RB.RigidBodyState(prop,nmcs,ro,Ro,ṙo,ωo)
    rb1 = RB.RigidBody(prop,state)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = RB.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
	apparatuses = (cables = ss,)
	hub = nothing
	connections = RB.connect(rbs,zeros(Int,0,0))
	jointedmembers = RB.unjoin()
	contacts = [RB.ID(1,i) for i = 1:2]
	cnt = RB.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers,contacts)
	st = RB.Structure(rbs,apparatuses,cnt)
    bot = RB.Robot(st,nothing)
end

rod = make_rod()

plot_traj!(rod;textsize=60)

function make_bar()
    contactable = true
    visible = true
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
    ri = [0.0,0.0,l/2/sin(θ)]
    rj = [l/2/sin(θ),0.0,0.0]
    ro = (ri + rj)./2
    ṙo = zero(ro)
    ωo = zero(ro)
    Ro = Matrix(RotY(θ))
    prop = RB.RigidBodyProperty(1,contactable,m,Ī,mass_locus,loci;visible)

    nmcs = RB.NCF.NC3D2P(ri,rj,ro,Ro)
    state = RB.RigidBodyState(prop,nmcs,ro,Ro,ṙo,ωo)
    rb1 = RB.RigidBody(prop,state)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = RB.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
	apparatuses = (cables = ss,)
	hub = nothing
	connections = RB.connect(rbs,zeros(Int,0,0))
	jointedmembers = RB.unjoin()
	contacts = [RB.ID(1,i) for i = 1:2]
	cnt = RB.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers,contacts)
	st = RB.Structure(rbs,apparatuses,cnt)
    bot = RB.Robot(st,nothing)
end

bar = make_bar()

plot_traj!(bar;textsize=60)
#
# 𝐌,𝚽,𝚽𝐪,𝐅,Jacobians,contact_funcs = contact_dynfuncs(rod)
# 𝐠,get_idx,get_FCs,get_D = contact_funcs
# q0,v0 = RB.get_coords(rod.st)
#
# 𝐠(q0)
# get_idx(q0)

tspan = (0.0,0.4)
h = 1e-3


prob = RB.DynamicsProblem(rod,dynfuncs)

contacts_traj = RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-7,maxiters=50,exception=true)

contacts_traj = RB.solve!(prob,RB.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(rod;textsize=60)

tspan = (0.0,0.4)

prob = RB.DynamicsProblem(bar,dynfuncs)

contacts_traj = RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-7,maxiters=50,exception=true)


contacts_traj_voa = VectorOfArray(contacts_traj)



plot_traj!(bar;textsize=60)

function find_θs(bot)
	(;st,traj) = bot
	rbs = RB.get_bodies(st)
    [
        begin
            RB.update_bodies!(st,q)
			RB.update_orientations!(st)
            R = RotXYZ(RotXYZ(rbs[1].state.R))
            R.theta2
        end
        for q in traj.q
    ]
end
θs = find_θs(rod)
plot(rod.traj.t,θs)
xlims!(0,0.4)
ylims!(0,0.8)
