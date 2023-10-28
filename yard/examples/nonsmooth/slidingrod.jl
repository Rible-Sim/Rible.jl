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
    movable = true
    constrained = false
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
    prop = RB.RigidBodyProperty(1,movable,m,Ī,mass_locus,loci;constrained)

    lncs, _ = RB.NCF.NC1P3V(ri,ro,Ro,ṙo,ωo)
    state = RB.RigidBodyState(prop,lncs,ro,Ro,ṙo,ωo)
    rb1 = RB.RigidBody(prop,state)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = RB.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
	tensiles = (cables = ss,)
	hub = nothing
	connections = RB.connect(rbs,zeros(Int,0,0))
	jointedmembers = RB.unjoin()
	contacts = [RB.ID(1,i) for i = 1:2]
	cnt = RB.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers,contacts)
	st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,nothing)
end

rod = make_rod()

plot_traj!(rod;textsize=60)

function make_bar()
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
    ri = [0.0,0.0,l/2/sin(θ)]
    rj = [l/2/sin(θ),0.0,0.0]
    ro = (ri + rj)./2
    ṙo = zero(ro)
    ωo = zero(ro)
    Ro = Matrix(RotY(θ))
    prop = RB.RigidBodyProperty(1,movable,m,Ī,mass_locus,loci;constrained)

    lncs, _ = RB.NCF.NC3D2P(ri,rj,ro,Ro,ṙo,ωo)
    state = RB.RigidBodyState(prop,lncs,ro,Ro,ṙo,ωo)
    rb1 = RB.RigidBody(prop,state)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = RB.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
	tensiles = (cables = ss,)
	hub = nothing
	connections = RB.connect(rbs,zeros(Int,0,0))
	jointedmembers = RB.unjoin()
	contacts = [RB.ID(1,i) for i = 1:2]
	cnt = RB.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers,contacts)
	st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,nothing)
end

bar = make_bar()

plot_traj!(bar;textsize=60)

function contact_dynamics(st)
	# bars = RB.get_rigidbars(st)
	rbs = RB.get_bodies(st)
	μ = 0.1
	e = 0.0
	contacts = [RB.Contact(i,μ,e) for i = 1:2]

    function prepare_contacts!(contacts, q)
		RB.update_rigids!(st,q)
		body = rbs[1]
		for pid = 1:2
			gap = body.state.loci_states[pid][3]
			contacts[pid].state.gap = gap
			contacts[pid].state.active = ifelse(gap<0,true,false)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
            state.frame = RB.spatial_frame(n)
			(;n,t1,t2) = state.frame
			bodyid = 1
            C = rbs[bodyid].state.cache.Cps[id]
			CT = C*RB.build_T(st,bodyid)
			Dn = Matrix(transpose(n)*CT)
			Dt1 = Matrix(transpose(t1)*CT)
			Dt2 = Matrix(transpose(t2)*CT)
			D[3(i-1)+1,:] = Dn
			D[3(i-1)+2,:] = Dt1
			D[3(i-1)+3,:] = Dt2
			inv_μ_vec[3(i-1)+1] = 1/ac.μ
		end
        restitution_coefficients = [ac.e for ac in active_contacts]
		gaps = [ac.state.gap for ac in active_contacts]
		H = Diagonal(inv_μ_vec)
        active_contacts, na, gaps, D, H, restitution_coefficients
    end

	function update_contacts!(active_contacts, v, Λ)
		vs = RB.split_by_lengths(v,3)
		Λs = RB.split_by_lengths(Λ,3)
		for (ac,va,Λa) in zip(active_contacts,vs,Λs)
			ac.state.v = SVector{3}(va)
			ac.state.Λ = SVector{3}(Λa)
		end
	end
	contacts, prepare_contacts!, update_contacts!
end

function dynfuncs(bot)
    (;st) = bot
    function F!(F,q,q̇,t)
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        RB.apply_gravity!(st)
        F .= RB.generate_forces!(st)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        RB.build_∂Q̌∂q̌!(∂F∂q̌,st)
        RB.build_∂Q̌∂q̌̇!(∂F∂q̌̇,st)
    end

    F!,Jac_F!,contact_dynamics(st)
end
#
# 𝐌,𝚽,𝚽𝐪,𝐅,Jacobians,contact_funcs = contact_dynfuncs(rod)
# 𝐠,get_indices,get_FCs,get_D = contact_funcs
# q0,v0 = RB.get_q(rod.st)
#
# 𝐠(q0)
# get_indices(q0)

tspan = (0.0,0.4)
h = 1e-3


prob = RB.SimProblem(rod,dynfuncs)

contacts_traj = RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-7,maxiters=50,exception=true)

contacts_traj = RB.solve!(prob,RB.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(rod;textsize=60)

tspan = (0.0,0.4)

prob = RB.SimProblem(bar,dynfuncs)

contacts_traj = RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-7,maxiters=50,exception=true)


contacts_traj_voa = VectorOfArray(contacts_traj)



plot_traj!(bar;textsize=60)

function find_θs(bot)
	(;st,traj) = bot
	rbs = RB.get_bodies(st)
    [
        begin
            RB.update_rigids!(st,q)
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
