using LinearAlgebra, Statistics
using StaticArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using StructArrays
using BenchmarkTools
using TypeSortedCollections
using CoordinateTransformations
using Rotations
using ForwardDiff
using Interpolations
import GeometryBasics as GB
using OffsetArrays
using Printf
using Unitful
using Match
using EponymTuples
using Integrals
using Revise
import Meshes
using Meshing
import TensegrityRobots as TR
cd("examples/nonsmooth")
include("../analysis.jl"); includet("../analysis.jl")
include("../vis.jl"); includet("../vis.jl")
include("../dyn.jl"); includet("../dyn.jl")
figdir::String = raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ns"

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


∂Q∂e |> issymmetric

@btime ForwardDiff.jacobian!($out,$Q,$e)
ne = length(e)
eT = eltype(e)
out = zeros(eT,ne)
∂Q∂e_forwarddiff!(out,ancs,0.5L,e)

cablemesh = sample(ancs,e,1000)

mesh(cablemesh,transparency=false)

function make_cube(id,r̄ijkl,ro,R,ri,rj=nothing,rk=nothing,rl=nothing;
		movable = true,
		constrained = false,
		pres_idx = Int[],
		Φ_mask = collect(1:6),
	)
	# free_idx = collect(1:6)
	m = 0.3
	Īg = SMatrix{3,3}(
		Matrix(Diagonal([1.0E-03,1.0E-03,1.0E-03]))
	)
	ū = (r̄ijkl[2] - r̄ijkl[1])
	v̄ = (r̄ijkl[3] - r̄ijkl[1])
	w̄ = (r̄ijkl[4] - r̄ijkl[1])
	r̄g = r̄ijkl[1] .+ ū./2 .+ v̄./2 .+ w̄./2
	r̄ps = deepcopy(r̄ijkl)
	push!(r̄ps,r̄ijkl[1] .+ ū .+ v̄ .+ w̄)
	@show m,diag(Īg),r̄g
	prop = TR.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
				)
	ṙo = zero(ro)
	ω = zero(ro)
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
	
	box = Meshes.Box(
		Meshes.Point3(r̄ps[1]),
		Meshes.Point3(r̄ps[5])
	) |> Meshes.boundary
	trimesh = Meshes.discretize(box,) |> simple2mesh
	state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,pres_idx,Φ_mask)
	rb = TR.RigidBody(prop,state,trimesh)
end

function cable_ancf(pres_idx, 𝐞, L = 1.0) 
    radius = 0.01
    # ancs = TR.ANCF.ANC3DRURU(8.96e3;E=110e6,L,radius)
    ancs = TR.ANCF.ANC3DRURU(1.11e3;E=1.3e6,L,radius)
    mass = TR.ANCF.build_mass(ancs)
	@show mass
    T = typeof(L)
    r̄g = SVector(L/2,T(0),T(0))
    r̄p1 = SVector(T(0),T(0),T(0))
    r̄p2 = SVector(L,T(0),T(0))
    r̄ps = [
        r̄p1,r̄p2
    ]
    prop = TR.FlexibleBodyProperty(
        1,
        :cable,
        mass,
        r̄g,
        length(r̄ps),
        r̄ps
    )
    # cache = TR.get_CoordinatesCache(prop,ancs,𝐞)
    state = TR.FlexibleBodyState(prop,ancs,𝐞;pres_idx)
    fb = TR.FlexibleBody(prop,state)
end

function make_flexcable(
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.5*√3, 0.0, 0.5-1e-2],
        rjx = [ 0.0,-1.0, 0.0];
		doDR=false,
        μ = 0.5,
        e = 0.9,
        L = 1.0,
		nx = 2,
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
	R = RotX(0.0)
	rb = make_cube(nx+1,r̄ijkl,rj,R,rj;pres_idx=rb_pres_idx,constrained=ifelse(!isempty(rb_pres_idx),true,false))
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
	cst1 = TR.FixedIndicesConstraint([1,2,3],ri)
	jointed = TR.join((cst1,),indexedcoords)
    cnt = TR.Connectivity(numberedpoints,indexedcoords,tensioned,jointed)
    contacts = [TR.Contact(i,μ,e) for i = [5]]
    tg = TR.TensegrityStructure(fbs,tensiles,cnt,contacts)
    bot = TR.TensegrityRobot(tg)
end

flexcable_DR = make_flexcable(;L=1.5,nx=5,doDR=true)
flexcable = make_flexcable(;L=1.5,nx=5)

TR.GDR!(flexcable_DR;β=1e-3,maxiters=1e5)

TR.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])

plot_traj!(flexcable,)

flexcable[1].state.cache.pres_idx

flexcable.tg.connectivity.indexed.mem2sysincst

θ = -30 |> deg2rad
inclined_plane = TR.Plane([-tan(θ),0,1],[0,0,-0.1*√3])

function flexcable_contact_dynfuncs(bot,ground_plane)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        # TR.update_tensiles!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        # TR.update_tensiles!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        # TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end

	rbs = TR.get_rigidbodies(tg)
	rblast = rbs[end]
	nb = length(rbs)
	(;n) = ground_plane
    function prepare_contacts!(contacts,q)
		TR.update_rigids!(tg,q)
		for (pres_idxd,pid) in enumerate([5])
			# gap = rblast.state.rps[pid][3] 
			gap = TR.signed_distance(
				rblast.state.rps[pid],
				inclined_plane
			)
			TR.activate!(contacts[pres_idxd],gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		inv_μ_vec = ones(eltype(q),3na)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			state.frame = TR.SpatialFrame(n)
			inv_μ_vec[3(i-1)+1] = 1/ac.μ
		end
        es = [ac.e for ac in active_contacts]
		gaps = [ac.state.gap for ac in active_contacts]
		H = Diagonal(inv_μ_vec)
        active_contacts, gaps, H, es
    end
	
	function get_directions_and_positions(active_contacts,q)
		na = length(active_contacts)
		TR.update_rigids!(tg,q)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		ŕ = Vector{eltype(q)}(undef,3na)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			(;n,t1,t2) = state.frame
			dm = hcat(n,t1,t2) |> transpose
			ŕ[3(i-1)+1:3(i-1)+3] = dm*rblast.state.rps[id]
            # S = rblast.state.cache.Sps[id]
			# ST = S*TR.build_T(tg,nb)
			# D[3(i-1)+1:3(i-1)+3,:] = dm*ST
			C = rblast.state.cache.Cps[id]
			CT = C*TR.build_T(tg,nb)
			D[3(i-1)+1:3(i-1)+3,:] = dm*CT
		end
		D,ŕ
	end

	function get_∂Dq̇∂q(active_contacts,q,q̇)
		na = length(active_contacts)
		TR.update_rigids!(tg,q)
		T = eltype(q)
		nq = length(q)
		∂Dq̇∂q = zeros(T,3na,nq)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			(;n,t1,t2) = state.frame
			r̄p = rblast.prop.r̄ps[id]
			if rblast.state.cache.funcs.nmcs isa TR.QBF.QC
				∂Cẋ∂x = TR.QBF.make_∂Cẋ∂x(r̄p)
				TI = TR.build_T(tg,nb)
				∂Cq̇∂q = ∂Cẋ∂x(TI*q,TI*q̇)*TI
				∂Dq̇∂q[3(i-1)+1,:] = transpose(n)*∂Cq̇∂q
				∂Dq̇∂q[3(i-1)+2,:] = transpose(t1)*∂Cq̇∂q
				∂Dq̇∂q[3(i-1)+3,:] = transpose(t2)*∂Cq̇∂q
			end
		end
		∂Dq̇∂q
	end

	function get_∂DᵀΛ∂q(active_contacts,q,Λ)
		na = length(active_contacts)
		TR.update_rigids!(tg,q)
		T = eltype(q)
		nq = length(q)
		∂DᵀΛ∂q = zeros(T,nq,nq)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			(;n,t1,t2) = state.frame
			r̄p = rblast.prop.r̄ps[id]
			if rblast.state.cache.funcs.nmcs isa TR.QBF.QC
				∂Cᵀf∂x = TR.QBF.make_∂Cᵀf∂x(r̄p)
				TI = TR.build_T(tg,nb)
				Λi = @view Λ[3(i-1)+1:3(i-1)+3]
				fi = hcat(n,t1,t2)*Λi
				∂DᵀΛ∂q .+= transpose(TI)*∂Cᵀf∂x(TI*q,fi)*TI
			end
		end
		∂DᵀΛ∂q
	end

    # @eponymtuple(F!,Jac_F!)
	@eponymtuple(F!,Jac_F!,prepare_contacts!,get_directions_and_positions,get_∂Dq̇∂q,get_∂DᵀΛ∂q)
end

tspan = (0.0,2.0)
h = 1e-4
TR.solve!(
	TR.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
	TR.Zhong06();
	tspan,dt=h,ftol=1e-14,maxiters=50,exception=true,verbose=false)

TR.solve!(
	TR.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
	TR.ZhongCCP();
	tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false)

with_theme(theme_pub;
		Axis3=(
			azimuth = 7.925530633326988,
			elevation = 0.11269908169872402
		)
	) do
	plot_traj!(
	flexcable;
	AxisType=Axis3,
	# dorecord=true,
	showpoints=false,
	zlims=[-0.5,1.51],
	doslide=true,
	showinfo=false,
	ground=inclined_plane,
	# figname="cable.mp4",
	)
end

plot_traj!(flexcable)

ME = TR.mechanical_energy!(flexcable)
lines(ME.E)

contacts_traj_voa = VectorOfArray(flexcable.contacts_traj)
csa = StructArray(contacts_traj_voa[1,1:5])

plotsave_contact_persistent(flexcable)
    
#dt
# dts = [1e-1,3e-2,1e-2,3e-3,1e-3,1e-4]
dts = [1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,1e-6]
flexcables_dt = [
	begin
		flexcable = make_flexcable(;L=1.5,nx=5)
		TR.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])
		TR.solve!(TR.SimProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
				TR.Zhong06();
				tspan=(0.0,1e-1),dt,ftol=1e-14,maxiters=50,exception=false)
	end
	for dt in dts
]

rp2 = get_trajectory!(flexcables_dt[begin],1,2)
lines(rp2[2,:])
GM.activate!(); plotsave_error(flexcables_dt,dts,bid=6,pid=1,di=1)

flexcables_dt[end].traj.t
