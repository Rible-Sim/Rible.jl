using LinearAlgebra
using StaticArrays
# using Parameters
using TypeSortedCollections
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using CoordinateTransformations
using BenchmarkTools
using Rotations
import GeometryBasics as GB
using Unitful, Match, Printf
import Meshes
# using MeshViz
# using Meshing
using EponymTuples
using OffsetArrays
using Revise
import TensegrityRobots as TR
# GLMakie.enable_SSAO[] = true
cd("examples/nonsmooth")
include("../analysis.jl"); includet("../analysis.jl")
include("../vis.jl"); includet("../vis.jl")
include("../dyn.jl"); includet("../dyn.jl")

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
    prop = TR.RigidBodyProperty(1,true,m,inertia,r̄g,r̄ps)
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

function dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end

	rbs = TR.get_rigidbodies(tg)
	rb1 = rbs[1]
	function prepare_contacts!(contacts, q)
		TR.update_rigids!(tg,q)
		for (cid,pid) in enumerate(collect(1:8))
			gap = rb1.state.rps[pid][3]
			TR.activate!(contacts[cid],gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
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
			rbid = 1
            C = rbs[rbid].state.cache.Cps[id]
			CT = C*TR.build_T(tg,rbid)
			dm = hcat(n,t1,t2) |> transpose
			D[3(i-1)+1:3(i-1)+3,:] = dm*CT
			ŕ[3(i-1)+1:3(i-1)+3] = dm*rbs[rbid].state.rps[id]
		end
		D,ŕ
	end

	@eponymtuple(F!,Jac_F!,prepare_contacts!,get_directions_and_positions)
    
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


prob = TR.SimProblem(cube,dynfuncs)
TR.solve!(prob,TR.ZhongCCP();tspan,dt,ftol=1e-14,maxiters=50,exception=false)

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
