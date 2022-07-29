using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
# using Parameters
using GeometryBasics
using OffsetArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using RecursiveArrayTools
using EponymTuples
using TypeSortedCollections
using Printf
using CoordinateTransformations
using Test
using Unitful
using Match
using FileIO
using Revise
import TensegrityRobots as TR
import Meshes
cd("examples\\quadruped")
includet("def.jl")
# includet("plotting.jl")
# includet("../analysis.jl")
includet("../vis.jl")


quadbot = quad(10.0)
plot_traj!(quadbot;
	zlims = (0,1)
)


function quad_dynfuncs(bot)
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
	bars = TR.get_rigidbars(tg)
	function prepare_contacts!(contacts, q)
		TR.update_rigids!(tg,q)
		foreach(bars) do rb
			(;id) = rb.prop
			gap = rb.state.rps[1][3]
			TR.activate!(contacts[id],gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
		for (i,ac) in enumerate(active_contacts)
			(;state) = ac
			state.frame = TR.SpatialFrame(n)
			(;n,t1,t2) = state.frame
			rbid = ac.id
			C = bars[rbid].state.cache.Cps[1]
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

tspan = (0.0,1.91)
tspan = (0.0,1.0)
h = 1e-3

prob = TR.SimProblem(quadbot,quad_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true)

TR.solve!(prob,TR.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(quadbot)

me = TR.mechanical_energy!(quadbot)

me.E |> lines
contacts_traj_voa = VectorOfArray(quadbot.contacts_traj)
