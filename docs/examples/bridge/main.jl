using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using TypeSortedCollections
using Rotations
using CoordinateTransformations
using Interpolations
using OffsetArrays
using RecursiveArrayTools
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using Unitful
using Printf
using CircularArrays
using StructArrays
using EponymTuples
using GeometryBasics, Meshing
# using FloatingTableView
using Match
using Revise
import TensegrityRobots as TR
cd("examples/bridge")
includet("bridge.jl")
includet("../vis.jl")
includet("../analysis.jl")
includet("../dyn.jl")


const CA = CircularArray
## deployable bridge
n = 5; m = 2
tgbridge1 = nbridge(n,m;θ=0.0,c=0.0,h=1.0)
plot_traj!(tgbridge1;showmesh=false,showground=false)
ℓ = TR.get_cables_len(tgbridge1.tg)
Y = TR.build_Y(tgbridge1)
ax = float.(Y)\ℓ
B,F̃ = TR.build_inverse_statics_for_actuation(tgbridge1,tgbridge1)
ap,an = TR.get_solution_set(B,F̃)
a1 = ap.-0.8*an[:,1]

tgbridge0 = nbridge(n,m;r=0.7714704084729095,c=25.0,h=0.45,l=0)

plot_traj!(tgbridge0;showmesh=false,showground=false)

ℓ = TR.get_cables_len(tgbridge0.tg)
Y = TR.build_Y(tgbridge0)
ax = float.(Y)\ℓ
B,F̃ = TR.build_inverse_statics_for_actuation(tgbridge0,tgbridge0)
ap,an = TR.get_solution_set(B,F̃)
a0 = ap.+0.4*an[:,1]
μ0,μ1 = Y*a0, Y*a1
TR.actuate!(tgbridge0,a0)
TR.update!(tgbridge0.tg)
TR.get_cables_tension(tgbridge0.tg)
_ = TR.check_static_equilibrium_output_multipliers(tgbridge0.tg)
l = 0
actor = make_pres_actor(repeat(μ0,l+1),repeat(μ1,l+1),1.0,4.0)
twobridges = nbridge(n,m;r=0.7714704084729095,c=25.0,h=0.45,l)

dpbridge = TR.TensegrityRobot(twobridges.tg,(actuators=[actor],))
TR.actuate!(dpbridge,[0.0])
TR.update!(dpbridge.tg)
_ = TR.check_static_equilibrium_output_multipliers(dpbridge.tg)

dt = 1e-3

TR.solve!(TR.SimProblem(dpbridge,
            (x)->dynfuncs(x;actuate=true,gravity=false)
            ),
            TR.Zhong06(),
            dt=dt,
            tspan=(0.0,5.0),
            ftol=1e-14)


plot_traj!(dpbridge;showmesh=false,showground=false,actuate=true)

ẍ,ẋ,x = random_vibration(;h=1e-5,tend=100.0)
lines(ẍ)
lines(ẋ)
lines(x)

dpbridge.tg.state.rigids[1].q̃

function pres_rand!(sysstate,tg)
    (;t,q̃,q̃̇,q̃̈) = sysstate
    (;rigidbodies,connectivity) = tg
    (;mem2syspres) = connectivity.indexed
    Random.seed!(06546)
    ẍ,ẋ,x = random_vibration(;tend=5.0,h=dt)
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
        if rbid in 1:3:13
            r0 = zeros(3)
            if rbid == 1
                r0 .= [
                    0.1,
                    0.0,
                    0.0
                ]
                q̃[mem2syspres[rbid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[mem2syspres[rbid]] .=       [amp*p*cos(p*t),0,0]
                q̃̈[mem2syspres[rbid]] .=      [-amp*p^2*sin(p*t),0,0]
            elseif rbid == 2
                r0 .= [
                    -0.05,
                     0.08660254037844387,
                    0.0
                ]
                q̃[mem2syspres[rbid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[mem2syspres[rbid]] .=       [amp*p*cos(p*t),0,0]
                q̃̈[mem2syspres[rbid]] .=      [-amp*p^2*sin(p*t),0,0]
            elseif rbid == 3
                r0 .= [
                    -0.05,
                    -0.08660254037844387,
                    0.0
                ]
                q̃[mem2syspres[rbid]] .= r0 .+ [amp*sin(p*t),0,0]
                q̃̇[mem2syspres[rbid]] .=       [amp*p*cos(p*t),0,0]
                q̃̈[mem2syspres[rbid]] .=      [-amp*p^2*sin(p*t),0,0]
            end
        end
    end
end

framerate = 30
1/dt/framerate
record(fig, "deploy_bridge.mp4", 1:33:length(dpbridge.traj);framerate) do this_step
    fig.content[2].value[] = this_step
end

ms = TR.mechanical_energy!(dpbridge;gravity=false)
lines(dpbridge.traj.t,ms.E)
