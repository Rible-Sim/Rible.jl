using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
#using DifferentialEquations
#using ForwardDiff
using DynamicPolynomials
using Printf
using CoordinateTransformations
# using Makie
# AbstractPlotting.__init__()
using Observables
using MakieLayout
using AbstractPlotting
using Makie
AbstractPlotting.__init__()
using Revise
using HomotopyContinuation
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

include("links_define.jl")
include("link_u_plotting.jl")

n = 2
h = 0.07
R = RotX(0.0)
linkn = links(n,h,R)
q0,q̇0,λ0 = TR.get_initial(linkn)


Y = Array(build_Y(linkn))

λ,Δu,a = TR.inverse(linkn,deepcopy(linkn),Y,gravity=true)
Δu
u0 = [s.original_restlen for s in linkn.strings]
rl = u0 + Δu
ℓ = [s.state.length for s in linkn.strings]
ℓ - rl
s = 1 ./ℓ



scene, layout = layoutscene( resolution = (1600, 1000))
ax1 = layout[1, 1] = LScene(scene, camera = cam3d!, raw = false)
bars,cables,_ = plotstructure!(ax1,linkn)
ls_g = labelslider!(scene, "Gravity:", 0:0.01:1; format = x -> "$(x)m/s²")
layout[2, 1] = ls_g.layout
ls_u = [labelslider!(scene, "u$i:", 0:(l/10):l; format = x -> begin @sprintf "%6fm" x end, sliderkw=Dict(:startvalue=>u0[i])) for (i,l) in enumerate(ℓ)]
slider_sublayout = GridLayout()
layout[1,2] = slider_sublayout
for i = 1:length(ℓ)
    slider_sublayout[i,1] = ls_u[i].layout
end
targets_ob = [ls_g.slider.value]
for ls in ls_u
    push!(targets_ob,ls.slider.value)
end

onany(targets_ob...) do targets...
    g1 = targets[1]
    u1 = [u for u in targets[2:end]]
    result =  TR.forward(linkn,(q0,s,λ,rl,1.0),(u1,g1))
    rsols = real_solutions(result)
    # @show length(rsols)
    rsol1 = rsols[1]
    qsol = rsol1[1:linkn.ncoords]
    update_scene!(linkn,bars,cables,qsol)
    # asol = rsol1[linkn.ncoords+1:linkn.ncoords+linkn.nstrings]
    # λsol = rsol1[linkn.ncoords+linkn.nstrings+1:end]
end
colsize!(layout, 1, Relative(2/3))
rowsize!(layout, 1, Relative(9/10))
scene
