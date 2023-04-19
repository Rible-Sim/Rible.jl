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
using HomotopyContinuation
using GLMakie
using Revise
import TensegrityRobots as TR

include("links_define.jl")
include("link_u_plotting.jl")

n = 2
h = 0.066
R = RotX(0.0)
linkn = links(n,h,R)
plotstructure(linkn)
q0,q̇0,λ0 = TR.get_initial(linkn)

Y = Array(build_Y(linkn))

λ,u,a= TR.inverse(linkn,deepcopy(linkn),Y,gravity=false)
u0 = [s.original_restlen for s in linkn.strings]
rl = u
ℓ = [s.state.length for s in linkn.strings]
s = 1 ./ℓ
TR.actuate!(linkn,a)
TR.check_static_equilibrium(linkn,q0,λ)

TR.forward(linkn,(q0,s,λ),(rl,0.0),(rl,1.0))

fig
fig = Figure(resolution = (1600, 1000))
ax1 = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false,resolution = (1200, 900)))
bars,cables,_ = plotstructure!(ax1,linkn)
ls_g = labelslider!(fig, "Gravity:", 0:0.01:1; format = x -> "$(x)×9.8 m/s²")
fig[2, 1] = ls_g.layout

nu = length(rl)
ulabels = ["u$i:" for i in 1:nu]
uranges = [0:(l/100):l for l in ℓ]
formats = Ref(x -> begin @sprintf "%6fm" x end)
startvalue = u0[1]
ls_u = labelslidergrid!(fig,ulabels,uranges;formats, sliderkw=Dict(:startvalue=>startvalue))
fig[:,2] = ls_u.layout

targets_ob = [ls_g.slider.value]
for sl in ls_u.sliders
    push!(targets_ob,sl.value)
end

onany(targets_ob...) do targets...
    g1 = targets[1]
    u1 = [u for u in targets[2:end]]
    result =  TR.forward(linkn,(q0,s,λ),(rl,0.0),(u1,g1))
    rsols = real_solutions(result)
    rsol1 = rsols[1]
    qsol = rsol1[1:linkn.ncoords]
    update_scene!(linkn,bars,cables,qsol)
end

colsize!(fig.layout, 1, Relative(3/4))
rowsize!(fig.layout, 1, Relative(10/10))
fig
