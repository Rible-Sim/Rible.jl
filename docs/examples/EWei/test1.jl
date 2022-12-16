using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
using HomotopyContinuation
using Printf
using GLMakie
using Unitful
using NLsolve
using Revise
using StructArrays
using EponymTuples
using LaTeXStrings
using StaticArrays
using TypeSortedCollections
using Match
using XLSX
import TensegrityRobots as TR
includet("define.jl")
includet("../vis.jl")

function rigidbody(i, rij, apss, rg)
    if i == 1
        movable = false
        constrained = true
        ci = collect(1:6)
        Φi = Int[]
    else
        movable = true
        if i == 2
            constrained = true
            ci = collect(1:2)
            Φi = [1]
        else
            constrained = false
            ci = Int[]
            if i == 4
                Φi = [1]
            else
                Φi = collect(1:3)
            end
        end
    end
    ri, rj = rij
    aps = apss
    u = rj - ri
    α = atan(u[2], u[1])
    r̄g = rg
    m = 1.0
    I = .5
    nrp = length(aps)
    ṙo = zeros(2); ω = 0.0
    prop = TR.RigidBodyProperty(
        i,
        movable,
        m,
        SMatrix{2, 2}([
            I 0
            0 I
        ]),
        SVector{2}(r̄g),
        [SVector{2}(aps[i]) for i in 1:nrp],
        constrained=constrained
    )
    ro = SVector{2}(ri)
    if isodd(i)
        lncs, _ = TR.NaturalCoordinates.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α, ṙo, ω)
    else
        lncs, _ = TR.NaturalCoordinates.NC2D2P(SVector{2}(ri), SVector{2}(rj), ro, α, ṙo, ω)
    end
    state = TR.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)
    rb = TR.RigidBody(prop, state)
end
rij = [[[0.0, 0.0], [1.0, 0.0]],
       [[0.0, 0.0], [0.0, 1.0]],
       [[0.0, 1.0], [1.0, 1.0]],
       [[0.0, 1.0], [0.0, 2.0]],
       [[0.0, 2.0], [1.0, 2.0]]]
apss = [[[1.0, 0.0], [-1.0, 0.0],[0.0, 0.0]],
        [[0.0, 0.0], [1.0, 0.0]],
        [[1.0, 0.0], [-1.0, 0.0],[0.0, 0.0]],
        [[0.0, 0.0], [1.0, 0.0]],
        [[1.0, 0.0], [-1.0, 0.0],[0.0, 0.0]]]
rg = [[0.0, 0.0], [0.5, 0.0],
      [0.0, 0.0], [0.5, 0.0],
      [0.0, 0.0]]

rbs = [rigidbody(i, rij[i], apss[i], rg[i]) for i in 1:5]
rigdibodies = TypeSortedCollection(rbs)
numberedpoints = TR.number(rigdibodies)
matrix_sharing_raw = Vector{Matrix{Int}}()
for i in 1:4
    s = zeros(2, 5)
    if !isodd(i)
        s[1:2, i] = 3:4
        s[1:2, i+1] = 1:2
    else
        s[1:2, i] = 1:2
        s[1:2, i+1] = 1:2
    end
    push!(matrix_sharing_raw, s)
end
matrix_sharing = reduce(vcat, matrix_sharing_raw)
indexedcoords = TR.index(rigdibodies, matrix_sharing)
nstrings = 4
ks = 1.0; restlens = 0.5
cables = [TR.Cable2D(i, restlens, ks, 0.0) for i in 1:nstrings]
tensiles = (cables=cables,)
acs = []
hub = (actuators=acs,)

s = zeros(4, 5)
s[1, 1] = 2; s[1, 3] = -2
s[2, 1] = 1; s[2, 3] = -1
s[3, 3] = 2; s[3, 5] = -2
s[4, 3] = 1; s[4, 5] = -1
connections = (cables=TR.connect(rigdibodies, s), )
cnt = TR.Connectivity(numberedpoints, indexedcoords, connections)
tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt)
bot = TR.TensegrityRobot(tg, hub)
plot_traj!(bot)

# bot.traj.q̇[begin][tail.tg.connectivity.indexed.mem2sysfull[end][1:4]] .= [0.1,0.0,0.1,0.0]
# function dynfuncs(bot)
#     (;tg) = bot
#     function F!(F,q,q̇,t)
#         TR.clear_forces!(tg)
#         TR.update_rigids!(tg,q,q̇)
#         TR.update_tensiles!(tg)
#         # TR.apply_gravity!(tg)
#         TR.generate_forces!(tg)
#         TR.get_force!(F,tg)
#         ## F .= 0
#     end
#     Jac_F! = nothing
#     @eponymtuple(F!)
# end

# prob = TR.SimProblem(bot,dynfuncs)
# TR.solve!(prob,TR.Zhong06();dt=0.01,tspan=(0.0,10.0),ftol=1e-13,verbose=true)

ω², δq̌ = TR.undamped_eigen(bot.tg)
@show [sqrt(ω²[i]) for i in 3:4]

_,λ0 = TR.check_static_equilibrium_output_multipliers(bot.tg,gravity=false)
bot_plot = deepcopy(bot)
q = TR.get_q(bot_plot.tg)
(;sysfree) = bot_plot.tg.connectivity.indexed
bot_plot.traj[1].q[sysfree] .+= δq̌[1]
plot_traj!(bot_plot)

display(Matrix(TR.build_M(bot.tg)))
display(TR.build_Ǩ(bot.tg,λ0))
