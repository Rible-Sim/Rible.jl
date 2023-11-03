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
import Rible as RB
includet("define.jl")
includet("../vis.jl")
nb = 24
r̄gC = [44.71/2, 0.0]
r̄gA = [50/2, 0.0]
r̄gB = [50/2, 0.0]
r̄g1 = [0.0, 0.0]
r̄g2 = [0.0, 0.0]
r̄g3 = [0.0, 0.0]
r̄g4⁺ = [50-7.788, 0.0]
r̄g4 = [35-3.835, 0.0]
r̄g5 = [35-1.841, 0.0]
r̄g6 = [25-0.079, 0.0]
r̄g7 = [20+0.688, 0.0]
r̄g8 = [20+0.766, 0.0]
r̄g9 = [20+0.865, 0.0]
r̄g10 = [20+1.074, 0.0]

mC = 1.; mA = 22.; mB = 26.; m1 = 60.; m2 = 59.
m3 = 54.; m4⁺ = 78.; m4 = 48.; m5 = 48.; m6 = 45.
m7 = 40.; m8 = 36.; m9 = 32.; m10 = 40.
iC = 1.; iA = 6.071; iB = 8.646; i1 = 60.660;
i2 = 54.204; i3 = 41.680; i4⁺ = 68.548; i4 = 39.220
i5 = 47.819; i6 = 53.351; i7 = 45.136; i8 = 36.216
i9 = 24.286; i10 = 13.179

# r̄lC = 44.71; r̄lA = 50; r̄lB = 50
r̄lC = [0.0, 0.0, 0.0, 0.0]; r̄lA = deepcopy(r̄gA); r̄lB = deepcopy(r̄gB)
r̄l1 = 40; r̄l2 = 40; r̄l3 = 40
r̄l4⁺ = [50, 32, 15]; r̄l4 = [35, 32, 15]
r̄l5 = [35, 28, 15]; r̄l6 = [25, 27, 15]
r̄l7 = [20, 21, 15]; r̄l8 = [20, 20, 15]
r̄l9 = [20, 15, 15]; r̄l10 = [20, 15, 15]

function Get_apsAB(r̄, mass_locus)
    ri = r̄
    rj = r̄ + 2* [mass_locus[2], mass_locus[1]]
    aps = [[0.0, 0.0], 2r̄g, [15.1, 0.0]]
    return ri, rj, aps
end

function Get_apsC()
    ri = [0.0, 0.0]
    rj = [0.0, 44.71]
    aps = [[44.71, 0.0], [0.0, -40.0], [0.0, 40.0], [15.1, 0.0]]
    return ri, rj, aps
end

function Get_aps123(r̄, r̄l)
    ri = r̄
    rj = r̄ + [r̄l, 0.0]
    aps = [[r̄l, 0.0], -[r̄l, 0.0], [0.0, 0.0]]
    return ri, rj, aps
end

function Get_apsOther(r̄, r̄ls)
    down, right, up = r̄ls
    ri = r̄ - [0.0, down]
    rj = r̄ + [0.0, up]
    aps = [[0.0, 0.0], [down, right], [down, -right], [up+down, 0.0], [down, 0.0]]
    return ri, rj, aps
end

function Get_aps(r̄, r̄l)
    return @match length(r̄l) begin
        2 => Get_apsAB(r̄, r̄l)
        1 => Get_aps123(r̄, r̄l)
        3 => Get_apsOther(r̄, r̄l)
        4 => Get_apsC()
    end
end

rds = [
    RigidData(r̄gC, mC, iC, r̄lC),
    RigidData(r̄g1, m1, i1, r̄l1),
    RigidData(r̄gB, mB, iB, r̄lB),
    RigidData(r̄g1, m1, i1, r̄l1),
    RigidData(r̄gA, mA, iA, r̄lA),
    RigidData(r̄g1, m1, i1, r̄l1),
    RigidData(r̄gB, mB, iB, r̄lB),
    RigidData(r̄g1, m1, i1, r̄l1),
    RigidData(r̄gA, mA, iA, r̄lA),
    RigidData(r̄g2, m2, i2, r̄l2),
    RigidData(r̄gB, mB, iB, r̄lB),
    RigidData(r̄g2, m2, i2, r̄l2),
    RigidData(r̄gA, mA, iA, r̄lA),
    RigidData(r̄g2, m2, i2, r̄l2),
    RigidData(r̄gB, mB, iB, r̄lB),
    RigidData(r̄g3, m3, i3, r̄l3),
    RigidData(r̄g4⁺, m4⁺, i4⁺, r̄l4⁺),
    RigidData(r̄g4, m4, i4, r̄l4),
    RigidData(r̄g5, m5, i5, r̄l5),
    RigidData(r̄g6, m6, i6, r̄l6),
    RigidData(r̄g7, m7, i7, r̄l7),
    RigidData(r̄g8, m8, i8, r̄l8),
    RigidData(r̄g9, m9, i9, r̄l9),
    RigidData(r̄g10, m10, i10, r̄l10)]

function Computer̄s(rds)
    n = length(rds)
    r̄s = [[0.0, 0.0] for _ in 1:n]
    for i in 1:n
        if length(rds[i].r̄l) == 1
            r̄s[i][2] += r̄s[i-1][2]
        elseif length(rds[i].r̄l) == 2
            len = rds[i].r̄l[1] * 2
            r̄s[i+1][2] += r̄s[i][2] + len
            if i != 1
                r̄s[i][2] += r̄s[i-1][2]
            end
        elseif length(rds[i].r̄l) == 3
            len = rds[i].r̄l[1]
            r̄s[i][2] += r̄s[i-1][2] + len
            if i != n
                r̄s[i+1][2] += rds[i].r̄l[3]
            end
        elseif length(rds[i].r̄l) == 4
            r̄s[i+1][2] += 44.71
        end
    end
    return r̄s
end

r̄s = Computer̄s(rds)
rigidIndex = vcat([1], [2i for i in 1:8], [i for i in 17:24])
barIndex = collect(3:2:13)
function rigidbody(i, rdsi, r̄si; rigidIndex=rigidIndex, barIndex=barIndex)
    if i == 1
        movable = false
        constrained = true
        ci = collect(1:6)
        # ci = collect(1:4)
        cstr_idx = Int[]
    else
        movable = true
        if i in 2:3
            constrained = true
            ci = collect(1:2)
            if i == 2
                cstr_idx = collect(1:3)
                # cstr_idx = [1]
            else
                # cstr_idx = Int[]
                cstr_idx = [1]
            end
        else
            constrained = false
            ci = Int[]
            if i in rigidIndex
                cstr_idx = collect(1:3)
                # cstr_idx = [1]
            else
                cstr_idx = [1]
                # cstr_idx = collect(1:3)
            end
        end
    end

    ri, rj, aps = Get_aps(r̄si, rdsi.r̄l)
    u = rj - ri
    α = atan(u[2], u[1])
    mass_locus = rdsi.mass_locus
    m = rdsi.m*1e-3
    I = rdsi.i
    nrp = length(aps)
    ṙo = zeros(2); ω = 0.0
    prop = RB.RigidBodyProperty(
        i,
        movable,
        m,
        SMatrix{2, 2}([
            I 0
            0 0
        ]),
        SVector{2}(mass_locus),
        [SVector{2}(aps[i]) for i in 1:nrp],
        constrained=constrained
    )
    ro = SVector{2}(ri)
    if i in rigidIndex
        nmcs = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α)
        # nmcs = RB.NCF.NC2D2P(SVector{2}(ri), SVector{2}(rj), ro, α)
    else
        nmcs = RB.NCF.NC2D2P(SVector{2}(ri), SVector{2}(rj), ro, α)
        # nmcs = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α)
    end
    state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)
    body = RB.RigidBody(prop, state)
end

rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]
# rigidbody(1, rds[1], r̄s[1])
rigdibodies = TypeSortedCollection(rbs)
numberedpoints = RB.number(rigdibodies)
matrix_sharing_raw = Vector{Matrix{Int}}()
for i in 1:nb-1
    s = zeros(2, nb)
    if i == 1
        s[1:2, 1] = 3:4
        s[1:2, 2] = 1:2
    elseif i in 2:2:16
        # if i in rigidIndex
        s[1:2, i] = 1:2
        s[1:2, i+1] = 1:2
    else
        s[1:2, i] = 3:4
        s[1:2, i+1] = 1:2
        # end
    end
    push!(matrix_sharing_raw, s)
end
matrix_sharing = reduce(vcat, matrix_sharing_raw)
# display(matrix_sharing)
indexedcoords = RB.index(rigdibodies, matrix_sharing)
nstrings = 8 * 4 + 8 * 2
ks = zeros(nstrings); restlens = zeros(nstrings)
for i in 1:4*8
    j = i % 4
    restlens[i] = ifelse(j in [1, 0], 25, 20)
    ks[i] = ifelse(j in [1, 0], 0.3, 0.5)
end
for i in 4*8+1:nstrings
    restlens[i] = 25
    ks[i] = 0.3
end

cables = [RB.Cable2D(i, restlens[i], ks[i], 0.0) for i in 1:nstrings]

ncsegs = 3
lens = [44.71, 50, 50]
pre = [10.0 for i in 1:ncsegs]
c_section = StructArray([RB.CableSegment(i, lens[i], .2, prestress=pre[i]) for i in 1:ncsegs])
cs1 = RB.ClusterCables(1, ncsegs-1, deepcopy(c_section); μ=0.02)
cs2 = RB.ClusterCables(2, ncsegs-1, deepcopy(c_section); μ=0.02)

tensiles = (cables=cables, clustercables=[cs1, cs2])
# tensiles = (cables=cables,)
acs = []
hub = (actuators = acs, )

matrix_cnt_raw = Vector{Matrix{Int}}()
s = zeros(4, nb)
s[1, 1] = 3; s[1, 2] = -2
s[2, 1] = 4; s[2, 2] = -2
s[3, 1] = 4; s[3, 2] = -1
s[4, 1] = 2; s[4, 2] = -1
push!(matrix_cnt_raw, s)
for i in 1:7
    s = zeros(4, nb)
    s[1, 2i] = 2; s[1, 2i+2] = -2
    s[2, 2i+1] = 3; s[2, 2i+2] = -2
    s[3, 2i+1] = 3; s[3, 2i+2] = -1
    s[4, 2i] = 1; s[4, 2i+2] = -1
    push!(matrix_cnt_raw, s)
end
s = zeros(2, nb)
s[1, 16] = 1; s[1, 17] = -3
s[2, 16] = 2; s[2, 17] = -2
push!(matrix_cnt_raw, s)
for i in 17:23
    s = zeros(2, nb)
    s[1, i] = 3; s[1, i+1] = -3
    s[2, i] = 2; s[2, i+1] = -2
    push!(matrix_cnt_raw, s)
end
matrix_cnt = reduce(vcat, matrix_cnt_raw)

matrix_cnt_raw2 = Vector{Matrix{Int}}()
s = zeros(2, nb)
s[1, 1] = 3; s[1, 2] = 2; s[1, 4] = 2; s[1, 6] = 2
s[2, 1] = 2; s[2, 2] = 1; s[2, 4] = 1; s[2, 6] = 1
# s[1, 1] = 3; s[1, 2] = 2; s[1, 4] = 2
# s[2, 1] = 2; s[2, 2] = 1; s[2, 4] = 1
# for i in 1:4
#     s[1, 2i] = 2
#     s[2, 2i] = 1
# end
matrix_cnt2 = s
# display(matrix_cnt)
connections = RB.connect(rigdibodies, matrix_cnt, matrix_cnt2)
# connections = (cables=RB.connect(rigdibodies, matrix_cnt),)

cnt = RB.Connectivity(numberedpoints, indexedcoords, connections)

st = RB.ClusterTensegrityStructure(rigdibodies, tensiles, cnt)
# st = RB.Structure(rigdibodies, tensiles, cnt)
bot = RB.Robot(st, hub)
# bot.traj.q̇[begin][bot.st.connectivity.indexed.bodyid2sys_full_coords[end][1:4]] .= [50,0.0,50,0.0]
# function dynfuncs(bot)
#     (;st) = bot
#     function F!(F,q,q̇,s,t)
#         RB.clear_forces!(st)
#         RB.update_bodies!(st,q,q̇)
#         RB.distribute_s̄!(st,s)
#         RB.update_tensiles!(st)
#         ## RB.apply_gravity!(st)
#         RB.assemble_force!(st)
#         RB.get_force!(F,st)
#         ## F .= 0
#     end
#     function apply_acu!(st, t; dt=1e-2)
#         function inner(t;dt=dt)
#             a = 10.0
#             if 0<t<1
#                 return -a*t*dt
#             elseif 1<=t<5
#                 return -a*dt
#             elseif 5<=t<6
#                 return a*t*dt - 6a*dt
#             else
#                 return 0
#             end
#         end
#         st.clustercables[1].segs[1].state.restlen += inner(t)
#     end
#     # apply_acu! = nothing
#     Jac_F! = true
#     @eponymtuple(F!, Jac_F!, apply_acu!)
# end
# prob = RB.SimProblem(bot,dynfuncs)
# RB.solve!(prob,RB.FBZhong06();dt=0.01,tspan=(0.0,5.0),ftol=1e-7,verbose=true)

plot_traj!(bot)
_,λ0 = RB.check_static_equilibrium_output_multipliers(bot.st,gravity=false)
ω², δq̌ = RB.undamped_eigen(bot.st)
# @show ω²
@show ω = [sqrt(ωi) for ωi in ω²[1:end]]

# bot_plot = deepcopy(bot)
# q = RB.get_coords(bot_plot.st)
# (;sys_free_coords_idx) = bot_plot.st.connectivity.indexed
# bot_plot.traj[1].q[sys_free_coords_idx] .+= 10 * δq̌[1]
# plot_traj!(bot_plot)
temp = Matrix(RB.build_M̌(bot.st))
@show size(temp),rank(temp)