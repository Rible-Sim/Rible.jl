using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
#plot(rand(10))
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
# using NLsolve
using Revise
using TensegrityRobots; const TR = TensegrityRobots
#cd("examples/tail")
includet("tail_define.jl")
includet("plotting.jl")
includet("../analysis.jl")

n=4
nver = n
nhor = n + 1
nb = nver + nhor
n_revolute = 2 + 2(n-1)

ver_lengths = zeros(nver)
hor_lengths = zeros(nhor)
ver_lengths .= 2
hor_lengths .= 4
O = zeros(2,nhor)
for i in 2:nhor
    O[2,i] = O[2,i-1] - ver_lengths[i-1]
end


P = zeros(2,nhor)
for i in 1:nhor
    P[:,i] = O[:,i] .+ [hor_lengths[i]/2,0.0]
end

m = fill(1.0,nb)
L = zeros(nb)

ver_index = 2:2:nb
hor_index = 1:2:nb
#ver_index = [2,4,6,8,10,12,14]
for (i,j) in enumerate(ver_index)
    L[j] = ver_lengths[i]
end

for (i,j) in enumerate(hor_index)
    L[j] = hor_lengths[i]
end

r̄g = [zeros(2) for i = 1:nb]
for (i,j) in enumerate(ver_index)
    r̄g[j] .= [ver_lengths[i]/2,0.0]
end
#这个东西没懂

inertia = zeros(nb)
for (i,j) in enumerate(ver_index)
    inertia[j] = m[j]*L[j]^2/3
end
for (i,j) in enumerate(hor_index)
    inertia[j] = m[j]*L[j]^2/12
end
ri = [zeros(2) for i = 1:nb]
rj = [zeros(2) for i = 1:nb]


for (i,j) in enumerate(ver_index)
    ri[j] .= O[:,i]
    rj[j] .= O[:,i+1]
end

for (i,j) in enumerate(hor_index)
    ri[j] .= O[:,i] # .=0
    rj[j] .= P[:,i]
end

p1 = [zeros(2) for i = 1:nb]
p2 = [zeros(2) for i = 1:nb]
for (i,j) in enumerate(ver_index)
    p1[j] .= [0.0,0.0]
    p2[j] .= [ver_lengths[i],0.0]
end


for (i,j) in enumerate(hor_index)
    p1[j] .= [-hor_lengths[i]/2,0.0]
    p2[j] .= [ hor_lengths[i]/2,0.0]
end
#推测p为约束

function rigidbody(i,r̄g,m,inertia,ri,rj,aps)
    if i == 1
        movable = true
        constrained = true
        pres_idx = [1,2,4]
    else
        movable = true
        constrained = false
        pres_idx = Int[]
    end
    nap = length(aps)
    aps = [SVector{2}(aps[i]) for i = 1:nap]
    prop = TR.RigidBodyProperty(i,movable,
                m,inertia,
                SVector{2}(r̄g),
                aps;constrained=constrained
                )
    state = TR.RigidBodyState(prop,ri,rj,pres_idx)
    rb = TR.RigidBody(prop,state)
end
rbs = [rigidbody(i,r̄g[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]

ncables = 4n
original_restlens = zeros(ncables)
ks = zeros(ncables)
for i = 1:ncables
    j = i % 4
    original_restlens[i] = ifelse(j∈[1,0],2,√2*2)
    ks[i] = ifelse(j∈[1,0],100,100.0)
end
ss = [TR.SString2D(i,0.5original_restlens[i],ks[i],0.0) for i = 1:ncables]

tensiles = (cables=ss,)
acs = [TR.ManualActuation(TR.SimpleActuator(4(i-1)+j,5original_restlens[4(i-1)+j])) for i = 1:n  for j = 1:4]
hub = (actuators=acs,)

rb2p = [
    ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
]
body2q_raw = [[2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p]
body2q = TR.filter_body2q(body2q_raw,rbs)
string2ap = Vector{Tuple{TR.ID,TR.ID}}()
for i = 1:n
    push!(string2ap,(TR.ID(2i+1,1),TR.ID(2i-1,1)))
    push!(string2ap,(TR.ID(2i+1,1),TR.ID(2i  ,1)))
    push!(string2ap,(TR.ID(2i+1,2),TR.ID(2i  ,1)))
    push!(string2ap,(TR.ID(2i+1,2),TR.ID(2i-1,2)))
end
cnt = TR.Connectivity(body2q,string2ap)
tg = TR.TensegrityStructure(rbs,tensiles,cnt)
TR.update_cables_apply_forces!(tg)
TR.jac_singularity_check(tg)
tr = TR.TensegrityRobot(tg,hub)

Figur=plotstructure(tr)
save("figure.png", Figur)
q,_ = TR.get_q(tg)
