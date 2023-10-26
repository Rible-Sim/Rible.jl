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
using Rible; const TR = Rible
includet("tail_define.jl")
includet("plotting.jl")
n = 4
nver = n
nhor = n + 1
nb = nver + nhor
n_revolute = 2 + 2(n-1) #意义暂不明

ver_lengths = zeros(nver)
hor_lengths = zeros(nhor)
ver_lengths .= 4
hor_lengths .= 4
#定义每根刚体杆的长度均为4
O = zeros(2,nhor)

for i in 2:nhor
    O[2,i] = O[2,i-1] - ver_lengths[i-1]
end

P = zeros(2,nhor)
for i in 1:nhor
    P[:,i]=O[:,i] .+ [hor_lengths[i]/2,0]
end

#上述过程在定义节点的坐标，以矩阵的形式

m = fill(1.3e-2,nb)
L = zeros(nb)

ver_index = 2:2:nb
hor_index = 1:2:nb

for (i,j) in enumerate(ver_index)
    L[j] = ver_lengths[i]
end

for (i,j) in enumerate(hor_index)
    L[j] = hor_lengths[i]
end

CoM = [zeros(2) for i=1:nb]
for (i,j) in enumerate(ver_index)
    CoM[j] .= [ver_lengths[i]/2,0]
end

inertia = zeros(nb)

for (i,j) in enumerate(ver_index)
    inertia[j] = m[j]*L[j]^2/3
end

#惯性矩
for (i,j) in enumerate(hor_index)
    inertia[j] = m[j]*L[j]^2/12
end
#依然是惯性矩

ri = [zeros(2) for i = 1:nb]
rj = [zeros(2) for i = 1:nb]

for (i,j) in enumerate(ver_index)
    ri[j] .=O[:,i]
    rj[j] .=O[:,i+1]
end

for (i,j) in enumerate(hor_index)
    ri[j] .= O[:,i]
    rj[j] .= P[:,i]
end
#定义自然坐标的ri和rj，其实就是自然坐标的向量

p1 = [zeros(2) for i = 1:nb]
p2 = [zeros(2) for i = 1:nb]
for (i,j) in  enumerate(ver_index)
    p1[j] .= [0,0]
    p2[j] .= [ver_lengths[i],0]
end

for (i,j) in enumerate(hor_index)
    p1[j] .= [-hor_lengths[i]/2,0]
    p2[j] .= [hor_lengths[i]/2,0]
end
#p1和p2的意义暂不明

function rigidbody(i,CoM,m,inertia,ri,rj,aps)
    if i==1
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

    prop = RB.RigidBodyProperty(i,movable,     #所以TR.RigidBodyProperty这个包怎么用啊？
                    m,inertia,
                    SVector{2}(CoM),
                    aps;constrained=constrained
                    )
    state = RB.RigidBodyState(prop,ri,rj,pres_idx)
    rb = RB.RigidBody(prop,state)

end

rbs=[rigidbody(i,CoM[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]
ncables = 4n
original_restlens = zeros(ncables)
ks = zeros(ncables)

for i = 1:ncables
    j = i % 4        #除法后取余数
    original_restlens[i] = ifelse(j∈[1,0],0.04,√3*0.02)    #if...else？？？
    ks[i] = ifelse(j∈[1,0],50.0,100.0)
end

ss = [RB.SString2D(0.5original_restlens[i],ks[i],10.0) for i = 1:ncables]

    # @code_warntype   RB.DString(k[i],original_restlen[i],
    #         restlen[i],actuallength[i],zeros(MVector{4}))

acs = [RB.Actuator(SVector{4}(ss[4(i-1)+1:4i])) for i = 1:n]

rb2p = [
    ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
]
body2q_raw = [[2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p]
body2q = RB.filter_body2q(body2q_raw,rbs)
string2ap = Vector{Tuple{RB.ID,RB.ID}}()
for i = 1:n
    push!(string2ap,(RB.ID(2i+1,1),RB.ID(2i-1,1)))        #专用的push函数来向vector string2ap里面增加元素
    push!(string2ap,(RB.ID(2i+1,1),RB.ID(2i  ,1)))
    push!(string2ap,(RB.ID(2i+1,2),RB.ID(2i  ,1)))
    push!(string2ap,(RB.ID(2i+1,2),RB.ID(2i-1,2)))
end
cnt = RB.Connectivity(body2q,string2ap)
st = RB.Structure(rbs,ss,acs,cnt)
RB.update_cables_apply_forces!(st)
st
