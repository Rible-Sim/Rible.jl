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
#cd("examples/tail")
includet("tail_define.jl")
includet("plotting.jl")
includet("../analysis.jl")
#张拉整体段为竖直4根刚体水平5根刚体，随动段为竖直4根刚体
n=4
nver = 8  #4+4
nhor = 5
nb = nver + nhor
n_revolute = 2 + 2(n-1)
ver_lengths = zeros(nver)
ver_lengths .= 2
hor_lengths = zeros(nhor)
hor_lengths .= 4
O = zeros(2,nhor)

for i in 2:nhor
    O[2,i] = O[2,i-1] - ver_lengths[i-1]
end

P = zeros(2,nhor)
for i in 1:nhor
    P[:,i] = O[:,i] .+ [hor_lengths[i]/2,0.0]
end

#矩阵P表示结构整体所有右端点的合集即：点（2,0）、（2，-2）、（2，-4）等

m = fill(1.0,nb)
L = zeros(nb)

#ver_index = 2:2:14
#hor_index = 1:2:15
ver_index = [2,4,6,8,10,11,12,13]
hor_index = [1,3,5,7,9]
for (i,j) in enumerate(ver_index)
    L[j] = ver_lengths[i]
end

for (i,j) in enumerate(hor_index)
    L[j] = hor_lengths[i]
end

#L为长度矩阵

mass_locus = [zeros(2) for i = 1:nb]
for (i,j) in enumerate(ver_index)
    mass_locus[j] .= [ver_lengths[i]/2,0.0]
end

inertia = zeros(nb)
for (i,j) in enumerate(ver_index)
    inertia[j] = m[j]*L[j]^2/3
end

for (i,j) in enumerate(hor_index)
    inertia[j] = m[j]*L[j]^2/12
end

ri = [zeros(2) for i = 1:nb]
rj = [zeros(2) for i = 1:nb]

O1=[[0,0] [0,-2] [0,-4] [0,-6] [0,-8] [0,-10] [0,-12] [0,-14] [0,-16]]
P1=[[2,0] [2,-2] [2,-4] [2,-6] [2,-8] [2,-10] [2,-12] [2,-14] [2,-16]]
for (i,j) in enumerate(ver_index)
    ri[j] .= O1[:,i]
    rj[j] .= O1[:,i+1]
end

for (i,j) in enumerate(hor_index)
    ri[j] .= O1[:,i]
    rj[j] .= P1[:,i]
end

#自然坐标的ri和rj
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
#推测p为自由度
######################################分割线####################################
function rigidbody(i,mass_locus,m,inertia,ri,rj,aps)
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
    prop = RB.RigidBodyProperty(i,movable,
                m,inertia,
                SVector{2}(mass_locus),
                aps;constrained=constrained
                )
    state = RB.RigidBodyState(prop,ri,rj,pres_idx)
    body = RB.RigidBody(prop,state)
end
rbs = [rigidbody(i,mass_locus[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]
#length(rbs)

ncables = 4n
original_restlens = zeros(ncables)
ks = zeros(ncables)
for i = 1:ncables
    j = i % 4
    original_restlens[i] = ifelse(j∈[1,0],2,√2*2)
    ks[i] = ifelse(j∈[1,0],100,100.0)
end
ss = [RB.SString2D(i,0.5original_restlens[i],ks[i],0.0) for i = 1:ncables]

tensiles = (cables=ss,)
acs = [RB.ManualActuation(RB.SimpleActuator(4(i-1)+j,5original_restlens[4(i-1)+j])) for i = 1:n  for j = 1:4]
hub = (actuators=acs,)

rb2p = [
    ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
]
bodyid2q_raw = [[2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p]
bodyid2q = RB.filter_bodyid2q(bodyid2q_raw,rbs)
string2ap = Vector{Tuple{RB.ID,RB.ID}}()
for i = 1:n
    push!(string2ap,(RB.ID(2i+1,1),RB.ID(2i-1,1)))
    push!(string2ap,(RB.ID(2i+1,1),RB.ID(2i  ,1)))
    push!(string2ap,(RB.ID(2i+1,2),RB.ID(2i  ,1)))
    push!(string2ap,(RB.ID(2i+1,2),RB.ID(2i-1,2)))
end
cnt = RB.Connectivity(bodyid2q,string2ap)
st = RB.Structure(rbs,tensiles,cnt)
RB.update_cables_apply_forces!(st)
function jac_singularity_check(st)
    q,_ = RB.get_coords(st)
    q[22] = -10
    A = RB.build_A(st)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    for (bodyid,rb) in enumerate(st.rigidbodies)
        if body.prop.movable && body.prop.constrained
            q_rb = body.state.coords.q
            Aq_rb = vcat(body.state.cache.cfuncs.Φq(q_rb),
                         body.state.cache.funcs.Φq(q_rb))
            rb_rank = rank(Aq_rb)
            intrinsic_Aq = body.state.cache.funcs.Φq(q_rb)
            # @show bodyid,lucompletepiv!(copy(intrinsic_Aq))
            # col_index = GECP(intrinsic_Aq)
            # @show bodyid,col_index
            # @show rank(intrinsic_Aq[:,col_index[1:6]])
            if rb_rank < minimum(size(Aq_rb))
                @warn "The $(bodyid)th rigid body's Jacobian is singular: rank(A(q))=$(rb_rank)<$(minimum(size(Aq_rb)))"
            end
        end
    end
end
jac_singularity_check(st)   #运行到这一步报warningSystem's Jacobian is singular: rank(A(q))=15<16
                                # @ Rible D:\Rible.jl\src\tensegrity.jl:475
tr = RB.Robot(st,hub)

Figur=plotstructure(tr)
save("figure.png", Figur)
q,_=RB.get_coords(st)
#q,_ = RB.get_coords(st)
#q = zeros(T,ncoords)
#ncoords = st.ncoords
#A = RB.build_A(st)
#rbs = st.rigidbodies
#@unpack bodyid2q = st.connectivity
#ncoords = st.ncoords
#T = RB.get_numbertype(st)
#q = zeros(T,ncoords)
q =Vec{28,Float64}(0,0,2,0,0,-2,2,-2,0,-4,2,-4,0,-6,2,-6,0,-8,2,-8,0,-10,0,-12,0,-14,0,-16)
A = RB.build_A(st)
Aq = A(q)
sys_rank = rank(Aq)#求矩阵的秩
minimum(size(Aq))
###############################################分割线##################################
rbs = st.rigidbodies
@unpack bodyid2q = st.connectivity
ncoords = st.ncoords
T = RB.get_numbertype(st)
q = zeros(T,ncoords)


@code_lowered for bodyid in st.mvbodyindex
    pindex = bodyid2q[bodyid]
    q[pindex] .= rbs[bodyid].state.coords.q
end

st.mvbodyindex
for bodyid in st.mvbodyindex
    pindex = bodyid2q[bodyid]
    println(pindex)
    q[pindex] .= rbs[bodyid].state.coords.q
    println(q)
    println("next step")
end

for bodyid in st.mvbodyindex
    pindex = bodyid2q[bodyid]
    q[pindex] .= rbs[bodyid].state.coords.q
    q[22] = -10
end
q


rbs[1].state.coords.q
rbs[2].state.coords.q
rbs[3].state.coords.q
rbs[4].state.coords.q
rbs[5].state.coords.q
rbs[6].state.coords.q
rbs[7].state.coords.q
rbs[8].state.coords.q
rbs[9].state.coords.q
rbs[10].state.coords.q
rbs[11].state.coords.q
rbs[12].state.coords.q
rbs[13].state.coords.q



st.mvbodyindex
