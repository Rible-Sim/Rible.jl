using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
# using NLsolve
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

function make_tail(n)
    nver = n
    nhor = n + 1
    nb = nver + nhor
    n_revolute = 2 + 2(n-1)
    ver_lengths = zeros(nver)
    hor_lengths = zeros(nhor)
    ver_lengths .= 0.04
    hor_lengths .= 0.04
    O = zeros(2,nhor)
    for i in 2:nhor
        O[2,i] = O[2,i-1] - ver_lengths[i-1]
    end
    P = zeros(2,nhor)
    for i in 1:nhor
        P[:,i] = O[:,i] .+ [hor_lengths[i]/2,0.0]
    end

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
    CoM = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        CoM[j] .= [ver_lengths[i]/2,0.0]
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
    for (i,j) in enumerate(ver_index)
        ri[j] .= O[:,i]
        rj[j] .= O[:,i+1]
    end
    for (i,j) in enumerate(hor_index)
        ri[j] .= O[:,i]
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
    function rigidbody(i,CoM,m,inertia,ri,rj,aps)
        if i == 1
            movable = false
        else
            movable = true
        end
        if i == 2
            constrained = true
            constrained_index = [1,2]
        else
            constrained = false
            constrained_index = Vector{Int}()
        end
        nap = length(aps)
        aps = [SVector{2}(aps[i]) for i = 1:nap]
        prop = TR.RigidBodyProperty(i,movable,
                    m,inertia,
                    SVector{2}(CoM),
                    aps;constrained=constrained
                    )
        state = TR.RigidBodyState(prop,ri,rj,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,CoM[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]

    nstring = 4n
    original_restlens = zeros(nstring)
    ks = zeros(nstring)
    for i = 1:nstring
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],0.04,√3*0.02)
        ks[i] = ifelse(j∈[1,0],50.0,100.0)
    end
    ss = [TR.SString2D(original_restlens[i],ks[i],1.0) for i = 1:nstring]
    # @code_warntype   TR.DString(k[i],original_restlen[i],
    #         restlen[i],actuallength[i],zeros(MVector{4}))

    acs = [TR.Actuator(SVector{4}(ss[4(i-1)+1:4i])) for i = 1:n]

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
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
n = 4
tail = make_tail(n)
@code_warntype make_tail(n)
q0,q̇0,λ0 = TR.get_initial(tail)
q̇0[end-3:end] .= 5*[0.1,0.0,0.1,0.0]

TR.get_nbodyconstraint(tail)
TR.get_nbodydof(tail)
TR.get_nbodycoords(tail)
TR.get_fixindex(tail)
TR.get_fixA(tail)

@code_warntype make_tail(n)

function dynfuncs(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    #Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,A,F!,Jacs = dynfuncs(tail,q0)
Φ(q0)
A(q0)
dt = 0.01
prob = TS.DyProblem(dynfuncs(tail,q0),q0,q̇0,λ0,(0.0,20.0))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
@code_warntype A(q0)
