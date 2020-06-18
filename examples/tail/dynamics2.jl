using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
using LaTeXStrings
# using NLsolve
using Revise
using SPARK
using Robot2D
const R2 = Robot2D

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
        nap = length(aps)
        anchorpoints = [SVector{2}(aps[i]) for i = 1:nap]
        prop = R2.RigidBody2DProperty(i,movable,
                    m,inertia,
                    SVector{2}(CoM),
                    nap,anchorpoints
                    )
        state = R2.RigidBody2DState(prop,ri,rj)
        rb = R2.RigidBody2D(prop,state)
    end
    rbs = [rigidbody(i,CoM[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]

    nstring = 4n
    original_restlengths = zeros(nstring)
    restlengths = zeros(nstring)
    actuallengths = zeros(nstring)
    for i = 1:nstring
        j = i % 4
        original_restlengths[i] =
                 restlengths[i] =
               actuallengths[i] = ifelse(j∈[1,0],0.04,√3*0.02)
    end
    ks = fill(100.0,nstring)
    ss = [R2.SString(ks[i],original_restlengths[i],
        R2.SStringState(restlengths[i],actuallengths[i],0.0)) for i = 1:nstring]
    # @code_warntype   R2.DString(k[i],original_restlength[i],
    #         restlength[i],actuallength[i],zeros(MVector{4}))

    acs = [R2.Actuator(ss[4(i-1)+1:4i]) for i = 1:n]

    rb2p = [
        ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
    ]
    body2q = [[2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p]

    string2ap = Vector{Tuple{R2.ID,R2.ID}}()
    for i = 1:n
        push!(string2ap,(R2.ID(2i+1,1),R2.ID(2i-1,1)))
        push!(string2ap,(R2.ID(2i+1,1),R2.ID(2i  ,1)))
        push!(string2ap,(R2.ID(2i+1,2),R2.ID(2i  ,1)))
        push!(string2ap,(R2.ID(2i+1,2),R2.ID(2i-1,2)))
    end
    cnt = R2.Connectivity(body2q,string2ap)
    R2.Structure2D(rbs,ss,acs,cnt)
end
n = 4
tail = make_tail(n)

function dynfuncs(st2d)

    M = R2.build_massmatrix(st2d)
    Φ = R2.build_Φ(st2d)
    A = R2.build_A(st2d)

    Q̃=R2.build_Q̃(st2d)

    function F!(F,q,q̇,t)
        R2.reset_forces!(st2d)
        R2.q2rbstate!(st2d,q,q̇)
        R2.update_forces!(st2d)
        # F .= Q̃*R2.fvector(st2d)
        F .= 0.0
        R2.assemble_forces!(F,st2d)
    end

    M,Φ,A,F!,nothing
end


M,Φ,A,F!,Jacs = dynfuncs(tail)
q0,q̇0,λ0 = R2.get_initial(tail)
q̇0[end-3:end] .= [0.1,0.0,0.1,0.0]

dt = 0.01
prob = SPARK.DyProblem(dynfuncs(tail),q0,q̇0,λ0,(0.0,20.0))
sol = SPARK.solve(prob,dt=dt,ftol=1e-14,verbose=true)
