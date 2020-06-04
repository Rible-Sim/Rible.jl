using SPARK
using LinearAlgebra
using Parameters
using StaticArrays
using Revise
using Robot2D
const R2 = Robot2D

function tail(n)
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
    # aps = [p1[i],p2[i]]
    # nap = length(aps)
    # ap = [R2.AnchorPoint2D(p) for p in aps]
    #
    # prop = R2.RigidBody2DProperty(movable,name,type,
    #             m[i],
    #             SMatrix{2,2}(inertia[i]),
    #             SVector{2}(CoM[i]),
    #             SVector{nap}(ap)
    #             )
    # @code_warntype R2.RigidBody2DProperty(movable,name,type,
    #             m[i],
    #             SMatrix{2,2}(inertia[i]),
    #             SVector{2}(CoM[i]),
    #             SVector{nap}(ap)
    #             )
    # state = R2.RigidBody2DState(prop,ri[i],rj[i])
    # @code_warntype R2.RigidBody2DState(prop,ri[i],rj[i])
    # rb = R2.RigidBody2D(prop,state)
    # @code_warntype R2.RigidBody2D(prop,state)
    function rigidbody(i,CoM,m,inertia,ri,rj,aps)
        name = Symbol("rb"*string(i))
        type = :generic
        if i == 1
            movable = false
        else
            movable = true
        end
        nap = length(aps)
        anchorpoints = [SVector{2}(aps[i]) for i = 1:nap]
        prop = R2.RigidBody2DProperty(movable,name,type,
                    m,inertia,
                    SVector{2}(CoM),
                    nap,
                    anchorpoints
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

    string2bp = Vector{Tuple{R2.ID,R2.ID}}()
    for i = 1:n
        push!(string2bp,(R2.ID(2i+1,1),R2.ID(2i-1,1)))
        push!(string2bp,(R2.ID(2i+1,1),R2.ID(2i  ,1)))
        push!(string2bp,(R2.ID(2i+1,2),R2.ID(2i  ,1)))
        push!(string2bp,(R2.ID(2i+1,2),R2.ID(2i-1,2)))
    end

    R2.Structure2D(rbs,ss,acs,R2.Connectivity(body2q,string2bp))
end
n = 4
tailstruct = tail(n)

function tail_spark(n,st2d)
    rbs = st2d.rigidbodies
    ss = st2d.strings
    cnt = st2d.connectivity

    @unpack body2q = cnt
    @unpack nbody,nfixbody,nstring = st2d
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
    n_revolute = 2 + 2(n-1)
    nq = 4nbody - 2n_revolute
    @show nq
    #total_mass_matrix = zeros(4nbody,4nbody)
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt.body2q[rbid]
        mass_matrix[pindex,pindex] .+= rbs[rbid].state.auxs.M
    end
    function M!(mm,q)
        mm .= mass_matrix
    end

    # System functions
    function ∂T∂q̇!(p,M,q,q̇)
        mul!(p,mass_matrix,q̇)
    end

    function actuate!(ac,t)
        Δ = sin(t)
        @unpack strings = ac
        for i = 1:2
            str = strings[i]
            str.state.restlength = str.original_restlength*(1+0.2Δ)
        end
        for i = 3:4
            str = strings[i]
            str.state.restlength = str.original_restlength*(1-0.2Δ)
        end
    end
    function F!(F,q,q̇,t)
        # for i = 1:4
        #     actuate!(st2d.actuators[i],t)
        # end
        R2.reset_forces!(rbs)
        R2.q2rbstate!(st2d,q,q̇)
        R2.update_forces!(st2d)
        R2.genforces!(rbs)
        F .= 0.0
        for (rbid,rb) in enumerate(rbs)
            pindex = cnt.body2q[rbid]
            F[pindex] .+= rb.state.auxs.Q
        end
    end

    function Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        for (rbid,rb) in enumerate(rbs)
            pindex = cnt.body2q[rbid]
            ret[3nfixbody+rbid] = rb.state.auxs.Φ(q[pindex])
        end
        xi,yi,xj,yj = q[1:4]
        ret[1] = xi - 0.0
        ret[2] = yi - 0.0
        ret[3] = yj - 0.0
        ret
    end

    function A(q)
        ret = zeros(eltype(q),nconstraint,nq)
        for (rbid,rb) in enumerate(rbs)
            pindex = cnt.body2q[rbid]
            ret[3nfixbody+rbid,pindex] .= rb.state.auxs.Φq(q[pindex])
        end
        ret[1,1:4] = [1.0,0.0,0.0,0.0]
        ret[2,1:4] = [0.0,1.0,0.0,0.0]
        ret[3,1:4] = [0.0,0.0,0.0,1.0]
        ret
    end

    A,Φ,∂T∂q̇!,F!,M!,nothing
end

A,Φ,∂T∂q̇!,F!,M!,jacs = tail_spark(n,tailstruct)

function initial(n,st2d)
    rbs = st2d.rigidbodies
    cnt = st2d.connectivity
    n_revolute = 2 + 2(n-1)
    @unpack nbody,nfixbody = st2d
    nq = 4nbody - 2n_revolute
    q0 = zeros(nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt.body2q[rbid]
        q0[pindex] .= rb.state.coords.q
    end
    q̇0 = zero(q0)
    #q̇0[end] = 0.01
    q̇0[end-3:end] .= [0.1,0.0,0.1,0.0]
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
    λ0 = zeros(nconstraint)
    q0,q̇0,λ0
end
q0,q̇0,λ0 = initial(n,tailstruct)

s = 1
tab = SPARKTableau(s)
tspan = (0.0,20.0)
dt = 0.01
cache = SPARKCache(size(A(q0))[2],size(A(q0))[1],dt,tspan,s,(A,Φ,∂T∂q̇!,F!,M!,jacs))

state = SPARKsolve!(q0,q̇0,λ0,cache,tab)
