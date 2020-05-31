using SPARK
using LinearAlgebra
using Parameters
using StaticArrays
using Revise
using Robot2D
const R2 = Robot2D
function man(n)
    nbody = 2n - 1
    n_revolute = 2(n-1)
    nbp = 2 + 2(n-1)
    a_lower = fill(0.1,n)
    a_upper = fill(0.08,n-1)
    m_lower = fill(0.02,n)
    m_upper = fill(0.013,n-1)
    Ic_lower = fill(11.3,n)
    Ic_upper = fill(5.4,n-1)
    lower_index = 1:2:nbody
    upper_index = 2:2:nbody
    a = zeros(nbody)
    m = zeros(nbody)
    Ia = zeros(nbody)
    for (i,j) in enumerate(lower_index)
        a[j] = a_lower[i]
        m[j] = m_lower[i]
        Ia[j] = Ic_lower[i] + m[j]*1/3*a[j]^2
    end
    for (i,k) in enumerate(upper_index)
        a[k] = a_upper[i]
        m[k] = m_lower[i]
        Ia[k] = Ic_upper[i] + m[k]*1/3*a[k]^2
    end
    A = zeros(2,nbp)
    for i in 2:nbp
        A[:,i] .= A[:,i-1] .+ [a[i-1],0.0]
    end

    function rigidbody(i,m,a,Ia,ri,rj)
        name = Symbol("rb"*string(i))
        type = :generic
        if i == 1
            movable = false
        else
            movable = true
        end
        @show ri[1],rj[1]
        CoM_x = a/2
        CoM_y = √3/6*a
        if isodd(i)
            CoM_y = -CoM_y
        end
        CoM = SVector{2}([CoM_x,CoM_y])

        nap = 3
        ap1 = SVector{2}([0.0,0.0])
        ap2 = SVector{2}([a,0.0])
        ap3_x = a/2
        ap3_y = √3/2*a
        if isodd(i)
            ap3_y = -ap3_y
        end
        ap3 = SVector{2}([ap3_x,ap3_y])
        anchorpoints = [ap1,ap2,ap3]

        prop = R2.RigidBody2DProperty(movable,name,type,
                    m,Ia,
                    CoM,
                    nap,
                    anchorpoints
                    )
        state = R2.RigidBody2DState(prop,ri,rj)
        rb = R2.RigidBody2D(prop,state)
    end
    rbs = [rigidbody(i,m[i],a[i],
            Ia[i],A[:,i],A[:,i+1]) for i = 1:nbody]

    nstring = 2(nbody-1)
    upstringlen = norm(rbs[2].state.p[3] - rbs[1].state.p[1])
    lostringlen = norm(rbs[2].state.p[2] - rbs[1].state.p[3])
    original_restlengths = zeros(nstring)
    restlengths = zeros(nstring)
    actuallengths = zeros(nstring)
    ks = zeros(nstring)
    for i = 1:nstring
        j = i % 4
        original_restlengths[i] =
                 restlengths[i] =
               actuallengths[i] =  ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],1.3e2,1.6e2)
    end
    vss = [R2.SString(ks[i],0.95original_restlengths[i],
        R2.SStringState(restlengths[i],actuallengths[i],0.0)) for i = 1:nstring]
    # @code_warntype   R2.DString(k[i],original_restlength[i],
    #         restlength[i],actuallength[i],zeros(MVector{4}))
    rbs,vss
end

n = 4
man(4)
#@code_warntype man(4)
function manstructure(n)
    rbs,vss = man(n)

    rb2p = [
        [i,i+1] for i = 1:length(rbs)
    ]
    body2q = [
        SVector(2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]) for pid in rb2p
    ]

    string2p_raw = [
        [zeros(Int,2),zeros(Int,2)] for i = 1:length(vss)
    ]
    string2p = Vector{Tuple{R2.ID,R2.ID}}()
    for i = 1:length(rbs)-1
        # string2p_raw[2i-1][1] .= [i  ,1]
        # string2p_raw[2i-1][2] .= [i+1,3]
        # string2p_raw[2i][1]   .= [i  ,3]
        # string2p_raw[2i][2]   .= [i+1,2]
        push!(string2p,(R2.ID(i,1),R2.ID(i+1,3)))
        push!(string2p,(R2.ID(i,3),R2.ID(i+1,2)))
    end
    cnt = R2.Connectivity(body2q,string2p)
    R2.Structure2D(rbs,vss,cnt)
end

manipulator = manstructure(4)
R2.reset_forces!(manipulator.rigidbodies)
R2.update_forces!(manipulator)
[ss.state.tension for ss in manipulator.strings]
manipulator.rigidbodies[2].state.Fanc[3]
function spark(n,st2d)
    rbs = st2d.rigidbodies
    vss = st2d.strings
    cnt = st2d.connectivity
    @unpack body2q = cnt
    @unpack nbody,nfixbody,nstring = st2d
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
    n_revolute = 2(n-1)
    nq = body2q[end][end]
    @assert nq == 4nbody - 2n_revolute
    #total_mass_matrix = zeros(4nbody,4nbody)
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= rbs[rbid].state.auxs.M
    end
    function M!(mm,q)
        mm .= mass_matrix
    end

    # System functions
    function ∂T∂q̇!(p,M,q,q̇)
        mul!(p,mass_matrix,q̇)
    end

    function F!(F,q,q̇,t)
        R2.reset_forces!(rbs)
        R2.q2rbstate!(st2d,q,q̇)
        R2.update_forces!(st2d)
        R2.genforces!(rbs)
        F .= 0.0
        for (rbid,rb) in enumerate(rbs)
            pindex = body2q[rbid]
            F[pindex] .+= rb.state.auxs.Q
        end
    end

    function Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        for (rbid,rb) in enumerate(rbs)
            pindex = body2q[rbid]
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
            pindex = body2q[rbid]
            ret[3nfixbody+rbid,pindex] .= rb.state.auxs.Φq(q[pindex])
        end
        ret[1,1:4] = [1.0,0.0,0.0,0.0]
        ret[2,1:4] = [0.0,1.0,0.0,0.0]
        ret[3,1:4] = [0.0,0.0,0.0,1.0]
        ret
    end

    A,Φ,∂T∂q̇!,F!,M!,nothing
end

A,Φ,∂T∂q̇!,F!,M!,jacs = spark(n,manipulator)

function initial(n,st2d)
    rbs = st2d.rigidbodies
    cnt = st2d.connectivity
    n_revolute = 2(n-1)
    @unpack nbody,nfixbody = st2d
    nq = cnt.body2q[end][end]
    q0 = zeros(nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt.body2q[rbid]
        q0[pindex] .= rb.state.coords.q
    end
    q̇0 = zero(q0)
    #q̇0[end] = 0.01
    q̇0[end-1:end] .= [0.0,0.001]
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
    λ0 = zeros(nconstraint)
    q0,q̇0,λ0
end
q0,q̇0,λ0 = initial(n,manipulator)

s = 1
tab = SPARKTableau(s)
tspan = (0.0,20.0)
cache = SPARKCache(size(A(q0))[2],size(A(q0))[1],0.01,tspan,s,(A,Φ,∂T∂q̇!,F!,M!,jacs))

state = SPARKsolve!(q0,q̇0,λ0,cache,tab)
