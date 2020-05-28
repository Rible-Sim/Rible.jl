using SPARK
using LinearAlgebra
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
    i = 1
    rb1 = rigidbody(i,CoM[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]])
    rbs = [rigidbody(i,CoM[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]

    original_restlength = fill(SVector(0.04,√3*0.02,√3*0.02,0.04),n)
    restlength = fill(MVector(0.04,√3*0.02,√3*0.02,0.04),n)
    actuallength = fill(MVector(0.04,√3*0.02,√3*0.02,0.04),n)
    k = fill(100.0,n)
    vss = [R2.DString(k[i],original_restlength[i],
        restlength[i],actuallength[i],zeros(MVector{4})) for i = 1:n]
    # @code_warntype   R2.DString(k[i],original_restlength[i],
    #         restlength[i],actuallength[i],zeros(MVector{4}))
    rbs,vss
end

n = 4
rbs,vss = tail(n)

@code_warntype tail(4)

function tailstructure(n)
    rbs,vss = tail(n)

    connectivity = [
        ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
    ]
    p_connectivity = [
        [2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in connectivity
    ]

    R2.Structure2D(rbs,vss,p_connectivity)
end
tailstruct = tailstructure(4)
@code_warntype tailstructure(4)

R2.reset_forces!(tailstruct.rigidbodies)
R2.update_forces!(tailstruct)
tailstruct.rigidbodies
function tail_spark(st2d)
    rbs = st2d.rigidbodies
    vss = st2d.strings
    cnt = st2d.connectivity

    ntotalbody = length(rbs)
    mvbodyindex = collect(2:ntotalbody)
    mvrigidbodies = rbs[mvbodyindex]
    nbody = length(mvrigidbodies)
    ninconstraint = ntotalbody
    nfix = ntotalbody - nbody
    nexconstraint = 3nfix
    nconstraint = ninconstraint + nexconstraint
    n_revolute = 2 + 2(n-1)
    nq = 4ntotalbody - 2n_revolute
    @show nq
    #total_mass_matrix = zeros(4nbody,4nbody)
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt[rbid]
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
        vs1 = vss[1]
        Δ = sin(t)
        rls = vs1.restlengths
        orls = vs1.original_restlengths
        rls[1] = orls[1]*(1+0.2Δ)
        rls[2] = orls[2]*(1+0.2Δ)
        rls[3] = orls[3]*(1-0.2Δ)
        rls[4] = orls[4]*(1-0.2Δ)
        R2.reset_forces!(rbs)
        R2.q2rbstate!(st2d,q,q̇)
        R2.update_forces!(st2d)
        R2.genforces!(rbs)
        F .= 0.0
        for (rbid,rb) in enumerate(rbs)
            pindex = cnt[rbid]
            F[pindex] .+= rb.state.auxs.Q
        end
    end

    function Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        for (rbid,rb) in enumerate(rbs)
            pindex = cnt[rbid]
            ret[3nfix+rbid] = rb.state.auxs.Φ(q[pindex])
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
            pindex = cnt[rbid]
            ret[3nfix+rbid,pindex] .= rb.state.auxs.Φq(q[pindex])
        end
        ret[1,1:4] = [1.0,0.0,0.0,0.0]
        ret[2,1:4] = [0.0,1.0,0.0,0.0]
        ret[3,1:4] = [0.0,0.0,0.0,1.0]
        ret
    end

    A,Φ,∂T∂q̇!,F!,M!,nothing
end

A,Φ,∂T∂q̇!,F!,M!,jacs = tail_spark(tailstruct)

function initial(n,st2d)
    rbs = st2d.rigidbodies
    cnt = st2d.connectivity
    n_revolute = 2 + 2(n-1)
    ntotalbody = length(rbs)
    nq = 4ntotalbody - 2n_revolute
    q0 = zeros(nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt[rbid]
        q0[pindex] .= rb.state.coords.q
    end
    q̇0 = zero(q0)
    #q̇0[end] = 0.01
    #q̇0[end-3:end] .= [0.01,0.0,0.01,0.0]
    ninconstraint = ntotalbody
    nfix = 1
    nexconstraint = 3nfix
    nconstraint = ninconstraint + nexconstraint
    λ0 = zeros(nconstraint)
    q0,q̇0,λ0
end
q0,q̇0,λ0 = initial(n,tailstruct)
Φ(q0)
A(q0)
@code_warntype Φ(q0)
@code_warntype A(q0)
F0 = similar(q0)
F!(F0,q0,q̇0,0.0)
@code_warntype F!(F0,q0,q̇0,0.0)

mm = zeros(length(q0),length(q0))
M!(mm,q0)
@code_warntype M!(mm,q0)

p0 = similar(q0)
∂T∂q̇!(p0,mm,q0,q̇0)
@code_warntype ∂T∂q̇!(p0,mm,q0,q̇0)
size(A(q0))
s = 1
tab = SPARKTableau(s)
tspan = (0.0,20.0)
cache = SPARKCache(size(A(q0))[2],size(A(q0))[1],0.01,tspan,s,(A,Φ,∂T∂q̇!,F!,M!,jacs))

state = SPARKsolve!(q0,q̇0,λ0,cache,tab)


@code_warntype SPARKsolve!(q0,q̇0,λ0,cache,tab)
te = [energy(state.qs[i],state.q̇s[i],rbs,vss) for i = 1:length(state.ts)]





@unpack mass,inertia = rbs[1].prop

Ic = inertia[2,2]
T = 1/2*ω*Ic*ω

1/2*transpose(q̇0)*M*q̇0
