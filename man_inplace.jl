using LinearAlgebra
using Parameters
using StaticArrays
using BenchmarkTools
using Revise
using SPARK
using Robot2D
const R2 = Robot2D
function man(n,θ=0.0)
    nbody = 2n - 1
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
        A[:,i] .= A[:,i-1] .+ a[i-1]*[cos(θ*(i-1)),sin(θ*(i-1))]
    end

    function rigidbody(i,m,a,Ia,ri,rj)
        name = Symbol("rb"*string(i))
        type = :generic
        if i == 1
            movable = false
        else
            movable = true
        end
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
    upstringlen = 0.95norm(rbs[2].state.p[3] - rbs[1].state.p[1])
    lostringlen = norm(rbs[2].state.p[2] - rbs[1].state.p[3])
    original_restlengths = zeros(nstring)
    restlengths = zeros(nstring)
    actuallengths = zeros(nstring)
    ks = zeros(nstring)
    cs = zeros(nstring)
    cs .= 1.0
    for i = 1:nstring
        j = i % 4
        original_restlengths[i] =
                 restlengths[i] =
               actuallengths[i] =  ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],1.0e2,1.0e2)
    end
    ss = [R2.SString(ks[i],cs[i],original_restlengths[i],
        R2.SStringState(restlengths[i],actuallengths[i],0.0)) for i = 1:nstring]
    # @code_warntype   R2.DString(k[i],original_restlength[i],
    #         restlength[i],actuallength[i],zeros(MVector{4}))
    rbs,ss

    acs = [
        ifelse(isodd(i),R2.Actuator(ss[2(i-1)+1:2i]),
                        R2.Actuator(ss[2i:-1:2(i-1)+1]))
        for i = 1:n
    ]

    rb2p = [
        [i,i+1] for i = 1:length(rbs)
    ]
    body2q = [
        SVector(2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]) for pid in rb2p
    ]

    string2p_raw = [
        [zeros(Int,2),zeros(Int,2)] for i = 1:length(ss)
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
    R2.Structure2D(rbs,ss,acs,cnt)
end
n = 4
manipulator = man(n,-π/12)
refman = man(n,-π/12)
function update_angles(st2d)
    rbs = st2d.rigidbodies
    angles = zeros(st2d.nbody)
    for (rbid,rb) in enumerate(rbs)
        if rbid > 1
            state0 = rbs[rbid-1].state
            V0 = state0.p[2]-state0.p[1]
            state1 = rbs[rbid].state
            V1 = state1.p[2]-state1.p[1]
            angles[rbid] = atan(V1[1]*V0[2]-V1[2]*V0[1],V0[1]*V1[1]+V0[2]*V1[2])
        end
    end
    angles
end
function const_massmatrix_for_spark(mass_matrix)
    function M!(mm,q)
        mm .= mass_matrix
    end
    function ∂T∂q̇!(p,M,q,q̇)
        mul!(p,mass_matrix,q̇)
    end
    M!, ∂T∂q̇!
end
function man_wend(n,st2d)
    rbs = st2d.rigidbodies
    vss = st2d.strings
    cnt = st2d.connectivity
    @unpack body2q = cnt
    @unpack nbody,nfixbody,nstring = st2d
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
    nq = body2q[end][end]
    q0,q̇0 = R2.get_q(st2d)

    M = R2.build_massmatrix(rbs,body2q)
    M!, ∂T∂q̇! = const_massmatrix_for_spark(M)
    ref_angles = update_angles(refman)
    pids = [R2.PIDController.PID(0.01,0.1,1.0,
            setpoint=ref_angles[i+1],dt=0.01) for i = 1:2]
    function actuate!(st2d,pids,t)
        acs = st2d.actuators
        angles = update_angles(st2d)
        for (id,ac) = enumerate(acs)
            pid = pids[id]
            input_angle = angles[id]
            output = R2.PIDController.update!(pid,input_angle,t)
            #@show id,output
            @unpack strings = ac
            for i in [1]
                str = strings[i]
                #str.state.restlength = str.original_restlength*(1+0.02Δ)
                str.state.restlength = str.original_restlength + output
                #@show i,str.state.restlength
            end
            for i in [2]
                str = strings[i]
                #str.state.restlength = str.original_restlength*(1-0.02Δ)
                str.state.restlength = str.original_restlength - output
                #@show i,str.state.restlength
            end
        end
    end

    function F!(F,q,q̇,t)
        R2.reset_forces!(rbs)
        R2.q2rbstate!(st2d,q,q̇)
        #actuate!(st2d,pids,t)
        R2.update_forces!(st2d)
        F .= 0.0
        R2.assemble_forces!(F,st2d)
    end

    Φ = R2.build_Φ(st2d)
    A = R2.build_A(st2d)

    #M,Φ,A,F!,nothing

    A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = man_wend(n,manipulator)
A,Φ,∂T∂q̇!,F!,M!,Jacs = man_wend(n,manipulator)
q0,q̇0,λ0 = R2.get_initial(manipulator,Φ)
@code_warntype Φ(q0)
@code_warntype A(q0)
F(q0,q̇0,0.0)
@code_warntype F(q0,q̇0,0.0)
dt = 0.01
prob = SPARK.DyProblem(man_wend(n,manipulator),q0,q̇0,λ0,(0.0,40.0))
state = SPARK.solve(prob,dt=dt,ftol=1e-13)

s = 1
tab = SPARKTableau(s)
tspan = (0.0,20.0)
cache = SPARKCache(size(A(q0))[2],size(A(q0))[1],dt,tspan,s,man_wend(n,manipulator))
state = SPARKsolve!(q0,q̇0,λ0,cache,tab,ftol=1e-14)

@btime SPARK.solve(prob,dt=dt)
@time SPARK.solve(prob,dt=dt)
SPARK.solve(prob,dt=dt)
@code_warntype SPARK.solve(prob,dt=dt)
