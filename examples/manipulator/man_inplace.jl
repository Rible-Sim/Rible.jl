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

function man_ndof(ndof,θ=0.0)
    nbody = ndof + 1
    nbp = 2nbody - ndof
    n_lower = count(isodd,1:nbody)
    n_upper = count(iseven,1:nbody)
    a_lower = fill(0.1,n_lower)
    a_upper = fill(0.08,n_upper)
    m_lower = fill(0.02,n_lower)
    m_upper = fill(0.013,n_upper)
    Ic_lower = fill(11.3,n_lower)
    Ic_upper = fill(5.4,n_upper)
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
    A[:,2] .= A[:,1] .+ a[1]*[1.0,0.0]
    for i in 3:nbp
        A[:,i] .= A[:,i-1] .+ a[i-1]*[cos(θ*(i-2)),sin(θ*(i-2))]
    end

    function rigidbody(i,m,a,Ia,ri,rj)
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
        aps = [ap1,ap2,ap3]

        prop = R2.RigidBodyProperty(i,movable,m,Ia,
                    CoM,aps
                    )
        state = R2.RigidBodyState(prop,ri,rj)
        rb = R2.RigidBody(prop,state)
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
    cs = zeros(nstring)
    cs .= 10.0
    for i = 1:nstring
        j = i % 4
        original_restlengths[i] =
                 restlengths[i] =
               actuallengths[i] =  ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],1.0e2,1.0e2)
    end
    ss = [R2.SString(ks[i],cs[i],original_restlengths[i],
        R2.SStringState(restlengths[i],actuallengths[i],0.0)) for i = 1:nstring]
    acs = [
        ifelse(isodd(i),R2.Actuator(ss[2(i-1)+1:2i]),
                        R2.Actuator(ss[2i:-1:2(i-1)+1]))
        for i = 1:nbody-1
    ]

    rb2p = [
        [i,i+1] for i = 1:length(rbs)
    ]
    body2q = [
        SVector(2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]) for pid in rb2p
    ]

    string2ap_raw = [
        [zeros(Int,2),zeros(Int,2)] for i = 1:length(ss)
    ]
    string2ap = Vector{Tuple{R2.ID,R2.ID}}()
    for i = 1:length(rbs)-1
        push!(string2ap,(R2.ID(i,1),R2.ID(i+1,3)))
        push!(string2ap,(R2.ID(i,3),R2.ID(i+1,2)))
    end
    cnt = R2.Connectivity(body2q,string2ap)
    R2.TensegrityStructure(rbs,ss,acs,cnt)
end
# ------------------Create Tensegrity Struture --------------------------
ndof = 6
refman = man_ndof(ndof,-π/12) # reference
manipulator = man_ndof(ndof,0.0)
# ------------------Create Tensegrity Struture\\-------------------------

q0,q̇0,λ0 = R2.get_initial(manipulator) #backup
# ----------------------Inverse Kinematics ------------------------------
function inverse(st2d,refst2d)
    function ikfuncs(st2d)

        A = R2.build_A(st2d)

        Q̃=R2.build_Q̃(st2d)

        function F!(F,u)
            R2.reset_forces!(st2d)
            R2.actuate!(st2d,u)
            R2.update_forces!(st2d)
            F .= Q̃*R2.fvector(st2d)
        end

        A,F!
    end
    q0,q̇0,λ0 = R2.get_initial(refst2d)
    R2.q2rbstate!(st2d,q0,q̇0)
    nu = length(st2d.actuators)
    u0 = zeros(nu)
    ikprob = SPARK.IKProblem(ikfuncs(st2d),q0,u0,λ0)
    R2.iksolve(ikprob)
end
u,_ = inverse(manipulator,refman) # PID setpoints

R2.q2rbstate!(manipulator,q0,q̇0)
R2.actuate!(manipulator,zero(u)) # reverse to initial
# ----------------------Inverse Kinematics\\-----------------------------

dt = 0.01 # Same dt used for PID AND Dynamics solver
# ---------------------Create Controllers --------------------------------
# use angles to define errors
function get_angles(st2d)
    rbs = st2d.rigidbodies
    angles = zeros(st2d.nbody-1)
    for (rbid,rb) in enumerate(rbs)
        if rbid > 1
            state0 = rbs[rbid-1].state
            V0 = state0.p[2]-state0.p[1]
            state1 = rbs[rbid].state
            V1 = state1.p[2]-state1.p[1]
            angles[rbid-1] = atan(V1[1]*V0[2]-V1[2]*V0[1],V0[1]*V1[1]+V0[2]*V1[2])
        end
    end
    angles
end
refangles = get_angles(refman)
refman
pids = [R2.PIDController.PID(0.01,0.0,0.01,
            setpoint=ui,dt=dt) for ui in refangles]
# ---------------------Create Controllers\\-------------------------------

# --------------------Create Robot ---------------------------------------
rob = R2.TGRobot2D(manipulator,R2.ControlHub(pids))
# --------------------Create Robot\\--------------------------------------

# --------------------Define Control Action-------------------------------
function make_control!(get_feedback)
    function inner_control!(robot2d,t)
        @unpack st2d, hub = robot2d
        @unpack ctrls,trajs = hub
        inputs = get_feedback(st2d)
        for (id,pid,traj,actuator) in zip(eachindex(ctrls),ctrls,trajs,st2d.actuators)
            input = inputs[id]
            output = R2.PIDController.update!(pid,input)
            R2.actuate!(actuator,output)
            R2.record!(traj,pid)
        end
    end
end

control! = make_control!(get_angles)
# --------------------Define Control Action\\-----------------------------

# ----------------------------Dynamics-----------------------------------
function dynfuncs(st2d)

    M = R2.build_massmatrix(st2d)
    M!, ∂T∂q̇! = SPARK.const_massmatrix_functions(M)
    Φ = R2.build_Φ(st2d)
    A = R2.build_A(st2d)

    Q̃=R2.build_Q̃(st2d)

    function F!(F,q,q̇,t)
        R2.reset_forces!(st2d)
        R2.q2rbstate!(st2d,q,q̇)
        R2.update_forces!(st2d)
        F .= Q̃*R2.fvector(st2d)
        #F .= 0.0
        #R2.assemble_forces!(F,st2d)
        #@show isapprox(F,)
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(manipulator)
prob = SPARK.DyProblem(dynfuncs(manipulator),q0,q̇0,λ0,(0.0,40.0))
# R2.actuate!(manipulator,u)
# sol = SPARK.solve(prob,dt=dt,ftol=1e-13,verbose=true)

function make_affect!(robot2d,control!)
    function inner_affect!(intor)
        R2.q2rbstate!(robot2d.st2d,intor.qprev,intor.q̇prev)
        R2.update_forces!(robot2d.st2d)
        control!(robot2d,intor.t)
    end
end
cb = SPARK.DiscreteCallback((x)->true,make_affect!(rob,make_control!(get_angles)))
R2.PIDController.tune!(rob.hub.ctrls[1],1.4,0.004,13)
R2.PIDController.tune!(rob.hub.ctrls[2],1.3,0.004,12)
R2.PIDController.tune!(rob.hub.ctrls[3],1.2,0.004,11)
R2.PIDController.tune!(rob.hub.ctrls[4],1.1,0.004,10)
R2.PIDController.tune!(rob.hub.ctrls[5],1.0,0.005,8.5)
R2.PIDController.tune!(rob.hub.ctrls[6],0.9,0.005,7.0)
R2.reset!.(rob.st2d.actuators)
R2.reset!.(rob.hub.trajs)
sol = SPARK.solve(prob,dt=dt,ftol=1e-13,callback=cb)
pltfig.clear(); pltfig = controlplot(rob.hub.trajs)
pltfig = controlplot(rob.hub.trajs)
sol = SPARK.solve(prob,dt=dt,ftol=1e-13,callback=cb,verbose=true)
# ----------------------------Dynamics-----------------------------------

function controlplot(trajs)
    ntraj = length(trajs)
    fig,axs_raw = plt.subplots(ntraj,1,num="PID",figsize=(5,15))

    if typeof(axs_raw)<:Array
        axs = axs_raw
    else
        axs = [axs_raw]
    end
    for (id,ax) in enumerate(axs)
        @unpack ts,es,us = trajs[id]
        bx = ax.twinx()
        ep = ax.plot(ts,es,label=latexstring("\\epsilon_$id"), lw = 3)
        up = bx.plot(ts,us,label=latexstring("u_$id"), lw = 3, color=:orange)
        ps = [ep[1],up[1]]
        bx.set_ylabel(L"u")
        bx.set_ylim([-0.1,0.4])
        ax.set_ylim(es[1].*[-0.1,1.1])
        ax.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(0.2es[1]))
        ax.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(xmax=es[1]))
        ax.set_ylabel(L"\epsilon(\%)")
        ax.set_xlabel(L"t")
        ax.legend(ps, [p_.get_label() for p_ in ps])
        ax.grid("on")
    end
    fig.savefig("manpid.png",dpi=300,bbox_inches="tight")
    plt.close(fig)
end
# @code_warntype man_wend(n,manipulator)
# A,Φ,∂T∂q̇!,F!,M!,Jacs = man_wend(n,manipulator)
# @code_warntype Φ(q0)
# @code_warntype A(q0)
# Fholder = similar(q0)
# @code_warntype F!(Fholder,q0,q̇0,0.0)
# F!(Fholder,q0,q̇0,0.0)

sol = SPARK.solve(prob,dt=dt,ftol=1e-13,verbose=true)
sol = SPARK.solve(prob,dt=dt,ftol=1e-13)






R2.q2rbstate!(manipulator,sol.qs[end],sol.q̇s[end])
endangles = get_angles(manipulator)

sol = SPARK.solve(prob,dt=dt,ftol=1e-13,verbose=true,callback=cb)

@time SPARK.solve(prob,dt=dt,ftol=1e-13)
s = 1
tab = SPARKTableau(s)
tspan = (0.0,20.0)
cache = SPARKCache(size(A(q0))[2],size(A(q0))[1],dt,tspan,s,man_wend(n,manipulator))
state = SPARKsolve!(q0,q̇0,λ0,cache,tab,ftol=1e-14)

@btime SPARK.solve(prob,dt=dt)
@time SPARK.solve(prob,dt=dt)
SPARK.solve(prob,dt=dt)
@code_warntype SPARK.solve(prob,dt=dt)

function showactuator(st2d)
    for (acid,actuator) in enumerate(st2d.actuators)
        @unpack strings = actuator
        str1 = strings[1]
        u1 = str1.state.restlength - str1.original_restlength
        str2 = strings[2]
        u2 = str2.state.restlength - str2.original_restlength
        @show acid,u1,u2
    end
end
showactuator(manipulator)
R2.actuate!(manipulator,0.1*ones(6))
