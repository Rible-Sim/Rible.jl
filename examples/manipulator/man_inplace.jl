using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
using LaTeXStrings
# using NLsolve
using Revise
# using TS
# using Robot2D
# const TR = Robot2D
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

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

        prop = TR.RigidBodyProperty(i,movable,m,Ia,
                    CoM,aps
                    )
        state = TR.RigidBodyState(prop,ri,rj)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,m[i],a[i],
            Ia[i],A[:,i],A[:,i+1]) for i = 1:nbody]

    nstring = 2(nbody-1)
    upstringlen = norm(rbs[2].state.p[3] - rbs[1].state.p[1])
    lostringlen = norm(rbs[2].state.p[2] - rbs[1].state.p[3])
    original_restlens = zeros(nstring)
    restlens = zeros(nstring)
    actuallengths = zeros(nstring)
    ks = zeros(nstring)
    cs = zeros(nstring)
    cs .= 10.0
    for i = 1:nstring
        j = i % 4
        original_restlens[i] =
                 restlens[i] =
               actuallengths[i] =  ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],1.0e2,1.0e2)
    end
    ss = [TR.SString2D(original_restlens[i],ks[i],cs[i]) for i = 1:nstring]
    acs = [
        ifelse(isodd(i),TR.Actuator(ss[2(i-1)+1:2i]),
                        TR.Actuator(ss[2i:-1:2(i-1)+1]))
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
    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:length(rbs)-1
        push!(string2ap,(TR.ID(i,1),TR.ID(i+1,3)))
        push!(string2ap,(TR.ID(i,3),TR.ID(i+1,2)))
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
# ------------------Create Tensegrity Struture --------------------------
ndof = 6
refman = man_ndof(ndof,-π/12) # reference
manipulator = man_ndof(ndof,0.0)
# ------------------Create Tensegrity Struture\\-------------------------

q0,q̇0,λ0 = TR.get_initial(manipulator) #backup
# ----------------------Inverse Kinematics ------------------------------
function inverse(tgstruct,refst2d)
    function ikfuncs(tgstruct)

        A = TR.build_A(tgstruct)

        Q̃=TR.build_Q̃(tgstruct)

        function F!(F,u)
            TR.reset_forces!(tgstruct)
            TR.actuate!(tgstruct,u)
            TR.update_strings_apply_forces!(tgstruct)
            F .= Q̃*TR.fvector(tgstruct)
        end

        A,F!
    end
    q0,q̇0,λ0 = TR.get_initial(refst2d)
    TR.distribute_q_to_rbs!(tgstruct,q0,q̇0)
    nu = length(tgstruct.actuators)
    u0 = zeros(nu)
    ikprob = TS.IKProblem(ikfuncs(tgstruct),q0,u0,λ0)
    TR.iksolve(ikprob)
end
u,_ = inverse(manipulator,refman) 

TR.distribute_q_to_rbs!(manipulator,q0,q̇0)
TR.actuate!(manipulator,zero(u)) # reverse to initial
# ----------------------Inverse Kinematics\\-----------------------------

dt = 0.01 # Same dt used for PID AND Dynamics solver
# ---------------------Create Controllers --------------------------------
# use angles to define errors
function get_angles(tgstruct)
    rbs = tgstruct.rigidbodies
    angles = zeros(tgstruct.nbody-1)
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
pids = [TR.PIDController.PID(0.01,0.0,0.01,
            setpoint=ui,dt=dt) for ui in refangles]
# ---------------------Create Controllers\\-------------------------------

# --------------------Create Robot ---------------------------------------
rob = TR.TGRobot2D(manipulator,TR.ControlHub(pids))
# --------------------Create Robot\\--------------------------------------

# --------------------Define Control Action-------------------------------
function make_control!(get_feedback)
    function inner_control!(robot2d,t)
        @unpack tgstruct, hub = robot2d
        @unpack ctrls,trajs = hub
        inputs = get_feedback(tgstruct)
        for (id,pid,traj,actuator) in zip(eachindex(ctrls),ctrls,trajs,tgstruct.actuators)
            input = inputs[id]
            output = TR.PIDController.update!(pid,input)
            TR.actuate!(actuator,output)
            TR.record!(traj,pid)
        end
    end
end

control! = make_control!(get_angles)
# --------------------Define Control Action\\-----------------------------

# ----------------------------Dynamics-----------------------------------
function dynfuncs(tgstruct)

    M = TR.build_massmatrix(tgstruct)
    #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
    Φ = TR.build_Φ(tgstruct)
    A = TR.build_A(tgstruct)

    Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        #F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
        #@show isapprox(F,)
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(manipulator)
prob = TS.DyProblem(dynfuncs(manipulator),q0,q̇0,λ0,(0.0,40.0))
# TR.actuate!(manipulator,u)
# sol = TS.solve(prob,dt=dt,ftol=1e-13,verbose=true)

function make_affect!(robot2d,control!)
    function inner_affect!(intor)
        TR.distribute_q_to_rbs!(robot2d.tgstruct,intor.qprev,intor.q̇prev)
        TR.update_strings_apply_forces!(robot2d.tgstruct)
        control!(robot2d,intor.t)
    end
end
cb = TS.DiscreteCallback((x)->true,make_affect!(rob,make_control!(get_angles)))
TR.PIDController.tune!(rob.hub.ctrls[1],1.4,0.004,13)
TR.PIDController.tune!(rob.hub.ctrls[2],1.3,0.004,12)
TR.PIDController.tune!(rob.hub.ctrls[3],1.2,0.004,11)
TR.PIDController.tune!(rob.hub.ctrls[4],1.1,0.004,10)
TR.PIDController.tune!(rob.hub.ctrls[5],1.0,0.005,8.5)
TR.PIDController.tune!(rob.hub.ctrls[6],0.9,0.005,7.0)
TR.reset!.(rob.tgstruct.actuators)
TR.reset!.(rob.hub.trajs)
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-13,callback=cb)
pltfig.clear(); pltfig = controlplot(rob.hub.trajs)
pltfig = controlplot(rob.hub.trajs)
sol = TS.solve(prob,dt=dt,ftol=1e-13,callback=cb,verbose=true)
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

sol = TS.solve(prob,dt=dt,ftol=1e-13,verbose=true)
sol = TS.solve(prob,dt=dt,ftol=1e-13)






TR.distribute_q_to_rbs!(manipulator,sol.qs[end],sol.q̇s[end])
endangles = get_angles(manipulator)

sol = TS.solve(prob,dt=dt,ftol=1e-13,verbose=true,callback=cb)

@time TS.solve(prob,dt=dt,ftol=1e-13)
s = 1
tab = SPARKTableau(s)
tspan = (0.0,20.0)
cache = SPARKCache(size(A(q0))[2],size(A(q0))[1],dt,tspan,s,man_wend(n,manipulator))
state = SPARKsolve!(q0,q̇0,λ0,cache,tab,ftol=1e-14)

@btime TS.solve(prob,dt=dt)
@time TS.solve(prob,dt=dt)
TS.solve(prob,dt=dt)
@code_warntype TS.solve(prob,dt=dt)

function showactuator(tgstruct)
    for (acid,actuator) in enumerate(tgstruct.actuators)
        @unpack strings = actuator
        str1 = strings[1]
        u1 = str1.state.restlen - str1.original_restlen
        str2 = strings[2]
        u2 = str2.state.restlen - str2.original_restlen
        @show acid,u1,u2
    end
end
showactuator(manipulator)
TR.actuate!(manipulator,0.1*ones(6))
