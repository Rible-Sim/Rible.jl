using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
using LaTeXStrings
# using NLsolve
using Revise

using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

include("man_define.jl")

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
u,refλ0 = inverse(manipulator,deepcopy(manipulator))
TR.distribute_q_to_rbs!(manipulator,q0,q̇0)
TR.actuate!(manipulator,zero(u)) # reverse to initial
# ----------------------Inverse Kinematics\\-----------------------------

dt = 0.01 # Same dt used for PID AND Dynamics solver
# ---------------------Create Controllers --------------------------------
# use angles to define errors
function get_angles(tgstruct)
    rbs = tgstruct.rigidbodies
    angles = zeros(tgstruct.nbodies-1)
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
function dynfuncs(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        F .= Q̃*TR.fvector(tgstruct)
        # F .= 0.0
        # TR.assemble_forces!(F,tgstruct)
        # @show isapprox(F,Q̃*TR.fvector(tgstruct))
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(manipulator,q0)
Φ(q0)
Φq = A(q0)

function perturbation!(q0,q̇0)
    q̇0[end] = 0.0001
end
perturbation!(q0,q̇0)
prob = TS.DyProblem(dynfuncs(manipulator,q0),q0,q̇0,λ0,(0.0,50.0))
# TR.actuate!(manipulator,u)
# sol = TS.solve(prob,dt=dt,ftol=1e-13,verbose=true)
q = TR.undamped_modal_solve!(manipulator,q0,q̇0,λ0,50.0,dt)
plot(transpose(q))

function make_affect!(robot2d,control!)
    function inner_affect!(intor)
        TR.distribute_q_to_rbs!(robot2d.tgstruct,intor.qprev,intor.q̇prev)
        TR.update_strings_apply_forces!(robot2d.tgstruct)
        control!(robot2d,intor.tprev)
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
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-13,callback=cb,verbose=true)
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-12,callback=cb,verbose=true)
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

sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-12,verbose=true)
δq = [q-q0 for q in sol.qs]
PyPlot.plot()
using Plots
plot(transpose(hcat(δq...)))
pyplot()

kes = [TR.kinetic_energy_coords(manipulator,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
pes = [TR.potential_energy(manipulator,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
es = kes .+ pes

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
