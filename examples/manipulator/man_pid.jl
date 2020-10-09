using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
# using NLsolve
using Makie
AbstractPlotting.__init__()
using Revise

using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")

# ------------------Create Tensegrity Struture --------------------------
ndof = 6
refman = man_ndof(ndof,θ=-π/12) # reference
manipulator = man_ndof(ndof,k=4.e3,c=1.e3,θ=0.0)
# ------------------Create Tensegrity Struture\\-------------------------
Y = build_Y(manipulator)
q0,q̇0,λ0 = TR.get_initial(manipulator) #backup
# ----------------------Inverse Kinematics ------------------------------
refλ0,_,a = TR.inverse(manipulator,refman,Y)
# ----------------------Inverse Kinematics\\-----------------------------

dt = 0.01 # Same dt used for Dynamics solver
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
pids = [TR.PIDController.PID(0.0,0.0,0.0,
            setpoint=ui,dt=0.1) for ui in refangles]
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
            output = TR.PIDController.update!(pid,input,t)
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
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        TR.assemble_forces!(F,tgstruct)
        # @show isapprox(F,Q̃*TR.fvector(tgstruct))
    end

    M,Φ,A,F!,nothing
end
M,Φ,A,F!,_ = dynfuncs(manipulator,q0)

prob = TS.DyProblem(dynfuncs(manipulator,q0),q0,q̇0,λ0,(0.0,25.01))

function make_affect!(robot2d,control!)
    function inner_affect!(intor)
        TR.distribute_q_to_rbs!(robot2d.tgstruct,intor.qprev,intor.q̇prev)
        TR.update_strings_apply_forces!(robot2d.tgstruct)
        control!(robot2d,intor.tprev)
    end
end
cb = TS.DiscreteCallback((x)->true,make_affect!(rob,make_control!(get_angles)))
# TR.PIDController.tune!(rob.hub.ctrls[1],0.2,0.8,3.0)
# TR.PIDController.tune!(rob.hub.ctrls[2],0.1,1.0,1.2)
# TR.PIDController.tune!(rob.hub.ctrls[3],0.1,1.4,1.0)
TR.PIDController.tune!(rob.hub.ctrls[1],5.0,0.0,10.0)
TR.PIDController.tune!(rob.hub.ctrls[2],5.0,0.0,15.0)
TR.PIDController.tune!(rob.hub.ctrls[3],5.0,0.0,10.0)
TR.PIDController.tune!(rob.hub.ctrls[4],5.0,0.0,6.0)
TR.PIDController.tune!(rob.hub.ctrls[5],5.0,0.0,4.0)
TR.PIDController.tune!(rob.hub.ctrls[6],5.0,0.0,2.0)
TR.reset!.(rob.tgstruct.actuators)
TR.reset!.(rob.hub.trajs)
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-9,callback=cb,verbose=true)
pltfig = controlplot(rob.hub.trajs,6)
pltfig.savefig("manpid_error.png",dpi=300,bbox_inches="tight")
plt.close(pltfig)
plotstructure(manipulator,sol,sliderplot)
# ----------------------------Dynamics-----------------------------------
function controlplot(trajs,n=6)
    ntraj = length(trajs)
    fig,axs_raw = plt.subplots(2,3,figsize=(15,6))

    if typeof(axs_raw)<:Array
        axs = permutedims(axs_raw,[2,1])
    else
        axs = [axs_raw]
    end
    for (id,ax) in enumerate(axs)
        if id <= n
            @unpack ts,es,us = trajs[id]
            bx = ax.twinx()
            ax.grid(which="major", axis="both")
            # bx.grid("on")
            ep = ax.plot(ts,es,label=latexstring("\\epsilon_$id"), lw = 3)
            up = bx.plot(ts,us,label=latexstring("u_$id"), lw = 0.02, color=:orange)
            ps = [ep[1],up[1]]
            bx.set_ylabel(L"u")
            bx.set_ylim([0.0,1.0])
            ax.set_ylim(es[1].*[-0.2,1.2])
            ax.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(0.2es[1]))
            ax.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(xmax=es[1]))
            ax.set_ylabel(L"\epsilon(\%)")
            ax.set_xlabel(L"t")
            ax.legend(ps, [p_.get_label() for p_ in ps])
            if id <= 3
                ax.set_xticklabels([])
                ax.xaxis.label.set_visible(false)
            end
            if !(id ∈[1,4])
                ax.set_yticklabels([])
                ax.yaxis.label.set_visible(false)
            end
            if !(id ∈ [3,6])
                bx.set_yticklabels([])
                bx.yaxis.label.set_visible(false)
            end

        end

    end
    fig
end

tstops = [0,10,20,30,40,50]
man_fig = pyplotstructure(manipulator,sol,tstops)
man_fig.savefig("manpid.png",dpi=300,bbox_inches="tight")
plt.close(man_fig)

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
