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
using Robot
const TR = Robot

cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")

# ------------------Create Struture --------------------------
num_of_dof = 6
refman = man_ndof(num_of_dof,θ=-π/12) # reference
manipulator = man_ndof(num_of_dof,k=4.e3,c=1.e3,θ=0.0)
# ------------------Create Struture\\-------------------------
Y = build_Y(manipulator)
q0,q̇0,λ0 = RB.get_initial(manipulator) #backup
# ----------------------Inverse Kinematics ------------------------------
refλ0,_,a = RB.inverse(manipulator,refman,Y)
# ----------------------Inverse Kinematics\\-----------------------------

dt = 0.01 # Same dt used for Dynamics solver
# ---------------------Create Controllers --------------------------------
# use angles to define errors

refangles = get_angles(refman)
refman
pids = [RB.PIDController.PID(0.0,0.0,0.0,
            setpoint=ui,dt=0.1) for ui in refangles]
# ---------------------Create Controllers\\-------------------------------

# --------------------Create Robot ---------------------------------------
rob = RB.TGRobot2D(manipulator,RB.ControlHub(pids))
# --------------------Create Robot\\--------------------------------------

# --------------------Define Control Action-------------------------------
function make_control!(get_feedback)
    function inner_control!(robot2d,t)
        @unpack tgstruct, hub = robot2d
        @unpack ctrls,trajs = hub
        inputs = get_feedback(tgstruct)
        for (id,pid,traj,actuator) in zip(eachindex(ctrls),ctrls,trajs,tgstruct.actuators)
            input = inputs[id]
            output = RB.PIDController.update!(pid,input,t)
            RB.actuate!(actuator,output)
            RB.record!(traj,pid)
        end
    end
end

control! = make_control!(get_angles)
# --------------------Define Control Action\\-----------------------------

# ----------------------------Dynamics-----------------------------------
function dynfuncs(tgstruct,q0)

    M = RB.build_massmatrix(tgstruct)
    Φ = RB.build_Φ(tgstruct,q0)
    A = RB.build_A(tgstruct)

    function F!(F,q,q̇,t)
        RB.reset_forces!(tgstruct)
        RB.distribute_q_to_rbs!(tgstruct,q,q̇)
        RB.update_cables_apply_forces!(tgstruct)
        RB.assemble_forces!(F,tgstruct)
        # @show isapprox(F,Q̃*RB.fvector(tgstruct))
    end

    M,Φ,A,F!,nothing
end
M,Φ,A,F!,_ = dynfuncs(manipulator,q0)

prob = TS.DyProblem(dynfuncs(manipulator,q0),q0,q̇0,λ0,(0.0,25.01))

function make_affect!(robot2d,control!)
    function inner_affect!(sim)
        RB.distribute_q_to_rbs!(robot2d.tgstruct,sim.qprev,sim.q̇prev)
        RB.update_cables_apply_forces!(robot2d.tgstruct)
        control!(robot2d,sim.tprev)
    end
end
cb = TS.DiscreteCallback((x)->true,make_affect!(rob,make_control!(get_angles)))
# RB.PIDController.tune!(rob.hub.ctrls[1],0.2,0.8,3.0)
# RB.PIDController.tune!(rob.hub.ctrls[2],0.1,1.0,1.2)
# RB.PIDController.tune!(rob.hub.ctrls[3],0.1,1.4,1.0)
RB.PIDController.tune!(rob.hub.ctrls[1],5.0,0.0,10.0)
RB.PIDController.tune!(rob.hub.ctrls[2],5.0,0.0,15.0)
RB.PIDController.tune!(rob.hub.ctrls[3],5.0,0.0,10.0)
RB.PIDController.tune!(rob.hub.ctrls[4],5.0,0.0,6.0)
RB.PIDController.tune!(rob.hub.ctrls[5],5.0,0.0,4.0)
RB.PIDController.tune!(rob.hub.ctrls[6],5.0,0.0,2.0)
RB.reset!.(rob.tgstruct.actuators)
RB.reset!.(rob.hub.trajs)
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
            @unpack ts,restitution_coefficients,us = trajs[id]
            bx = ax.twinx()
            ax.grid(which="major", axis="both")
            # bx.grid("on")
            ep = ax.plot(ts,restitution_coefficients,label=latexstring("\\epsilon_$id"), lw = 3)
            up = bx.plot(ts,us,label=latexstring("u_$id"), lw = 0.02, color=:orange)
            ps = [ep[1],up[1]]
            bx.set_ylabel(L"u")
            bx.set_ylim([0.0,1.0])
            ax.set_ylim(restitution_coefficients[1].*[-0.2,1.2])
            ax.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(0.2es[1]))
            ax.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(xmax=restitution_coefficients[1]))
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

kes = [RB.kinetic_energy_coords(manipulator,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
pes = [RB.potential_energy(manipulator,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
restitution_coefficients = kes .+ pes

function showactuator(tgstruct)
    for (acid,actuator) in enumerate(tgstruct.actuators)
        @unpack cables = actuator
        str1 = cables[1]
        u1 = str1.state.restlen - str1.original_restlen
        str2 = cables[2]
        u2 = str2.state.restlen - str2.original_restlen
        @show acid,u1,u2
    end
end
showactuator(manipulator)
RB.actuate!(manipulator,0.1*ones(6))
