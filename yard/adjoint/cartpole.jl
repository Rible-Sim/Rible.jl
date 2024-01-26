using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/nonsmooth.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\IMSD 2025\LaTex_Abstract" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/IMSD 2024/LaTex_Abstract" #src
end #src
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

tw = 455.8843 #pt |> pt2px
scalefactor = 4

#--  slider crank 
include(joinpath(pathof(RB),"../../examples/robots/cartpole.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/cartpole.jl"))#jl

coordsType=RB.QCF.QC
coordsType=RB.NCF.NC

cp_terminal = cart_pole(;coordsType)
cp_sim = cart_pole(;θ0=π/4,y0 = -1.0,coordsType)

plot_traj!(cp_terminal;showmesh=false,showground=false)
plot_traj!(cp_sim;showmesh=false,showground=false)

tspan = (0.0,0.1)
dt = 1e-2
prob = RB.DynamicsProblem(cp_sim)
dsolver = RB.DynamicsSolver(
    RB.Zhong06()
)
RB.solve!(
    prob,
    dsolver;
    tspan,
    dt
)
plot_traj!(cp_sim;showmesh=false,showground=false)

# gauges terminal
m1 = RB.measure(cp_terminal.structure,cp_terminal.hub.gauges.data[1][1])
## m2 = RB.measure(cp_terminal,cp_terminal.hub.gauges.data[2][1])

# gauges
m1 = RB.measure(cp_sim.structure,cp_sim.hub.gauges.data[1][1])
## m2 = RB.measure(cp_sim,cp_sim.hub.gauges.data[2][1])

# cost in errors
RB.cost!(cp_sim,cp_sim.traj.q[end],cp_sim.traj.q̇[end],cp_sim.traj.t[end])

RB.error_jacobian(cp_sim)
## RB.measure_jacobian(cp_sim,cp_sim.hub.gauges.data[2][1])

RB.cost_jacobian!(cp_sim,)

# cost in actuators
RB.get_num_of_actions(cp_sim.hub.actuators.data[1][1])
cp_sim.hub.actuators.data[1][1]
RB.generalized_force(cp_sim,cp_sim.hub.actuators.data[1][1])
RB.generalized_force(cp_sim,cp_sim.hub.actuators.data[2][1])
RB.actions_jacobian(cp_sim,cp_sim.hub.actuators.data[1][1])
RB.actions_jacobian(cp_sim,cp_sim.hub.actuators.data[2][1])

RB.actions_jacobian(cp_terminal,cp_terminal.hub.actuators.data[2][1])

## path_pos_vel_cost(cp_sim) # interpolate reference trajectory

## terminal_pos_vel_cost(cp_sim)

## cost = (
##     path_cost(x,ẋ,u,t),
##     terminal_cost(x,ẋ,u,tend),
## )
cp_sim = cart_pole(;θ0=π/4,y0 = -1.0,coordsType)
dsprob = RB.DynamicsSensitivityProblem(cp_sim)
adsolver = RB.DiscreteAdjointDynamicsSolver(
    RB.DynamicsSolver(
        RB.Zhong06()
    )
)
dssolver = RB.AdjointDynamicsSensitivitySolver(
    dsolver,
    adsolver
)
adprob  = RB.AdjointDynamicsProblem(cp_sim,nothing)

_,solvercache = RB.solve!(
    dsprob,
    dssolver;
    tspan,
    dt
)

solvercache.cache.∂J∂uᵀ[1,:]

RB.cost!(cp_sim,cp_sim.traj.q[end],cp_sim.traj.q̇[end],cp_sim.control_traj.u[end],cp_sim.traj.t[end])
interpolation
Interpolations

function J(u)
    tspan = (0.0,0.1)
    dt = 1e-2
    xs = tspan[1]+dt/2:dt:tspan[2]
    scaled_itp = extrapolate(scale(interpolate(u, BSpline(Constant())), xs), Line())
    bot = cart_pole(;θ0=π/4,y0 = -1.0,coordsType,f=(t)->[scaled_itp(t)])
    RB.solve!(
        RB.DynamicsProblem(bot),
        RB.DynamicsSolver(
            RB.Zhong06()
        );
        tspan,
        dt
    )
    RB.cost!(bot,bot.traj.q[end],bot.traj.q̇[end],bot.traj.t[end])
end
J(ones(10))

FiniteDiff.finite_difference_gradient(J,ones(10))

import FiniteDiff
plot_traj!(cp_sim;showmesh=false,showground=false)
