using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/adjoint.jl"))
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
using OptimizationOptimJL
tw = 455.8843 #pt |> pt2px
scalefactor = 4

#--  slider crank 
include(joinpath(pathof(RB),"../../examples/robots/cartpole.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/cartpole.jl"))#jl

coordsType=RB.QCF.QC
coordsType=RB.NCF.NC

cp_terminal = cart_pole(;coordsType)
f(t) = [1.0]
@time cp_sim = cart_pole(;θ0=π/4,y0 = -1.0,coordsType,f)


function forward_dyn(u = ones(10),params=nothing;
        tspan = (0.0,0.1),
        dt = 1e-2
    )
    times = tspan[1]+dt/2:dt:tspan[2]
    f(t) = [extrapolate(scale(interpolate(u, BSpline(Constant())), times), Line())(t)]
    bot = cart_pole(;θ0=π/4,y0 = -1.0,coordsType,f)
    RB.solve!(
        RB.DynamicsProblem(bot),
        RB.DynamicsSolver(
            RB.Zhong06()
        );
        tspan,
        dt
    )
    bot
end

function J(u = ones(10),params=nothing;
        tspan = (0.0,0.1),
        dt = 1e-2
    )
    bot = forward_dyn(u;tspan,dt)    
    RB.cost!(bot,bot.traj.q[end],bot.traj.q̇[end],bot.traj.t[end])
end

function dJ!(G, u = ones(10),params=nothing;
        tspan = (0.0,0.1),
        dt = 1e-2
    )
    times = tspan[1]+dt/2:dt:tspan[2]
    f(t) = [extrapolate(scale(interpolate(u, BSpline(Constant())), times), Line())(t)]
    bot = cart_pole(;θ0=π/4,y0 = -1.0,coordsType,f)
    dsprob = RB.DynamicsSensitivityProblem(bot)
    adsolver = RB.DiscreteAdjointDynamicsSolver(
        RB.DynamicsSolver(
            RB.Zhong06()
        )
    )
    dssolver = RB.AdjointDynamicsSensitivitySolver(
        dsolver,
        adsolver
    )
    ## adprob  = RB.AdjointDynamicsProblem(bot,nothing)  )
    _,solvercache = RB.solve!(
        dsprob,
        dssolver;
        tspan,
        dt
    )
    G .= solvercache.cache.∂J∂uᵀ[1,:]
end

cp_sim = forward_dyn()
J()
dJ()

J(ones(10))

FiniteDiff.finite_difference_gradient(J,ones(10))

optprob = OptimizationFunction(J, Optimization.AutoFiniteDiff())
prob = Optimization.OptimizationProblem(optprob, ones(10), )
sol = solve(prob, Optim.LBFGS())
sol.u

optprob = OptimizationFunction(J, grad = dJ!)
prob = Optimization.OptimizationProblem(optprob, ones(10),)
sol = solve(prob, Optim.LBFGS())
sol.u
cp_sim = forward_dyn(sol.u)

plot_traj!(cp_terminal;showmesh=false,showground=false)
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



plot_traj!(bot;showmesh=false,showground=false)