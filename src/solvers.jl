struct Wendlandt end
struct Zhong06 end

struct SimProblem{BotType,FuncsType}
    bot::BotType
    dynfuncs::FuncsType
    function SimProblem(bot,make_dynfuncs)
        dynfuncs = make_dynfuncs(bot)
        new{typeof(bot),typeof(dynfuncs)}(bot,dynfuncs)
    end
end

struct DiscreteCallback{condType,aType}
    condition::condType
    affect!::aType
end

const FALSE_CALLBACK_CONDITION = (intor) -> false
const TRUE_CALLBACK_CONDITION  = (intor) -> true
const DEFAULT_CALLBACK = DiscreteCallback(
                        FALSE_CALLBACK_CONDITION,
                        (integrator)->nothing)
const NO_CONTROL = (intor,cache) -> nothing

mutable struct IntegratorState{StateT}
    now::StateT
    prv::StateT
    convergence::Bool
end

struct Integrator{ProbType,StateType,CtrlType,T}
    prob::ProbType
    state::StateType
    control!::CtrlType
    tspan::Tuple{T,T}
    restart::Bool
    totalstep::Int
end

function solve!(prob::SimProblem,solver,
                control! = nothing;
                tspan=(0.0,1.0),restart=true,dt,karg...)
    (;bot,dynfuncs) = prob
    (;tg,traj) = bot
    if restart
        resize!(traj,1)
    end
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = ceil(Int,totaltime/dt)
    for istep = 1:totalstep
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = tstart + dt*istep
        if !isa(control!, Nothing)
            control!(traj[end],tg)
        end
    end
    now = deepcopy(traj[end])
    prv = deepcopy(traj[end])
    convergence = false
    state = IntegratorState(now,prv,convergence)
    intor = Integrator(prob,state,control!,tspan,restart,totalstep)
    solve!(intor,solver;dt,karg...)
end

function solve!(intor,solver;karg...)
    solvercache = generate_cache(solver,intor;karg...)
    solve!(intor,solvercache;karg...)
    # retrieve!(intor,solvercache)
    # intor,solvercache
    intor.prob.bot
end

include("solvers/Wendlandt.jl")
include("solvers/Zhong06.jl")
include("solvers/nonsmooth.jl")
include("solvers/Zhong06NSNH.jl")
