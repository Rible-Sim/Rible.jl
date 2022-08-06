abstract type AbstractSolver end

struct SimProblem{BotType,FuncsType}
    bot::BotType
    dynfuncs::FuncsType
    function SimProblem(bot,make_dynfuncs)
        dynfuncs_raw = make_dynfuncs(bot)
        if dynfuncs_raw isa NamedTuple{(:F!,)}
            dynfuncs = (F! = dynfuncs_raw.F!, Jac_F! = nothing)
        elseif dynfuncs_raw isa NamedTuple{(:F!,:Jac_F!)}
            dynfuncs = dynfuncs_raw
        elseif dynfuncs_raw isa NamedTuple{(:F!,:Jac_F!,:prepare_contacts!)}
            dynfuncs = dynfuncs_raw
        elseif dynfuncs_raw isa NamedTuple{(:F!, :apply_acu!)}
            dynfuncs = (F! = dynfuncs_raw.F!, Jac_F! = nothing, apply_acu! = dynfuncs_raw.apply_acu!)
        elseif dynfuncs_raw isa NamedTuple{(:F!, :Jac_F!, :apply_acu!)}
            dynfuncs = dynfuncs_raw
        else
            @warn "`dynfuncs` not recognized, but proceed anyway."
            dynfuncs = dynfuncs_raw
        end
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

struct Integrator{ProbType,CtrlType,T}
    prob::ProbType
    controller::CtrlType
    tspan::Tuple{T,T}
    restart::Bool
    totalstep::Int
end

function prepare_traj!(traj;tspan,dt,restart=true)
    if restart
        resize!(traj,1)
    end
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = ceil(Int,totaltime/dt)
    for istep = 1:totalstep
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = tstart + dt*istep
    end
    totalstep
end

function solve!(prob::SimProblem,solver::AbstractSolver,
                controller = (prescribe! = nothing, actuate! = nothing);
                tspan,dt,restart=true,karg...)
    (;bot,dynfuncs) = prob
    (;tg,traj,contacts_traj) = bot
    if restart; resize!(contacts_traj,1); end
    totalstep = prepare_traj!(traj;tspan,dt,restart)
    if !isa(controller.prescribe!, Nothing)
        for i in enumerate(traj)
            prescribe!(traj[i],tg)
        end
    end
    intor = Integrator(prob,controller,tspan,restart,totalstep)
    solve!(intor,solver;dt,karg...)
end

function solve!(intor,solver;karg...)
    solvercache = generate_cache(solver,intor;karg...)
    solve!(intor,solvercache;karg...)
    # retrieve!(intor,solvercache)
    # intor,solvercache
    # intor.prob.bot
end

# include("solvers/Wendlandt.jl")
include("solvers/Zhong06.jl")
include("solvers/FBZhong06.jl")
include("solvers/Alpha.jl")
include("solvers/CCP/CCPsolvers.jl")
include("solvers/CCP/ZhongCCP.jl")
include("solvers/CCP/AlphaCCP.jl")

# include("solvers/nonsmooth.jl")
# include("solvers/Zhong06NSNH.jl")
