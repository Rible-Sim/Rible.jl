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

function prepare_traj!(traj,contacts_traj;tspan,dt,restart=true)
    if restart
        resize!(traj,1)
        resize!(contacts_traj,1)
    end
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = ceil(Int,totaltime/dt)
    for istep = 1:totalstep
        push!(traj,deepcopy(traj[end]))        
        push!(contacts_traj,deepcopy(contacts_traj[end]))
        for c in contacts_traj[end]
            c.state.active = false
        end
        traj.t[end] = tstart + dt*istep
    end
    totalstep
end

function Integrator(prob::SimProblem,solver::AbstractSolver,
        controller = (prescribe! = nothing, actuate! = nothing);
        tspan,dt,restart=true,karg...
    )
    (;bot,dynfuncs) = prob
    (;structure,traj,contacts_traj) = bot
    totalstep = prepare_traj!(traj,contacts_traj;tspan,dt,restart)
    (;prescribe!, actuate!) = controller
    if !isa(prescribe!, Nothing)
        for i in eachindex(traj)
            prescribe!(traj[i],structure)
        end
    end
    intor = Integrator(prob,controller,tspan,restart,totalstep)
end

function solve!(prob::SimProblem,solver::AbstractSolver,
                controller = (prescribe! = nothing, actuate! = nothing);
                tspan,dt,restart=true,karg...)
    intor = Integrator(prob,solver,controller;tspan,dt,restart)
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
include("solvers/Zhong06Q.jl")
include("solvers/FBZhong06.jl")
include("solvers/Alpha.jl")
include("solvers/CCP/CCPsolvers.jl")
include("solvers/CCP/ZhongCCP.jl")
include("solvers/CCP/ZhongQCCP.jl")
include("solvers/CCP/AlphaCCP.jl")

# include("solvers/nonsmooth.jl")
# include("solvers/Zhong06NSNH.jl")
