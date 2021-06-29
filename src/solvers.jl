
struct Wendlandt end
struct Zhong06 end

struct SimProblem{BotType,FuncsType,ControlType,T}
    bot::BotType
    dyfuncs::FuncsType
    control!::ControlType
    tspan::Tuple{T,T}
    restart::Bool
end

mutable struct IntegratorState{T,qT}
    t::T
    q::qT
    q̇::qT
    tprev::T
    qprev::qT
    q̇prev::qT
end

struct Integrator{ProbType,StateType}
    prob::ProbType
    state::StateType
    convergence::Bool
    nx::Int
    nq::Int
    nλ::Int
end

const FALSE_CALLBACK_CONDITION = (intor) -> false
const TRUE_CALLBACK_CONDITION  = (intor) -> true
struct DiscreteCallback{condType,aType}
    condition::condType
    affect!::aType
end
const DEFAULT_CALLBACK = DiscreteCallback(
                        FALSE_CALLBACK_CONDITION,
                        (integrator)->nothing)

const NO_CONTROL = (intor,cache) -> nothing

function SimProblem(bot,make_dyfuncs,tspan::Tuple{T,T};restart=true) where T
    SimProblem(bot,make_dyfuncs(bot),NO_CONTROL,tspan,restart)
end

function SimProblem(bot,make_dyfuncs,control!,tspan::Tuple{T,T};restart=true) where T
    SimProblem(bot,make_dyfuncs(bot),control!,tspan,restart)
end

function solve(prob::SimProblem,solver;karg...)
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack tg,traj = bot
    M,Φ,A,F!,Jacs = dyfuncs
    if restart
        reset!(traj)
        q0 = traj.qs[begin]
        q̇0 = traj.q̇s[begin]
        λ0 = traj.λs[begin]
    else
        q0 = traj.qs[end]
        q̇0 = traj.q̇s[end]
        λ0 = traj.λs[end]
    end
    ts = [tspan[1]]
    Asize = size(A(q0))
    nq = Asize[2]
    @assert nq == length(q0)
    nλ = Asize[1]
    @assert nλ == length(λ0)
    nx = nq + nλ
    state = IntegratorState(ts[end],copy(q0),copy(q̇0),ts[end],copy(q0),copy(q̇0))
    convergence = true
    intor = Integrator(prob,state,convergence,nx,nq,nλ)
    cache = generate_cache(solver,intor;karg...)
    solve!(intor,cache;karg...)
end

function solve!(prob::SimProblem,solver;karg...)
    @unpack bot = prob
    intor,cache = solve(prob,solver;karg...)
    append!(bot.traj.ts,cache.ts[2:end])
    append!(bot.traj.qs,cache.qs[2:end])
    append!(bot.traj.q̇s,cache.q̇s[2:end])
    append!(bot.traj.λs,cache.λs[2:end])
    bot
end

include("solvers/Wendlandt.jl")
include("solvers/Zhong06.jl")
include("solvers/nonsmooth.jl")
include("solvers/Zhong06NSNH.jl")
