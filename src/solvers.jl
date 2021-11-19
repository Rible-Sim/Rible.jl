abstract type ConstrainedSolver end
abstract type SlidingConstrainedSolver end
struct Wendlandt <: ConstrainedSolver end
struct Zhong06 <: ConstrainedSolver end
struct Newmark{T} <: ConstrainedSolver
    γ::T
    β::T
end

Newmark(γ=0.5,β=0.25) = Newmark(γ,β)

struct SlidingZhong06 <: SlidingConstrainedSolver end
struct SlidingNewmark{T} <: SlidingConstrainedSolver
    newmark::Newmark{T}
end

SlidingNewmark() = SlidingNewmark(Newmark())

struct SimProblem{BotType,FuncsType,ControlType,T}
    bot::BotType
    dyfuncs::FuncsType
    control!::ControlType
    tspan::Tuple{T,T}
    restart::Bool
end

# mutable struct IntegratorState{T,qT}
#     t::T
#     q::qT
#     q̇::qT
#     tprev::T
#     qprev::qT
#     q̇prev::qT
# end

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

function solve!(prob::SimProblem,solver::ConstrainedSolver;karg...)
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack tg,traj = bot
    @unpack A = dyfuncs
    if restart
        reset!(bot)
        q0 = traj.qs[begin]
        q̇0 = traj.q̇s[begin]
        λ0 = traj.λs[begin]
    else
        q0 = traj.qs[end]
        q̇0 = traj.q̇s[end]
        λ0 = traj.λs[end]
    end
    ts = [tspan[1]]
    nλ,nq = size(A(q0))
    @assert nλ == length(λ0)
    @assert nq == length(q0)
    nx = nq + nλ
    current = (t=[ts[end]],q=copy(q0),q̇=copy(q̇0))
    lasttime = deepcopy(current)
    state = (current=current,lasttime=lasttime)
    convergence = true
    intor = Integrator(prob,state,convergence,nx,nq,nλ)
    cache = generate_cache(solver,intor;karg...)
    solve!(intor,cache;karg...)
    append!(bot.traj.ts,cache.ts[2:end])
    append!(bot.traj.qs,cache.qs[2:end])
    append!(bot.traj.q̇s,cache.q̇s[2:end])
    append!(bot.traj.λs,cache.λs[2:end])
    bot
end

function solve!(prob::SimProblem,solver::SlidingConstrainedSolver;karg...)
    @unpack bot = prob
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack tg,traj = bot
    @unpack A = dyfuncs
    if restart
        reset!(bot)
        q0 = traj.qs[begin]
        q̇0 = traj.q̇s[begin]
        λ0 = traj.λs[begin]
        s̄0 = traj.s̄s[begin]
    else
        q0 = traj.qs[end]
        q̇0 = traj.q̇s[end]
        λ0 = traj.λs[end]
        s̄0 = traj.s̄s[end]
    end
    ts = [tspan[1]]
    nλ,nq = size(A(q0))
    @assert nλ == length(λ0)
    @assert nq == length(q0)
    ns̄ = length(s̄0)
    nx = nq + nλ + ns̄
    current = (t=copy(ts),q=copy(q0),q̇=copy(q̇0),s̄=copy(s̄0))
    lasttime = deepcopy(current)
    state = (current=current,lasttime=lasttime)
    convergence = true
    intor = Integrator(prob,state,convergence,nx,nq,nλ)
    cache = generate_cache(solver,intor;karg...)
    solve!(intor,cache;karg...)
    append!(bot.traj.ts,cache.ts[2:end])
    append!(bot.traj.qs,cache.qs[2:end])
    append!(bot.traj.q̇s,cache.q̇s[2:end])
    append!(bot.traj.λs,cache.λs[2:end])
    append!(bot.traj.s̄s,cache.s̄s[2:end])
    bot
end

include("solvers/Wendlandt.jl")
include("solvers/Zhong06.jl")
include("solvers/Newmark.jl")
include("solvers/nonsmooth.jl")
include("solvers/Zhong06NSNH.jl")
