
struct Wendlandt end
struct Zhong06 end

struct DyProblem{FuncsType,qT,λT,T}
    funcs::FuncsType
    tspan::Tuple{T,T}
    q0::qT
    q̇0::qT
    λ0::λT
    nx::Int
    nq::Int
    nλ::Int
end

struct Solution{tT,qT,λT}
    ts::tT
    qs::qT
    q̇s::qT
    ps::qT
    λs::λT
end

mutable struct Integrator{ProbType,T,qT,solType}
    prob::ProbType
    t::T
    q::qT
    q̇::qT
    sol::solType
    tprev::T
    qprev::qT
    q̇prev::qT
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

function DyProblem(funcs,q0,q̇0,λ0,tspan)
    M,Φ,A,F!,Jacs = funcs
    Asize = size(A(q0))
    nq = Asize[2]
    @assert nq == length(q0)
    nλ = Asize[1]
    @assert nλ == length(λ0)
    nx = nq + nλ
    DyProblem(funcs,tspan,q0,q̇0,λ0,nx,nq,nλ)
end

function solve(prob::DyProblem,solver;karg...)
    @unpack tspan,q0,q̇0,λ0,funcs = prob
    M,Φ,A,F!,Jacs = funcs
    ts = [tspan[1]]
    qs = [copy(q0)]
    q̇s = [copy(q̇0)]
    ps = [M*copy(q̇0)]
    λs = [copy(λ0)]
    sol = Solution(ts,qs,q̇s,ps,λs)
    intor = Integrator(prob,tspan[1],copy(q0),copy(q̇0),sol,
                            tspan[1],copy(q0),copy(q̇0))
    cache = generate_cache(solver,intor)
    solve(intor,cache;karg...)
end

include("solvers/Wendlandt.jl")
include("solvers/Zhong06.jl")
