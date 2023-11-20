abstract type AbstractContactModel end
struct Contactless <: AbstractContactModel end
struct FrictionRestitutionCombined{frictionType,restitutionType} <: AbstractContactModel 
    restitution::restitutionType
    friction::frictionType
end

abstract type AbstractFrictionModel end
struct Frictionless <: AbstractFrictionModel end
struct CoulombFriction <: AbstractFrictionModel end
struct PolyhedralCoulombFriction <: AbstractFrictionModel end
struct MaximumDissipation <: AbstractFrictionModel end

abstract type AbstractRestitutionModel end
struct Inelastic <: AbstractRestitutionModel end
struct NewtonRestitution <: AbstractRestitutionModel end
struct PoissonRestitution <: AbstractRestitutionModel end
struct StrangeRestitution <: AbstractRestitutionModel end

abstract type AbstractTendonModel end
struct NoTendon <: AbstractTendonModel end
struct SlidingTendon <: AbstractTendonModel end

abstract type AbstractSolver end

abstract type AbstractComplementaritySolver <: AbstractSolver end
struct InteriorPointMethod <: AbstractComplementaritySolver end
const IPM = InteriorPointMethod
struct AcceleratedProjectedGradientDescent <: AbstractComplementaritySolver end
const APGD = AcceleratedProjectedGradientDescent
struct AlternatingDirectionMethodofMultipliers <: AbstractComplementaritySolver end
const ADMM = AlternatingDirectionMethodofMultipliers
struct ProjectedGaussSeidel <: AbstractComplementaritySolver end
const PGS = ProjectedGaussSeidel
struct ProjectedGaussJacobi <: AbstractComplementaritySolver end
const PGJ = ProjectedGaussJacobi
struct SmoothedFischerBurmeister  <: AbstractComplementaritySolver end
struct SemismoothNewton  <: AbstractComplementaritySolver end

abstract type AbstractContactSolver <: AbstractSolver end

struct MonolithicContactSolver{complementarity_solverType} <: AbstractContactSolver
    complementarity_solver::complementarity_solverType
end
struct InnerLayerContactSolver{complementarity_solverType} <: AbstractContactSolver
    complementarity_solver::complementarity_solverType
end

abstract type AbstractTendonSolver <: AbstractSolver end

struct MonolithicTendonSolver{complementarity_solverType} <: AbstractTendonSolver
    complementarity_solver::complementarity_solverType
end
struct InnerLayerTendonSolver{complementarity_solverType} <: AbstractTendonSolver
    complementarity_solver::complementarity_solverType
end

abstract type AbstractIntegrator end
struct Zhong06 <: AbstractIntegrator end
struct GeneralizedAlpha{T} <: AbstractIntegrator
    αm::T
    αf::T
    γ::T
    β::T
end

struct DynamicsSolver{integratorType,contact_solverType,tendon_solverType} <: AbstractSolver 
    integrator::integratorType
    contact_solver::contact_solverType
    tendon_solver::tendon_solverType
end

function DynamicsSolver(integrator)
    DynamicsSolver(integrator,nothing,nothing)
end

struct DynamicsProblem{RobotType,contact_modelType,tendon_modelType}
    bot::RobotType
    contact_model::contact_modelType
    tendon_model::tendon_modelType
end

function DynamicsProblem(bot)
    DynamicsProblem(bot,Contactless(),NoTendon())
end
struct Simulator{ProbType,CtrlType,T}
    prob::ProbType
    controller::CtrlType
    tspan::Tuple{T,T}
    restart::Bool
    totalstep::Int
end


function Simulator(prob::DynamicsProblem,solver::DynamicsSolver,
        controller = (prescribe! = nothing, actuate! = nothing);
        tspan,dt,restart=true,karg...
    )
    (;bot,) = prob
    (;structure,traj,contacts_traj) = bot
    totalstep = prepare_traj!(traj,contacts_traj;tspan,dt,restart)
    (;prescribe!, actuate!) = controller
    if !isa(prescribe!, Nothing)
        for i in eachindex(traj)
            prescribe!(traj[i],structure)
        end
    end
    Simulator(prob,controller,tspan,restart,totalstep)
end

function GeneralizedAlpha(ρ∞)
    αm = (2ρ∞-1)/(ρ∞+1)
    αf = ρ∞/(ρ∞+1)
    γ = 1/2 + αf - αm
    β = 1/4*(γ+1/2)^2
    GeneralizedAlpha(αm,αf,γ,β)
end

function GeneralizedAlpha(ρ∞,h)
    αm = (2ρ∞-1)/(ρ∞+1)
    αf = ρ∞/(ρ∞+1)
    γ = 1/2 + αf - αm
    β = 1/4*(γ+1/2)^2
    γₜ = (1-αm)/(1-αf)/(γ*h)
    βₜ = h*β/γ - h/2
    GeneralizedAlpha(αm,αf,γ,β,γₜ,βₜ)
end

function Newmark()
    αf = αm = 0.0
    γ = 1/2
    β = 1/4
    GeneralizedAlpha(αm,αf,γ,β)
end

struct DiscreteCallback{condType,aType}
    condition::condType
    affect!::aType
end

const FALSE_CALLBACK_CONDITION = (sim) -> false
const TRUE_CALLBACK_CONDITION  = (sim) -> true
const DEFAULT_CALLBACK = DiscreteCallback(
                        FALSE_CALLBACK_CONDITION,
                        (integrator)->nothing)
const NO_CONTROL = (sim,cache) -> nothing


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

function solve!(prob::DynamicsProblem,solver::DynamicsSolver,
                controller = (prescribe! = nothing, actuate! = nothing);
                tspan,dt,restart=true,karg...)
    simulator = Simulator(prob,solver,controller;tspan,dt,restart)
    solve!(simulator,solver;dt,karg...)
end

function solve!(simulator::Simulator,solver::DynamicsSolver;karg...)
    solvercache = generate_cache(simulator,solver;karg...)
    solve!(simulator,solvercache;karg...)
    # retrieve!(simulator,solvercache)
    # simulator,solvercache
    # simulator.prob.bot
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
include("solvers/CCP/ZhongQCCPN.jl")
include("solvers/CCP/ZhongQCCPNMono.jl")

# include("solvers/nonsmooth.jl")
# include("solvers/Zhong06NSNH.jl")
