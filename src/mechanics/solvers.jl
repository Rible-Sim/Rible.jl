abstract type AbstractContactModel end
struct Contactless <: AbstractContactModel end
struct RestitutionFrictionCombined{restitutionType,frictionType} <: AbstractContactModel 
    restitution::restitutionType
    friction::frictionType
end

abstract type AbstractFrictionModel <: AbstractContactModel end
struct Frictionless <: AbstractFrictionModel end
struct CoulombFriction <: AbstractFrictionModel end
struct PolyhedralCoulombFriction <: AbstractFrictionModel end
struct MaximumDissipation <: AbstractFrictionModel end

abstract type AbstractRestitutionModel <: AbstractContactModel end
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

function DynamicsSolver(integrator::AbstractIntegrator)
    DynamicsSolver(integrator,nothing,nothing)
end

function DynamicsSolver(integrator::AbstractIntegrator,contact_solver::AbstractContactSolver)
    DynamicsSolver(integrator,contact_solver,nothing)
end

struct DynamicsProblem{RobotType,envType,contact_modelType,tendon_modelType}
    bot::RobotType
    env::envType
    contact_model::contact_modelType
    tendon_model::tendon_modelType
end

function DynamicsProblem(bot::Robot)
    DynamicsProblem(bot::Robot,EmptySpace(),Contactless(),NoTendon())
end

function DynamicsProblem(bot::Robot,contact_model::AbstractContactModel)
    DynamicsProblem(bot::Robot,EmptySpace(),contact_model,NoTendon())
end

function DynamicsProblem(bot::Robot,env::AbstractContactEnvironment,contact_model::AbstractContactModel)
    DynamicsProblem(bot::Robot,env,contact_model,NoTendon())
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

# include("dynamics_solvers/Wendlandt.jl")
include("dynamics_solvers/complementarity_solvers.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_constant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_CCP_constant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_CCP_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_frictionless_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_frictionless_nonconstant_mass_mono.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_sliding_cable_FB.jl")
# include("dynamics_solvers/Zhong06_family/Zhong06_nonholonomic_nonsmooth.jl")

include("dynamics_solvers/Alpha_family/Alpha.jl")
include("dynamics_solvers/Alpha_family/AlphaCCP.jl")
# include("dynamics_solvers/Alpha_family/nonsmooth.jl")
# include("dynamics_solvers/Alpha_family/NSSFC.jl")

# include("dynamics_solvers/nonsmooth.jl")
