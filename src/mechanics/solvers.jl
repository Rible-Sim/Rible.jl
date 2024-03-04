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

abstract type AbstractApparatusModel end
struct Naive <: AbstractApparatusModel end
struct EulerEytelwein <: AbstractApparatusModel end

abstract type AbstractSolver end
abstract type AbstractDynamicsSolver <: AbstractSolver end

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

abstract type AbstractApparatusSolver <: AbstractSolver end

struct MonolithicApparatusSolver{complementarity_solverType} <: AbstractApparatusSolver
    complementarity_solver::complementarity_solverType
end

struct InnerLayerApparatusSolver{complementarity_solverType} <: AbstractApparatusSolver
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

struct Moreau{T} <: AbstractIntegrator
    θ::T
end

struct DynamicsSolver{integratorType,contact_solverType,apparatus_solverType,optionsType} <: AbstractDynamicsSolver 
    integrator::integratorType
    contact_solver::contact_solverType
    apparatus_solver::apparatus_solverType
    options::optionsType
end

function DynamicsSolver(integrator::AbstractIntegrator,options::NamedTuple)
    DynamicsSolver(integrator,nothing,nothing,options)
end

function DynamicsSolver(integrator::AbstractIntegrator;options...)
    DynamicsSolver(integrator,nothing,nothing,values(options))
end

function DynamicsSolver(integrator::AbstractIntegrator,contact_solver::AbstractContactSolver,options::NamedTuple)
    DynamicsSolver(integrator,contact_solver,nothing,options)
end
function DynamicsSolver(integrator::AbstractIntegrator,contact_solver::AbstractContactSolver;options...)
    DynamicsSolver(integrator,contact_solver,nothing,values(options))
end

function DynamicsSolver(integrator::AbstractIntegrator,apparatus_solver::AbstractApparatusSolver;options...)
    DynamicsSolver(integrator,nothing,apparatus_solver,values(options))
end

function DynamicsSolver(integrator::AbstractIntegrator,contact_solver::AbstractContactSolver,apparatus_solver::AbstractApparatusSolver;options...)
    DynamicsSolver(integrator,contact_solver,apparatus_solver,values(options))
end

struct DiscreteAdjointDynamicsSolver{integratorType,contact_solverType,apparatus_solverType,optionsType} <: AbstractDynamicsSolver 
    forward_solver::DynamicsSolver{integratorType,contact_solverType,apparatus_solverType,optionsType}
end

abstract type AbstractMechanicsProblem end

abstract type AbstractDynamicsProblem <: AbstractMechanicsProblem end

abstract type AbstractPolicy end
struct NoPolicy <: AbstractPolicy end
struct ActorPolicy{T} <: AbstractPolicy
    nt::T
end
struct TimePolicy{T} <: AbstractPolicy
    nt::T
end

struct DynamicsProblem{RobotType,policyType,envType,contact_modelType,apparatus_modelType,optionsType} <: AbstractDynamicsProblem
    bot::RobotType
    policy::policyType
    env::envType
    contact_model::contact_modelType
    apparatus_model::apparatus_modelType
    options::optionsType
end

function DynamicsProblem(bot::Robot;options...)
    DynamicsProblem(bot::Robot,NoPolicy(),GravitySpace(true),Contactless(),Naive(),values(options))
end

function DynamicsProblem(bot::Robot,policy::AbstractPolicy;options...)
    DynamicsProblem(bot::Robot,policy,GravitySpace(true),Contactless(),Naive(),values(options))
end

function DynamicsProblem(bot::Robot,policy::AbstractPolicy,appar_model::AbstractApparatusModel;options...)
    DynamicsProblem(bot::Robot,policy,GravitySpace(true),Contactless(),appar_model,values(options))
end

function DynamicsProblem(bot::Robot,appar_model::AbstractApparatusModel;options...)
    DynamicsProblem(bot::Robot,NoPolicy(),GravitySpace(true),Contactless(),appar_model,values(options))
end

function DynamicsProblem(bot::Robot,contact_model::AbstractContactModel;options...)
    DynamicsProblem(bot::Robot,NoPolicy(),GravitySpace(true),contact_model,Naive(),values(options))
end

function DynamicsProblem(bot::Robot,env::AbstractContactEnvironment,contact_model::AbstractContactModel;options...)
    DynamicsProblem(bot::Robot,NoPolicy(),env,contact_model,Naive(),values(options))
end

function DynamicsProblem(bot::Robot,policy::AbstractPolicy,env::AbstractContactEnvironment,contact_model::AbstractContactModel,appar_model::AbstractApparatusModel;options...)
    DynamicsProblem(bot::Robot,policy,env,contact_model,appar_model,values(options))
end

struct DynamicsSensitivityProblem{RobotType,policyType,envType,contact_modelType,apparatus_modelType,optionsType} <: AbstractDynamicsProblem
    bot::RobotType
    policy::policyType
    env::envType
    contact_model::contact_modelType
    apparatus_model::apparatus_modelType
    options::optionsType
end

function DynamicsSensitivityProblem(bot::Robot,policy=NoPolicy();options...)
    DynamicsSensitivityProblem(bot::Robot,policy,GravitySpace(),Contactless(),Naive(),values(options))
end

function DynamicsSensitivityProblem(bot::Robot,contact_model::AbstractContactModel;options...)
    DynamicsSensitivityProblem(bot::Robot,NoPolicy(),GravitySpace(),contact_model,Naive(),values(options))
end

function DynamicsSensitivityProblem(bot::Robot,env::AbstractContactEnvironment,contact_model::AbstractContactModel;options...)
    DynamicsSensitivityProblem(bot::Robot,NoPolicy(),env,contact_model,Naive(),values(options))
end

struct AdjointDynamicsSensitivitySolver{forward_solverType,adjoint_solverType} <: AbstractDynamicsSolver
    forward_solver::forward_solverType
    adjoint_solver::adjoint_solverType
end

struct AdjointDynamicsProblem{RobotType,policyType,envType,contact_modelType,apparatus_modelType,optionsType,costType} <: AbstractDynamicsProblem
    bot::RobotType
    policy::policyType
    env::envType
    contact_model::contact_modelType
    apparatus_model::apparatus_modelType
    cost::costType
    options::optionsType
end

function AdjointDynamicsProblem(bot::Robot,cost;options...)
    AdjointDynamicsProblem(bot::Robot,NoPolicy(),GravitySpace(),Contactless(),Naive(),cost,values(options))
end

struct SolverHistory{dataType}
    data::dataType
end

function SolverHistory(data::NamedTuple,n)
    SolverHistory(
        StructArray{typeof(data)}(undef,n)
    )
end

function SolverHistory(bot::Robot,solver::AbstractDynamicsSolver,n)
    T = get_numbertype(bot)
    ini_record = (
        residual=zero(T),
        iteration=1,
        walltime = 1.0,
        num_of_contacts = 2
    )
    SolverHistory(
        ini_record,
        n
    )
end

function SolverHistory(
        bot::Robot,
        solver::DynamicsSolver{IntorType,<:InnerLayerContactSolver},
        n
    ) where {IntorType}
    T = get_numbertype(bot)
    ini_record = (
        residual=zero(T),
        inner_iteration=1,
        outer_iteration=2,
        walltime = 1.0,
        num_of_contacts = 2,
        outer_condition_number = typemax(T),
        inner_condition_number = typemax(T)
    )
    SolverHistory(
        ini_record,
        n
    )
end

function SolverHistory(
        bot::Robot,
        solver::DynamicsSolver{IntorType,<:MonolithicContactSolver},
        n
    ) where {IntorType}
    T = get_numbertype(bot)
    ini_record = (
        residual=zero(T),
        iteration = 2,
        walltime = 1.0,
        num_of_contacts = 2,
        stepsizes = [one(T)],
        condition_number = typemax(T)
    )
    SolverHistory(
        ini_record,
        n
    )
end

function record!(sh::SolverHistory,data,k)
    sh.data[k] = data
end


struct Simulator{ProbType,CtrlType,T,dataType}
    prob::ProbType
    controller::CtrlType
    tspan::Tuple{T,T}
    restart::Bool
    totalstep::Int
    solver_history::SolverHistory{dataType}
end

function Simulator(prob::AbstractDynamicsProblem,solver::AbstractDynamicsSolver,
        controller = (prescribe! = nothing, actuate! = nothing);
        tspan,dt,restart=true,karg...
    )
    (;bot,) = prob
    (;structure,traj,contacts_traj,control_traj) = bot
    totalstep = prepare_traj!(bot;tspan,dt,restart)
    (;prescribe!, actuate!) = controller
    if !isa(prescribe!, Nothing)
        for i in eachindex(traj)
            prescribe!(traj[i],structure)
        end
    end
    solver_history = SolverHistory(
        bot,
        solver,
        totalstep
    )
    Simulator(prob,controller,tspan,restart,totalstep,solver_history)
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


function prepare_traj!(bot::Robot;tspan,dt,restart=true)
    (;traj,contacts_traj,control_traj) = bot
    if restart
        resize!(traj,1)
        resize!(contacts_traj,1)
        resize!(control_traj,1)
    end
    nlaststep = length(traj)
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = ceil(Int,totaltime/dt)
    resize!(traj,nlaststep+totalstep)
    resize!(contacts_traj,nlaststep+totalstep)
    resize!(control_traj,nlaststep+totalstep)
    for istep = nlaststep+1:nlaststep+totalstep
        traj[istep] = deepcopy(traj[nlaststep])
        traj.t[istep] = tstart + dt*(istep-1)
        contacts_traj[istep] = deepcopy(contacts_traj[nlaststep])
        for c in contacts_traj[istep]
            c.state.active = false
        end
        control_traj[istep] = deepcopy(control_traj[nlaststep])
        ## control_traj[istep].u .= get_actions!(bot,traj.t[istep])
    end
    totalstep
end

function solve!(prob::AbstractDynamicsProblem,solver::AbstractDynamicsSolver,
                controller = (prescribe! = nothing, actuate! = nothing);
                tspan,dt,restart=true,karg...)
    simulator = Simulator(prob,solver,controller;tspan,dt,restart)
    solvercache = solve!(simulator,solver;dt,karg...)
    simulator, solvercache
end

function solve!(simulator::Simulator,solver::AbstractDynamicsSolver;karg...)
    solvercache = generate_cache(simulator,solver;karg...)
    solve!(simulator,solvercache;karg...)
    # retrieve!(simulator,solvercache)
    solvercache
    # simulator.prob.bot
end

function generate_cache(
        simulator::Simulator,
        solver::AbstractDynamicsSolver;
        dt,kargs...
    ) 

    generate_cache(
        simulator,
        solver,
        has_constant_mass_matrix(simulator.prob.bot);
        dt,kargs...
    )

end



# include("dynamics_solvers/Wendlandt.jl")
include("dynamics_solvers/complementarity_solvers.jl")

include("dynamics_solvers/Moreau_family/Moreau_constant_mass.jl")
include("dynamics_solvers/Moreau_family/Moreau_CCP_constant_mass.jl")
# include("dynamics_solvers/Moreau_family/Moreau_CCP_nonconstant_mass.jl")

include("dynamics_solvers/Zhong06_family/Zhong06_momentum.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_constant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_CCP_constant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_CCP_constant_mass_mono.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_CCP_constant_mass_cluster_cables.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_CCP_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_frictionless_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_frictionless_nonconstant_mass_mono.jl")
include("dynamics_solvers/Zhong06_family/Zhong06_constant_mass_cluster_cables.jl")
# include("dynamics_solvers/Zhong06_family/Zhong06_nonholonomic_nonsmooth.jl")

include("dynamics_solvers/Zhong06_family/Adjoint_Zhong06_nonconstant_mass.jl")
include("dynamics_solvers/Zhong06_family/Adjoint_Zhong06_constant_mass.jl")


include("dynamics_solvers/Alpha_family/Alpha.jl")
include("dynamics_solvers/Alpha_family/AlphaCCP.jl")
# include("dynamics_solvers/Alpha_family/nonsmooth.jl")
# include("dynamics_solvers/Alpha_family/NSSFC.jl")

# include("dynamics_solvers/nonsmooth.jl")
