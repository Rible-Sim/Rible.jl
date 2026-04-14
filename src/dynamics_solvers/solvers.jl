

"""
Standard Dynamics Solver.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct DynamicsSolver{integratorType,body_solverType,apparatus_solverType,contact_solverType,optionsType} <: AbstractDynamicsSolver
    integrator::integratorType
    body_solver::body_solverType
    apparatus_solver::apparatus_solverType
    contact_solver::contact_solverType
    options::optionsType
end

function DynamicsSolver(integrator::AbstractIntegrator,options::NamedTuple)
    DynamicsSolver(integrator, NoBodySolver(), NoApparatusSolver(), NoContactSolver(), options)
end

function DynamicsSolver(integrator::AbstractIntegrator, body_solver=NoBodySolver();options...)
    DynamicsSolver(integrator, body_solver, NoApparatusSolver(), NoContactSolver(), values(options))
end

function DynamicsSolver(integrator::AbstractIntegrator,contact_solver::AbstractContactSolver,options::NamedTuple)
    DynamicsSolver(integrator, NoBodySolver(), NoApparatusSolver(), contact_solver, options)
end
function DynamicsSolver(integrator::AbstractIntegrator,contact_solver::AbstractContactSolver;options...)
    DynamicsSolver(integrator, NoBodySolver(), NoApparatusSolver(), contact_solver, values(options))
end

function DynamicsSolver(integrator::AbstractIntegrator,apparatus_solver::AbstractApparatusSolver;options...)
    DynamicsSolver(integrator, NoBodySolver(), apparatus_solver, NoContactSolver(), values(options))
end

function DynamicsSolver(integrator::AbstractIntegrator, apparatus_solver::AbstractApparatusSolver, contact_solver::AbstractContactSolver; options...)
    DynamicsSolver(integrator, NoBodySolver(), apparatus_solver, contact_solver, values(options))
end

abstract type AbstractObjective end

"""
Objective function parameters for dynamics problems.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Objective{weightsType,ηtype<:Function} <: AbstractObjective
    trajectory_error_gauges_weights::weightsType
    trajectory_actuators_weights::weightsType
    terminal_error_gauges_weights::weightsType
    terminal_actuators_weights::weightsType
    "weight of the trajectory cost as a function of time"
    η::ηtype 
end


function get_trajectory_cost_weights(objt::Objective, t, dt)
    return dt .* objt.η.(t)
end

function default_objective(bot::Robot)
    (;error_gauges,actuators) = bot.hub
    trajectory_error_gauges_weights = zeros(Int,length(error_gauges))
    trajectory_actuators_weights = zeros(Int,length(actuators))
    terminal_error_gauges_weights = zeros(Int,length(error_gauges))
    terminal_actuators_weights = zeros(Int,length(actuators))
    Objective(
        trajectory_error_gauges_weights,
        trajectory_actuators_weights,
        terminal_error_gauges_weights,
        terminal_actuators_weights,
        (x)->1.0
    )
end

struct DiscreteAdjointDynamicsSolver{integratorType,body_solverType,apparatus_solverType,contact_solverType,optionsType,objtType<:AbstractObjective} <: AbstractDynamicsSolver 
    forward_solver::DynamicsSolver{integratorType,body_solverType,apparatus_solverType,contact_solverType,optionsType}
    objt::objtType
end


"""
Dynamics Problem definition.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct DynamicsProblem{RobotType,policyType,envType,objtType,contact_modelType,apparatus_modelType,optionsType} <: AbstractDynamicsProblem
    bot::RobotType
    policy::policyType
    env::envType
    objt::objtType
    contact_model::contact_modelType
    apparatus_model::apparatus_modelType
    options::optionsType
end

"""
    DynamicsProblem(bot::Robot;
        env=GravityEnv(),
        policy=NoPolicy(),
        contact_model=Contactless(),
        apparatus_model=Naive(),
        kwargs...
    )

Construct a dynamics problem for the given robot.

# Arguments
- `bot::Robot`: the robot model (required, positional)
- `env`: environment model (default: `GravityEnv()`)
- `policy`: control policy (default: `NoPolicy()`)
- `contact_model`: contact model (default: `Contactless()`)
- `apparatus_model`: apparatus model (default: `Naive()`)
- `kwargs...`: additional options forwarded to the solver

# Examples
```julia
prob = DynamicsProblem(bot)
prob = DynamicsProblem(bot; env=GravityEnv(), contact_model=Zhong06())
prob = DynamicsProblem(bot; policy=my_policy, apparatus_model=my_apparatus)
```
"""
function DynamicsProblem(bot::Robot;
        env=GravityEnv(),
        policy=NoPolicy(),
        contact_model=Contactless(),
        apparatus_model=Naive(),
        kwargs...
    )
    DynamicsProblem(bot, policy, env, default_objective(bot), contact_model, apparatus_model, values(kwargs))
end

struct AdjointDynamicsSensitivitySolver{forward_solverType,adjoint_solverType} <: AbstractDynamicsSensitivitySolver
    forward_solver::forward_solverType
    adjoint_solver::adjoint_solverType
end

struct DirectDynamicsSensitivitySolver{forward_solverType,objtType<:AbstractObjective} <: AbstractDynamicsSensitivitySolver
    forward_solver::forward_solverType
    objt::objtType
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

"""
Simulator for running dynamics problems.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct Simulator{probType,ctrlType,T,dataType}
    prob::probType
    controller::ctrlType
    tspan::Tuple{T,T}
    restart::Bool
    totalstep::Int
    solver_history::SolverHistory{dataType}
    convergence::Bool
end

function Simulator(prob::AbstractDynamicsProblem,solver::AbstractDynamicsSolver,
        prescribe! = (state, structure) -> nothing;
        tspan,dt,restart=true,kwargs...
    )
    (;bot,) = prob
    (;structure,traj,contacts_traj,control_traj) = bot
    totalstep = prepare_traj!(bot;tspan,dt,restart)
    for i in eachindex(traj)
        prescribe!(traj[i],structure)
    end
    reset!(bot)
    solver_history = SolverHistory(bot,solver,totalstep)
    Simulator(prob,prescribe!,tspan,restart,totalstep,solver_history,true)
end

function prepare_traj!(bot::Robot;tspan,dt,restart=true)
    (;traj,contacts_traj,contact_caches_traj,control_traj) = bot
    if restart
        resize!(traj,1)
        resize!(contacts_traj,1)
        resize!(contact_caches_traj,1)
        resize!(control_traj,1)
    end
    nlaststep = length(traj)
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = ceil(Int,totaltime/dt)
    new_step_len = nlaststep+totalstep
    resize!(traj,new_step_len)
    resize!(contacts_traj,new_step_len)
    resize!(contact_caches_traj,new_step_len)
    resize!(control_traj,new_step_len)
    for istep = nlaststep+1:new_step_len
        traj[istep] = deepcopy(traj[nlaststep])
        contacts_traj[istep] = deepcopy(contacts_traj[nlaststep])
        for c in contacts_traj[istep]
            c.state.active = false
        end
        contact_caches_traj[istep] = deepcopy(contact_caches_traj[nlaststep])
        control_traj[istep] = deepcopy(control_traj[nlaststep])
        control_traj[istep].u .= get_initial_actions!(bot,traj.t[istep])
    end
    traj.t .= range(tstart,tend,step=dt)
    totalstep
end

"""
Result of a simulation.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct SimulationResult{simType,cacheType}
    simulator::simType
    solver_cache::cacheType
end

"""
    solve!(prob::AbstractDynamicsProblem,solver::AbstractDynamicsSolver;
        prescribe! = (state, structure) -> nothing,
        tspan,dt,restart=true,kwargs...
    )

Solve the dynamics problem with the solver, return the simulation result.
"""
function solve!(prob::AbstractDynamicsProblem,solver::AbstractDynamicsSolver;
        prescribe! = (state, structure) -> nothing,
        tspan,dt,restart=true,kwargs...
    )
    simulator = Simulator(prob,solver,prescribe!;tspan,dt,restart)
    solvercache = solve!(simulator,solver;dt,kwargs...)
    SimulationResult(simulator, solvercache)
end

"""
    solve!(prob::AbstractDynamicsProblem,solver::AdjointDynamicsSensitivitySolver;
        prescribe! = (state, structure) -> nothing,
        tspan,dt,restart=true,kwargs...
    )

Solve the dynamics problem with the adjoint solver, return the simulation result.
"""
function solve!(prob::AbstractDynamicsProblem,solver::AbstractDynamicsSensitivitySolver;
        prescribe! = (state, structure) -> nothing,
        tspan,dt,restart=true,kwargs...
    )
    simulator = Simulator(prob,solver,prescribe!;tspan,dt,restart)
    # forward
    forward_cache = solve!(simulator,solver.forward_solver;dt,kwargs...)
    # adjoint sensitivity
    #todo merge the forward and the direct solver
    (;totalstep) = simulator
    sen_cache = generate_cache(prob,solver;dt,totalstep,kwargs...)
    solve!(simulator,forward_cache,sen_cache;dt,kwargs...)
    SimulationResult(simulator, sen_cache)
end

"""
    solve!(simulator::Simulator,solver::AbstractDynamicsSolver;karg...)

Solve the simulator with the solver, return the solver cache.
"""
function solve!(simulator::Simulator,solver::AbstractDynamicsSolver;karg...)
    (;prob) = simulator
    solver_cache = generate_cache(prob,solver;karg...)
    solve!(simulator,solver_cache;karg...)
    solver_cache
end


"""
    generate_cache(
        simulator::Simulator,
        solver::AbstractDynamicsSolver;
        dt,kwargs...
    )
    
Generate cache for the solver, dispatch based on constant mass matrix trait.
"""
function generate_cache(
        prob::AbstractDynamicsProblem,
        solver::AbstractDynamicsSolver;
        kwargs...
    )

    generate_cache(
        prob,
        solver,
        has_constant_mass_matrix(prob.bot);
        kwargs...
    )
end

include("workspaces.jl")

include("complementarity_solvers.jl")

include("zhong06_solvers.jl")

include("RungKutta_family.jl")


