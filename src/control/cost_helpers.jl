struct InitialCost end
struct StructureCost end
struct ControlCost end

function dyn_solve!(bot,policy,env::AbstractEnvironment,integrator::AbstractIntegrator; 
        kwargs...
    )
    if iscontact(env)
        dyn_prob = DynamicsProblem(bot;
            policy, env,
            contact_model=RestitutionFrictionCombined(
                NewtonRestitution(),
                CoulombFriction(),
            )
        )
        dyn_solver = DynamicsSolver(
            integrator,
            MonolithicContactSolver(
                InteriorPointMethod()
            )
        )
    else
        dyn_prob = DynamicsProblem(bot; policy, env)
        dyn_solver = DynamicsSolver(
            integrator
        )
    end
    solve!(
        dyn_prob,dyn_solver;kwargs...
    )
end

function adj_sen_solve!(bot,policy,env::AbstractEnvironment,integrator::AbstractIntegrator,objt,;
        kwargs...
    )

    if iscontact(env)
        dyn_prob = DynamicsProblem(bot;
            policy, env,
            contact_model=RestitutionFrictionCombined(
                NewtonRestitution(),
                CoulombFriction())
        )
        dyn_solver = DynamicsSolver(
            integrator,
            MonolithicContactSolver(
                InteriorPointMethod()
            )
        )
    else
        dyn_prob = DynamicsProblem(bot; policy, env)
        dyn_solver = DynamicsSolver(
            integrator
        )
    end
    adj_solver = DiscreteAdjointDynamicsSolver(
        dyn_solver,
        objt,
    )
    adj_sen_solver = AdjointDynamicsSensitivitySolver(
        dyn_solver,
        adj_solver
    )
    sim = solve!(
        dyn_prob,
        adj_sen_solver;
        kwargs...
    )
end

function drc_sen_solve!(bot,policy,env::AbstractEnvironment,integrator, objt,;
        kwargs...
    )

    if iscontact(env)
        dyn_prob = DynamicsProblem(bot;
            policy, env,
            contact_model=RestitutionFrictionCombined(
                NewtonRestitution(),
                CoulombFriction())
        )
        dyn_solver = DynamicsSolver(
            integrator,
            MonolithicContactSolver(
                InteriorPointMethod()
            )
        )
    else
        dyn_prob = DynamicsProblem(bot; policy, env)
        dyn_solver = DynamicsSolver(
            integrator
        )
    end
    drc_sen_solver = DirectDynamicsSensitivitySolver(
        dyn_solver,
        objt,
    )
    sim = solve!(
        dyn_prob,
        drc_sen_solver;
        kwargs...
    )
end

function make_cost(bot::Robot,policy::AbstractPolicy,env::AbstractEnvironment,integrator::AbstractIntegrator,objt::AbstractObjective,mode;
        idx=:,
        params_func = (x)->x,
        dt,kwargs...
    ) 
    # Define the cost function here (e.g., energy consumption, trajectory tracking, etc.
    function inner_cost(x,opt_params=nothing)
        if mode isa InitialCost
            set_initial!(bot,params_func(x);idx)
        elseif mode isa StructureCost
            set_structure_params!(bot,params_func(x);idx)
        elseif mode isa ControlCost
            set_control_params!(bot,policy,params_func(x);idx)
        end
        dyn_solve!(bot,policy,env,integrator;dt,kwargs...)
        J = cost!(bot,objt,integrator,dt;)
        J
    end
end

function make_cost_gradient(bot::Robot,policy::AbstractPolicy,env::AbstractEnvironment,integrator::AbstractIntegrator,objt::AbstractObjective,mode;
        idx=:,
        ∂q0∂y0 = I,
        ∂q̇0∂ẏ0 = I,
        params_func = (x)->x,
        kwargs...
    )
    nq̌ = get_num_of_free_coords(bot.structure)
    (;M) = build_mass_matrices(bot)
    ∂p0∂q̇0 = M
    ∂p0∂ẏ0 = ∂p0∂q̇0*∂q̇0∂ẏ0
    function inner_cost_gradient(G, x,opt_params=nothing)
        if mode isa InitialCost
            set_initial!(bot,params_func(x);idx)
        elseif mode isa StructureCost
            set_structure_params!(bot,params_func(x);idx)
        elseif mode isa ControlCost
            set_control_params!(bot,policy,params_func(x);idx)
        end
        sim = adj_sen_solve!(
            bot,policy, env,integrator,objt;
            kwargs...
        )
        if mode isa InitialCost
            Gy0 = ∂q0∂y0'*sim.solver_cache.∂J∂x₀ᵀ[   1: nq̌]
            Gẏ0 = ∂p0∂ẏ0'*sim.solver_cache.∂J∂x₀ᵀ[nq̌+1:2nq̌]
            G .= vcat(Gy0,Gẏ0)
        elseif mode isa StructureCost
            G .= sum(sim.solver_cache.∂J∂cᵀ)[idx]
        elseif mode isa ControlCost
            G .= sum(sim.solver_cache.∂J∂θᵀ)[idx]
        end
    end
end

function make_cost_and_gradient(bot::Robot,policy::AbstractPolicy,env::AbstractEnvironment,integrator::AbstractIntegrator,objt::AbstractObjective,mode;
        idx=:,
        ∂q0∂y0 = I,
        ∂q̇0∂ẏ0 = I,
        params_func = (x)->x,
        tspan,dt,kwargs...
    )
    nq̌ = get_num_of_free_coords(bot.structure)
    (;M) = build_mass_matrices(bot)
    ∂p0∂q̇0 = M
    ∂p0∂ẏ0 = ∂p0∂q̇0*∂q̇0∂ẏ0
    function inner_cost_and_gradient(F, G, x,opt_params=nothing)
        if mode isa InitialCost
            set_initial!(bot,params_func(x);idx)
        elseif mode isa StructureCost
            set_structure_params!(bot,params_func(x);idx)
        elseif mode isa ControlCost
            set_control_params!(bot,policy,params_func(x);idx)
        end
        
        if iscontact(env)
            dyn_prob = DynamicsProblem(bot;
                policy, env,
                contact_model=RestitutionFrictionCombined(
                    NewtonRestitution(),
                    CoulombFriction())
            )
            dyn_solver = DynamicsSolver(
                integrator,
                MonolithicContactSolver(
                    InteriorPointMethod()
                )
            )
        else
            dyn_prob = DynamicsProblem(bot; policy, env)
            dyn_solver = DynamicsSolver(
                integrator
            )
        end
        adj_solver = DiscreteAdjointDynamicsSolver(
            dyn_solver,
            objt,
        )
        adj_sen_solver = AdjointDynamicsSensitivitySolver(
            dyn_solver,
            adj_solver
        )
        prescribe! = (state, structure) -> nothing
        restart = true
        simulator = Simulator(dyn_prob,adj_sen_solver,prescribe!;tspan,dt,restart)

        # forward
        forward_cache = solve!(simulator,dyn_solver;dt,kwargs...)
        if G !== nothing
            # adjoint sensitivity 
            (;totalstep) = simulator
            sen_cache = generate_cache(dyn_prob,adj_sen_solver;dt,totalstep,kwargs...)
            solve!(simulator,forward_cache,sen_cache;dt,kwargs...)
        
            if mode isa InitialCost
                Gy0 = ∂q0∂y0'*sen_cache.∂J∂x₀ᵀ[   1: nq̌]
                Gẏ0 = ∂p0∂ẏ0'*sen_cache.∂J∂x₀ᵀ[nq̌+1:2nq̌]
                G .= vcat(Gy0,Gẏ0)
            elseif mode isa StructureCost
                G .= sum(sen_cache.∂J∂cᵀ)[idx]
            elseif mode isa ControlCost
                G .= sum(sen_cache.∂J∂θᵀ)[idx]
            end
        end
        if F !== nothing
            return cost!(bot,objt,integrator,dt;)
        end
        nothing
    end
end
