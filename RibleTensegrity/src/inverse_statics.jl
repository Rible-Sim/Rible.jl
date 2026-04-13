

function build_Ci(body)
    Ci = reduce(hcat,[
        transpose(Cpi) for Cpi in body.state.cache.Cp
    ])
end

function build_Q(structure::AbstractStructure, num_of_full_coords, bodyid2sys_full_coords)
    T = get_numbertype(structure)
    num_of_dims = get_num_of_dims(structure)
    cables = get_cables(structure)
    ids, num_of_cables = get_ids(cables)
    ids_sorted = sort(ids)
    Q = zeros(T,num_of_full_coords,num_of_dims*num_of_cables)
    foreach(cables) do cable
        j = findfirst(x->x==cable.id,ids_sorted)
        (;hen,egg) = cable.joint.hen2egg
        C_hen = hen.body.cache.Cps[hen.pid]
        C_egg = egg.body.cache.Cps[egg.pid]
        sys_full_hen = bodyid2sys_full_coords[hen.body.prop.id]
        sys_full_egg = bodyid2sys_full_coords[egg.body.prop.id]
        js = (j-1)*num_of_dims
        Q[sys_full_egg,js+1:js+num_of_dims] .-= transpose(C_egg)
        Q[sys_full_hen,js+1:js+num_of_dims] .+= transpose(C_hen)
    end
    Q
end

function build_Q(structure::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = structure.connectivity
    build_Q(structure,num_of_full_coords,bodyid2sys_full_coords)
end

function build_Q(structure::AbstractStructure{bType,aType,<:PresFreeConnectivity}) where {bType, aType}
    (;num_of_full_coords,bodyid2sys_full_coords,sys_free_coords_idx) = structure.connectivity
    @view build_Q(structure,num_of_full_coords,bodyid2sys_full_coords)[sys_free_coords_idx,:]
end

function build_L̂(structure)
    (;connectivity,num_of_dim) = structure
    cables = get_cables(structure)
    ids, num_of_cables = get_ids(cables)
    ids_sorted = sort(ids)
    T = get_numbertype(structure)
    L̂ = spzeros(T, num_of_cables*num_of_dim, num_of_cables)
    foreach(cables) do cable
        j = findfirst(x->x==cable.id,ids_sorted)
        js = (j-1)*num_of_dim
        L̂[js+1:js+num_of_dim,j] = cable.force.state.direction
    end
    L̂
end

function build_K̂(structure)
    (;connectivity,num_of_dim) = structure
    (;apparatuses) = structure
    cables = get_cables(apparatuses)
    ncables = length(cables)
    T = get_numbertype(structure)
    K̂ = spzeros(T, ncables*num_of_dim, ncables)
    foreach(cables) do appar
        j = appar.id
        (;joint,force) = appar
        js = (j-1)*num_of_dim
        K̂[js+1:js+num_of_dim,j] = force.state.direction*force.k
    end
    K̂
end

function build_L(structure)
    (;connectivity,num_of_dim) = structure
    (;connected) = connectivity.
    (;cables) = structure.apparatuses
    ncables = length(cables)
    T = get_numbertype(structure)
    L = spzeros(T, ncables*num_of_dim, ncables)
    foreach(connected) do scnt
        j = scnt.id
        scable = cables[j]
        js = (j-1)*num_of_dim
        L[js+1:js+num_of_dim,j] = scable.state.direction*scable.state.length
    end
    L
end

function build_Γ(structure)
    function inner_Γ(q)
        clear_forces!(structure)
        update_bodies!(structure,q)
        update_cables_apply_forces!(structure)
        forces = [s.state.tension.*s.state.direction for s in structure.cables]
        reduce(vcat,forces)
    end
end

function build_G!(structure::AbstractStructure,field;)
    structure = deepcopy(structure)
    (;bodies) = structure
    (;members) = structure.state
    clear_forces!(structure)
    apply_field!(structure, field)
    foreach(bodies) do body
        (;prop,state,cache) = body
        (;F) = members[prop.id]
        mass_center_force!(F,state,cache)
    end
end

function build_G(structure::AbstractStructure,field)
    build_G!(structure,field)
    structure.state.system.F
end

function build_G(structure::AbstractStructure{bType,aType,<:PresFreeConnectivity},field) where{bType, aType}
    build_G!(structure,field)
    structure.state.system.F̌
end

function make_U(structure)
    (;num_of_dim) = structure
    (;num_of_full_coords) = structure.connectivity
    (;cables) = structure.apparatuses
    ncables = length(cables)
    function inner_U(s,u)
        ret = zeros(eltype(s),ncables*num_of_dim,num_of_full_coords)
        for i = 1:ncables
            k = cables[i].k
            Ji = Array(build_Ji(structure,i))
            ret[(i-1)*num_of_dim+1:i*num_of_dim,:] = k*Ji*(1-s[i]*u[i])
        end
        ret
    end
    function inner_U(s,u,k)
        ret = zeros(eltype(s),ncables*num_of_dim,num_of_full_coords)
        for i = 1:ncables
            Ji = Array(build_Ji(structure,i))
            ret[(i-1)*num_of_dim+1:i*num_of_dim,:] = k[i]*Ji*(1-s[i]*u[i])
        end
        ret
    end
    inner_U
end

function make_Q(structure,q0)
    (;num_of_dim) = structure
    cnt = structure.connectivity
    (;num_of_full_coords,sys_pres_coords_idx,sys_full_coords_idx,bodyid2sys_full_coords) = cnt
    (;connected) = 
    (;cables) = structure.apparatuses
    (;bodyid2sys_locus_id,sys_locus_id2coords_idx) = cnt
    function inner_Q̌(q̌,s,u)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_full_coords_idx] .= q̌
        ret = zeros(eltype(q̌),num_of_full_coords)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(connected) do scnt
            j = scnt.id
            (;k) = cables[j]
            rb1 = scnt.hen.body
            rb2 = scnt.egg.body
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Ūj = @view (transpose(Jj)*Jj)[sys_full_coords_idx,:]
            ret .+= k*(u[j]*s[j]-1)*Ūj*q
        end
        ret
    end
    function inner_Q̌(q̌,s,μ,k,c)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_full_coords_idx] .= q̌
        ret = zeros(eltype(q̌),num_of_full_coords)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.hen.body
            rb2 = scnt.egg.body
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
            c1 = c[sys_locus_id2coords_idx[bodyid2sys_locus_id[rb1id][ap1id]]]
            c2 = c[sys_locus_id2coords_idx[bodyid2sys_locus_id[rb2id][ap2id]]]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            # C1 = rb1.state.cache.Cps[ap1id]
            # C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Ūj = @view (transpose(Jj)*Jj)[sys_full_coords_idx,:]
            ret .+= k[j]*(μ[j]*s[j]-1)*Ūj*q
        end
        ret
    end
    # function inner_Q̌(q̌,γ,c)
    #     ret = zeros(eltype(γ),nfullcoords)
    #     foreach(string2ap) do scnt
    #         j = scnt.id
    #         rb1 = scnt.hen.body
    #         rb2 = scnt.egg.body
    #         rb1id = rb1.prop.id
    #         rb2id = rb2.prop.id
    #         ap1id = scnt.hen.pid
    #         ap2id = scnt.egg.pid
    #         is1 = (apnb[rb1id][ap1id]-1)*num_of_dim
    #         is2 = (apnb[rb2id][ap2id]-1)*num_of_dim
    #         c1 = c[is1+1:is1+num_of_dim]
    #         c2 = c[is2+1:is2+num_of_dim]
    #         C1 = rb1.state.cache.funcs.C(c1)
    #         C2 = rb2.state.cache.funcs.C(c2)
    #         T1 = build_Ti(structure,rb1id)
    #         T2 = build_Ti(structure,rb2id)
    #         Jj = C2*T2-C1*T1
    #         Uj = transpose(Jj)*Jj
    #         ret .+= γ[j]*Uj*q̌
    #     end
    #     ret
    # end
    # inner_Q̌
end


function check_inverse_sanity(B)
    n1,n2 = size(B)
    ret = false
    if n1>n2
        @warn "$n1(Eqns)>$n2(Unks). Inverse statics is generally unfeasible."
    elseif n1==n2
        rankB = rank(B)
        if rankB<n2
            @warn "$n1(Eqns)=$n2(Unks). But rank deficiency $(n2-rankB)."
        else #rankB==n2
            @info "$n1(Eqns)=$n2(Unks). Inverse statics is determinate."
            ret = true
        end
    else # n1<n2
        @info "$n1(Eqns)<$n2(Unks). Inverse statics is generally indeterminate by $(n2-n1)."
    end
    ret
end

function build_inverse_statics_core(structure_input,structure_ref::TensegrityStructure,field)
    structure = deepcopy(structure_input)
    clear_forces!(structure)
    update_bodies!(structure,structure_ref.state.system)
    update_apparatuses!(structure)
    G = build_G(structure, field,)
    build_inverse_statics_core(structure, structure_ref.state.system, G)
end

function build_inverse_statics_core(structure,inst_state,F)
    A = cstr_jacobian(structure,inst_state)
    Nᵀ = transpose(nullspace(A))
    Q̃ = build_Q(structure)
    structure,Nᵀ,Q̃,F
end

function build_inverse_statics_for_density(structure_input, structure_ref, field; scale=true)
    structure, Nᵀ, Q̃, F = build_inverse_statics_core(structure_input, structure_ref, field;)
    L = build_L(structure)
    # Left hand side
    Q̃L = Q̃*L

    B = -Nᵀ*Q̃L

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_stiffness(structure_input,structure_ref,field,;scale=true)
    structure, Nᵀ, Q̃, F = build_inverse_statics_core(structure_input, structure_ref, field;)

    L̂ = build_L̂(structure)
    ℓ = get_cables_len(structure)
    u = get_cables_restlen(structure)

    # Left hand side
    Q̃L̄ = Q̃*L̂*Diagonal(ℓ-u)

    B = -Nᵀ*Q̃L̄

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_restlength(structure_input, structure_ref, field; scale=true)
    structure, Nᵀ, Q̃, F = build_inverse_statics_core(structure_input, structure_ref, field;)
    ℓ = get_cables_len(structure)
    K̂ = build_K̂(structure)
    # Left hand side
    Q̃K̂ = Q̃*K̂

    B = Nᵀ*Q̃K̂

    # Right hand side
    F̃ = Nᵀ*(Q̃K̂*ℓ + F)
    B,F̃
end

function build_inverse_statics_for_force(structure_input, structure_ref, field; scale=true)
    structure, Nᵀ, Q̃, F = build_inverse_statics_core(structure_input, structure_ref, field; )

    L̂ = build_L̂(structure)

    # Left hand side
    Q̃L̂ = Q̃*L̂

    B = -Nᵀ*Q̃L̂

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_actuation(botinput,botref::Robot,field;Y=build_Y(botinput),scale=true)
    build_inverse_statics_for_actuation(botinput,botref.structure,field;Y,scale)
end

function build_inverse_statics_for_actuation(botinput,structure_ref::TensegrityStructure,field;Y=build_Y(botinput),scale=true)
    structure,Nᵀ,Q̃,F = build_inverse_statics_core(botinput.structure,structure_ref,field;)

    L̂ = build_L̂(structure)

    # Left hand side
    Q̃L̂Y = Q̃*L̂*Y

    B = -Nᵀ*Q̃L̂Y

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_actuation_force(botinput,botref::Robot,field;Y=build_Y(botinput),scale=true)
    build_inverse_statics_for_actuation_force(botinput,botref.structure,field;Y,scale)
end

function build_inverse_statics_for_actuation_force(botinput,structure_ref::TensegrityStructure,field;Y=build_Y(botinput),scale=true)
    structure,Nᵀ,Q̃,F = build_inverse_statics_core(botinput.structure,structure_ref,field;)
    L̂ = build_L̂(structure)
    # Left hand side
    Q̃L̂Y = Q̃*L̂*Y

    B = Nᵀ*Q̃L̂Y

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function get_solution_set(B,F̃)
    # A Particular solution
    x = pinv(B)*F̃
    # A basis for the null space of B
    nb = nullspace(B)
    return x,nb
end

function set_restlen!(structure,u)
    (;apparatuses,connectivity,state) = structure
    (;apparid2params_idx) = connectivity
    cables = get_cables(apparatuses)
    ids, num_of_cables = get_ids(cables)
    ids_sorted = sort(ids)
    foreach(cables) do cable
        (;force,id) = cable
        j = findfirst(==(id),ids_sorted)
        cable_params = @view state.system.c[apparid2params_idx[id]]
        cable_params .= (k = force.k, c = force.c, restlen = u[j]) |> collect
    end
    stretch!(structure)
end

function set_restlen!(bot::Robot, u)
    set_restlen!(bot.structure, u)
    bot.traj.c[begin] .= bot.structure.state.system.c
end


function check_static_equilibrium_output_multipliers(structure_input::TensegrityStructure, field::AbstractField = NoField();)
    q = get_coords(structure_input)
    check_static_equilibrium_output_multipliers(structure_input, field, q;)
end

function check_static_equilibrium_output_multipliers(structure_input::TensegrityStructure, field::AbstractField, q::AbstractVector;)
    structure = deepcopy(structure_input)
    check_static_equilibrium_output_multipliers!(structure, field, q;)
end

function gen_force!(inst_state, structure::AbstractStructure, field::AbstractField,)
    clear_forces!(structure)
    lazy_update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state)
    apply_field!(structure, field; )
    assemble_forces!(inst_state, structure)
end

function compute_multipliers(structure::AbstractStructure, )
    gen_forces = structure.state.system.F
    A = cstr_jacobian(structure,structure.state.system)
    λ = inv(A*transpose(A))*A*(-gen_forces)
    cstr_forces = transpose(A)*λ
    static_equilibrium = cstr_forces ≈ -gen_forces
    @debug "Res. forces = $(gen_forces+cstr_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(gen_forces+cstr_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium, λ
end

function compute_multipliers(structure::AbstractStructure{bType,aType,<:PresFreeConnectivity}, ) where {bType, aType}
    (;sys_free_coords_idx) = structure.connectivity
    gen_forces = @view structure.state.system.F[sys_free_coords_idx]
    A = cstr_jacobian(structure,structure.state.system)
    λ = inv(A*transpose(A))*A*(-gen_forces)
    cstr_forces = transpose(A)*λ
    static_equilibrium = cstr_forces ≈ -gen_forces
    @debug "Res. forces = $(gen_forces+cstr_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(gen_forces+cstr_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium, λ
end

function check_static_equilibrium_output_multipliers!(structure::TensegrityStructure, field::AbstractField, q::AbstractVector;)
    (;state) = structure
    state.system.q .= q
    gen_force!(state.system, structure, field)
    static_equilibrium, λ = compute_multipliers(structure)
end

function inverse_for_restlength(botinput, botref::Robot, field;
        scale=true,recheck=true,kwarg...
    )
    inverse_for_restlength(botinput.structure,botref.structure, field;scale,recheck,kwarg...)
end

function inverse_for_restlength(structure_input, structure_ref::TensegrityStructure, field;
        scale=true,recheck=true,
        eps_abs=1e-6,eps_rel=1e-6,verbose=false,
        fmin=0.0,
    )
    #todo dispatch on presfree connectivity type
    B,F̃ = build_inverse_statics_for_restlength(structure_input,structure_ref, field;scale)
    cables = get_cables(structure_input)
    nμ = cables |> length
    κ = get_cables_stiffness(structure_input)
    ℓ = get_cables_len(structure_ref)
    h = -ℓ.*κ
    if check_inverse_sanity(B)
        x0 = B\F̃
    else
        @info "Using Quadratic Programming."
        model = JuMP.Model(COSMO.Optimizer)
        JuMP.set_optimizer_attribute(model, "verbose", verbose)
        JuMP.set_optimizer_attribute(model, "eps_abs", eps_abs)
        JuMP.set_optimizer_attribute(model, "eps_rel", eps_rel)
        # JuMP.set_optimizer_attribute(model, "rho", 1e-5)
        # JuMP.set_optimizer_attribute(model, "adaptive_rho", false)
        # JuMP.set_optimizer_attribute(model, "alpha", 1.0)
        # JuMP.set_optimizer_attribute(model, "scaling", 0)
        JuMP.set_optimizer_attribute(model, "max_iter", 10_000)
        JuMP.@variable(model, x[1:nμ])
        # JuMP.@objective(model, Min, sum(x.^2))
        JuMP.@objective(model, Min, 
            0.5*transpose(x)*Diagonal(κ)*x + transpose(h)*x
        )
        _,xn = get_solution_set(B,F̃)
        @show size(xn)
        JuMP.@constraint(model, static, B*x .== F̃)
        # ϵ = 1e-3minimum(ℓ)
        JuMP.@constraint(model, positive_u, x .>= 0.0)
        JuMP.@constraint(model, positive_f, κ.*(ℓ .- x) .>= fmin)
        # JuMP.print(model)
        JuMP.optimize!(model)
        if JuMP.termination_status(model) == JuMP.OPTIMAL
            @info "Optimal Solution found."
        else
            @show JuMP.termination_status(model)
            error("Inverse statics optimization failed.")
        end
        # @show JuMP.primal_status(model)
        # @show JuMP.dual_status(model)
        # @show JuMP.objective_value(model)
        x0 = JuMP.value.(x)
        # @show abs.(B*x0 .- F̃)
        # x0 = pinv(B)*F̃
    end
    if recheck
        tgcheck = deepcopy(structure_input)
        q = get_coords(structure_ref)
        set_restlen!(tgcheck,x0)
        _, λ = check_static_equilibrium_output_multipliers!(tgcheck, field, q; )
    end
    x0
end

function inverse_for_multipliers(botinput::Robot, field, botref::Robot=botinput; scale=true, recheck=true)
    inverse_for_multipliers(botinput.structure,botref.structure, field;scale,recheck)
end

function inverse_for_multipliers(structure_input::TensegrityStructure, field, structure_ref::TensegrityStructure=structure_input; scale=true, recheck=true)
    q = get_coords(structure_ref)
    _,λ = check_static_equilibrium_output_multipliers(structure_input, field,q;)
    λ
end

function inverse_for_actuation(botinput,botref, field;Y=build_Y(botinput),
                    scale=true,recheck=true)
    structure_ref = botref.structure
    B,F̃ = build_inverse_statics_for_actuation(botinput,structure_ref, field;Y,scale)
    if check_inverse_sanity(B)
        y0 = B\F̃
    else
        @info "Using Moore-Penrose pseudoinverse"
        y0 = pinv(B)*F̃
    end
    na = size(Y,2)
    a = y0[1:na]
    if recheck
        botcheck = deepcopy(botinput)
        q = get_coords(structure_ref)
        actuate!(botcheck,a)
        _, λ = check_static_equilibrium_output_multipliers(botcheck.structure, field, q;)
    end
    λ,a
end

function optimize_for_stiffness_and_restlen(
        bot, field;
        fmin = 0.0,
        
    )

    B,F̃ = build_inverse_statics_for_actuation_force(
        bot, bot, field;
    )

    f0,Z = get_solution_set(B,F̃)
    @show size(Z)
    ncables = bot.structure.apparatuses.cables |> length
    Y = build_Y(bot)
    l0 = get_cables_len(bot.structure)

    # decision variables [k,x_ini,x]
    O_YZ = zero(Y*Z)
    O_Z = zero(Z)
    O_l = zeros(size(Z,1),ncables)
    A = [
        Diagonal(l0) -Y*Z ;
        O_l             Z ;
    ]
    b_ = [
        -Y*f0;
           f0 .- fmin;
    ]
    cstr1 = COSMO.Constraint(A,b_,COSMO.Nonnegatives)
    # P = Diagonal(vcat(ones(ncables),zeros(size(Z,2)),zeros(size(Z,2))))
    P = Diagonal(vcat(ones(ncables),zeros(size(Z,2))))
    q = zeros(size(P,2))

    model = COSMO.Model()
    custom_settings = COSMO.Settings(eps_abs = 1e-6)
    COSMO.assemble!(model, P, q, cstr1, settings = custom_settings)
    results = COSMO.optimize!(model)
    results.x
    results.obj_val
    if results.status != :Solved
        error("Status Code: $(results.status)")
    end
    nz = size(Z,2)
    k = results.x[1:ncables]
    z = results.x[ncables+1:end]

    f = f0+Z*z
    # extrema(f)

    μ = l0-inv(Diagonal(k))*Y*f
    l0 -inv(Diagonal(k))*Y*(Z*z+f0) |> minimum

    # minimum(μ)
    k,μ
end

