function build_T(st,i)
    (;indexed) = st.connectivity
    (;num_of_full_coords,bodyid2sys_full_coords) = indexed
    T = spzeros(Int,length(bodyid2sys_full_coords[i]),num_of_full_coords)
    # Ť = zeros(Int,length(bodyid2sys_full_coords[i]),num_of_free_coords)
    # T̃ = zeros(Int,length(bodyid2sys_full_coords[i]),num_of_pres_coords)
    for (i,j) in enumerate(bodyid2sys_full_coords[i])
        T[i,j] = 1
    end
    T
end

function build_Ci(body)
    Ci = reduce(hcat,[
        transpose(Cpi) for Cpi in body.state.cache.Cp
    ])
end

function build_Q̃(structure)
    (;indexed) = structure.connectivity
    (;num_of_free_coords,bodyid2sys_free_coords,num_of_force_apparatuses) = indexed
    T = get_numbertype(structure)
    num_of_dim = get_num_of_dims(structure)
    Q̃ = zeros(T,num_of_free_coords,num_of_dim*num_of_force_apparatuses)
    foreach(structure.apparatuses) do appar
        j = appar.id
        if appar.joint isa CableJoint
            (;hen,egg) = appar.joint
            C1 = hen.bodysig.cache.Cps[hen.pid]
            C2 = egg.bodysig.cache.Cps[egg.pid]
            uci1 = hen.bodysig.coords.free_idx
            uci2 = egg.bodysig.coords.free_idx
            m2sf1 = bodyid2sys_free_coords[hen.bodysig.prop.id]
            m2sf2 = bodyid2sys_free_coords[egg.bodysig.prop.id]
            Q̃[m2sf2,(j-1)*num_of_dim+1:j*num_of_dim] .-= transpose(C2)[uci2,:]
            Q̃[m2sf1,(j-1)*num_of_dim+1:j*num_of_dim] .+= transpose(C1)[uci1,:]
        end
    end
    Q̃
end


function iksolve(prob;ftol=1e-14)
    (;funcs,q0,u0,λ0,nq,nu,nλ) = prob
    A,F! = funcs
    F = zeros(nq)
    function IK_R!(R,x)
        u = @view x[   1:nu]
        λ = @view x[nu+1:nu+nλ]
        F!(F,u)
        R .= transpose(A(q0))*λ - F
    end
    initial_x = vcat(u0,λ0)
    ik_result = nlsolve(IK_R!,initial_x,ftol=ftol)
    @info "Convergence: $(converged(ik_result)), Iterations: $(ik_result.iterations), ftol=$ftol"
    u_result = ik_result.zero[   1:nu]
    λ_result = ik_result.zero[nu+1:nu+nλ]
    u_result,λ_result
end



function build_L̂(st)
    (;connectivity,num_of_dim) = st
    (;connected) = connectivity.tensioned
    (;cables) = st.apparatuses
    ncables = length(cables)
    T = get_numbertype(st)
    L̂ = spzeros(T, ncables*num_of_dim, ncables)
    foreach(connected) do scnt
        j = scnt.id
        scable = cables[j]
        js = (j-1)*num_of_dim
        L̂[js+1:js+num_of_dim,j] = scable.state.direction
    end
    L̂
end

function build_K̂(st)
    (;connectivity,num_of_dim) = st
    (;connected) = connectivity.tensioned
    (;cables) = st.apparatuses
    ncables = length(cables)
    T = get_numbertype(st)
    K̂ = spzeros(T, ncables*num_of_dim, ncables)
    foreach(connected) do scnt
        j = scnt.id
        scable = cables[j]
        js = (j-1)*num_of_dim
        K̂[js+1:js+num_of_dim,j] = scable.state.direction*scable.k
    end
    K̂
end

function build_L(st)
    (;connectivity,num_of_dim) = st
    (;connected) = connectivity.tensioned
    (;cables) = st.apparatuses
    ncables = length(cables)
    T = get_numbertype(st)
    L = spzeros(T, ncables*num_of_dim, ncables)
    foreach(connected) do scnt
        j = scnt.id
        scable = cables[j]
        js = (j-1)*num_of_dim
        L[js+1:js+num_of_dim,j] = scable.state.direction*scable.state.length
    end
    L
end

function build_Γ(st)
    function inner_Γ(q)
        clear_forces!(st)
        update_bodies!(st,q)
        update_cables_apply_forces!(st)
        forces = [s.state.tension.*s.state.direction for s in st.cables]
        reduce(vcat,forces)
    end
end

function build_Ǧ(tginput;factor=1.0)
    st = deepcopy(tginput)
    clear_forces!(st)
    apply_gravity!(st;factor)
    Ǧ = assemble_forces!(st)
end

function make_U(st)
    (;num_of_dim) = st
    (;num_of_full_coords) = st.connectivity.indexed
    (;cables) = st.apparatuses
    ncables = length(cables)
    function inner_U(s,u)
        ret = zeros(eltype(s),ncables*num_of_dim,num_of_full_coords)
        for i = 1:ncables
            k = cables[i].k
            Ji = Array(build_Ji(st,i))
            ret[(i-1)*num_of_dim+1:i*num_of_dim,:] = k*Ji*(1-s[i]*u[i])
        end
        ret
    end
    function inner_U(s,u,k)
        ret = zeros(eltype(s),ncables*num_of_dim,num_of_full_coords)
        for i = 1:ncables
            Ji = Array(build_Ji(st,i))
            ret[(i-1)*num_of_dim+1:i*num_of_dim,:] = k[i]*Ji*(1-s[i]*u[i])
        end
        ret
    end
    inner_U
end

function make_Q̌(st,q0)
    (;num_of_dim) = st
    (;numbered,indexed,tensioned) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_pres_coords_idx,sys_free_coords_idx,bodyid2sys_full_coords) = indexed
    (;connected) = tensioned
    (;cables) = st.apparatuses
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = numbered
    function inner_Q̌(q̌,s,u)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        ret = zeros(eltype(q̌),num_of_free_coords)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(connected) do scnt
            j = scnt.id
            (;k) = cables[j]
            rb1 = scnt.hen.bodysig
            rb2 = scnt.egg.bodysig
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
            Ūj = @view (transpose(Jj)*Jj)[sys_free_coords_idx,:]
            ret .+= k*(u[j]*s[j]-1)*Ūj*q
        end
        ret
    end
    function inner_Q̌(q̌,s,μ,k,c)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        ret = zeros(eltype(q̌),num_of_free_coords)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.hen.bodysig
            rb2 = scnt.egg.bodysig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
            c1 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb1id][ap1id]]]
            c2 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb2id][ap2id]]]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            # C1 = rb1.state.cache.Cps[ap1id]
            # C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Ūj = @view (transpose(Jj)*Jj)[sys_free_coords_idx,:]
            ret .+= k[j]*(μ[j]*s[j]-1)*Ūj*q
        end
        ret
    end
    # function inner_Q̌(q̌,γ,c)
    #     ret = zeros(eltype(γ),nfullcoords)
    #     foreach(string2ap) do scnt
    #         j = scnt.id
    #         rb1 = scnt.hen.bodysig
    #         rb2 = scnt.egg.bodysig
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
    #         T1 = build_Ti(st,rb1id)
    #         T2 = build_Ti(st,rb2id)
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

function build_inverse_statics_core(tginput,tgref::Structure,Fˣ=nothing;gravity=false)
    q = get_coords(tgref)
    q̌ = get_free_coords(tgref)
    st = deepcopy(tginput)
    clear_forces!(st)
    update_bodies!(st,q)
    update_apparatuses!(st)
    if gravity
        Ǧ = build_Ǧ(st)
    else
        Ǧ = zero(q̌)
    end
    if isnothing(Fˣ)
        F̌ = Ǧ
    else
        F̌ = F̌ˣ + Ǧ
    end

    build_inverse_statics_core(st,q,F̌)
end

function build_inverse_statics_core(st,q::AbstractVector,F)
    A = make_cstr_jacobian(st)
    Nᵀ = transpose(nullspace(A(q)))
    Q̃ = build_Q̃(st)
    st,Nᵀ,Q̃,F
end

function build_inverse_statics_for_density(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    st,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    L = build_L(st)
    # Left hand side
    Q̃L = Q̃*L

    B = -Nᵀ*Q̃L

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_stiffness(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    st,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)

    L̂ = build_L̂(st)
    ℓ = get_cables_len(st)
    u = get_cables_restlen(st)

    # Left hand side
    Q̃L̄ = Q̃*L̂*Diagonal(ℓ-u)

    B = -Nᵀ*Q̃L̄

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_restlength(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    st,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    ℓ = get_cables_len(st)
    K̂ = build_K̂(st)
    # Left hand side
    Q̃K̂ = Q̃*K̂

    B = Nᵀ*Q̃K̂

    # Right hand side
    F̃ = Nᵀ*(Q̃K̂*ℓ + F)
    B,F̃
end

function build_inverse_statics_for_force(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    st,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)

    L̂ = build_L̂(st)

    # Left hand side
    Q̃L̂ = Q̃*L̂

    B = -Nᵀ*Q̃L̂

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_actuation(botinput,botref::Robot,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    build_inverse_statics_for_actuation(botinput,botref.st,Fˣ;Y,gravity,scale)
end

function build_inverse_statics_for_actuation(botinput,tgref::Structure,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    st,Nᵀ,Q̃,F = build_inverse_statics_core(botinput.st,tgref,Fˣ;gravity)

    L̂ = build_L̂(st)

    # Left hand side
    Q̃L̂Y = Q̃*L̂*Y

    B = -Nᵀ*Q̃L̂Y

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_actuation_force(botinput,botref::Robot,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    build_inverse_statics_for_actuation_force(botinput,botref.st,Fˣ;Y,gravity,scale)
end

function build_inverse_statics_for_actuation_force(botinput,tgref::Structure,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    st,Nᵀ,Q̃,F = build_inverse_statics_core(botinput.st,tgref,Fˣ;gravity)
    L̂ = build_L̂(st)
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

function check_restlen(st)
    u = get_cables_restlen(st)
    check_restlen(st,u)
end

function check_restlen(st,u)
    ℓ = get_cables_len(st)
    if any((x)->x<0,u)
        @warn "Negative rest lengths"
    end
    if any((x)->x<0,ℓ-u)
        @warn "Nonpositive tension"
    end
end

function check_actuation(bot,Y,a)
    (;st) = bot
    Δu = Y*a
    u0 = get_original_restlen(bot)
    u = u0 + Δu
    check_restlen(bot.structure,u)
    u
end

function get_inverse_func(tg_input,reftg,Y;gravity=false,recheck=true,scale=true)
    # We only use the $q$ of the reference tgure.
    acttg,lhs,rhs,c = statics_equation_for_actuation(tg_input,reftg,Y;gravity,scale)
    ncoeffs = maximum(size(lhs))- rank(lhs)
    xp,nb = get_solution_set(lhs,rhs)
    @info "Number of coefficients = $ncoeffs"
    function inner_inverse_func(ξ)
        x = xp + nb*ξ
        λ = x[1:acttg.num_of_cstr].*c
        a = x[acttg.num_of_cstr+1:end]
        # check_actuation(acttg,Y,a)
        # refq = get_coords(reftg)
        # actuate!(acttg,a)
        # check_static_equilibrium(acttg,refq,λ;gravity)
        u0 = get_cables_restlen(tg_input)
        rl = u0 + Y*a
        λ,rl,a
    end
end

function check_static_equilibrium(tg_input,q,λ,F=nothing;gravity=false)
    st = deepcopy(tg_input)
    clear_forces!(st)
    distribute_q_to_rbs!(st,q)
    update_cables_apply_forces!(st)
    check_restlen(st,get_cables_restlen(st))
    if gravity
        apply_gravity!(st)
    end
    generalized_forces = assemble_forces(st)
    if !isnothing(F)
        generalized_forces .+= F[:]
    end
    cstr_forces = transpose(make_cstr_jacobian(st)(q))*λ
    static_equilibrium = cstr_forces ≈ generalized_forces
    @debug "Res. forces = $(generalized_forces-cstr_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(generalized_forces-cstr_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium
end

function check_static_equilibrium_output_multipliers(tg_input;F=nothing,gravity=false)
    q = get_coords(tg_input)
    check_static_equilibrium_output_multipliers(tg_input,q,F;gravity)
end

function check_static_equilibrium_output_multipliers(tg_input,q,F=nothing;gravity=false)
    st = deepcopy(tg_input)
    check_static_equilibrium_output_multipliers!(st,q,F;gravity)
end

function check_static_equilibrium_output_multipliers!(st,q,F=nothing;
        gravity=false,
        # stpt = nothing
    )
    clear_forces!(st)
    update_bodies!(st)
    update_apparatuses!(st)
    # check_restlen(st,get_cables_restlen(st))
    if gravity
        apply_gravity!(st)
    end
    generalized_forces = assemble_forces!(st)
    if !isnothing(F)
        generalized_forces .+= F[:]
    end
    q = get_coords(st)
    c = get_local_coords(st)
    q̌ = get_free_coords(st)
    A = make_cstr_jacobian(st,q)(q̌,c)
    # # @show A
    # s = get_s(st)
    # u = get_cables_restlen(st)
    # k = get_cables_stiffness(st)
    # if !(stpt isa Nothing) 
    #     @show stpt.q - q |> norm
    #     @show stpt.q̌ - q̌ |> norm
    #     @show stpt.s - s |> norm
    #     @show stpt.u - u |> norm
    #     @show stpt.k - k |> norm
    #     @show stpt.c - c |> norm
    # end

    # ǧeneralized_forces = make_Q̌(st,q)(q̌,s,u,k,c)
    # # @show s 
    # @show generalized_forces - ǧeneralized_forces
    λ = inv(A*transpose(A))*A*(-generalized_forces)
    cstr_forces = transpose(A)*λ
    static_equilibrium = cstr_forces ≈ -generalized_forces
    @debug "Res. forces = $(generalized_forces+cstr_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(generalized_forces+cstr_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium, λ
end

function inverse_for_restlength(botinput,botref::Robot,Fˣ=nothing;
        gravity=false,scale=true,recheck=true,kwarg...
    )
    inverse_for_restlength(botinput.st,botref.st,Fˣ;gravity,scale,recheck,kwarg...)
end

function inverse_for_restlength(tginput,tgref::Structure,Fˣ=nothing;
        gravity=false,scale=true,recheck=true,
        eps_abs=1e-6,eps_rel=1e-6,verbose=false,
        fmin=0.0,
    )
    B,F̃ = build_inverse_statics_for_restlength(tginput,tgref,Fˣ;gravity,scale)
    nμ = tgref.apparatuses.cables |> length
    κ = get_cables_stiffness(tginput)
    ℓ = get_cables_len(tgref)
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
        tgcheck = deepcopy(tginput)
        q = get_coords(tgref)
        set_restlen!(tgcheck,x0)
        _,λ = check_static_equilibrium_output_multipliers!(tgcheck,q;gravity)
    end
    x0
end

function inverse_for_multipliers(botinput::Robot,botref::Robot=botinput,F=nothing;gravity=false,scale=true,recheck=true)
    inverse_for_multipliers(botinput.st,botref.st,F;gravity,scale,recheck)
end

function inverse_for_multipliers(tginput::Structure,tgref::Structure=tginput,F=nothing;gravity=false,scale=true,recheck=true)
    q = get_coords(tgref)
    _,λ = check_static_equilibrium_output_multipliers(tginput,q;gravity)
    λ
end

function inverse_for_actuation(botinput,botref,Fˣ=nothing;Y=build_Y(botinput),
                    gravity=false,scale=true,recheck=true)
    tgref = botref.st
    B,F̃ = build_inverse_statics_for_actuation(botinput,tgref,Fˣ;Y,gravity,scale)
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
        q = get_coords(tgref)
        actuate!(botcheck,a)
        _,λ = check_static_equilibrium_output_multipliers(botcheck.st,q;gravity)
    end
    λ,a
end

function optimize_for_stiffness_and_restlen(
        bot;
        fmin = 0.0,
        gravity=true
    )

    B,F̃ = build_inverse_statics_for_actuation_force(
        bot,bot;gravity
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

