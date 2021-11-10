function build_Ti(tg::TensegrityStructure,i::Integer)
    nbodycoords = get_nbodycoords(tg)
    body2q = tg.connectivity.body2q
    ncoords = tg.ncoords
    build_Ti(nbodycoords,ncoords,body2q[i])
end

function build_Ti(nbodycoords,nq,q_index)
    Ti = spzeros(Int,nbodycoords,nq)
    for (body_q_id,q_id) in enumerate(q_index)
        Ti[body_q_id,q_id] = 1
    end
    Ti
end

function build_T(tg)
    nbodycoords = get_nbodycoords(tg)
    @unpack nmvbodies, ncoords, connectivity = tg
    body2q = connectivity.body2q
    T = spzeros(Int,ncoords,nbodycoords*tg.nmvbodies)
    for (mvrbid,rbid) in enumerate(tg.mvbodyindex)
        q_index = body2q[rbid]
        Ti = build_Ti(nbodycoords,ncoords,q_index)
        T[:,(mvrbid-1)*nbodycoords+1:mvrbid*nbodycoords] = transpose(Ti)
    end
    T
end
# T = build_T(manipulator)
# Array(T)
# @code_warntype build_T(manipulator)

function build_Ci(rb)
    Ci = reduce(hcat,[
        transpose(Cpi) for Cpi in rb.state.cache.Cp
    ])
end

function build_C(tgstruct)
    nbodycoords = get_nbodycoords(tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack nmvbodies,mvbodyindex,ndim = tgstruct
    block_size1 = repeat([nbodycoords],nmvbodies)
    block_size2 = [rb.prop.naps*ndim for rb in tgstruct.rigidbodies[mvbodyindex]]
    block_size1,block_size2
    T = get_numbertype(tgstruct)
    C = BlockArray{T}(undef,block_size1,block_size2)
    C .= 0
    for (mvrbid,rbid) in enumerate(mvbodyindex)
        C[Block(mvrbid,mvrbid)] = build_Ci(rbs[rbid])
    end
    C
end
# C = build_C(manipulator)
# Array(C)
# @code_warntype build_C(manipulator)
# manipulator.connectivity.string2ap[1]
function rbid2mvrbid(mvbodyindex)

end

function build_D(tgstruct)
    @unpack nstrings,nmvpoints,ndim,mvbodyindex = tgstruct
    @unpack body2q,string2ap = tgstruct.connectivity
    D = spzeros(Int,nmvpoints*ndim,nstrings*ndim)
    D_raw = spzeros(Int,nmvpoints,nstrings)
    iss = [0]
    for rbid in tgstruct.mvbodyindex
        rb = tgstruct.rigidbodies[rbid]
        push!(iss,iss[end]+rb.prop.naps)
    end

    for (sid,ap) in enumerate(string2ap)
        mvrbid1 = findfirst((x)->x==ap[1].rbid, mvbodyindex)
        if mvrbid1 != nothing
            D_raw[iss[mvrbid1]+ap[1].apid,sid] = 1
        end
        mvrbid2 = findfirst((x)->x==ap[2].rbid, mvbodyindex)
        if mvrbid2 != nothing
            D_raw[iss[mvrbid2]+ap[2].apid,sid] = -1
        end
    end
    D .= kron(D_raw,Matrix(1I,ndim,ndim))
end

function build_Q̃(tg)
    Q̃=build_T(tg)*build_C(tg)*build_D(tg)
end

function fvector(tgstruct)
    @unpack ndim, nstrings, strings = tgstruct
    ret = zeros(get_numbertype(tgstruct),ndim*nstrings)
    for (i,stri) in enumerate(strings)
        ret[(i-1)*ndim+1:i*ndim] = stri.state.tension*stri.state.direction
    end
    ret
end

function iksolve(prob;ftol=1e-14)
    @unpack funcs,q0,u0,λ0,nq,nu,nλ = prob
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
# D = build_D(manipulator)
# Array(D)

# Q̃=build_T(manipulator)*build_C(manipulator)*build_D(manipulator)
# Array(Q̃)


function build_L̂(tg)
    @unpack nstrings, ndim, strings = tg
    T = get_numbertype(tg)
    L̂ = spzeros(T, nstrings*ndim, nstrings)
    for (i,ss) in enumerate(strings)
        is = (i-1)*ndim
        L̂[is+1:is+ndim,i] = ss.state.direction
    end
    L̂
end

function build_L(tg)
    @unpack nstrings, ndim, strings = tg
    T = get_numbertype(tg)
    L = spzeros(T, nstrings*ndim, nstrings)
    for (i,ss) in enumerate(strings)
        is = (i-1)*ndim
        L[is+1:is+ndim,i] = ss.state.direction*ss.state.length
    end
    L
end

function build_Γ(tg)
    reset_forces!(tg)
    update_strings_apply_forces!(tg)
    forces = [s.state.tension.*s.state.direction for s in tg.strings]
    reduce(vcat,forces)
end

function build_K̂(tg)
     k = get_strings_stiffness(tg)
     K̂ = build_L̂(tg)*Diagonal(k)
end

function build_W(tgstruct)
    q,_ = get_q(tgstruct)
    A = build_A(tgstruct)
    Aq = A(q)
    W = transpose(Aq)*inv(Aq*transpose(Aq))*Aq
end

# function build_G(tginput;factor=1.0)
#     tg = deepcopy(tginput)
#     reset_forces!(tg)
#     apply_gravity!(tg)
#     G = assemble_forces(tg;factor=factor)
#     @show G
#     G
# end

function build_G(tg)
    ret = zeros(get_numbertype(tg),tg.ncoords)
    for (rbid,rb) in enumerate(tg.rigidbodies)
        @unpack prop, state = rb
        @unpack Cg = state.cache
        Tiᵀ = transpose(build_Ti(tg,rbid))
        Cgᵀ = transpose(Cg)
        f = prop.mass*get_gravity(rb)
        ret .+= Tiᵀ*Cgᵀ*f
    end
    ret
end

# Not ready
function build_Ji(tgstruct,i)
    rbs = tgstruct.rigidbodies
    cnt = tgstruct.connectivity
    @unpack string2ap = cnt
    ap = string2ap[i]
    C1 = rbs[ap[1].rbid].state.cache.Cp[ap[1].apid]
    C2 = rbs[ap[2].rbid].state.cache.Cp[ap[2].apid]
    T1 = build_Ti(tgstruct,ap[1].rbid)
    T2 = build_Ti(tgstruct,ap[2].rbid)
    Ji = C2*T2-C1*T1
end

function build_U(tgstruct)
    @unpack ncoords,nstrings,ndim,strings = tgstruct
    function inner_U(s,u)
        ret = zeros(eltype(s),nstrings*ndim,ncoords)
        for i = 1:nstrings
            k = strings[i].k
            Ji = Array(build_Ji(tgstruct,i))
            ret[(i-1)*ndim+1:i*ndim,:] = k*Ji*(1-s[i]*u[i])
        end
        ret
    end
    function inner_U(s,u,k)
        ret = zeros(eltype(s),nstrings*ndim,ncoords)
        for i = 1:nstrings
            Ji = Array(build_Ji(tgstruct,i))
            ret[(i-1)*ndim+1:i*ndim,:] = k[i]*Ji*(1-s[i]*u[i])
        end
        ret
    end
    inner_U
end

function build_S(tgstruct)
    @unpack nstrings = tgstruct
    function inner_S(q,s)
        ret = zeros(eltype(s),nstrings)
        for i = 1:nstrings
            Ji = build_Ji(tgstruct,i)
            Ui = Array(transpose(Ji)*Ji)
            ret[i] = transpose(q)*Ui*q*s[i]^2 - 1
        end
        ret
    end
end

function compensate_gravity_funcs(tgstruct)

    A = build_A(tgstruct)

    Q̃=build_Q̃(tgstruct)

    function F!(F,u)
        reset_forces!(tgstruct)
        actuate!(tgstruct,u)
        update_strings_apply_forces!(tgstruct)
        apply_gravity!(tgstruct)
        F .= 0
        #F .= Q̃*fvector(tgstruct)
        assemble_forces!(F,tgstruct)
    end

    A,F!
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

function build_inverse_statics_core(tginput,tgref::TensegrityStructure,Fˣ=nothing;gravity=false)
    q,_ = get_q(tgref)
    tg = deepcopy(tginput)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,q)
    update_strings_apply_forces!(tg)
    if gravity
        G = build_G(tg)
    else
        G = zero(q)
    end
    if isnothing(Fˣ)
        F = G
    else
        F = Fˣ + G
    end

    build_inverse_statics_core(tg,q,F)
end

function build_inverse_statics_core(tg,q::AbstractVector,F)
    A = build_A(tg)
    Aᵀ = transpose(A(q))
    Q̃ = build_Q̃(tg)
    tg,Aᵀ,Q̃,F
end

function build_inverse_statics_for_density(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Aᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    L = build_L(tg)
    # Left hand side
    Q̃L = Q̃*L
    if scale
        c = maximum(abs.(Q̃L))
    else
        c = one(eltype(Q̃L))
    end
    B = hcat(c*Aᵀ,-Q̃L)

    # Right hand side
    F̃ = F

    tg,B,F̃,c
end

function build_inverse_statics_for_stiffness(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Aᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)

    L̂ = build_L̂(tg)
    ℓ = get_strings_len(tg)
    u = get_strings_restlen(tg)

    # Left hand side
    Q̃L̄ = Q̃*L̂*Diagonal(ℓ-u)
    if scale
        c = maximum(abs.(Q̃L̄))
    else
        c = one(eltype(Q̃L̄))
    end
    B = hcat(c*Aᵀ,-Q̃L̄)

    # Right hand side
    F̃ = F

    tg,B,F̃,c
end

function build_inverse_statics_for_restlength(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Aᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    ℓ = get_strings_len(tg)
    K̂ = build_K̂(tg)
    # Left hand side
    Q̃K̂ = Q̃*K̂
    if scale
        c = maximum(abs.(Q̃K̂))
    else
        c = one(eltype(Q̃K̂))
    end
    B = hcat(c*Aᵀ,Q̃K̂)

    # Right hand side
    F̃ = Q̃K̂*ℓ + F

    tg,B,F̃,c
end

function build_inverse_statics_for_multipliers(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Aᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    ℓ = get_strings_len(tg)
    u = get_strings_restlen(tg)
    K̂ = build_K̂(tg)
    # Right hand side
    F̃ = Q̃*K̂*(ℓ-u) + F

    # Left hand side
    if scale
        c = maximum(abs.(F̃))
    else
        c = one(eltype(F̃))
    end
    B = c*Aᵀ

    tg,B,F̃,c
end

function build_inverse_statics_for_actuation(botinput,botref::TensegrityRobot,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    build_inverse_statics_for_actuation(botinput,botref.tg,Fˣ;Y,gravity,scale)
end

function build_inverse_statics_for_actuation(botinput,tgref::TensegrityStructure,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    tg,Aᵀ,Q̃,F = build_inverse_statics_core(botinput.tg,tgref,Fˣ;gravity)
    ℓ = get_strings_len(tg)
    K̂ = build_K̂(tg)
    u0 = get_original_restlen(botinput)
    # Left hand side
    Q̃K̂ = Q̃*K̂
    Q̃K̂Y = Q̃K̂*Y
    if scale
        c = maximum(abs.(Q̃K̂Y))
    else
        c = one(eltype(Q̃K̂Y))
    end
    B = hcat(c*Aᵀ,Q̃K̂Y)

    # Right hand side
    F̃ = Q̃K̂*(ℓ-u0) + F

    tg,B,F̃,c
end

function get_solution_set(B,F̃)
    # A Particular solution
    x = pinv(B)*F̃
    # A basis for the null space of B
    nb = nullspace(B)
    return x,nb
end

function check_restlen(tg,u)
    ℓ = get_strings_len(tg)
    if any((x)->x<0,u)
        @warn "Negative rest lengths"
    end
    if any((x)->x<=0,ℓ-u)
        @warn "Nonpositive tension"
    end
end

function check_actuation(bot,Y,a)
    @unpack tg = bot
    Δu = Y*a
    u0 = get_original_restlen(bot)
    u = u0 + Δu
    check_restlen(bot.tg,u)
    u
end

function get_inverse_func(tgstruct_input,refstruct,Y;gravity=false,recheck=true,scale=true)
    # We only use the $q$ of the reference structure.
    actstruct,lhs,rhs,c = statics_equation_for_actuation(tgstruct_input,refstruct,Y;gravity,scale)
    ncoeffs = maximum(size(lhs))- rank(lhs)
    xp,nb = get_solution_set(lhs,rhs)
    @info "Number of coefficients = $ncoeffs"
    function inner_inverse_func(ξ)
        x = xp + nb*ξ
        λ = x[1:actstruct.nconstraint].*c
        a = x[actstruct.nconstraint+1:end]
        # check_actuation(actstruct,Y,a)
        # refq,_ = get_q(refstruct)
        # actuate!(actstruct,a)
        # check_static_equilibrium(actstruct,refq,λ;gravity)
        u0 = get_strings_restlen(tgstruct_input)
        rl = u0 + Y*a
        λ,rl,a
    end
end

function check_static_equilibrium(tg_input,q,λ,F=nothing;gravity=false)
    tg = deepcopy(tg_input)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,q)
    update_strings_apply_forces!(tg)
    check_restlen(tg,get_strings_restlen(tg))
    if gravity
        apply_gravity!(tg)
    end
    generalized_forces = assemble_forces(tg)
    if !isnothing(F)
        generalized_forces .+= F[:]
    end
    constraint_forces = transpose(build_A(tg)(q))*λ
    static_equilibrium = constraint_forces ≈ generalized_forces
    @debug "Res. forces = $(generalized_forces-constraint_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(generalized_forces-constraint_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium
end

function inverse_for_restlength(botinput,botref::TensegrityRobot,Fˣ=nothing;gravity=false,scale=true,recheck=true)
    inverse_for_restlength(botinput.tg,botref.tg,Fˣ;gravity,scale,recheck)
end

function inverse_for_restlength(tginput,tgref::TensegrityStructure,Fˣ=nothing;gravity=false,scale=true,recheck=true)
    tg,B,F̃,c = build_inverse_statics_for_restlength(tginput,tgref,Fˣ;gravity,scale)
    nλ = get_nconstraint(tg)
    nu = tg.nstrings
    ny = nλ+nu
    if check_inverse_sanity(B)
        y0 = B\F̃
    else
        @info "Using Quadratic Programming."
        model = JuMP.Model(COSMO.Optimizer)
        JuMP.set_optimizer_attribute(model, "verbose", false)
        JuMP.set_optimizer_attribute(model, "eps_abs", 1e-15)
        JuMP.set_optimizer_attribute(model, "eps_rel", 1e-11)
        JuMP.@variable(model, y[1:ny])
        JuMP.@objective(model, Max, sum(y[nλ+1:nλ+nu].^2))
        JuMP.@constraint(model, static, B*y .== F̃)
        ϵ = 1e-10
        JuMP.@constraint(model, positive_u, y[nλ+1:nλ+nu].+ ϵ .>= 0)
        JuMP.@constraint(model, positive_f, y[nλ+1:nλ+nu].+ ϵ .<= get_strings_len(tg))
        # JuMP.print(model)
        JuMP.optimize!(model)
        if JuMP.termination_status(model) != JuMP.MathOptInterface.OPTIMAL
            error("Inverse statics optimization failed.")
        end
        # @show JuMP.primal_status(model)
        # @show JuMP.dual_status(model)
        # @show JuMP.objective_value(model)
        y0 = JuMP.value.(y)
        # @show abs.(B*y0 .- F̃)
        # y0 = pinv(B)*F̃
    end
    λ = y0[1:nλ].*c
    u = y0[nλ+1:nλ+nu]
    if recheck
        tgcheck = deepcopy(tginput)
        q,_ = get_q(tgref)
        set_restlen!(tgcheck,u)
        check_static_equilibrium(tgcheck,q,λ;gravity)
    end
    λ,u
end

function inverse_for_multipliers(botinput::TensegrityRobot,botref::TensegrityRobot=botinput,F=nothing;gravity=false,scale=true,recheck=true)
    inverse_for_multipliers(botinput.tg,botref.tg,F;gravity,scale,recheck)
end

function inverse_for_multipliers(tginput::TensegrityStructure,tgref::TensegrityStructure=tginput,F=nothing;gravity=false,scale=true,recheck=true)
    tg,B,F̃,c = build_inverse_statics_for_multipliers(tginput,tgref,F;gravity,scale)
    nλ = get_nconstraint(tg)
    y0 = B\F̃
    @debug "Res. : $(norm(B*y0 - F̃))"
    λ = y0.*c
    if recheck
        tgcheck = deepcopy(tginput)
        q,_ = get_q(tgref)
        check_static_equilibrium(tgcheck,q,λ,F;gravity)
    end
    λ
end

function inverse_for_actuation(botinput,botref,Fˣ=nothing;Y=build_Y(botinput),
                    gravity=false,scale=true,recheck=true)
    tgref = botref.tg
    tg,B,F̃,c = build_inverse_statics_for_actuation(botinput,tgref,Fˣ;Y,gravity,scale)
    if check_inverse_sanity(B)
        y0 = B\F̃
    else
        @info "Using Moore-Penrose pseudoinverse"
        y0 = pinv(B)*F̃
    end
    nλ = get_nconstraint(tg)
    na = size(Y,2)
    λ = y0[1:nλ].*c
    a = y0[nλ+1:nλ+na]
    u = check_actuation(botinput,Y,a)
    if recheck
        botcheck = deepcopy(botinput)
        q,_ = get_q(tgref)
        actuate!(botcheck,a)
        check_static_equilibrium(botcheck.tg,q,λ;gravity)
    end
    λ,u,a
end
