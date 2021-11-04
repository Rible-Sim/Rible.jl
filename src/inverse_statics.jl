function build_Ti(tgstruct::TensegrityStructure,i::Int)
    nbodycoords = get_nbodycoords(tgstruct)
    body2q = tgstruct.connectivity.body2q
    ncoords = tgstruct.ncoords
    build_Ti(nbodycoords,ncoords,body2q[i])
end

function build_Ti(nbodycoords,nq,q_index)
    Ti = spzeros(Int,nbodycoords,nq)
    for (body_q_id,q_id) in enumerate(q_index)
        Ti[body_q_id,q_id] = 1
    end
    Ti
end

function build_T(tgstruct)
    nbodycoords = get_nbodycoords(tgstruct)
    body2q = tgstruct.connectivity.body2q
    nmvbodies = tgstruct.nmvbodies
    ncoords = tgstruct.ncoords
    T = spzeros(Int,ncoords,nbodycoords*tgstruct.nmvbodies)
    for (mvrbid,rbid) in enumerate(tgstruct.mvbodyindex)
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
    Ci = hcat([
        transpose(Cpi)
        for Cpi in rb.state.cache.Cp
    ]...)
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

function build_Q̃(tgstruct)
    Q̃=build_T(tgstruct)*build_C(tgstruct)*build_D(tgstruct)
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


function build_L̂(tgstruct)
    @unpack nstrings, ndim, strings = tgstruct
    reset_forces!(tgstruct)
    update_strings_apply_forces!(tgstruct)
    T = get_numbertype(tgstruct)
    L̂ = spzeros(T, nstrings*ndim, nstrings)
    for (i,ss) in enumerate(strings)
        is = (i-1)*ndim
        L̂[is+1:is+ndim,i] = ss.state.direction
    end
    L̂
end

function build_L(tgstruct)
    @unpack nstrings, ndim, strings = tgstruct
    reset_forces!(tgstruct)
    update_strings_apply_forces!(tgstruct)
    T = get_numbertype(tgstruct)
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

function build_ℓ(tgstruct)
    reset_forces!(tgstruct)
    update_strings_apply_forces!(tgstruct)
    ℓ = [s.state.length for s in tgstruct.strings]
end

function build_K̂(tgstruct)
     ks = [s.k for s in tgstruct.strings]
     K̂ = build_L̂(tgstruct)*Diagonal(ks)
end

function build_W(tgstruct)
    q,_ = get_q(tgstruct)
    A = build_A(tgstruct)
    Aq = A(q)
    W = transpose(Aq)*inv(Aq*transpose(Aq))*Aq
end

# function build_K(tgstruct)
#     q,_ = get_q(tgstruct)
#     A = build_A(tgstruct)
#     Q̃ = build_Q̃(tgstruct)
#     K̂ = build_K̂(tgstruct)
#     # ℓ = build_ℓ(tgstruct)
#     K = hcat(
#         transpose(A(q)),
#         Q̃*K̂
#         )
# end

function build_G!(tgstruct;factor=1.0)
    reset_forces!(tgstruct)
    apply_gravity!(tgstruct)
    G = assemble_forces(tgstruct;factor=factor)
end

function build_Q̂(tgstruct)
    q,_ = get_q(tgstruct)
    A = build_A(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    L̂ = build_L̂(tgstruct)
    Q̂ = hcat(
        transpose(A(q)),
        -Q̃*L̂
    )
end
function build_RHS(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    K̂ = build_K̂(tgstruct)
    ℓ = build_ℓ(tgstruct)
    G = build_G!(tgstruct)
    RHS = Q̃*K̂*ℓ + G
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

function statics_equation_for_tension(tg_input,reftg;gravity=false,scale=true)
    tg = deepcopy(tg_input)
    refq0,_ = get_q(reftg)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,refq0)
    update_strings_apply_forces!(tg)

    L̂ = build_L̂(tg)
    Q̃ = build_Q̃(tg)
    A = build_A(tg)
    Aq = A(refq0)

    # Left hand side
    Q̃L̂ = Q̃*L̂
    if scale
        c = maximum(abs.(Q̃L̂))
    else
        c = one(eltype(Q̃L̂))
    end
    B = hcat(c*transpose(Aq),-Q̃L̂)

    # Right hand side
    F̃ = zero(refq0)
    if gravity
        G = build_G!(tg)
        F̃ .+= G
    end

    tg,B,F̃,c
end

function statics_equation_for_density(tg_input,reftg;gravity=false,scale=true)
    tg = deepcopy(tg_input)
    refq0,_ = get_q(reftg)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,refq0)
    update_strings_apply_forces!(tg)

    L = build_L(tg)
    Q̃ = build_Q̃(tg)
    A = build_A(tg)
    Aq = A(refq0)

    # Left hand side
    Q̃L = Q̃*L
    if scale
        c = maximum(abs.(Q̃L))
    else
        c = one(eltype(Q̃L))
    end
    B = hcat(c*transpose(Aq),-Q̃L)

    # Right hand side
    F̃ = zero(refq0)
    if gravity
        G = build_G!(tg)
        F̃ .+= G
    end

    tg,B,F̃,c
end

function statics_equation_for_restlength(tg_input,reftg;gravity=false,scale=true)
    tg = deepcopy(tg_input)
    refq0,_ = get_q(reftg)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,refq0)
    update_strings_apply_forces!(tg)

    K̂ = build_K̂(tg)
    ℓ = build_ℓ(tg)
    Q̃ = build_Q̃(tg)
    A = build_A(tg)
    Aq = A(refq0)

    # Left hand side
    Q̃K̂ = Q̃*K̂
    if scale
        c = maximum(abs.(Q̃K̂))
    else
        c = one(eltype(Q̃K̂))
    end
    B = hcat(c*transpose(Aq),Q̃K̂)

    # Right hand side
    F̃ = zero(refq0)
    F̃ .+= Q̃K̂*ℓ
    if gravity
        G = build_G!(tg)
        F̃ .+= G
    end

    tg,B,F̃,c
end

function statics_equation_for_stiffness(tg_input,reftg;gravity=false,scale=true)
    tg = deepcopy(tg_input)
    refq0,_ = get_q(reftg)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,refq0)
    update_strings_apply_forces!(tg)

    L̄ = build_L̄(tg)
    Q̃ = build_Q̃(tg)
    A = build_A(tg)
    Aq = A(refq0)

    # Left hand side
    Q̃L̄ = Q̃*L̄
    if scale
        c = maximum(abs.(Q̃L̄))
    else
        c = one(eltype(Q̃L̄))
    end
    B = hcat(c*transpose(Aq),-Q̃L̄)

    # Right hand side
    F̃ = zero(refq0)
    if gravity
        G = build_G!(tg)
        F̃ .+= G
    end

    tg,B,F̃,c
end

function get_solution_set(B,F̃)
    # A Particular solution
    x = pinv(B)*F̃
    # A basis for the null space of B
    nb = nullspace(B)
    return x,nb
end



function actuation_check(bot,Y,a)
    @unpack tg = bot
    Δu = Y*a
    original_restlens = get_original_restlen(bot)
    restlens = original_restlens + Δu
    if any((x)->x<0,restlens)
        @warn "Negative rest lengths"
    end
    lengths = get_strings_len(tg)
    if any((x)->x<=0,lengths-restlens)
        @warn "Nonpositive tension"
    end
end

function statics_equation_for_actuation(tg,reftg,Y;gravity=false,scale=true)
    acttg,B_old,F̃_old,c = statics_equation_for_restlength(tg,reftg;gravity,scale)
    nλ = get_nconstraint(acttg)
    nY1,nY2 = size(Y)
    IY = zeros(eltype(c),nλ+nY1,nλ+nY2)
    IY[1:nλ,1:nλ] .= Matrix(one(c)*I,nλ,nλ)
    IY[nλ+1:nλ+nY1,nλ+1:nλ+nY2] = Y
    B = B_old*IY
    u0 = [s.state.restlen for s in acttg.strings]
    F̃ = F̃_old - B_old*vcat(zeros(eltype(c),nλ),u0)
    acttg,B,F̃,c
end

function inverse_for_actuation(tg,reftg,Y;gravity=false,recheck=true,scale=true)
    # We only use the $q$ of the reference structure.
    acttg,B,F̃,c = statics_equation_for_actuation(tg,reftg,Y;gravity,scale)
    if rank(B) < minimum(size(B))
        @warn "LHS is singular: rank(B)=$(rank(B)) < $(minimum(size(B)))"
        @info "Using Moore-Penrose pseudoinverse"
        y0 = pinv(B)*F̃
    else
        y0 = B\F̃
    end
    nλ = get_nconstraint(acttg)
    _,na = size(Y)
    λ = y0[1:nλ].*c
    a = y0[nλ+1:nλ+na]
    actuation_check(acttg,Y,a)
    refq,_ = get_q(reftg)
    actuate!(acttg,a)
    check_static_equilibrium(acttg,refq,λ;gravity)
    u0 = [s.state.restlen for s in tg.strings]
    rl = u0 + Y*a
    λ,rl,a
end

function inverse(tr,reftr,Y;gravity=false,recheck=true,scale=true)
    # We only use the $q$ of the reference structure.
    acttg,lhs,rhs,c = statics_equation_for_actuation(tr.tg,reftr.tg,Y;gravity,scale)
    acttr = TensegrityRobot(acttg,deepcopy(tr.hub))
    if rank(lhs) < minimum(size(lhs))
        @warn "LHS is singular: rank(lhs)=$(rank(lhs)) < $(minimum(size(lhs)))"
        @info "Using Moore-Penrose pseudoinverse"
        x = pinv(lhs)*rhs
        @debug "Inv. Res. $(lhs*x-rhs)"
    else
        x = lhs\rhs
    end
    λ = x[1:acttg.nconstraint].*c
    a = x[acttg.nconstraint+1:end]
    actuation_check(acttr,Y,a)
    refq,_ = get_q(reftr.tg)
    actuate!(acttr,a)
    check_static_equilibrium(acttg,refq,λ;gravity)
    u0 = get_original_restlen(acttr)
    rl = u0 + Y*a
    λ,rl,a
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
        # actuation_check(actstruct,Y,a)
        # refq,_ = get_q(refstruct)
        # actuate!(actstruct,a)
        # check_static_equilibrium(actstruct,refq,λ;gravity)
        u0 = get_strings_restlen(tgstruct_input)
        rl = u0 + Y*a
        λ,rl,a
    end
end

function check_static_equilibrium(tgstruct_input,q,λ;gravity=false)
    tgstruct = deepcopy(tgstruct_input)
    reset_forces!(tgstruct)
    distribute_q_to_rbs!(tgstruct,q)
    update_strings_apply_forces!(tgstruct)
    if gravity
        apply_gravity!(tgstruct)
    end
    generalized_forces = assemble_forces(tgstruct)
    constraint_forces = transpose(build_A(tgstruct)(q))*λ
    static_equilibrium = constraint_forces ≈ generalized_forces
    @debug "Res. forces = $(generalized_forces-constraint_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(generalized_forces-constraint_forces))"
    end
    static_equilibrium
end
