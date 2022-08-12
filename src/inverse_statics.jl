function build_T(tg,i)
    (;indexed) = tg.connectivity
    (;nfull,mem2sysfull) = indexed
    T = zeros(Int,length(mem2sysfull[i]),nfull)
    # Ť = zeros(Int,length(mem2sysfull[i]),nfree)
    # T̃ = zeros(Int,length(mem2sysfull[i]),npres)
    for (i,j) in enumerate(mem2sysfull[i])
        T[i,j] = 1
    end
    T
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

function build_D(tg)
    @unpack ncables,nmvpoints,ndim,mvbodyindex = tg
    @unpack body2q,string2ap = tg.connectivity
    D = spzeros(Int,nmvpoints*ndim,ncables*ndim)
    D_raw = spzeros(Int,nmvpoints,ncables)
    iss = [0]

    foreach(tg.rigidbodies) do rb
        if rb.prop.id in tg.mvbodyindex
            push!(iss,iss[end]+rb.prop.naps)
        end
    end

    foreach(string2ap) do scnt
        id = scnt.id
        mvrbid1 = findfirst((x)->x==scnt.end1.rbsig.prop.id, mvbodyindex)
        if mvrbid1 != nothing
            D_raw[iss[mvrbid1]+scnt.end1.pid,id] = 1
        end
        mvrbid2 = findfirst((x)->x==scnt.end2.rbsig.prop.id, mvbodyindex)
        if mvrbid2 != nothing
            D_raw[iss[mvrbid2]+scnt.end2.pid,id] = -1
        end
    end
    D .= kron(D_raw,Matrix(1I,ndim,ndim))
end

function build_Q̃(tg)
    (;tensioned,indexed) = tg.connectivity
    (;connected) = tensioned
    (;cables) = tg.tensiles
    ncables = length(cables)
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    Q̃ = zeros(T,nfree,ndim*ncables)

    foreach(connected) do cc
        j = cc.id
        cable = cables[j]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        uci1 = rb1.state.cache.unconstrained_index
        uci2 = rb2.state.cache.unconstrained_index
        m2sf1 = mem2sysfree[rb1.prop.id]
        m2sf2 = mem2sysfree[rb2.prop.id]
        Q̃[m2sf2,(j-1)*ndim+1:j*ndim] .-= transpose(C2)[uci2,:]
        Q̃[m2sf1,(j-1)*ndim+1:j*ndim] .+= transpose(C1)[uci1,:]
    end
    Q̃
end

function fvector(tgstruct)
    @unpack ndim, ncables, cables = tgstruct
    ret = zeros(get_numbertype(tgstruct),ndim*ncables)
    for (i,stri) in enumerate(cables)
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



function build_L̂(tg)
    (;ndim) = tg
    (;cables) = tg.tensiles
    ncables = length(cables)
    T = get_numbertype(tg)
    L̂ = spzeros(T, ncables*ndim, ncables)
    for (i,ss) in enumerate(cables)
        is = (i-1)*ndim
        L̂[is+1:is+ndim,i] = ss.state.direction
    end
    L̂
end

function build_L(tg)
    (;ncables, ndim, cables) = tg
    T = get_numbertype(tg)
    L = spzeros(T, ncables*ndim, ncables)
    for (i,ss) in enumerate(cables)
        is = (i-1)*ndim
        L[is+1:is+ndim,i] = ss.state.direction*ss.state.length
    end
    L
end

function build_Γ(tg)
    function inner_Γ(q)
        clear_forces!(tg)
        update_rigids!(tg,q)
        update_cables_apply_forces!(tg)
        forces = [s.state.tension.*s.state.direction for s in tg.cables]
        reduce(vcat,forces)
    end
end

function build_K̂(tg)
     k = get_cables_stiffness(tg)
     K̂ = build_L̂(tg)*Diagonal(k)
end

function build_W(tgstruct)
    q = get_q(tgstruct)
    A = make_A(tgstruct)
    Aq = A(q)
    W = transpose(Aq)*inv(Aq*transpose(Aq))*Aq
end

function build_Ǧ(tginput;factor=1.0)
    tg = deepcopy(tginput)
    clear_forces!(tg)
    apply_gravity!(tg;factor)
    Ǧ = generate_forces!(tg)
end

# Not ready
# function build_Ji(tg::AbstractTensegrityStructure,i)
#     rbs = tg.rigidbodies
#     cnt = tg.connectivity
#     (;tensioned) = cnt.connectivity
#     ap = tensioned[1][i]
#     C1 = rbs[ap[1].rbid].state.cache.Cp[ap[1].apid]
#     C2 = rbs[ap[2].rbid].state.cache.Cp[ap[2].apid]
#     T1 = build_Ti(tg,ap[1].rbid)
#     T2 = build_Ti(tg,ap[2].rbid)
#     Ji = C2*T2-C1*T1
# end

function make_U(tg)
    (;ndim) = tg
    (;nfull) = tg.connectivity.indexed
    (;cables) = tg.tensiles
    ncables = length(cables)
    function inner_U(s,u)
        ret = zeros(eltype(s),ncables*ndim,nfull)
        for i = 1:ncables
            k = cables[i].k
            Ji = Array(build_Ji(tg,i))
            ret[(i-1)*ndim+1:i*ndim,:] = k*Ji*(1-s[i]*u[i])
        end
        ret
    end
    function inner_U(s,u,k)
        ret = zeros(eltype(s),ncables*ndim,nfull)
        for i = 1:ncables
            Ji = Array(build_Ji(tg,i))
            ret[(i-1)*ndim+1:i*ndim,:] = k[i]*Ji*(1-s[i]*u[i])
        end
        ret
    end
    inner_U
end

function make_Q̌(tg,q0)
    (;ndim) = tg
    (;numbered,indexed,tensioned) = tg.connectivity
    (;nfull,nfree,syspres,sysfree,mem2sysfull) = indexed
    (;connected) = tensioned
    (;cables) = tg.tensiles
    (;mem2num,num2sys) = numbered
    function inner_Q̌(q̌,s,u)
		q = Vector{eltype(q̌)}(undef,nfull)
		q[syspres] .= q0[syspres]
		q[sysfree] .= q̌
        ret = zeros(eltype(q̌),nfree)
        Jj = zeros(eltype(q̌),ndim,nfull)
        foreach(connected) do scnt
            j = scnt.id
            (;k) = cables[j]
            rb1 = scnt.end1.rbsig
            rb2 = scnt.end2.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.end1.pid
            ap2id = scnt.end2.pid
            # c1 = c[num2sys[mem2num[rb1id][ap1id]]]
            # c2 = c[num2sys[mem2num[rb2id][ap2id]]]
            # C1 = rb1.state.cache.funcs.C(c1)
            # C2 = rb2.state.cache.funcs.C(c2)
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = mem2sysfull[rb1.prop.id]
            mfull2 = mem2sysfull[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Ūj = @view (transpose(Jj)*Jj)[sysfree,:]
            # k = k[j]
            ret .+= k*(u[j]*s[j]-1)*Ūj*q
        end
        ret
    end
    # function inner_Q̌(q̌,γ,c)
    #     ret = zeros(eltype(γ),nfullcoords)
    #     foreach(string2ap) do scnt
    #         j = scnt.id
    #         rb1 = scnt.end1.rbsig
    #         rb2 = scnt.end2.rbsig
    #         rb1id = rb1.prop.id
    #         rb2id = rb2.prop.id
    #         ap1id = scnt.end1.pid
    #         ap2id = scnt.end2.pid
    #         is1 = (apnb[rb1id][ap1id]-1)*ndim
    #         is2 = (apnb[rb2id][ap2id]-1)*ndim
    #         c1 = c[is1+1:is1+ndim]
    #         c2 = c[is2+1:is2+ndim]
    #         C1 = rb1.state.cache.funcs.C(c1)
    #         C2 = rb2.state.cache.funcs.C(c2)
    #         T1 = build_Ti(tg,rb1id)
    #         T2 = build_Ti(tg,rb2id)
    #         Jj = C2*T2-C1*T1
    #         Uj = transpose(Jj)*Jj
    #         ret .+= γ[j]*Uj*q̌
    #     end
    #     ret
    # end
    # inner_Q̌
end

function build_KE(tg)
    (;ncables,ndim,cables,connectivity) = tg
    (;numbered,indexed,tensioned) = connectivity
    (;nfull,mem2sysfull) = indexed
    function inner_KE(q,s,k,c)
        ret = zeros(eltype(s),nfull,nfull)
        foreach(tensioned) do scnt
            j = scnt.id
            rb1 = scnt.end1.rbsig
            rb2 = scnt.end2.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.end1.pid
            ap2id = scnt.end2.pid
            is1 = (apnb[rb1id][ap1id]-1)*ndim
            is2 = (apnb[rb2id][ap2id]-1)*ndim
            c1 = c[is1+1:is1+ndim]
            c2 = c[is2+1:is2+ndim]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            T1 = build_Ti(tg,rb1id)
            T2 = build_Ti(tg,rb2id)
            Jj = C2*T2-C1*T1
            Uj = transpose(Jj)*Jj
            Ujq = Uj*q
            ret .+= k[j]*s[j]^2*(Ujq*transpose(Ujq))
        end
        ret
    end
end

function build_KG(tg)
    @unpack nfullcoords,ncables,ndim,cables = tg
    cnt = tg.connectivity
    @unpack string2ap,apnb = cnt
    function inner_KG(q,s,μ,k,c)
        ret = zeros(eltype(s),nfullcoords,nfullcoords)
        foreach(string2ap) do scnt
            j = scnt.id
            rb1 = scnt.end1.rbsig
            rb2 = scnt.end2.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.end1.pid
            ap2id = scnt.end2.pid
            is1 = (apnb[rb1id][ap1id]-1)*ndim
            is2 = (apnb[rb2id][ap2id]-1)*ndim
            c1 = c[is1+1:is1+ndim]
            c2 = c[is2+1:is2+ndim]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            T1 = build_Ti(tg,rb1id)
            T2 = build_Ti(tg,rb2id)
            Jj = C2*T2-C1*T1
            Uj = transpose(Jj)*Jj
            Ujq = Uj*q
            ret .+= k[j]*(1-μ[j]*s[j])*(Uj-s[j]^2*Ujq*transpose(Ujq))
        end
        ret
    end
end

function make_S(tg,q0)
    (;ndim) = tg
    (;numbered,indexed,tensioned) = tg.connectivity
    (;syspres,sysfree,nfull,mem2sysfull) = indexed
    (;mem2num,num2sys) = numbered
    (;connected) = tensioned
    (;cables) = tg.tensiles
    ncables = length(cables)
    function inner_S(q̌,s)
		q = Vector{eltype(q̌)}(undef,nfull)
		q[syspres] .= q0[syspres]
		q[sysfree] .= q̌
        ret = zeros(eltype(q̌),ncables)
        Jj = zeros(eltype(q̌),ndim,nfull)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.end1.rbsig
            rb2 = scnt.end2.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.end1.pid
            ap2id = scnt.end2.pid
            # c1 = c[num2sys[mem2num[rb1id][ap1id]]]
            # c2 = c[num2sys[mem2num[rb2id][ap2id]]]
            # C1 = rb1.state.cache.funcs.C(c1)
            # C2 = rb2.state.cache.funcs.C(c2)
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = mem2sysfull[rb1.prop.id]
            mfull2 = mem2sysfull[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            ret[j] = transpose(q)*Uj*q*s[j]^2 - 1
        end
        ret
    end
end

function compensate_gravity_funcs(tgstruct)

    A = make_A(tgstruct)

    Q̃=build_Q̃(tgstruct)

    function F!(F,u)
        clear_forces!(tgstruct)
        actuate!(tgstruct,u)
        update_cables_apply_forces!(tgstruct)
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
    q = get_q(tgref)
    q̌ = get_q̌(tgref)
    tg = deepcopy(tginput)
    clear_forces!(tg)
    update_rigids!(tg,q)
    update_tensiles!(tg)
    if gravity
        Ǧ = build_Ǧ(tg)
    else
        Ǧ = zero(q̌)
    end
    if isnothing(Fˣ)
        F̌ = Ǧ
    else
        F̌ = F̌ˣ + Ǧ
    end

    build_inverse_statics_core(tg,q,F̌)
end

function build_inverse_statics_core(tg,q::AbstractVector,F)
    A = make_A(tg)
    Nᵀ = transpose(nullspace(A(q)))
    Q̃ = build_Q̃(tg)
    tg,Nᵀ,Q̃,F
end

function build_inverse_statics_for_density(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    L = build_L(tg)
    # Left hand side
    Q̃L = Q̃*L

    B = -Nᵀ*Q̃L

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_stiffness(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)

    L̂ = build_L̂(tg)
    ℓ = get_cables_len(tg)
    u = get_cables_restlen(tg)

    # Left hand side
    Q̃L̄ = Q̃*L̂*Diagonal(ℓ-u)

    B = -Nᵀ*Q̃L̄

    # Right hand side
    F̃ = Nᵀ*F

    B,F̃
end

function build_inverse_statics_for_restlength(tginput,tgref,Fˣ=nothing;gravity=false,scale=true)
    tg,Nᵀ,Q̃,F = build_inverse_statics_core(tginput,tgref,Fˣ;gravity)
    ℓ = get_cables_len(tg)
    K̂ = build_K̂(tg)
    # Left hand side
    Q̃K̂ = Q̃*K̂

    B = Nᵀ*Q̃K̂

    # Right hand side
    F̃ = Nᵀ*(Q̃K̂*ℓ + F)
    B,F̃
end

function build_inverse_statics_for_actuation(botinput,botref::TensegrityRobot,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    build_inverse_statics_for_actuation(botinput,botref.tg,Fˣ;Y,gravity,scale)
end

function build_inverse_statics_for_actuation(botinput,tgref::TensegrityStructure,
                                        Fˣ=nothing;Y=build_Y(botinput),gravity=false,scale=true)
    tg,Nᵀ,Q̃,F = build_inverse_statics_core(botinput.tg,tgref,Fˣ;gravity)
    ℓ = get_cables_len(tg)
    K̂ = build_K̂(tg)
    u0 = get_original_restlen(botinput)
    # Left hand side
    Q̃K̂ = Q̃*K̂
    Q̃K̂Y = Q̃K̂*Y

    B = Nᵀ*Q̃K̂Y

    # Right hand side
    F̃ = Nᵀ*(Q̃K̂*(ℓ-u0) + F)

    B,F̃
end

function get_solution_set(B,F̃)
    # A Particular solution
    x = pinv(B)*F̃
    # A basis for the null space of B
    nb = nullspace(B)
    return x,nb
end

function check_restlen(tg)
    u = get_cables_restlen(tg)
    check_restlen(tg,u)
end

function check_restlen(tg,u)
    ℓ = get_cables_len(tg)
    if any((x)->x<0,u)
        @warn "Negative rest lengths"
    end
    if any((x)->x<0,ℓ-u)
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

function get_inverse_func(tg_input,reftg,Y;gravity=false,recheck=true,scale=true)
    # We only use the $q$ of the reference tgure.
    acttg,lhs,rhs,c = statics_equation_for_actuation(tg_input,reftg,Y;gravity,scale)
    ncoeffs = maximum(size(lhs))- rank(lhs)
    xp,nb = get_solution_set(lhs,rhs)
    @info "Number of coefficients = $ncoeffs"
    function inner_inverse_func(ξ)
        x = xp + nb*ξ
        λ = x[1:acttg.nconstraint].*c
        a = x[acttg.nconstraint+1:end]
        # check_actuation(acttg,Y,a)
        # refq = get_q(reftg)
        # actuate!(acttg,a)
        # check_static_equilibrium(acttg,refq,λ;gravity)
        u0 = get_cables_restlen(tg_input)
        rl = u0 + Y*a
        λ,rl,a
    end
end

function check_static_equilibrium(tg_input,q,λ,F=nothing;gravity=false)
    tg = deepcopy(tg_input)
    clear_forces!(tg)
    distribute_q_to_rbs!(tg,q)
    update_cables_apply_forces!(tg)
    check_restlen(tg,get_cables_restlen(tg))
    if gravity
        apply_gravity!(tg)
    end
    generalized_forces = assemble_forces(tg)
    if !isnothing(F)
        generalized_forces .+= F[:]
    end
    constraint_forces = transpose(make_A(tg)(q))*λ
    static_equilibrium = constraint_forces ≈ generalized_forces
    @debug "Res. forces = $(generalized_forces-constraint_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(generalized_forces-constraint_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium
end

function check_static_equilibrium_output_multipliers(tg_input;F=nothing,gravity=false)
    q = get_q(tg_input)
    check_static_equilibrium_output_multipliers(tg_input,q,F;gravity)
end

function check_static_equilibrium_output_multipliers(tg_input,q,F=nothing;gravity=false)
    tg = deepcopy(tg_input)
    clear_forces!(tg)
    update_rigids!(tg)
    update_tensiles!(tg)
    # check_restlen(tg,get_cables_restlen(tg))
    if gravity
        apply_gravity!(tg)
    end
    generalized_forces = generate_forces!(tg)
    if !isnothing(F)
        generalized_forces .+= F[:]
    end
    A = make_A(tg)(q)
    λ = inv(A*transpose(A))*A*(generalized_forces)
    constraint_forces = transpose(A)*λ
    static_equilibrium = constraint_forces ≈ generalized_forces
    @debug "Res. forces = $(generalized_forces-constraint_forces)"
    if !static_equilibrium
        @error "System not in static equilibrium. Err = $(norm(generalized_forces-constraint_forces))"
        @info "This error could be harmless, if the error is sufficiently small, or nonpositive tension occurs."
    end
    static_equilibrium, λ
end

function inverse_for_restlength(botinput,botref::TensegrityRobot,Fˣ=nothing;gravity=false,scale=true,recheck=true)
    inverse_for_restlength(botinput.tg,botref.tg,Fˣ;gravity,scale,recheck)
end

function inverse_for_restlength(tginput,tgref::TensegrityStructure,Fˣ=nothing;gravity=false,scale=true,recheck=true)
    B,F̃ = build_inverse_statics_for_restlength(tginput,tgref,Fˣ;gravity,scale)
    nu = tgref.tensiles.cables |> length
    if check_inverse_sanity(B)
        y0 = B\F̃
    else
        @info "Using Quadratic Programming."
        # COSMO.Settings(verbose = false, eps_abs = 1e-7, eps_rel = 1e-16)
        model = JuMP.Model(COSMO.Optimizer)
        JuMP.set_optimizer_attribute(model, "verbose", false)
        JuMP.set_optimizer_attribute(model, "eps_abs", 1e-7)
        JuMP.set_optimizer_attribute(model, "eps_rel", 1e-16)
        JuMP.@variable(model, y[1:nu])
        JuMP.@objective(model, Min, sum(y[1:nu].^2))
        # JuMP.@objective(model, Max, 0.0)
        JuMP.@constraint(model, static, B*y .== F̃)
        ϵ = 1e-2
        JuMP.@constraint(model, positive_u, y[1:nu].- ϵ .>= 0)
        JuMP.@constraint(model, positive_f, y[1:nu].+ ϵ .<= get_cables_len(tgref))
        # JuMP.print(model)
        JuMP.optimize!(model)
        if JuMP.termination_status(model) == JuMP.MathOptInterface.OPTIMAL
            @info "Optimal Solution found."
        else
            @show JuMP.termination_status(model)
            error("Inverse statics optimization failed.")
        end
        # @show JuMP.primal_status(model)
        # @show JuMP.dual_status(model)
        # @show JuMP.objective_value(model)
        y0 = JuMP.value.(y)
        # @show abs.(B*y0 .- F̃)
        # y0 = pinv(B)*F̃
    end
    u = y0[1:nu]
    if recheck
        tgcheck = deepcopy(tginput)
        q = get_q(tgref)
        set_restlen!(tgcheck,u)
        _,λ = check_static_equilibrium_output_multipliers(tgcheck,q;gravity)
    end
    λ,u
end

function inverse_for_multipliers(botinput::TensegrityRobot,botref::TensegrityRobot=botinput,F=nothing;gravity=false,scale=true,recheck=true)
    inverse_for_multipliers(botinput.tg,botref.tg,F;gravity,scale,recheck)
end

function inverse_for_multipliers(tginput::TensegrityStructure,tgref::TensegrityStructure=tginput,F=nothing;gravity=false,scale=true,recheck=true)
    q = get_q(tgref)
    _,λ = check_static_equilibrium_output_multipliers(tginput,q;gravity)
    λ
end

function inverse_for_actuation(botinput,botref,Fˣ=nothing;Y=build_Y(botinput),
                    gravity=false,scale=true,recheck=true)
    tgref = botref.tg
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
        q = get_q(tgref)
        actuate!(botcheck,a)
        _,λ = check_static_equilibrium_output_multipliers(botcheck.tg,q;gravity)
    end
    λ,a
end
