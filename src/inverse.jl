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
        setblock!(C,build_Ci(rbs[rbid]),mvrbid,mvrbid)
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

function build_K(tgstruct)
    q,_ = get_q(tgstruct)
    A = build_A(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    K̂ = build_K̂(tgstruct)
    # ℓ = build_ℓ(tgstruct)
    K = hcat(
        transpose(A(q)),
        Q̃*K̂
        )
end

function build_G(tgstruct;factor=1.0)
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
    G = build_G(tgstruct)
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

function build_U(tgstruct,a,u)
    @unpack ncoords,nstrings,ndim,strings = tgstruct
    ret = zeros(eltype(a),nstrings*ndim,ncoords)
    for i = 1:nstrings
        k = strings[i].k
        Ji = Array(build_Ji(tgstruct,i))
        ret[(i-1)*ndim+1:i*ndim,:] = k*Ji*(1-a[i]*u[i])
    end
    ret
end

function build_S(tgstruct,a,q)
    @unpack nstrings = tgstruct
    ret = zeros(eltype(a),nstrings)
    for i = 1:nstrings
        Ji = build_Ji(tgstruct,i)
        Wi = Array(transpose(Ji)*Ji)
        ret[i] = transpose(q)*Wi*q*a[i]^2 - 1
    end
    ret
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



function inverse(tgstruct_input,refstruct,Y;gravity=false,recheck=true)
    tgstruct = deepcopy(tgstruct_input)
    refq0,_,_ = get_initial(refstruct)

    reset_forces!(tgstruct)
    distribute_q_to_rbs!(tgstruct,refq0)
    update_strings_apply_forces!(tgstruct)

    L̂ = build_L̂(tgstruct)
    K̂ = build_K̂(tgstruct)
    ℓ = build_ℓ(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    A = build_A(tgstruct)
    Aq = A(refq0)
    u0 = [s.original_restlen for s in tgstruct.strings]
    rhs = Q̃*K̂*(ℓ-u0)
    if gravity == true
        G = build_G(tgstruct)
        rhs += G
    end
    lhs = hcat(transpose(Aq),Q̃*K̂*Y)
    # @show rank(lhs),size(lhs)
    x = lhs\rhs
    λ = x[1:tgstruct.nconstraint]
    a = x[tgstruct.nconstraint+1:end]
    Δu = Y*a

    original_restlens = [s.original_restlen for s in tgstruct.strings]
    restlens = original_restlens + Δu
    if any((x)->x<0,restlens)
        @error "Negative rest lengths"
    end
    lengths = [s.state.length for s in tgstruct.strings]
    if any((x)->x<0,lengths-restlens)
        @error "Zero tension"
    end
    # s = 1 ./ lengths
    actuate!(tgstruct,a)
    reset_forces!(tgstruct)
    update_strings_apply_forces!(tgstruct)
    appproximated = transpose(build_A(tgstruct)(refq0))*λ ≈ build_Q̃(tgstruct)*fvector(tgstruct)
    @info "Recheck result: $appproximated"
    λ,Δu,a
end

function forward(tgstruct,starts,targets)
    @var q[1:tgstruct.ncoords]
    @var a[1:tgstruct.nstrings]
    @var λ[1:tgstruct.nconstraint]
    @var u[1:tgstruct.nstrings]
    @var g
    polyq = 1.0q
    polya = 1.0a
    polyλ = 1.0λ
    polyu = 1.0u

    q0,a0,λ0,u0,g0 = starts
    u1,g1 = targets

    Φ = build_Φ(tgstruct,q0)
    A = build_A(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    U = build_U(tgstruct,polya,polyu)
    S = build_S(tgstruct,polya,polyq)

    G = build_G(tgstruct)
    F = [-transpose(A(polyq))*polyλ  + g*G + Q̃*U*polyq;
        S;
        Φ(polyq)]

    # reset_forces!(tgstruct)
    # update_strings_apply_forces!(tgstruct)
    # l = [s.state.length for s = tgstruct.strings]
    # u1 = u0
    # g0 = 1.0
    # g1 = 0.0

    startparas = vcat(u0,g0)
    endparas = vcat(u1,g1)
    startsols = [[q0; a0; λ0]]
    vg = [q,a,λ]
    # FFF = [subs(f, u=>l,g=>0.0) for f in F]
    Fsys = System(F; variable_groups=vg, parameters = [u;g])

    ParameterHomotopy(Fsys,startparas,endparas),startsols
    # [subs(f, q=>q0,λ=>λ0,a=>a0,u=>u0,g=>1.0) for f in F]
    result = solve(Fsys, startsols; start_parameters=startparas, target_parameters=endparas)
end
