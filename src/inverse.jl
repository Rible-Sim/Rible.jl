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

function build_G(tgstruct)
    reset_forces!(tgstruct)
    apply_gravity!(tgstruct)
    G = assemble_forces(tgstruct)
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
