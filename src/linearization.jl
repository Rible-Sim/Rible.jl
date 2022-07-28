function swapcols!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = 1:size(X,1)
        X[k,i], X[k,j] = X[k,j], X[k,i]
    end
end
function swaprows!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = 1:size(X,2)
        X[i,k], X[j,k] = X[j,k], X[i,k]
    end
end

function GECP(A_input)
    A = float(copy(A_input))
    n,m = size(A)
    col_index = collect(1:m)
    for k = 1:n
        Akrow2 = (@view A[k:end,k]).^2
        Akcol2 = (@view A[k,k:end]).^2
        ck_i = argmax(Akrow2)
        ck_j = argmax(Akcol2)
        # @show k,"before"
        # display(A)
        if Akrow2[ck_i] < Akcol2[ck_j]
            # swap columns
            swapcols!(A,k,k-1+ck_j)
            col_index[k],col_index[k-1+ck_j] = col_index[k-1+ck_j], col_index[k]
        else
            # swap rows
            swaprows!(A,k,k-1+ck_i)
        end
        # @show k,"after"
        # display(A)
        for i in k+1:n
            m_ik = A[i,k]/A[k,k]
            A[i,:] .-= m_ik*A[k,:]
        end
    end
    col_index
end

function find_full_constrained_index(lncs,q)
    cf = NaturalCoordinates.CoordinateFunctions(lncs)
    Aq = cf.Φq(q)
    col_index = GECP(Aq)
    col_index[size(Aq,1)+1:end]
end

<<<<<<< Updated upstream
function ∂Aᵀλ∂q(tg,λ)
=======
function ∂Aᵀλ∂q̌(tg::AbstractTensegrityStructure,λ)
    (;nfree) = tg.connectivity.indexed
    ret = zeros(eltype(λ),nfree,nfree)
    (;rigidbodies,nconstraints) = tg
    (;indexed,jointed) = tg.connectivity
    (;ninconstraints,mem2sysfree,mem2sysincst) = indexed
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
        memfree = mem2sysfree[rbid]
        memincst = mem2sysincst[rbid]
        uci = rb.state.cache.unconstrained_index
        if !isempty(memincst)
            ret[memfree,memfree] .+= rb.state.cache.funcs.∂Aᵀλ∂q(λ[memincst])[:,uci]
        end
    end
    ret
end

function ∂Aq̇∂q(tg,q̇)
>>>>>>> Stashed changes
    body2q = tg.connectivity.body2q
    ncoords = tg.ncoords
    nbodyc = get_nbodyconstraint(tg)
    ret = zeros(eltype(λ),ncoords,ncoords)
    is = 0
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        rb = tg.rigidbodies[rbid]
        # nc = rb.state.cache.nc
        # if nc > 0
        #     is += nc
        # end
        ret[pindex,pindex] .+= rb.state.cache.cfuncs.∂Aᵀλ∂q(λ[is+1:is+nbodyc])
        is += nbodyc
    end
    ret
end

function ∂Aq̇∂q(tgstruct,q̇)
    body2q = tgstruct.connectivity.body2q
    @unpack ncoords, nconstraint = tgstruct
    nbodyc = get_nbodyconstraint(tgstruct)
    ret = zeros(get_numbertype(tgstruct),nconstraint,ncoords)
    is = 0
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        rb = tgstruct.rigidbodies[rbid]
        q̇_rb = q̇[pindex]
        nc = rb.state.cache.nc
        if nc > 0
            is += nc
        end
        ret[is+1:is+nbodyc,pindex] .+= rb.state.cache.cfuncs.∂Aq̇∂q(q̇_rb)
        is += nbodyc
    end
    ret
end

function test_fvector(tgstruct,q0)
    function L(q)
        reset_forces!(tgstruct)
        distribute_q_to_rbs!(tgstruct,q,zero(q))
        update_strings_apply_forces!(tgstruct)
        fvector(tgstruct)
        [tgstruct.strings[i].state.length for i = 1:2]
    end
    FiniteDiff.finite_difference_jacobian(L,q0)
end

function build_K(tg,λ)
    Q̃ = build_Q̃(tg)
    q,_ = get_q(tg)
    ∂L∂q,_ = build_tangent(tg)
    K = -Q̃*∂L∂q .+ ∂Aᵀλ∂q(tg,λ)
end

function linearize(tginput,λ,u,q,q̇=zero(q))
    tg = deepcopy(tginput)
    set_restlen!(tg,u)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,q,q̇)
    update_strings!(tg)
    M = build_massmatrix(tg)
    A = build_A(tg)
    Q̃ = build_Q̃(tg)
    Jac_Γ = build_Jac_Γ(tg)
    ∂L∂q,∂L∂q̇ = Jac_Γ(q,q̇)
    @unpack ncoords,nconstraint = tg
    nz = ncoords + nconstraint
    M̂ = zeros(eltype(q),nz,nz)
    Ĉ  = zeros(eltype(q),nz,nz)
    K̂ = zeros(eltype(q),nz,nz)
    M̂[1:ncoords,1:ncoords] .= M
    Ĉ[1:ncoords,1:ncoords] .= -Q̃*∂L∂q̇

    # fjac = test_fvector(tgstruct,q)
    K̂[1:ncoords,1:ncoords] .= -Q̃*∂L∂q .+ ∂Aᵀλ∂q(tg,λ)
    Aq = A(q)
    c = maximum(abs.(K̂[1:ncoords,1:ncoords]))
    K̂[1:ncoords,ncoords+1:nz] .= c.*transpose(Aq)
    K̂[ncoords+1:nz,1:ncoords] .= c.*Aq
    M̂,Ĉ,K̂
end


function frequencyshift(M̂,Ĉ,K̂,α::Real)
    M̄ = M̂
    C̄ = 2α*M̂ + Ĉ
    K̄ = α^2*M̂ + α*Ĉ + K̂
    M̄,C̄,K̄
end

function frequencyshift(M̂,K̂,α::Real)
    M̄ = M̂
    K̄ = α*M̂ + K̂
    M̄,K̄
end


function enlarge(M̄,C̄,K̄)
    T = eltype(M̄)
    nz = size(M̄)[1]
    M̃ = zeros(T,2nz,2nz)
    M̃[1:nz,1:nz] .= -C̄
    M̃[1:nz,nz+1:2nz] .= -M̄
    M̃[nz+1:2nz,1:nz] .= Matrix(one(T)*I,nz,nz)
    K̃ = zeros(T,2nz,2nz)
    K̃[1:nz,1:nz] .= K̄
    K̃[nz+1:2nz,nz+1:2nz] .= Matrix(one(T)*I,nz,nz)
    M̃,K̃
end

function find_finite(ω2,Z,ndof)
    first_frequency_index = findfirst((x)->x>0,ω2)
    finite_ω2 = ω2[first_frequency_index:first_frequency_index+ndof-1]
    finite_Z = Z[:,first_frequency_index:first_frequency_index+ndof-1]
    finite_ω2,finite_Z
end

<<<<<<< Updated upstream
function normalize_wrt_mass!(Z,M)
    n = size(Z)[2]
=======
function build_Ǩ(tg)
    _,λ = check_static_equilibrium_output_multipliers(tg)
    build_Ǩ(tg,λ)
end

function build_∂Q̌∂q̌!(∂Q̌∂q̌,tg)
    (;cables,connectivity) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    # ∂Q̌∂q̌ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected.cables) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        uci1 = rb1.state.cache.unconstrained_index
        uci2 = rb2.state.cache.unconstrained_index
        mfree1 = mem2sysfree[rb1.prop.id]
        mfree2 = mem2sysfree[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension,length,lengthdot) = state
        if slack && (tension==0)
            ∂Q̌∂q̌ .-= 0
        else
            D .= direction*transpose(direction)
            density = tension/length
            β = c*lengthdot/length + density
            D .*= k-β
            D .+= β.*Im
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,uci2]
            J̌[:,mfree1] .-= C1[:,uci1]
            ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌
end

function build_∂Q̌∂q̌(tg,@eponymargs(cables,))
    (;cables,connectivity) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected.cables) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        uci1 = rb1.state.cache.unconstrained_index
        uci2 = rb2.state.cache.unconstrained_index
        mfree1 = mem2sysfree[rb1.prop.id]
        mfree2 = mem2sysfree[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension,length,lengthdot) = state
        if slack && (tension==0)
            ∂Q̌∂q̌ .-= 0
        else
            D .= direction*transpose(direction)
            density = tension/length
            β = c*lengthdot/length + density
            D .*= k-β
            D .+= β.*Im
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,uci2]
            J̌[:,mfree1] .-= C1[:,uci1]
            ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌
end

function build_∂Q̌∂q̌(tg,@eponymargs(clustercables))
    (;clustercables,connectivity) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    i = 0
    foreach(connected.clustercables) do clustercable
        i += 1
        foreach(clustercable) do cc
            cable = tg.clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            uci1 = rb1.state.cache.unconstrained_index
            uci2 = rb2.state.cache.unconstrained_index
            mfree1 = mem2sysfree[rb1.prop.id]
            mfree2 = mem2sysfree[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension,length,lengthdot) = state
            if tension==0
                ∂Q̌∂q̌ .-= 0
            else
                D .= direction*transpose(direction)
                density = tension/length
                β = c*lengthdot/length + density
                D .*= k-β
                D .+= β.*Im
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,uci2]
                J̌[:,mfree1] .-= C1[:,uci1]
                ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
            end
        end
    end
    return ∂Q̌∂q̌
end

function build_∂Q̌∂q̌(tg)
    build_∂Q̌∂q̌(tg, tg.tensiles)
end

function build_∂Q̌∂q̌(tg, @eponymargs(cables, clustercables))
    ∂Q̌∂q̌1 = build_∂Q̌∂q̌(tg, @eponymtuple(cables))
    ∂Q̌∂q̌2 = build_∂Q̌∂q̌(tg, @eponymtuple(clustercables))
    return ∂Q̌∂q̌1 + ∂Q̌∂q̌2
end

function build_∂Q̌∂q̌̇!(∂Q̌∂q̌̇,tg)
    (;cables,connectivity) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    # ∂Q̌∂q̌̇ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected.cables) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        uci1 = rb1.state.cache.unconstrained_index
        uci2 = rb2.state.cache.unconstrained_index
        mfree1 = mem2sysfree[rb1.prop.id]
        mfree2 = mem2sysfree[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension) = state
        if slack && (tension == 0)
            ∂Q̌∂q̌̇ .-= 0
        else
            D .= direction*transpose(direction)
            D .*= c
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,uci2]
            J̌[:,mfree1] .-= C1[:,uci1]

            ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
end

function build_∂Q̌∂q̌̇(tg, @eponymargs(cables, ))
    (;cables,connectivity) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌̇ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected.cables) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        uci1 = rb1.state.cache.unconstrained_index
        uci2 = rb2.state.cache.unconstrained_index
        mfree1 = mem2sysfree[rb1.prop.id]
        mfree2 = mem2sysfree[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension) = state
        if slack && (tension == 0)
            ∂Q̌∂q̌̇ .-= 0
        else
            D .= direction*transpose(direction)
            D .*= c
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,uci2]
            J̌[:,mfree1] .-= C1[:,uci1]

            ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfrjiexee2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌̇
end

function build_∂Q̌∂q̌̇(tg, @eponymargs(clustercables, ))
    (;cables,connectivity) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌̇ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    i = 0
    foreach(connected.clustercables) do clustercable
        i += 1
        foreach(clustercable) do cc
            cable = clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            uci1 = rb1.state.cache.unconstrained_index
            uci2 = rb2.state.cache.unconstrained_index
            mfree1 = mem2sysfree[rb1.prop.id]
            mfree2 = mem2sysfree[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension) = state
            if tension == 0
                ∂Q̌∂q̌̇ .-= 0
            else
                D .= direction*transpose(direction)
                D .*= c
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,uci2]
                J̌[:,mfree1] .-= C1[:,uci1]

                ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
            end
        end
    end
    return ∂Q̌∂q̌̇
end

function build_∂Q̌∂q̌̇(tg)
    return build_∂Q̌∂q̌̇(tg, tg.tensiles)
end

function build_∂Q̌∂q̌̇(tg, @eponymargs(cables, clustercables))
    ∂Q̌∂q̌̇1 = build_∂Q̌∂q̌̇(tg, @eponymtuple(cables))
    ∂Q̌∂q̌̇2 = build_∂Q̌∂q̌̇(tg, @eponymtuple(clustercables))
    return ∂Q̌∂q̌̇1 + ∂Q̌∂q̌̇2
end

function build_∂Q̌∂s̄(tg)
    (;cables,clustercables,connectivity, nclustercables) = tg
    (;connected,indexed) = connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂s̄ = zeros(T,2ns,nfree)
    D = zeros(T, ndim)
    lkn = zeros(T, 2ns, ndim)
    # Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)

    N_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    kc = Vector{Float64}()
    for (cid,clustercable) in enumerate(clustercables)
        @unpack segs = clustercable
        nsegs = length(segs)
        for (sid, seg) in enumerate(segs)
            push!(kc, seg.k)
        end
        N = sparse(zeros(Float64, nsegs, 2nsegs-2))
        N[1,1:2]=[1 -1]; N[end,end-1:end]=[-1 1]
        for i in 2:nsegs-1
            N[i, 2i-3:2i] = [-1 1 1 -1]
        end
        push!(N_list, N)
    end
    N = reduce(blockdiag,N_list)
    i = 0; j = 0
    foreach(connected.clustercables) do clustercable
        i += 1
        foreach(clustercable) do cc
            j += 1
            cable = clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            uci1 = rb1.state.cache.unconstrained_index
            uci2 = rb2.state.cache.unconstrained_index
            mfree1 = mem2sysfree[rb1.prop.id]
            mfree2 = mem2sysfree[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension) = state
            if tension == 0
                ∂Q̌∂s̄ .-= 0
            else
                D .= direction
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,uci2]
                J̌[:,mfree1] .-= C1[:,uci1]
                kN = kc[j] .* N[j,:]
                @tullio lkn[k, l] = D[l] * kN[k]
                ∂Q̌∂s̄ .-= lkn * J̌
            end
        end
    end
    return ∂Q̌∂s̄'
end

function build_Ǩ(tg,λ)
    (;nfree) = tg.connectivity.indexed
    T = get_numbertype(tg)
    # Ǩ = zeros(T,nfree,nfree)
    Ǩ = build_∂Q̌∂q̌(tg)
    Ǩ .= ∂Aᵀλ∂q̌(tg,λ) .-Ǩ
    # Ǩ .= Ǩ
    Ǩ
end

function norm_wrt!(Z,M)
    n = size(Z,2)
>>>>>>> Stashed changes
    for i = 1:n
        zmz = transpose(Z[:,i])*M*Z[:,i]
        Z[:,i] ./= sqrt(zmz)
    end
    Z
end
function undamped_eigen(tg,q0,λ0)
    if !check_static_equilibrium(tg,q0,λ0;gravity=false)
        @warn "Statics check failed, but proceed anyway."
    end
    u0 = get_strings_restlen(tg)
    M̂,Ĉ,K̂ = linearize(tg,λ0,u0,q0)
    α = 10
    M̄,K̄ = frequencyshift(M̂,K̂,α)
    # @show size(K̄),rank(K̄),cond(K̄),rank(M̄)
    d,aug_Z = eigen(K̄,M̄)
    aug_ω2 = d .- α
    @unpack ncoords, ndof = tg
    # @show aug_ω2
    ω2,Z = find_finite(aug_ω2,aug_Z,ndof)
    ω = sqrt.(ω2)
    Zq = Z[1:ncoords,:]
    M = build_massmatrix(tg)
    normalize_wrt_mass!(Zq,M)
    ω, Zq#, Z
end

function undamped_modal_solve!(tgstruct,q0,q̇0,λ0,tf,dt)
    M̂,Ĉ,K̂ = linearize(tgstruct,q0,λ0)
    # show(stdout,"text/plain",K̂)
    # showtable(K̂)
    # M̄,C̄,K̄ = TR.frequencyshift(M̂,Ĉ,K̂,0.1)
    # M̃,K̃ = TR.enlarge(M̄,C̄,K̄)
    aug_ω2,aug_Z = eigen(K̂,M̂)
    ω2,Z = find_finite(aug_ω2,aug_Z)
    # @show aug_ω2,ω2
    normalize_wrt_mass!(Z,M̂)
    # @show transpose(Z)*M̂*Z
    # @show transpose(Z)*K̂*Z
    ω = sqrt.(ω2)
    z0 = vcat(zero(q0),λ0)
    ż0 = vcat(q̇0,zero(λ0))
    ζ0 = transpose(Z)*M̂*z0
    ζd0 = transpose(Z)*M̂*ż0

    d = length(ζ0)
    step = Integer(tf/dt)
    ζ = Matrix{eltype(q0)}(undef,d,step+1)
    for it in 0:step
        t = dt*it
        ζ[:,it+1] .= ζ0.*cos.(ω.*t) .+ ζd0./ω.*sin.(ω.*t)
    end
    z = Z*ζ
    q = z[1:length(q0),:]
end

function find_nullspace(c)
    Nc,Nu = size(c)
    P = VectorOfArray([ begin
                            p = zeros(eltype(c),Nu)
                            p[i] = 1
                            p
                        end for i = 1:Nu])
    for i = 1:Nc
        # @show i
        cPⁱ⁻¹ = c*P
        rowi = cPⁱ⁻¹[i,:]
        _ℓⁱ = findall((x)->!iszero(x),rowi)
        ℓⁱ = _ℓⁱ[sortperm(rowi[_ℓⁱ],order=Base.Order.Reverse)]
        if isempty(ℓⁱ)
            continue
        end
        # display(cPⁱ⁻¹)
        # @show ℓⁱ
        # for k = 1:(length(ℓ)-1)
        for k = 1:(length(ℓⁱ)-1)
            αⁱₖ = -cPⁱ⁻¹[i,ℓⁱ[k]]./cPⁱ⁻¹[i,ℓⁱ[k+1]]
            P[:,ℓⁱ[k]] .= P[:,ℓⁱ[k]] + αⁱₖ*P[:,ℓⁱ[k+1]]
        end
        deleteat!(P.u,ℓⁱ[end])
    end
    Array(P)
end

function check_stability(tg;verbose=false)
    λ = inverse_for_multipliers(tg,tg)
    check_stability(tg,λ;verbose)
end

function check_stability(tg,λ;verbose=false)
    K = build_K(tg,λ)
    A = build_A(tg)
    q,_ = get_q(tg)
    N = nullspace(A(q))
    # N = find_nullspace(A(q))
    Ǩ = transpose(N)*K*N
    eigen_result = eigen(Ǩ)
    eigen_min = eigen_result.values[1]
    if verbose
        @show eigen_min
    end
    if eigen_min<0
        @warn "Instability detected!"
        return false
    else
        return true
    end
end
