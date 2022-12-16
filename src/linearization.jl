"""
交换列。
$(TYPEDSIGNATURES)
"""
function swapcols!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = 1:size(X,1)
        X[k,i], X[k,j] = X[k,j], X[k,i]
    end
end
"""
交换行。
$(TYPEDSIGNATURES)
"""
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

function find_full_pres_indices(lncs,q)
    cf = NaturalCoordinates.CoordinateFunctions(lncs,
            collect(1:NaturalCoordinates.get_ncoords(lncs)),
            collect(1:NaturalCoordinates.get_nconstraints(lncs))
        )
    Aq = cf.Φq(q)
    col_index = GECP(Aq)
    col_index[size(Aq,1)+1:end] |> sort
end

function ∂Aᵀλ∂q̌(tg::AbstractTensegrityStructure,λ)
    (;nfree) = tg.connectivity.indexed
    ret = zeros(eltype(λ),nfree,nfree)
    (;bodies,nconstraints) = tg
    (;indexed,jointed) = tg.connectivity
    (;ninconstraints,mem2sysfree,mem2sysincst) = indexed
    foreach(bodies) do rb
        rbid = rb.prop.id
        memfree = mem2sysfree[rbid]
        memincst = mem2sysincst[rbid]
        free_idx = rb.state.cache.free_idx
        if !isempty(memincst)
            ret[memfree,memfree] .+= rb.state.cache.funcs.∂Aᵀλ∂q(λ[memincst])[:,free_idx]
        end
    end
    ret
end

function ∂Aq̇∂q(tg,q̇)
    (;nfree) = tg.connectivity.indexed
    (;bodies,nconstraints) = tg
    (;indexed,jointed) = tg.connectivity
    (;ninconstraints,mem2sysfree,mem2sysincst) = indexed
    ret = zeros(eltype(q̇),nconstraints,nfree)
    foreach(bodies) do rb
        rbid = rb.prop.id
        memfree = mem2sysfree[rbid]
        memincst = mem2sysincst[rbid]
        free_idx = rb.state.cache.free_idx
        Φi = rb.state.cache.Φi
        if !isempty(memincst)
            ret[memincst,memfree] .+= rb.state.cache.funcs.∂Aq̇∂q(q̇[memfree])[Φi,free_idx]
        end
    end
    ret
end

function test_fvector(tg,q0)
    function L(q)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg,q,zero(q))
        update_cables_apply_forces!(tg)
        fvector(tg)
        [tg.cables[i].state.length for i = 1:2]
    end
    FiniteDiff.finite_difference_jacobian(L,q0)
end

"""
线性化。
$(TYPEDSIGNATURES)
"""
function linearize(tginput,λ,u,q,q̇=zero(q))
    tg = deepcopy(tginput)
    set_restlen!(tg,u)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,q,q̇)
    update_cables_apply_forces!(tg)
    M = build_massmatrix(tg)
    A = build_A(tg)
    Q̃ = build_Q̃(tg)
    ∂L∂q,∂L∂q̇ = build_tangent(tg)
    @unpack ncoords,nconstraint = tg
    nz = ncoords + nconstraint
    M̂ = zeros(eltype(q),nz,nz)
    Ĉ  = zeros(eltype(q),nz,nz)
    K̂ = zeros(eltype(q),nz,nz)
    M̂[1:ncoords,1:ncoords] .= M
    Ĉ[1:ncoords,1:ncoords] .= -Q̃*∂L∂q̇

    # fjac = test_fvector(tg,q)
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

function build_Ǩ(tg)
    _,λ = check_static_equilibrium_output_multipliers(tg)
    build_Ǩ(tg,λ)
end

function make_Ǩm_Ǩg(tg,q0)
    (;ndim) = tg
    (;numbered,indexed,tensioned) = tg.connectivity
    (;nfull,nfree,syspres,sysfree,mem2sysfull) = indexed
    (;connected) = tensioned
    (;cables) = tg.tensiles
    (;mem2num,num2sys) = numbered
    function inner_Ǩm_Ǩg(q̌,s,μ,k,c)
		q = Vector{eltype(q̌)}(undef,nfull)
		q[syspres] .= q0[syspres]
		q[sysfree] .= q̌
        Jj = zeros(eltype(q̌),ndim,nfull)
        retǨm = zeros(eltype(q̌),nfree,nfree)
        retǨg = zeros(eltype(q̌),nfree,nfree)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.end1.rbsig
            rb2 = scnt.end2.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.end1.pid
            ap2id = scnt.end2.pid
            c1 = c[num2sys[mem2num[rb1id][ap1id]]]
            c2 = c[num2sys[mem2num[rb2id][ap2id]]]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = mem2sysfull[rb1.prop.id]
            mfull2 = mem2sysfull[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            Ǔj = @view Uj[sysfree,sysfree]
            Ūjq = Uj[sysfree,:]*q
            retǨm .+= k[j]*s[j]^2*(Ūjq*transpose(Ūjq))
            retǨg .+= k[j]*(1-μ[j]*s[j])*(Ǔj-s[j]^2*Ūjq*transpose(Ūjq))
        end
        retǨm,retǨg
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
    function inner_S(q̌,s,c)
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
            c1 = c[num2sys[mem2num[rb1id][ap1id]]]
            c2 = c[num2sys[mem2num[rb2id][ap2id]]]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
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
    inner_S
end

# Out-of-place ∂Q̌∂q̌
function build_∂Q̌∂q̌(tg)
    build_∂Q̌∂q̌(tg, tg.connectivity.tensioned)
end

# Out-of-place ∂Q̌∂q̌ for cables and clustered cables
function build_∂Q̌∂q̌(tg, @eponymargs(connected, clustered))
    ∂Q̌∂q̌1 = build_∂Q̌∂q̌(tg, @eponymtuple(connected))
    ∂Q̌∂q̌2 = build_∂Q̌∂q̌(tg, @eponymtuple(clustered))
    return ∂Q̌∂q̌1 + ∂Q̌∂q̌2
end

# Out-of-place ∂Q̌∂q̌ for cables
function build_∂Q̌∂q̌(tg,@eponymargs(connected,))
    (;cables) = tg.tensiles
    (;indexed) = tg.connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        free_idx1 = rb1.state.cache.free_idx
        free_idx2 = rb2.state.cache.free_idx
        mfree1 = mem2sysfree[rb1.prop.id]
        mfree2 = mem2sysfree[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension,length,lengthdot) = state
        if slack && (tension==0)
            ∂Q̌∂q̌ .= 0
        else
            D .= direction*transpose(direction)
            density = tension/length
            β = c*lengthdot/length + density
            D .*= k-β
            D .+= β.*Im
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]
            ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌
end

# Out-of-place ∂Q̌∂q̌ for cluster cables
function build_∂Q̌∂q̌(tg,@eponymargs(clustered))
    (;clustercables) = tg.tensiles
    (;indexed) = tg.connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    i = 0
    foreach(clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            cable = clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            free_idx1 = rb1.state.cache.free_idx
            free_idx2 = rb2.state.cache.free_idx
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
                J̌[:,mfree2] .+= C2[:,free_idx2]
                J̌[:,mfree1] .-= C1[:,free_idx1]
                ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
            end
        end
    end
    return ∂Q̌∂q̌
end

# In-place ∂Q̌∂q̌ for cables and flexible bodies
function build_∂Q̌∂q̌!(∂Q̌∂q̌,tg)
    (;bodies,connectivity) = tg
    (;tensioned,indexed) = connectivity
    (;cables) = tg.tensiles
    (;connected) = tensioned
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    # ∂Q̌∂q̌ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        free_idx1 = rb1.state.cache.free_idx
        free_idx2 = rb2.state.cache.free_idx
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
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]
            ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end

    foreach(bodies) do body
        if body isa FlexibleBody
            (;cache) = body.state
            (;e,funcs) = cache
            (;ancs) = funcs
            ∂Q∂e = ANCF.make_∂Q∂e(ancs)(e)
            mfree = mem2sysfree[body.prop.id]
            free_idx = body.state.cache.free_idx
            ∂Q̌∂q̌[mfree,mfree] .-= ∂Q∂e[free_idx,free_idx]
        end
    end

    return ∂Q̌∂q̌
end

# Out-of-place ∂Q̌∂q̌̇
function build_∂Q̌∂q̌̇(tg)
    build_∂Q̌∂q̌̇(tg, tg.connectivity.tensioned)
end

# Out-of-place ∂Q̌∂q̌̇ for cables and clustered cables
function build_∂Q̌∂q̌̇(tg, @eponymargs(connected, clustered))
    ∂Q̌∂q̌̇1 = build_∂Q̌∂q̌̇(tg, @eponymtuple(connected))
    ∂Q̌∂q̌̇2 = build_∂Q̌∂q̌̇(tg, @eponymtuple(clustered))
    return ∂Q̌∂q̌̇1 + ∂Q̌∂q̌̇2
end

# Out-of-place ∂Q̌∂q̌̇ for cables
function build_∂Q̌∂q̌̇(tg, @eponymargs(connected, ))
    (;cables) = tg.tensiles
    (;indexed) = tg.connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌̇ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        free_idx1 = rb1.state.cache.free_idx
        free_idx2 = rb2.state.cache.free_idx
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
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]

            ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfrjiexee2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌̇
end

# Out-of-place ∂Q̌∂q̌̇ for clustered cables
function build_∂Q̌∂q̌̇(tg, @eponymargs(clustered, ))
    (;clustercables) = tg.tensiles
    (;indexed) = tg.connectivity
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    ∂Q̌∂q̌̇ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    i = 0
    foreach(clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            cable = clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            free_idx1 = rb1.state.cache.free_idx
            free_idx2 = rb2.state.cache.free_idx
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
                J̌[:,mfree2] .+= C2[:,free_idx2]
                J̌[:,mfree1] .-= C1[:,free_idx1]

                ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
            end
        end
    end
    return ∂Q̌∂q̌̇
end

# In-place ∂Q̌∂q̌̇ for cables
function build_∂Q̌∂q̌̇!(∂Q̌∂q̌̇,tg)
    (;tensioned,indexed) = tg.connectivity
    (;connected) = tensioned
    (;cables) = tg.tensiles
    (;nfull,nfree,sysfree,mem2sysfree,mem2sysfull) = indexed
    T = get_numbertype(tg)
    ndim = get_ndim(tg)
    # ∂Q̌∂q̌̇ = zeros(T,nfree,nfree)
    D = @MMatrix zeros(T,ndim,ndim)
    Im = Symmetric(SMatrix{ndim,ndim}(one(T)*I))
    J̌ = zeros(T,ndim,nfree)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;end1,end2) = cc
        rb1 = end1.rbsig
        rb2 = end2.rbsig
        C1 = rb1.state.cache.Cps[end1.pid]
        C2 = rb2.state.cache.Cps[end2.pid]
        free_idx1 = rb1.state.cache.free_idx
        free_idx2 = rb2.state.cache.free_idx
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
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]

            ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
end

function build_∂Q̌∂s̄(tg)
    (;connectivity) = tg
    (;cables,clustercables) = tg.tensiles
    nclustercables = length(clustercables)
    (;tensioned,indexed) = connectivity
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
    foreach(tensioned.clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            j += 1
            cable = clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            free_idx1 = rb1.state.cache.free_idx
            free_idx2 = rb2.state.cache.free_idx
            mfree1 = mem2sysfree[rb1.prop.id]
            mfree2 = mem2sysfree[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension) = state
            if tension == 0
                ∂Q̌∂s̄ .-= 0
            else
                D .= direction
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,free_idx2]
                J̌[:,mfree1] .-= C1[:,free_idx1]
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
    Ǩ = -build_∂Q̌∂q̌(tg) .- ∂Aᵀλ∂q̌(tg,λ)
    # Ǩ .= Ǩ
    Ǩ
end

function norm_wrt!(Z,M)
    n = size(Z,2)
    for i = 1:n
        z = @view Z[:,i]
        zmz = transpose(z)*M*z
        z ./= sqrt(zmz)
    end
    Z
end

function undamped_eigen(tg;gravity=false)
    _,λ = check_static_equilibrium_output_multipliers(tg;gravity)
    q = get_q(tg)
    q̌ = get_q̌(tg)
    M̌ = build_M̌(tg)
    Ǩ = build_Ǩ(tg,λ)
    Ǎ = make_A(tg)(q)
    Ň = nullspace(Ǎ)
    ℳ = transpose(Ň)*M̌*Ň
    𝒦 = transpose(Ň)*Ǩ*Ň
    # @show ℳ, 𝒦
    ω²,ξ = eigen(𝒦,ℳ)
    # @show transpose(ξ)*ℳ*ξ
    Ňξ = Ň*ξ
    # @show transpose(Ňξ)*M̌*Ňξ
    norm_wrt!(Ňξ,M̌)
    δq̌ = [v for v in eachcol(Ňξ)]
    ω²,δq̌
    # nq = length(q̌)
    # nλ = length(λ)
    # nx = nq + nλ
    # M̂ = zeros(eltype(q),nx,nx)
    # K̂ = zeros(eltype(q),nx,nx)
    # M̂[1:nq,1:nq] .= M̌
    # K̂[1:nq,1:nq] .= Ǩ
    # c = maximum(abs.(K̂[1:nq,1:nq]))
    # K̂[1:nq,nq+1:nx] .= c.*transpose(Ǎ)
    # K̂[nq+1:nx,1:nq] .= c.*Ǎ
    #
    # eigen(K̂,M̂)
end

function undamped_eigen!(bot::TensegrityRobot;gravity=false,scaling=0.01)
    (;tg,traj) = bot
    q̌ = get_q̌(tg)
    ω²,δq̌ = undamped_eigen(tg;gravity)
    neg_indices = ω².<=0
    if !isempty(neg_indices)
        @warn "Negative ω² occurs. zeroing."
        ω²[neg_indices] .= 0
    end
    ω = sqrt.(ω²)
    resize!(traj,1)
    nω = length(ω)
    for i = 1:nω
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = ω[i]
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    bot
end

function old_undamped_eigen(tg)
    λ0 = check_static_equilibrium_output_multipliers(tg)
    M̂,Ĉ,K̂ = linearize(tg,q0,λ0)
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

function undamped_modal_solve!(tg,q0,q̇0,λ0,tf,dt)
    M̂,Ĉ,K̂ = linearize(tg,q0,λ0)
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

"""
返回零空间。
$(TYPEDSIGNATURES)
"""
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

"""
校核稳定性。
$(TYPEDSIGNATURES)
"""
function check_stability(tg::TensegrityStructure;F̌=nothing,verbose=false)
    static_equilibrium,λ = check_static_equilibrium_output_multipliers(tg;F=F̌)
    @assert static_equilibrium
    check_stability(tg,λ;verbose)
end

function check_stability(tg::TensegrityStructure,λ;verbose=false)
    q = get_q(tg)
    c = get_c(tg)
    A = make_A(tg,q)
    Ň(q̌,c) = nullspace(A(q̌))
    check_stability(tg,λ,Ň;verbose)
end

function check_stability(tg::TensegrityStructure,λ,Ň;verbose=false)
    q̌ = get_q̌(tg)
    c = get_c(tg)
    Ǩ0 = build_Ǩ(tg,λ)
    Ň0 = Ň(q̌,c)
    𝒦0 = transpose(Ň0)*Ǩ0*Ň0
    eigen_result = eigen(𝒦0)
    eigen_min = eigen_result.values[1]
    if eigen_min<0
        @warn "Instability detected! Minimum eigen value: $eigen_min"
        isstable = false
    else
        isstable = true
    end
    isstable, Ň0, eigen_result
end

function check_stability!(bot::TensegrityRobot,Ň;
        gravity=false,
        scaling=0.01,
        scalings=nothing
    )
    (;tg,traj) = bot
    static_equilibrium,λ = check_static_equilibrium_output_multipliers(tg)
    @assert static_equilibrium
    q̌ = get_q̌(tg)
    _, Ň0, er = check_stability(bot.tg,λ,Ň;verbose=true)
    resize!(traj,1)
    for i in 1:length(er.values)
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = er.values[i]
        δq̌i = Ň0*er.vectors[:,i]
        # @show δq̌i, er.vectors[:,i]
        if scalings isa Nothing
            si = scaling
        else
            si = scalings[i]
        end
        ratio = norm(δq̌i) / norm(q̌) 
        traj.q̌[end] .= q̌ .+ si.*δq̌i/ratio
    end
    bot
end

function make_N(tg::TensegrityStructure,q0::AbstractVector)
	(;bodies,connectivity) = tg
    (;nfree,nfull,syspres,sysfree,mem2sysfree,mem2sysincst,ninconstraints) = connectivity.indexed
    function inner_N(q̌)
        T = eltype(q̌)
		q = Vector{T}(undef,nfull)
		q[syspres] .= q0[syspres]
		q[sysfree] .= q̌
        ret = zeros(T,nfree,nfree-ninconstraints)
        foreach(bodies) do rb
            rbid = rb.prop.id
            (;lncs) = rb.state.cache.funcs
			memfree = mem2sysfree[rbid]
            if !isempty(mem2sysincst[rbid])
                if lncs isa NaturalCoordinates.LNC3D12C
                        u,v,w = NaturalCoordinates.get_uvw(lncs,q̌[memfree])
                        N = @view ret[mem2sysfree[rbid],mem2sysincst[rbid]]
                        N[1:3,1:3]   .= Matrix(1I,3,3)
                        N[4:6,4:6]   .= -NaturalCoordinates.skew(u)
                        N[7:9,4:6]   .= -NaturalCoordinates.skew(v)
                        N[10:12,4:6] .= -NaturalCoordinates.skew(w)
                elseif lncs isa NaturalCoordinates.LNC2D6C                    
                        u,v = NaturalCoordinates.get_uv(lncs,q̌[memfree])
                        N = @view ret[mem2sysfree[rbid],mem2sysincst[rbid]]
                        N[1:2,1:2] .= Matrix(1I,2,2)
                        N[3:4,3] .= -NaturalCoordinates.skew(u)
                        N[5:6,3] .= -NaturalCoordinates.skew(v)
                end
            end
        end
        ret
    end
end

function get_poly(bot_input;
        Ň
    )
    bot = deepcopy(bot_input)
    (;tg) = bot
    # (;ndof,nconstraints,connectivity) = bot.tg
    # (;cables) = tg.tensiles
    # (;nfull,nfree) = connectivity.indexed
    # ncables = length(cables)
    # nλ = nconstraints
    gue = get_initial(tg)
    Φ = make_Φ(tg,gue.q)
    A = make_A(tg,gue.q)
    Q̌ = make_Q̌(tg,gue.q)
    S = make_S(tg,gue.q)
    Ǩm_Ǩg = make_Ǩm_Ǩg(tg,gue.q)

    pv = get_polyvar(tg)

    pnq̌ = 1.0pv.q̌ .+ 0.0
    pns = 1.0pv.s .+ 0.0
    pnλ = 1.0pv.λ .+ 0.0
    pnd = 1.0pv.d .+ 0.0
    pnc = 1.0pv.c .+ 0.0
    pnk = 1.0pv.k .+ 0.0
    pnμ = 1.0pv.μ .+ 0.0
    polyΦ = Φ(pnq̌,pnd,pnc)
    polyA = A(pnq̌,pnc)
    polyQ̌ = Q̌(pnq̌,pns,pnμ,pnk,pnc)
    polyS = S(pnq̌,pns,pnc)
    polyQ̌a = transpose(polyA)*pnλ
    polyǨa = reduce(hcat,differentiate.(-polyQ̌a,Ref(pv.q̌))) |> transpose
    polyǨm, polyǨg = Ǩm_Ǩg(pnq̌,pns,pnμ,pnk,pnc)
    polyǨ = polyǨm .+ polyǨg .+ polyǨa
    polyŇ = Ň(pnq̌,pnc)
    poly𝒦 = transpose(polyŇ)*polyǨ*polyŇ

    polyP = [
        - polyQ̌ .- transpose(polyA)*pnλ ;
        polyS;
        polyΦ;
        # poly𝒦*pnξ.-pnζ.*pnξ;
        # transpose(pnξ)*pnξ-1;
    ]

    # Ǩ0 = TR.build_Ǩ(bot.tg,gue.λ)
    # Ǩx = map(polyǨ) do z
    # 		z(
    # 			pv.q̌=>gue.q̌,
    # 			pv.s=>gue.s,
    # 			pv.λ=>gue.λ,
    # 			pv.μ=>gue.μ,
    # 			pv.k=>gue.k,
    # 			pv.d=>gue.d,
    # 			pv.c=>gue.c
    # 		)
    # 	end
    # # @show Ǩ0
    # @show Ǩ0.- Ǩx |> norm

    # P0 = map(polyP) do z
    # 	z(
    # 		pvq̌=>q̌0,
    # 		pvs=>s0,
    # 		pvλ=>λ0,
    # 		# pvξ=>ξ0,
    # 		pvμ=>μ0,
    # 		pvk=>k0,
    # 	    pvd=>d0,
    # 		pvc=>c0,
    # 		# pv.ζ=>ζ0
    # 	)
    # end
    # @show P0[                 1:nfree] |> norm
    # @show P0[           nfree+1:nfree+ncables] |> norm
    # @show P0[   nfree+ncables+1:nfree+ncables+nλ] |> norm
    # @show P0[nfree+ncables+nλ+1:nfree+ncables+nλ+ndof]
    # @show P0[end]
    polyP,poly𝒦,gue,pv
end

function pinpoint(bot_input;
        Ň
    )
	polyP, poly𝒦, gue, pv = get_poly(bot_input;Ň)
	ň = length(pv.q̌)
	ns = length(pv.s)
	nλ = length(pv.λ)
    function make_bf()
        function inner_pp!(f,x)
            q̌x = @view x[        1:ň]
            sx = @view x[      ň+1:ň+ns]
            λx = @view x[   ň+ns+1:ň+ns+nλ]
			Px = map(polyP) do z
                z(
                    pv.q̌=>q̌x,
                    pv.s=>sx,
                    pv.λ=>λx,
                    pv.μ=>gue.μ,
                    pv.k=>gue.k,
                    pv.d=>gue.d,
                    pv.c=>gue.c,
                )
			end

			f .= Px
        end
    end
    f_holder = zeros(ň+ns+nλ)
    x_initial = vcat(gue.q̌,gue.s,gue.λ)
    pp! = make_bf()

    pp = nlsolve(pp!,x_initial,ftol=1e-10,iterations=100,method=:newton)
    # @show
    pp!(f_holder,pp.zero)
	# @show f_holder |> norm
    # @show f_holder[                 1:ň+ns+nλ] |> norm
    # @show f_holder[ň+ns+nλ+1:ň+ns+nλ+ndof] |> norm
    # @show f_holder[end]
    # @show  pp.zero[ň+ns+nλ+1:ň+ns+nλ+ndof]
    # @show  pp.zero[end]
    q̌ = pp.zero[        1:ň]
    s = pp.zero[      ň+1:ň+ns]
    λ = pp.zero[   ň+ns+1:ň+ns+nλ]
    ini = @eponymtuple(
			q̌,s,λ,
			isconverged=converged(pp),
			d=gue.d, c=gue.c, μ=gue.μ, k=gue.k
	)
	polyP, poly𝒦, ini, pv
end

function path_follow(bot_input;Ň)
	polyP, poly𝒦, ini, pv = pinpoint(bot_input;Ň)
	variable_groups = [pv.q̌,pv.s,pv.λ]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ]]
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
	Psys = System(polyP;parameters)
	result = HomotopyContinuation.solve(
			Psys,
			startsols;
			start_parameters,
			target_parameters,
			threading = false
	)
	path_results = results(result)
	if length(path_results) != 1
		@show failed(result)
		error("Tracking failed.")
	end
	path_result1 = path_results[1]
	sol = real(solution(path_result1))
	q̌,s,λ = split_by_lengths(sol,length.(variable_groups))
	@eponymtuple(q̌,s,λ)
end

function path_follow_critical(bot_input)
	polyP, ini, pv = pinpoint_critical(bot_input)
	variable_groups = [pv.q̌,pv.s,pv.λ,pv.ξ,[pv.ζ]]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ;ini.ξ;ini.ζ]]
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
	Psys = System(polyP;parameters)
	result = HomotopyContinuation.solve(
			Psys,
			startsols;
			start_parameters,
			target_parameters,
			threading = false
	)
	path_results = results(result)
	if length(path_results) != 1
		@show failed(result)
		error("Tracking failed.")
	end
	path_result1 = path_results[1]
	sol = real(solution(path_result1))
	q̌,s,λ,ξ,ζ = TR.split_by_lengths(sol,length.(variable_groups))
	@eponymtuple(q̌,s,λ,ξ,ζ)
end


