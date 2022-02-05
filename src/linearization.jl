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

function ∂Aᵀλ∂q(tg,λ)
    body2q = tg.connectivity.body2q
    ncoords = tg.ncoords
    nbodyc = get_nbodyconstraint(tg)
    ret = zeros(eltype(λ),ncoords,ncoords)
    is = Ref(0)
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        if rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            # nc = rb.state.cache.nc
            # if nc > 0
            #     is += nc
            # end
            ret[pindex,pindex] .+= rb.state.cache.cfuncs.∂Aᵀλ∂q(λ[is[]+1:is[]+nbodyc])
            is[] += nbodyc
        end
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
    update_strings_apply_forces!(tg)
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

function normalize_wrt_mass!(Z,M)
    n = size(Z)[2]
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
