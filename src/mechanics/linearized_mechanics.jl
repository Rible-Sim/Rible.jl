
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

function find_finite(ω2,Z,num_of_dof)
    first_frequency_index = findfirst((x)->x>0,ω2)
    finite_ω2 = ω2[first_frequency_index:first_frequency_index+num_of_dof-1]
    finite_Z = Z[:,first_frequency_index:first_frequency_index+num_of_dof-1]
    finite_ω2,finite_Z
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

function  check_static_equilibrium_output_multipliers(st::AbstractStructure, field::AbstractField;)
    @error "Not implemented for AbstractStructure"
end

function undamped_eigen(st, field;)
    _, λ = check_static_equilibrium_output_multipliers(st, field;)
    q = get_coords(st)
    M̌ = assemble_M̌(st)
    Ǩ = build_Ǩ(st,q,λ)
    Ǎ = cstr_jacobian(st,st.state.system)
    Ň = nullspace(Ǎ)
    ℳ = transpose(Ň)*M̌*Ň
    𝒦 = transpose(Ň)*Ǩ*Ň
    # @show ℳ, 𝒦
    ω²,ξ = eigen(Symmetric(𝒦),Symmetric(ℳ))
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

function old_undamped_eigen(st)
    λ0 = check_static_equilibrium_output_multipliers(st)
    M̂,Ĉ,K̂ = linearize(st,q0,λ0)
    α = 10
    M̄,K̄ = frequencyshift(M̂,K̂,α)
    # @show size(K̄),rank(K̄),cond(K̄),rank(M̄)
    d,aug_Z = eigen(K̄,M̄)
    aug_ω2 = d .- α
    (;ncoords, num_of_dof) = st
    # @show aug_ω2
    ω2,Z = find_finite(aug_ω2,aug_Z,num_of_dof)
    ω = sqrt.(ω2)
    Zq = Z[1:ncoords,:]
    M = build_massmatrix(st)
    normalize_wrt_mass!(Zq,M)
    ω, Zq#, Z
end

function undamped_modal_solve!(st,q0,q̇0,λ0,tf,dt)
    M̂,Ĉ,K̂ = linearize(st,q0,λ0)
    # show(stdout,"text/plain",K̂)
    # showtable(K̂)
    # M̄,C̄,K̄ = RB.frequencyshift(M̂,Ĉ,K̂,0.1)
    # M̃,K̃ = RB.enlarge(M̄,C̄,K̄)
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

function undamped_eigen!(bot::Robot, field; scaling=0.01)
    (;structure,traj) = bot
    q̌ = get_free_coords(structure)
    ω², δq̌ = undamped_eigen(structure, field;)
    neg_idx = findall(ω².<=0)
    if !isempty(neg_idx)
        @warn "Negative ω² occurs, idx $neg_idx, zeroing."
        ω²[neg_idx] .= 0
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


"""
$(TYPEDSIGNATURES)
"""
function check_stability(st::AbstractStructure;F̌=nothing,verbose=false)
    static_equilibrium,λ = check_static_equilibrium_output_multipliers(st;F=F̌)
    @assert static_equilibrium
    check_stability(st,λ;verbose)
end

function check_stability(st::AbstractStructure,λ;verbose=false)
    q = get_coords(st)
    c = get_params(st)
    A = make_cstr_jacobian(st,q)
    Ň(q̌,c) = nullspace(A(q̌))
    check_stability(st,λ,Ň;verbose)
end

function check_stability(st::AbstractStructure,λ,Ň;verbose=false)
    q̌ = get_free_coords(st)
    c = get_params(st)
    Ǩ0 = build_Ǩ(st,λ)
    Ň0 = Ň(q̌,c)
    𝒦0 = transpose(Ň0)*Ǩ0*Ň0
    eigen_result = eigen(𝒦0)
    nn = count(x -> x < 0, eigen_result.values)
    if nn > 1
        @warn "Instability detected! Number of negative eigenvalues: $nn"
        isstable = false
    else
        isstable = true
    end
    isstable, Ň0, eigen_result
end