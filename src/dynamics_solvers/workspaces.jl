
# ============================================================================
# Workspace Structs for Efficient Argument Passing
# ============================================================================


struct GaugeWorkspace{T}
    d::Vector{T}
    Jq::Matrix{T}
    Jq̇::Matrix{T}
    Js::Matrix{T}
    tmp_grad::Vector{T}
    tmp_grad_s::Vector{T}
end

function GaugeWorkspace(T::Type{<:Number}, len::Int, nq::Int, ns::Int)
    return GaugeWorkspace{T}(
        zeros(T, len),
        zeros(T, len, nq),
        zeros(T, len, nq),
        zeros(T, len, ns),
        zeros(T, nq),
        zeros(T, ns)
    )
end

struct CostGradient{VT}
    ∂ϕ∂qᵀ::VT
    ∂ϕ∂q̇ᵀ::VT
    ∂ϕ∂pᵀ::VT
    ∂ϕ∂uᵀ::VT
    ∂ϕ∂sᵀ::VT
    ∂ϕ∂θᵀ::VT
    ∂ϕ∂cᵀ::VT
end

struct CostHessian{MT}
    ∂ϕ∂qᵀ∂q::MT
    ∂ϕ∂q̇ᵀ∂q̇::MT
    ∂ϕ∂qᵀ∂p::MT
    ∂ϕ∂pᵀ∂p::MT
    ∂ϕ∂qᵀ∂u::MT
    ∂ϕ∂q̇ᵀ∂u::MT
    ∂ϕ∂pᵀ∂u::MT
    ∂ϕ∂uᵀ∂u::MT
end

"""
Workspace for Newton iteration variables.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct NewtonWorkspace{VT,MT}
    x::VT           # Solution vector
    Res::VT         # Residual vector
    Jac::MT         # Jacobian matrix
    Δx::VT          # Newton step
    JacΔx::VT       # Jacobian * Δx workspace
    xₖ::VT          # Temporary solution vector
    𝐰::VT           # Work vector for contact
    ∂Γ∂x::MT        # Partial derivative matrix
    lu_tmp::MT      # Factorization scratch for Jac
    ipiv::Vector{LinearAlgebra.BlasInt} # Pivot workspace for LU
end

function NewtonWorkspace(x::VT, Res::VT, Jac::MT, Δx::VT, xₖ::VT, 𝐰::VT, ∂Γ∂x::MT) where {VT,MT}
    JacΔx = similar(Δx)
    lu_tmp = similar(Jac)
    ipiv = Vector{LinearAlgebra.BlasInt}(undef, size(Jac, 1))
    return NewtonWorkspace(x, Res, Jac, Δx, JacΔx, xₖ, 𝐰, ∂Γ∂x, lu_tmp, ipiv)
end

function NewtonWorkspace(::Type{T}, nx::Int, nΛ::Int, n2::Int) where {T}
    x = zeros(T, nx)
    Res = zeros(T, nx)
    Jac = zeros(T, nx, nx)
    Δx = zeros(T, nx)
    JacΔx = zeros(T, nx)
    xₖ = zeros(T, nx)
    𝐰 = zeros(T, nΛ)
    ∂Γ∂x = zeros(T, nΛ, n2)
    lu_tmp = zeros(T, nx, nx)
    ipiv = Vector{LinearAlgebra.BlasInt}(undef, nx)
    
    return NewtonWorkspace(x, Res, Jac, Δx, JacΔx, xₖ, 𝐰, ∂Γ∂x, lu_tmp, ipiv)
end


struct InnerContactWorkspace{T}
    Λ::Vector{T}
    Γ::Vector{T}
    Λʳ::Vector{T}
    ΔΛ::Vector{T}
    𝐁::Matrix{T}
    𝐁t::Matrix{T}
    𝐛::Vector{T}
    𝐜ᵀ::Matrix{T}
    𝐍::Matrix{T}
    𝐲::Vector{T}
end

function InnerContactWorkspace(T, nx, nΛ)
    return InnerContactWorkspace(
        zeros(T, nΛ), # Λ
        zeros(T, nΛ), # Γ
        zeros(T, nΛ), # Λʳ
        zeros(T, nΛ), # ΔΛ
        zeros(T, nx, nΛ), # 𝐁
        zeros(T, nx, nΛ), # 𝐁t
        zeros(T, nΛ), # 𝐛
        zeros(T, nΛ, nx), # 𝐜ᵀ
        zeros(T, nΛ, nΛ), # 𝐍
        zeros(T, nΛ) # 𝐲
    )
end


"""
Contact-specific variables for predictor-corrector method.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct ContactVariables{T,VT,SVT,SAT}
    Λ::SAT          # Contact force multipliers (SubArray)
    Γ::SAT          # Dual variables (SubArray)
    Λ_split::SVT    # Split view of Λ
    Γ_split::SVT    # Split view of Γ
    Λp::VT          # Predicted Λ
    Γp::VT          # Predicted Γ
    Λp_split::SVT   # Split view of Λp
    Γp_split::SVT   # Split view of Γp
    Δxp::VT         # Predictor step
    ΔΛp::SAT        # SubArray view of ΔΛp
    ΔΓp::SAT        # SubArray view of ΔΓp
    ΔΛp_split::SVT  # Split view of ΔΛp
    ΔΓp_split::SVT  # Split view of ΔΓp
    Δxc::VT         # Corrector step
    ΔΛc::SAT        # SubArray view of ΔΛc
    ΔΓc::SAT        # SubArray view of ΔΓc
    ΔΛc_split::SVT  # Split view of ΔΛc
    ΔΓc_split::SVT  # Split view of ΔΓc
    𝐞_split::Vector{SVector{3, T}}    # Unit vectors for centering
    J::Diagonal{T,SVector{3,T}}  # Cone metric
end


function ContactVariables(x::AbstractVector{T}, n2::Int, na::Int, nx::Int) where {T}
    nΛ = 3 * na
    
    # Create views into x for Λ and Γ
    Λ = @view x[(n2+1):n2+nΛ]
    Γ = @view x[n2+nΛ+1:n2+2nΛ]
    Λ_split = split_by_lengths(Λ, 3)
    Γ_split = split_by_lengths(Γ, 3)
    
    # Create predicted variables
    Λp = zero(Λ)
    Γp = zero(Γ)
    Λp_split = split_by_lengths(Λp, 3)
    Γp_split = split_by_lengths(Γp, 3)
    
    # Create predictor step variables
    Δxp = zeros(T, nx)
    ΔΛp = @view Δxp[(n2+1):n2+nΛ]
    ΔΓp = @view Δxp[n2+nΛ+1:n2+2nΛ]
    ΔΛp_split = split_by_lengths(ΔΛp, 3)
    ΔΓp_split = split_by_lengths(ΔΓp, 3)
    
    # Create corrector step variables
    Δxc = zeros(T, nx)
    ΔΛc = @view Δxc[(n2+1):n2+nΛ]
    ΔΓc = @view Δxc[n2+nΛ+1:n2+2nΛ]
    ΔΛc_split = split_by_lengths(ΔΛc, 3)
    ΔΓc_split = split_by_lengths(ΔΓc, 3)
    
    # Create centering unit vectors and cone metric
    𝐞_split = SVector{3,T}[SVector(one(T), zero(T), zero(T)) for i = 1:na]
    J = Diagonal(SVector(one(T), -one(T), -one(T)))
    
    return ContactVariables(
        Λ, Γ, Λ_split, Γ_split,
        Λp, Γp, Λp_split, Γp_split,
        Δxp, ΔΛp, ΔΓp, ΔΛp_split, ΔΓp_split,
        Δxc, ΔΛc, ΔΓc, ΔΛc_split, ΔΓc_split,
        𝐞_split, J
    )
end


"""
Workspace for Jacobian computations - contains all intermediate matrices needed.
Unified workspace used by Primal, Adjoint, and Direct sensitivity solvers.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Zhong06JacobianWorkspace{T}
    # Force Jacobians
    Fₘ::Vector{T}
    ∂F∂q::Matrix{T}
    ∂F∂q̇::Matrix{T}
    ∂Fₘ∂u::Matrix{T}
    ∂Fₘ∂c::Matrix{T}
    ∂F∂s::Matrix{T}
    
    # Control Jacobians
    ∂C∂qₖ::Matrix{T}
    ∂C∂pₖ::Matrix{T}
    ∂C∂sₖ::Matrix{T}
    ∂C∂qₖ₊₁::Matrix{T}
    ∂C∂pₖ₊₁::Matrix{T}
    
    # Mass matrices
    Mₘ::SparseMatrixCSC{T,Int}
    M⁻¹ₘ::SparseMatrixCSC{T,Int}
    ∂Mₘhq̇ₘ∂qₘ::SparseMatrixCSC{T,Int}
    M::SparseMatrixCSC{T,Int}
    Ḿ::SparseMatrixCSC{T,Int}
    M̌::SparseMatrixCSC{T,Int}
    M̄::SparseMatrixCSC{T,Int}
    M̌⁻¹::SparseMatrixCSC{T,Int}
    
    # Constraint Jacobians
    Aₖ₊₁::Matrix{T}
    Aₖ::Matrix{T}
    ∂Aᵀλ∂q::Matrix{T}
    ϕbuf::Vector{T}
    # Auxiliary Jacobians
    ∂S∂q::Matrix{T}
    ∂S∂s::Matrix{T}
    
    # Adjoint-specific gradients (for cost functions)
    ∂ϕ∂qᵀ::Vector{T}
    ∂ϕ∂q̇ᵀ::Vector{T}
    ∂ϕ∂pᵀ::Vector{T}
    ∂ϕ∂uᵀ::Vector{T}
    ∂ϕ∂sᵀ::Vector{T}
    ∂ϕf∂qᵀ::Vector{T}
    ∂ϕf∂q̇ᵀ::Vector{T}
    ∂ϕf∂pᵀ::Vector{T}
    ∂ϕf∂uᵀ::Vector{T}
    ∂ϕf∂sᵀ::Vector{T}
    cost_∂g∂q::Matrix{T}
    cost_∂g∂q̇::Matrix{T}
    cost_∂g∂s::Matrix{T}
    cost_∂g∂u::Matrix{T}
    cost_tmp_vec::Vector{T}
    gauge_workspaces::Vector{GaugeWorkspace{T}}
end

function Zhong06JacobianWorkspace(bot::Robot)
    (;structure, hub) = bot
    T = get_numbertype(structure)
    strip_sym(mat) = mat isa Symmetric ? mat.data : mat
    nq = get_num_of_free_coords(structure)
    nλ = get_num_of_cstr(structure)
    nu = get_num_of_actions(bot)
    nc = get_num_of_params(structure)
    ns = get_num_of_aux_var(structure)

    mass_mats = build_mass_matrices(structure)
    M = strip_sym(mass_mats.M)
    M⁻¹ = strip_sym(getfield(mass_mats, Symbol("M⁻¹")))
    M̌ = strip_sym(mass_mats.M̌)
    M̌⁻¹ = strip_sym(mass_mats.M̌⁻¹)
    Ḿ = strip_sym(mass_mats.Ḿ)
    M̄ = strip_sym(mass_mats.M̄)

    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)

    gauge_lengths = zeros(Int, hub.coalition.num_of_error_gauges)
    foreach(hub.error_gauges) do gauge
        gauge_lengths[gauge.id] = get_num_of_capta(gauge)
    end

    Zhong06JacobianWorkspace(
        T, nq, nλ, nu, nc, ns,
        M, M⁻¹, ∂Mₘhq̇ₘ∂qₘ,
        M, Ḿ, M̌, M̄, M̌⁻¹,
        gauge_lengths;
        num_error_gauges=hub.coalition.num_of_error_gauges,
        num_actuators=hub.coalition.num_of_actuators,
        num_actions=hub.coalition.num_of_actions,
    )
end

function Zhong06JacobianWorkspace(T::Type{<:Number}, 
        nq::Int, nλ::Int, nu::Int, nc::Int, ns::Int,
        Mₘ::SparseMatrixCSC{NumType,Int}, 
        M⁻¹ₘ::SparseMatrixCSC{NumType,Int}, 
        ∂Mₘhq̇ₘ∂qₘ::SparseMatrixCSC{NumType,Int}, 
        M::SparseMatrixCSC{NumType,Int}=spzeros(T,0,0),
        Ḿ::SparseMatrixCSC{NumType,Int}=spzeros(T,0,0),
        M̌::SparseMatrixCSC{NumType,Int}=spzeros(T,0,0),
        M̄::SparseMatrixCSC{NumType,Int}=spzeros(T,0,0),
        M̌⁻¹::SparseMatrixCSC{NumType,Int}=spzeros(T,0,0),
        gauge_lengths::Vector{Int}=Int[];

        num_error_gauges::Int=0,
        num_actuators::Int=nu,
        num_actions::Int=nu,
    ) where NumType <: Number
    return Zhong06JacobianWorkspace(
        zeros(T, nq),                                    # Fₘ 
        zeros(T, nq, nq),                                # ∂F∂q
        zeros(T, nq, nq),                                # ∂F∂q̇
        zeros(T, nq, nu),    # ∂Fₘ∂u
        zeros(T, nq, nc),    # ∂Fₘ∂c
        zeros(T, nq, ns),    # ∂F∂s
        zeros(T, nq, nq),                                # ∂C∂qₖ
        zeros(T, nq, nq),                                # ∂C∂pₖ
        zeros(T, nq, ns),                                # ∂C∂sₖ
        zeros(T, nq, nq),                                # ∂C∂qₖ₊₁
        zeros(T, nq, nq),                                # ∂C∂pₖ₊₁
        Mₘ,                                              # Mₘ
        M⁻¹ₘ,                                            # M⁻¹ₘ
        ∂Mₘhq̇ₘ∂qₘ,                                       # ∂Mₘhq̇ₘ∂qₘ
        M, Ḿ, M̌, M̄, M̌⁻¹,                                 # Mass variants
        zeros(T, nλ, nq),                                # Aₖ₊₁ (A is nλ x nq usually)
        zeros(T, nλ, nq),                                # Aₖ
        zeros(T, nq, nq),                                # ∂Aᵀλ∂q
        zeros(T, nλ),                                    # ϕbuf
        zeros(T, ns, nq),                                # ∂S∂q
        zeros(T, ns, ns),                                # ∂S∂s
        zeros(T, nq),                                    # ∂ϕ∂qᵀ
        zeros(T, nq),                                    # ∂ϕ∂q̇ᵀ
        zeros(T, nq),                                    # ∂ϕ∂pᵀ
        zeros(T, nu),                                    # ∂ϕ∂uᵀ
        zeros(T, ns),                                    # ∂ϕ∂sᵀ
        zeros(T, nq),                                    # ∂ϕf∂qᵀ
        zeros(T, nq),                                    # ∂ϕf∂q̇ᵀ
        zeros(T, nq),                                    # ∂ϕf∂pᵀ
        zeros(T, nu),                                    # ∂ϕf∂uᵀ
        zeros(T, ns),                                    # ∂ϕf∂sᵀ
        zeros(T, num_error_gauges, nq),                  # cost_∂g∂q
        zeros(T, num_error_gauges, nq),                  # cost_∂g∂q̇
        zeros(T, num_error_gauges, ns),                  # cost_∂g∂s
        zeros(T, num_actuators, num_actions),            # cost_∂g∂u
        zeros(T, nq),                                    # cost_tmp_vec
        GaugeWorkspace{T}[GaugeWorkspace(T, len, nq, ns) for len in gauge_lengths] # gauge_workspaces
    )
end


"""
Unified solver state containing all trajectory information needed for Jacobian computation.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Zhong06SolverState{InstStateType, T}
    # States
    state_k::InstStateType
    state_kp1::InstStateType
    state_mid::InstStateType

    # Time parameters (dt is useful to keep for quick access, but others are in states)
    dt::T
end

function Base.getproperty(solver_state::Zhong06SolverState, sym::Symbol)
    # Forward common properties to the relevant state object or internal field
    if sym == :qₖ
        return solver_state.state_k.q
    elseif sym == :q̇ₖ || sym == :vₖ
        return solver_state.state_k.q̇
    elseif sym == :pₖ
        return solver_state.state_k.p
    elseif sym == :sₖ
        return solver_state.state_k.s
    elseif sym == :tₖ
        return solver_state.state_k.t
    elseif sym == :q̌ₖ
        return solver_state.state_k.q̌
    elseif sym == :p̌ₖ
        return solver_state.state_k.p̌
    elseif sym == :q̌̇ₖ
        return solver_state.state_k.q̌̇
        
    elseif sym == :qₖ₊₁
        return solver_state.state_kp1.q
    elseif sym == :q̇ₖ₊₁ || sym == :vₖ₊₁
        return solver_state.state_kp1.q̇
    elseif sym == :pₖ₊₁
        return solver_state.state_kp1.p
    elseif sym == :sₖ₊₁
        return solver_state.state_kp1.s
    elseif sym == :tₖ₊₁
        return solver_state.state_kp1.t
    elseif sym == :q̌ₖ₊₁
        return solver_state.state_kp1.q̌
    elseif sym == :p̌ₖ₊₁
        return solver_state.state_kp1.p̌
    elseif sym == :Λₖ₊₁
        return solver_state.state_kp1.Λ
    elseif sym == :Γₖ₊₁
        return solver_state.state_kp1.Γ

    # state_mid λ is not used in Zhong06 
    elseif sym == :λₘ || sym == :λₖ₊₁
        return solver_state.state_kp1.λ
    elseif sym == :qₘ
        return solver_state.state_mid.q
    elseif sym == :q̇ₘ
        return solver_state.state_mid.q̇
    elseif sym == :sₘ
        return solver_state.state_mid.s
    elseif sym == :tₘ
        return solver_state.state_mid.t
    elseif sym == :F̌
        return solver_state.state_mid.F̌

    # Forward free coordinates for PresFreeCoordinatesState if available
    elseif sym == :q̌ₖ
        return solver_state.state_k.q̌
    elseif sym == :q̌ₖ₊₁
        return solver_state.state_kp1.q̌
    elseif sym == :q̌̇ₖ₊₁
        return solver_state.state_kp1.q̌̇
    elseif sym == :q̌ₘ
        return solver_state.state_mid.q̌
    else
        return getfield(solver_state, sym)
    end
end
"""
Constants for Jacobian computation.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Zhong06Constants{T}
    h::T
    mass_norm::T
    nq::Int
    nq̌::Int
    nλ::Int
    nu::Int
    nc::Int
    nθ::Int
    ns::Int
    n1::Int # nq̌
    n2::Int # 2nq̌
    n3::Int # 2nq̌ + nλ
end

function Zhong06Constants(bot::Robot,policy::AbstractPolicy,structure::AbstractStructure,mass_norm::Number,h::Number)
    (;structure,) = bot
    nq = get_num_of_full_coords(structure)
    nq̌ = get_num_of_free_coords(structure)
    nλ = get_num_of_cstr(structure)
    nu = get_num_of_actions(bot)
    nc = get_num_of_params(structure)
    nθ = get_num_of_params(policy)
    ns = get_num_of_aux_var(structure)
    n1 = nq̌
    n2 = 2nq̌
    n3 = n2 + nλ
    return Zhong06Constants(h, mass_norm, nq, nq̌, nλ, nu, nc, nθ, ns, n1, n2, n3)
end

"""
Jacobian blocks for the system at time k+1.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Zhong06JacobianBlocks{MatType,BackupMatType,CMatType}
    Jacᵏ⁺¹ₖ₊₁::MatType
    Jacᵏ⁺¹ₖ::MatType
    Jacᵏ⁺¹ₖ_backup::BackupMatType
    Jacᵏ⁺¹ₘu::MatType
    Jacᵏ⁺¹ₘc::CMatType
end

function Zhong06JacobianBlocks(T::Type{<:Number},nx::Int, nu::Int, nc::Int, Jacᵏ⁺¹ₖ_backup)
    return Zhong06JacobianBlocks(
        zeros(T, nx, nx),  # Jacᵏ⁺¹ₖ₊₁
        zeros(T, nx, nx),  # Jacᵏ⁺¹ₖ
        Jacᵏ⁺¹ₖ_backup,    # Jacᵏ⁺¹ₖ_backup
        zeros(T, nx, nu),     # Jacᵏ⁺¹ₘu
        zeros(T, nx, nc)      # Jacᵏ⁺¹ₘc
    )
end
"""
Workspace for storing cost gradients and Hessians.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct CostGradientHessianWorkspace{VT,MT}
    # Gradients
    ∂ϕ∂xᵀ::VT
    gradient::CostGradient{VT}
    
    # Hessians
    ∂ϕ∂xᵀ∂x::MT
    hessian::CostHessian{MT}
end

function CostGradientHessianWorkspace(T::Type{<:Number}, n3::Int, nq::Int, ns::Int, nu::Int, nθ::Int, nc::Int)
    hessian = CostHessian(
        zeros(T, nq, nq),    # ∂ϕ∂qᵀ∂q
        zeros(T, nq, nq),    # ∂ϕ∂q̇ᵀ∂q̇
        zeros(T, nq, nq),    # ∂ϕ∂qᵀ∂p
        zeros(T, nq, nq),    # ∂ϕ∂pᵀ∂p
        zeros(T, nq, nu),    # ∂ϕ∂qᵀ∂u
        zeros(T, nq, nu),    # ∂ϕ∂q̇ᵀ∂u
        zeros(T, nq, nu),    # ∂ϕ∂pᵀ∂u
        zeros(T, nu, nu),    # ∂ϕ∂uᵀ∂u
    )
    gradient = CostGradient(
        zeros(T, nq),        # ∂ϕ∂qᵀ
        zeros(T, nq),        # ∂ϕ∂q̇ᵀ
        zeros(T, nq),        # ∂ϕ∂pᵀ
        zeros(T, nu),        # ∂ϕ∂uᵀ
        zeros(T, ns),        # ∂ϕ∂sᵀ
        zeros(T, nθ),        # ∂ϕ∂θᵀ
        zeros(T, nc),        # ∂ϕ∂cᵀ
    )
    return CostGradientHessianWorkspace(
        zeros(T, n3),        # ∂ϕ∂xᵀ
        gradient,
        zeros(T, n3, n3),    # ∂ϕ∂xᵀ∂x
        hessian
    )
end

function clear!(cost_workspace::CostGradientHessianWorkspace)
    cost_workspace.∂ϕ∂xᵀ .= 0
    cost_workspace.∂ϕ∂xᵀ∂x .= 0
    cost_workspace.gradient.∂ϕ∂qᵀ .= 0
    cost_workspace.gradient.∂ϕ∂q̇ᵀ .= 0
    cost_workspace.gradient.∂ϕ∂pᵀ .= 0
    cost_workspace.gradient.∂ϕ∂uᵀ .= 0
    cost_workspace.gradient.∂ϕ∂sᵀ .= 0
    cost_workspace.hessian.∂ϕ∂qᵀ∂q .= 0
    cost_workspace.hessian.∂ϕ∂q̇ᵀ∂q̇ .= 0
    cost_workspace.hessian.∂ϕ∂pᵀ∂p .= 0
    cost_workspace.hessian.∂ϕ∂qᵀ∂p .= 0
    cost_workspace.hessian.∂ϕ∂qᵀ∂u .= 0
    cost_workspace.hessian.∂ϕ∂q̇ᵀ∂u .= 0
    cost_workspace.hessian.∂ϕ∂pᵀ∂u .= 0
    cost_workspace.hessian.∂ϕ∂uᵀ∂u .= 0
end

struct DirectSensitivityWorkspace{T}
    Jac_state::Vector{Matrix{T}}
    Jac_action::Vector{Matrix{T}}
    Jac_control_params::Vector{Matrix{T}}
    traj_cost_gradients_wrt_state::Vector{Vector{T}}
    traj_cost_hessians_wrt_state::Vector{Matrix{T}}
    traj_cost_gradients_wrt_action::Vector{Vector{T}}
    traj_cost_hessians_wrt_action::Vector{Matrix{T}}
    term_cost_gradient_wrt_state::Vector{T}
    term_cost_hessian_wrt_state::Matrix{T}
    term_cost_gradient_wrt_action::Vector{T}
    term_cost_hessian_wrt_action::Matrix{T}
end

function DirectSensitivityWorkspace(T::Type{<:Number}, n3::Int, nu::Int)
    return DirectSensitivityWorkspace(
        Matrix{T}[], # Jac_state
        Matrix{T}[], # Jac_action
        Matrix{T}[], # Jac_control_params
        Vector{T}[], # traj_cost_gradients_wrt_state
        Matrix{T}[], # traj_cost_hessians_wrt_state
        Vector{T}[], # traj_cost_gradients_wrt_action
        Matrix{T}[], # traj_cost_hessians_wrt_action
        zeros(T, n3), # term_cost_gradient_wrt_state
        zeros(T, n3, n3), # term_cost_hessian_wrt_state
        zeros(T, nu), # term_cost_gradient_wrt_action
        zeros(T, nu, nu) # term_cost_hessian_wrt_action
    )
end



@inline function _lu_solve_from_jacobian!(ws::NewtonWorkspace)
    lu_tmp = ws.lu_tmp
    ipiv = ws.ipiv
    copyto!(lu_tmp, ws.Jac)
    @. ws.Δx = -ws.Res
    (_, _, info) = LinearAlgebra.LAPACK.getrf!(lu_tmp, ipiv; check=false)
    @assert info == 0 "LU factorization failed with info=$info"
    LinearAlgebra.LAPACK.getrs!('N', lu_tmp, ipiv, ws.Δx)
    nothing
end

"""
Linear solver that uses the analytic Jacobian.

Expected signature for custom solvers:
    linear_solver!(ws::NewtonWorkspace, solver_state, solver_cache, bot, policy, env)
"""
function default_linear_solver!(ws::NewtonWorkspace,
        solver_state, solver_cache, bot, policy, env)
    # Assumes ws.Res is already the residual at ws.x
    compute_constant_mass_jacobian!(ws.Jac, ws.x, solver_state, solver_cache, bot, policy, env)
    _lu_solve_from_jacobian!(ws)
end

"""
Linear solver that builds the Jacobian with finite differences of the residual.

`residual_func!` must populate `ws.Res` using the values in `ws.x`.
"""
function finite_diff_linear_solver!(ws::NewtonWorkspace,
        solver_state, solver_cache, bot, policy, env;
        fdtype::Type{<:Val}=Val{:central})
    # Reuse xₖ as scratch to avoid an extra allocation for the backup.
    copyto!(ws.xₖ, ws.x)

    fd_residual! = function(res_out, xvec)
        ws.x .= xvec
        compute_constant_mass_residual!(res_out, ws.x, solver_state, solver_cache, bot, policy, env)
        nothing
    end

    FiniteDiff.finite_difference_jacobian!(ws.Jac, fd_residual!, ws.x, fdtype)

    # Restore the base iterate and residual before solving
    copyto!(ws.x, ws.xₖ)
    compute_constant_mass_residual!(ws.Res, ws.x, solver_state, solver_cache, bot, policy, env)
    _lu_solve_from_jacobian!(ws)
end