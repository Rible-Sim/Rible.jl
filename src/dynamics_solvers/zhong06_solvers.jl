

# ============================================================================
# Data Structures
# ============================================================================

struct Zhong06 <: AbstractIntegrator end


# only bilateral
function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,mass_norm::Real,h::Real)
    pₖ = -pₖ₋₁ .+
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        mass_norm/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ
end

function Momentum_k!(p̌ₖ,p̌ₖ₋₁,qₖ,qₖ₋₁,λₖ,Ḿ,A::Function,Aᵀₖ₋₁,mass_norm::Real,h::Real)
    p̌ₖ .= -p̌ₖ₋₁.+2/h.*Ḿ*(qₖ.-qₖ₋₁) .+
        mass_norm/h.*(transpose(A(qₖ))-Aᵀₖ₋₁)*λₖ
end

function Momentum_k!(p̌ₖ,p̌ₖ₋₁,qₖ,qₖ₋₁,λₖ,Ḿ,structure::AbstractStructure,Aᵀₖ₋₁,mass_norm::Real,h::Real)
    # p̌ₖ .= -p̌ₖ₋₁.+2/h.*Ḿ*(qₖ.-qₖ₋₁) .+
    #     mass_norm/h.*(transpose(cstr_jacobian(structure, qₖ))-Aᵀₖ₋₁)*λₖ
end


# + unilateral
function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,mass_norm::Real,h::Real)
    pₖ = -pₖ₋₁ .+
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        mass_norm/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ .+
        mass_norm.*(transpose(Dₖ)-transpose(Dₖ₋₁))*H*Λₘ
end

# + unilateral scaling
function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,mass_norm::Real,scalingΛ::Real,h::Real)
    pₖ = -pₖ₋₁ .+
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        mass_norm/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ .+
        mass_norm*scalingΛ/h.*(transpose(Dₖ)-transpose(Dₖ₋₁))*H*Λₘ
end

# nonholonomic 
function Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A::Function,B::Function,h::Real)
    pᵏ = -pᵏ⁻¹ .+ 2/h.*M*(qᵏ.-qᵏ⁻¹) .+
        1/(2h).*(transpose(A(qᵏ))-transpose(A(qᵏ⁻¹)))*λᵏ .+
        1/(2h).*(transpose(B(qᵏ))-transpose(B(qᵏ⁻¹)))*μᵏ
end

abstract type AbstractZhong06Cache end

struct Zhong06_Constant_Mass_Cache{solT,T,opiontsType,state_midType} <: AbstractZhong06Cache
    solver::solT
    jacobian_workspace::Zhong06JacobianWorkspace{T}
    consts::Zhong06Constants{T}
    options::opiontsType
    state_mid::state_midType
end

struct Zhong06_CCP_Constant_Mass_Inner_Cache{
        T,
        RobotType,
        PolicyType,
        FieldType,
        EnvType,
        OptionsType,
    }  <: AbstractZhong06Cache
    bot::RobotType
    policy::PolicyType
    field::FieldType
    env::EnvType
    jacobian_workspace::Zhong06JacobianWorkspace{T}
    consts::Zhong06Constants{T}
    options::OptionsType
end

struct Zhong06_Nonconstant_Mass_Cache{
        T,
        RobotType,
        PolicyType,
        EnvType,
        OptionsType,
    } <: AbstractZhong06Cache
    bot::RobotType
    policy::PolicyType
    env::EnvType
    jacobian_workspace::Zhong06JacobianWorkspace{T}
    consts::Zhong06Constants{T}
    options::OptionsType
end

include("zhong06_res_and_jac.jl")

include("Zhong06_Constant_Mass_family/base/primal_solver.jl")
include("Zhong06_Constant_Mass_family/base/direct_solver.jl")
include("Zhong06_Constant_Mass_family/base/adjoint_solver.jl")

include("Zhong06_Constant_Mass_family/CCP_mono/primal_solver.jl")
include("Zhong06_Constant_Mass_family/CCP_mono/direct_solver.jl")
include("Zhong06_Constant_Mass_family/CCP_mono/adjoint_solver.jl")

include("Zhong06_Constant_Mass_family/CCP_inner/primal_solver.jl")

include("Zhong06_Nonconstant_Mass_family/base/primal_solver.jl")