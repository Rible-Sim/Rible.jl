
# only bilateral
function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,mass_norm::Real,h::Real)
    pₖ = -pₖ₋₁ .+ 
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        mass_norm/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ
end

function Momentum_k!(p̌ₖ,p̌ₖ₋₁,qₖ,qₖ₋₁,λₖ,Ḿ,A,Aᵀₖ₋₁,mass_norm::Real,h::Real)
    p̌ₖ .= -p̌ₖ₋₁.+2/h.*Ḿ*(qₖ.-qₖ₋₁) .-
        mass_norm/h.*(transpose(A(qₖ))-Aᵀₖ₋₁)*λₖ
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

function get_state(bot::Robot,::DynamicsSolver{<:Zhong06})
    (;q,p,λ) = bot.structure.state.system
    ComponentArray(
        @eponymtuple(
            q,p,λ
        )
    )
end