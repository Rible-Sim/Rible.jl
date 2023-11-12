
function mass_matrix(J,q::AbstractVector)
    4Lᵀmat(q)*J*Lmat(q)
end

function get_Mγ(Jγ,γ,q::AbstractVector)
    4(Lᵀmat(q)*Jγ*Lmat(q) + γ*I)
end

function get_M⁻¹γ(J⁻¹γ,γ⁻¹,q::AbstractVector)
    (1//4)*(Lᵀmat(q)*J⁻¹γ*Lmat(q) + γ⁻¹*I)
end

function get_∂Mγq̇∂q(Jγ,q::AbstractVector,q̇::AbstractVector)
    L = Lmat(q)
    η = Jγ*L*q̇
    -4(transpose(L)*Jγ*Lmat(q̇) .+ ∂Lᵀη∂q(η))
end

function get_∂M⁻¹γp∂q(J⁻¹γ,q::AbstractVector,p::AbstractVector)
    L = Lmat(q)
    η = J⁻¹γ*L*p
    -1//4*(transpose(L)*J⁻¹γ*Lmat(p) .+ ∂Lᵀη∂q(η))
end

function get_∂Tγ∂qᵀ(Jγ,q::AbstractVector,q̇::AbstractVector)
    4Lᵀmat(q̇)*Jγ*Lmat(q̇)*q
end

function get_∂Tγ∂qᵀ∂q(Jγ,q̇::AbstractVector)
    4Lᵀmat(q̇)*Jγ*Lmat(q̇)
end


function make_M(qcs::QC)
    (;m,γ,Jγ) = qcs
    function inner_M(x)
        q = @view x[4:7]
        T = eltype(x)
        O34 = @SMatrix zeros(T,3,4)
        hcat(
            vcat(
                SMatrix{3,3}(m*I),
                transpose(O34)
            ),
            vcat(
                O34,
                get_Mγ(Jγ,γ,q)
            )
        )
    end
end

function make_M⁻¹(qcs::QC)
    (;m⁻¹,γ⁻¹,J⁻¹γ) = qcs
    function inner_M(x)
        q = @view x[4:7]
        T = eltype(q)
        O34 = @SMatrix zeros(T,3,4)
        hcat(
            vcat(
                SMatrix{3,3}(m⁻¹*I),
                transpose(O34)
            ),
            vcat(
                O34,
                get_M⁻¹γ(J⁻¹γ,γ⁻¹,q)
            )
        )
    end
end

function make_∂Mẋ∂x(qcs::QC)
    (;Jγ) = qcs
    function inner_∂Mẋ∂x(x,ẋ)
        q = @view x[4:7]
        q̇ = @view ẋ[4:7]
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                get_∂Mγq̇∂q(Jγ,q,q̇)
            )
        )
    end
end

function make_∂M⁻¹y∂x(qcs::QC)
    (;J⁻¹γ) = qcs
    function inner_∂M⁻¹y∂x(x,y)
        q = @view x[4:7]
        p = @view y[4:7]
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                get_∂M⁻¹γp∂q(J⁻¹γ,q,p)
            )
        )
    end
end


cstr_forces_jacobian(λ::AbstractVector) = cstr_forces_jacobian(first(λ))

function cstr_forces_jacobian(λ)
    o = zero(λ)    
    Diagonal(SA[o,o,o,λ,λ,λ,λ])
end


function make_∂T∂xᵀ(qcs::QC)
    (;Jγ) = qcs
    function inner_∂T∂xᵀ(x,ẋ)
        q = @view x[4:7]
        q̇ = @view ẋ[4:7]
        T = eltype(q)
        o3 = @SVector zeros(T,3)
        vcat(
            o3,
            get_∂Tγ∂qᵀ(Jγ,q,q̇)
        )            
    end
end

function make_∂T∂xᵀ∂x(qcs::QC)
    (;Jγ) = qcs
    function inner_∂T∂xᵀ∂x(ẋ)
        q̇ = @view ẋ[4:7]
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                get_∂Tγ∂qᵀ∂q(Jγ,q̇)
            )
        )
    end
end