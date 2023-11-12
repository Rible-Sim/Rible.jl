
vec(q::Quaternion) = SA[q.s, q.v1, q.v2, q.v3]


struct QC{T,JT}
    m::T
    m⁻¹::T
    γ::T
    Jγ::JT
    γ⁻¹::T
    J⁻¹γ::JT
end

function QC(m::T,J::AbstractMatrix{T};γ=maximum(diag(J))) where {T}
    m⁻¹ = inv(m)
    Jγ = Diagonal(J) - γ*I
    γ⁻¹ = inv(γ)
    J⁻¹γ = Diagonal(inv(J)) - γ⁻¹*I
    QC(m,m⁻¹,γ,Jγ,γ⁻¹,J⁻¹γ)
end

get_num_of_cstr(::QC) = 1
get_num_of_coords(::QC) = 7
get_num_of_dof(::QC) = 6
get_num_of_local_dims(::QC) = 3