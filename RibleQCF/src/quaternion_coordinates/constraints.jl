
function nullspace_mat(::QC,x::AbstractVector{T}) where{T}
    q = @view x[4:7]
    O43 = @SMatrix zeros(T,4,3)
    O33 = @SMatrix zeros(T,3,3)
    hcat(
        vcat(
            SMatrix{3,3,T}(I(3)),
            O43
        ),
        vcat(
            O33,
            Lᵀmat(q)
        )
    )

end

function get_deform(::QC{T}) where T
    one(T)
end

function cstr_function!(ret, qcs::QC, x::AbstractVector, d=get_deform(qcs))
    q = @view x[4:7]
    ret[begin] = (transpose(q)*q - d)/2
end

function cstr_jacobian!(ret, ::QC, jac, x::AbstractVector)
    q = @view x[4:7]
    o = zero(eltype(q))
    ret .=  SA[
        o o o q[1] q[2] q[3] q[4];
    ]
end

function cstr_forces_jacobian!(ret, ::QC, x::AbstractVector, λ::AbstractVector)
    fill!(ret, zero(eltype(ret)))
    λval = λ[begin]
    @inbounds for i in 4:7
        ret[i, i] += λval
    end
end

function add_cstr_forces_jacobian!(ret, coords::QC, λ)
    λval = λ[begin]
    @inbounds for i in 4:7
        ret[i, i] += λval
    end
end

function find_independent_free_idx(::QC,q)
    [1,2,3,5,6,7]
end
