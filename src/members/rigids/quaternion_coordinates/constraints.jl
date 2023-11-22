function get_deform(::QC{T}) where T
    one(T)
end

function cstr_function(qcs::QC,cstr_idx,x,d=get_deform(qcs))
    q = @view x[4:7]
    (transpose(q)*q - d)/2
end

function cstr_jacobian(::QC,free_idx,cstr_idx,x)
    q = @view x[4:7]
    o = zero(eltype(q))
    SA[
        o o o q[1] q[2] q[3] q[4];
    ]
end

cstr_forces_jacobian(qcs::QC,free_idx,cstr_idx,λ::AbstractVector) = cstr_forces_jacobian(qcs,free_idx,cstr_idx,first(λ))

function cstr_forces_jacobian(qcs::QC,free_idx,cstr_idx,λ::Number)
    o = zero(eltype(λ))    
    [Diagonal(SA[o,o,o,λ,λ,λ,λ])[free_idx,free_idx]][cstr_idx] |> sum
end

function find_independent_free_idx(::QC,q)
    [1,2,3,5,6,7]
end