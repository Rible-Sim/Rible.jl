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

cstr_forces_jacobian(λ::AbstractVector) = cstr_forces_jacobian(first(λ))

function cstr_forces_jacobian(λ)
    o = zero(λ)    
    Diagonal(SA[o,o,o,λ,λ,λ,λ])
end

# wrong
function make_cstr_hessians(::QC{T}) where T
    [
        @SMatrix zeros(T,7,7)
    ]
end

function make_cstr_forces_jacobian(::QC,free_idx,cstr_idx,cstr_hessians)
    function cstr_forces_jacobian(λ)
        ret = [
            begin
                a = -λ[i] .* cstr_hessians[j][free_idx,free_idx]
                # display(a)
                a 
            end
            for (i,j) in enumerate(cstr_idx)
        ]
        sum(ret)
    end
end

function find_independent_idx(::QC,q)
    [1,2,3,5,6,7]
end