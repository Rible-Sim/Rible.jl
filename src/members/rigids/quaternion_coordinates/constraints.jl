
function make_cstr_function(::QC,cstr_idx)
    function inner_cstr_function(x::AbstractVector)
        q = @view x[4:7]
        (transpose(q)*q - 1)/2
    end
end

function make_cstr_jacobian(::QC,free_idx,cstr_idx)
    function inner_cstr_jacobian(x::AbstractVector)
        q = @view x[4:7]
        o = zero(eltype(q))
        SA[
            o o o q[1] q[2] q[3] q[4];
        ]
    end
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