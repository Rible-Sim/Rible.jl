
function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end


function lucompletepiv!(A)
    n=size(A, 1)
    rowpiv=zeros(Int, n)
    colpiv=zeros(Int, n)
    for k=1:n
        Asub = abs.(A[k:n, k:n])#Search for next pivot
        _, index_max = findmax(Asub)
        μ,λ = index_max.I
        μ += k-1; λ += k-1
        rowpiv[k] = μ
        A[[k,μ], 1:n] = A[[μ,k], 1:n]
        colpiv[k] = λ
        A[1:n, [k,λ]] = A[1:n, [λ,k]]
        if A[k,k]≠0
            ρ = k+1:n
            A[ρ,k] = A[ρ,k]./A[k,k]
            A[ρ,ρ] = A[ρ,ρ] - A[ρ,k:k]*A[k:k,ρ]
        end
    end
    return (rowpiv, colpiv)
end

rotation_matrix(θ::Real) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]


"""
Return skew matrix
$(TYPEDSIGNATURES)
"""
function skew(w::AbstractVector)
    w1,w2,w3 = w
    o = zero(w1)
    @SMatrix [o -w3 w2;
              w3 o -w1;
             -w2 w1 o]
end

function skew(a::SVector{2})
	[-a[2],a[1]]
end

@inline @inbounds function LinearAlgebra.cross(a::Number,b::AbstractVector)
    ret = similar(b)
    ret[1] = -a*b[2]
    ret[2] =  a*b[1]
    ret
end

@inline @inbounds function LinearAlgebra.cross(a::StaticArray{Tuple{1}, T, 1},b::StaticArray{Tuple{2}, T, 1}) where T
    ret = similar(b)
    ret[1] = -a[1]*b[2]
    ret[2] =  a[1]*b[1]
    ret
end


@inline @inbounds function HouseholderOrthogonalization(n)
    n̂ = norm(n)
    _, i = findmax(n.+n̂)
    ei = 1:3 .== i
    h = n + n̂*ei
    H = I - 2h*transpose(h)/(transpose(h)*h)
    u, v  = [H[:,j] for j = 1:3 if j != i]
    if i == 2
        u, v = v, u
    end
    u, v
end

function get_orthonormal_axes(normal::AbstractVector)
    normal /= norm(normal)
    tangent, bitangent = HouseholderOrthogonalization(normal)
    SMatrix{3,3}(
        normal[1], normal[2], normal[3],
        tangent[1], tangent[2], tangent[3],
        bitangent[1], bitangent[2], bitangent[3],
    )
end
