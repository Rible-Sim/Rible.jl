
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

function split_by_lengths(x::AbstractVector, n::AbstractVector{<:Int})
    result = Vector{Vector{eltype(x)}}()
    start = firstindex(x)
    for len in n
      push!(result, x[start:(start + len - 1)])
      start += len
    end
    result
end

function split_by_lengths(x::AbstractVector{T}, len::Int) where {T}
    num_of_segments, remainder = divrem(length(x),len)
    istart = firstindex(x)
    [
        (@view x[istart+len*(i-1):istart+len*i-1])
        for i = 1:num_of_segments
    ]
end

"""
交换列。
$(TYPEDSIGNATURES)
"""
function swapcols!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k in axes(X,1)
        X[k,i], X[k,j] = X[k,j], X[k,i]
    end
end
"""
交换行。
$(TYPEDSIGNATURES)
"""
function swaprows!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k in axes(X,2)
        X[i,k], X[j,k] = X[j,k], X[i,k]
    end
end

function GECP(A_input)
    A = float(copy(A_input))
    n,m = size(A)
    col_index = collect(1:m)
    for k = 1:n
        Akrow2 = (@view A[k:end,k]).^2
        Akcol2 = (@view A[k,k:end]).^2
        ck_i = argmax(Akrow2)
        ck_j = argmax(Akcol2)
        # @show k,"before"
        # display(A)
        if Akrow2[ck_i] < Akcol2[ck_j]
            # swap columns
            swapcols!(A,k,k-1+ck_j)
            col_index[k],col_index[k-1+ck_j] = col_index[k-1+ck_j], col_index[k]
        else
            # swap rows
            swaprows!(A,k,k-1+ck_i)
        end
        # @show k,"after"
        # display(A)
        for i in k+1:n
            m_ik = A[i,k]/A[k,k]
            A[i,:] .-= m_ik*A[k,:]
        end
    end
    col_index
end


"""
ID
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct ID{sigType,pidType,aidType}
    "Signifier of body"
    bodysig::sigType
    "Index of the anchor point"
    pid::pidType
    "Index of the translational axis"
    trlid::aidType
    "Index of the rotational axis"
    rotid::aidType
end

function ID(bodysig,pid,trlid)
    ID(bodysig,pid,trlid,trlid)
end

function ID(bodysig,pid)
    ID(bodysig,pid,1,1)
end

"""
Hen 2 Egg
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct Hen2Egg{henType<:ID,eggType<:ID}
    "hen/parent/predecessor"
    hen::henType
    "egg/child/successor"
    egg::eggType
end


function get_ids(things)
    ids = mapreduce(get_id,vcat,things;init=Int[])
    nb = length(ids)
    ids,nb
end

function check_id_sanity(things)
    ids,nb = get_ids(things)
    @assert minimum(ids) == 1
    @assert maximum(ids) == nb
    @assert allunique(ids)
    ids,nb
end