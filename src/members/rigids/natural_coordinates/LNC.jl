"""
Data local natural coords。
"""
struct LNCData{N,M,T,L}
    r̄i::SArray{Tuple{N},T,1,N}
    X̄::SArray{Tuple{N,M},T,2,L}
    invX̄::SArray{Tuple{M,N},T,2,L}
end

function LNCData(r̄i::SVector,X̄::SMatrix)
    LNCData(r̄i,X̄,(pinv(X̄)))
end

const I2_Bool = IMatrix(2)
const I3_Bool = IMatrix(3)

"""
2D local natural coords of a point mass.
$(TYPEDEF)
"""
struct NC{N,M,T,L}
    np::Int64
    nv::Int64
    data::LNCData{N,M,T,L}
    conversion_to_std::SparseMatrixCSC{Int64,Int64}
    conversion_to_X::SparseMatrixCSC{Int64,Int64}
end


"""
Return the dimension of space.
$(TYPEDSIGNATURES)
"""
get_num_of_dims(::LNCData{N,M,T,L}) where {N,M,T,L} = N
get_num_of_dims(nmcs::NC) = get_num_of_dims(nmcs.data)

"""
Return local dimension of natural coodinates.
$(TYPEDSIGNATURES)
"""
get_num_of_local_dims(::LNCData{N,M,T,L}) where {N,M,T,L} = M
get_num_of_local_dims(nmcs::NC) = get_num_of_local_dims(nmcs.data)

"""
Return the number of coordinates.
$(TYPEDSIGNATURES)
"""
get_num_of_coords(::LNCData{N,M,T,L}) where {N,M,T,L} = N+L
get_num_of_coords(nmcs::NC) = get_num_of_coords(nmcs.data)

"""
Return the number of degrees of freedom
$(TYPEDSIGNATURES)
"""
get_num_of_dof(nmcs::NC) =  get_num_of_coords(nmcs) - get_num_of_cstr(nmcs)

function Base.getproperty(nmcs::NC,p::Symbol)
    if (p == :r̄i) 
        return nmcs.data.r̄i
    elseif (p == :ū) 
        return nmcs.data.X̄[:,1]
    elseif (p == :v̄)
        return nmcs.data.X̄[:,2]
    elseif (p == :w̄)
        return nmcs.data.X̄[:,3]
    elseif (p == :X̄)
        return nmcs.data.X̄
    elseif (p == :invX̄)
        return nmcs.data.invX̄
    elseif (p == :r̄j)
        return nmcs.data.r̄i + nmcs.data.ū
    elseif  (p == :r̄k)
        return nmcs.data.r̄i + nmcs.data.v̄
    elseif  (p == :r̄l)
        return nmcs.data.r̄i + nmcs.data.w̄
    else # fallback to getfield
        return getfield(nmcs, p)
    end
end

#  UnionAll types
const NC2D{M,T,L} = NC{2,M,T,L}
const NC3D{M,T,L} = NC{3,M,T,L}

const NC2D2C{T} = NC{2,0,T,0}
const NC2D4C{T} = NC{2,1,T,2}
const NC2D6C{T} = NC{2,2,T,4}
const NC3D3C{T} = NC{3,0,T,0}
const NC3D6C{T} = NC{3,1,T,3}
const NC3D12C{T} = NC{3,3,T,9}

"""
Return the number of constraints.
$(TYPEDSIGNATURES)
"""
get_num_of_cstr(::NC2D2C) = 0
get_num_of_cstr(::NC3D3C) = 0
get_num_of_cstr(::NC2D4C) = 1
get_num_of_cstr(::NC2D6C) = 3
get_num_of_cstr(::NC3D6C) = 1
get_num_of_cstr(::NC3D12C) = 6

"""
！！！！！！！！！
$(TYPEDSIGNATURES)
"""
function get_conversion_to_std(ndim,np,nv)
    nb = np+nv
    B = Matrix(1I,nb,nb)
    for i = 2:np
        B[i,1] = -1
    end
    kron(
        sparse(B),
        IMatrix(ndim)
    )
end

function get_conversion_to_X(ndim,nld,np,nv)
    nb = np+nv
    B = Matrix(1I,nb,nb)
    B[1] = 0
    for i = 2:np
        B[i,1] = -1
    end
    kron(
        sparse(B),
        IMatrix(ndim)
    )
end

function NC(np,nv,ndim,nld,r̄i::AbstractVector{T},X̄) where T
    NC(
        np,
        nv,
        LNCData(
            SVector{ndim}(r̄i),
            SMatrix{ndim,nld,T}(X̄)
        ),
        get_conversion_to_std(ndim,np,nv),
        get_conversion_to_X(ndim,nld,np,nv)
    )
end

"""
Return 2D point mass natural coodinates.
$(TYPEDSIGNATURES)
"""
function NC2D1P(ri::AbstractVector{T}) where T
    r̄i = @SVector zeros(T,2)
    np = 1
    nv = 0
    ndim = 2
    nld = 0
    NC(
        np,
        nv,
        LNCData(
            SVector{ndim}(r̄i),
            SMatrix{ndim,nld,T}()
        ),
        get_conversion_to_std(ndim,np,nv),
        get_conversion_to_X(ndim,nld,np,nv)
    )
end

"""
Return 2D rigid bar natural coodinates.
$(TYPEDSIGNATURES)
"""
function NC2D1P1V(ri::AbstractVector{T},u::AbstractVector{T},
               origin_position=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
               ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    ū = invR*u
    X̄ = ū
    np = 1
    nv = 1
    ndim = 2
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄)
end

function NC2D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
    ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    r̄j = invR*(rj-origin_position) #
    ū = r̄j-r̄i
    X̄ = ū
    np = 2
    nv = 0
    ndim = 2
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 2D rigid bodies natural coodinates, using 1 basic point and 2 base vectors.
$(TYPEDSIGNATURES)
"""
function NC1P2V(ri::AbstractVector{T},
               origin_position=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
               ) where T
    o = zero(T)
    i = one(T)
    ū = SVector(i,o)
    v̄ = SVector(o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    u = R*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    np = 1
    nv = 2
    ndim = 2
    nld = 2
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 2D rigid bodies natural coodinates, using 2 basic points and 1 base vector.
$(TYPEDSIGNATURES)
"""
function NC2P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                origin_position=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    ū = r̄j-r̄i
    v̄ = rotation_matrix(π/2)*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    np = 2
    nv = 1
    ndim = 2
    nld = 2
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 2D rigid bodies natural coodinates, using 3 basic points
$(TYPEDSIGNATURES)
"""
function NC3P(ri::AbstractVector{T},rj::AbstractVector{T},rk::AbstractVector{T},
              origin_position=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    r̄k = invR*(rk-origin_position)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    X̄ = hcat(ū,v̄)
    np = 3
    nv = 0
    ndim = 2
    nld = 2
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 3D point mass natural coodinates.
$(TYPEDSIGNATURES)
"""
function NC3D1P(ri::AbstractVector{T}) where T
    r̄i = @SVector zeros(T,3)
    np = 1
    nv = 0
    ndim = 3
    nld = 0
    NC(
        np,
        nv,
        LNCData(
            SVector{ndim}(r̄i),
            SMatrix{ndim,nld,T}()
        ),
        get_conversion_to_std(ndim,np,nv),
        get_conversion_to_X(ndim,nld,np,nv)
    )
end

"""
Return 3D rigid bar natural coodinates.
$(TYPEDSIGNATURES)
"""
function NC3D1P1V(ri::AbstractVector{T},u::AbstractVector{T},
               origin_position=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
               ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) 
    ū = invR*u
    X̄ = ū
    np = 1
    nv = 1
    ndim = 3
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄)
end

function NC3D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
    ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    r̄j = invR*(rj-origin_position)
    ū = r̄j-r̄i
    X̄ = ū
    np = 1
    nv = 0
    ndim = 3
    nld = 0
    NC(np,nv,ndim,nld,r̄i,X̄)
end

## 3D Rigid body
"""
Return 3D rigid bodies natural coodinates, using 1 basic point and 3 base vectors
$(TYPEDSIGNATURES)
"""
function NC1P3V(ri::AbstractVector{T},
                origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                ) where T
    o = zero(T)
    i = one(T)
    ū = SVector(i,o,o)
    v̄ = SVector(o,i,o)
    w̄ = SVector(o,o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    u = R*ū
    v = R*v̄
    w = R*w̄
    X̄ = hcat(ū,v̄,w̄)
    np = 1
    nv = 3
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄)
end

function NC1P3V(ri::AbstractVector{T},u::AbstractVector{T},v::AbstractVector{T},w::AbstractVector{T},
                origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                ) where T
    o = zero(T)
    i = one(T)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    np = 1
    nv = 3
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 3D rigid bodies natural coodinates, using 2 basic points and 2 base vectors.
$(TYPEDSIGNATURES)
"""
function NC2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    u = rj-ri
    ū = invR*u
    v̄,w̄ = HouseholderOrthogonalization(ū./norm(ū)).*norm(ū)
    v = R*v̄
    w = R*w̄
    X̄ = hcat(ū,v̄,w̄)
    np = 2
    nv = 2
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄)
end

function NC2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                v::AbstractVector{T},w::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    u = rj-ri
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    np = 2
    nv = 2
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 3D rigid bodies natural coodinates, using 3 basic points and 1 base vector.
$(TYPEDSIGNATURES)
"""
function NC3P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    u = rj-ri
    v = rk-ri
    w = u×v
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    r̄k = invR*(rk-origin_position)
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    np = 3
    nv = 1
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄)
end

"""
Return 3D rigid bodies natural coodinates, using 4 basic points
$(TYPEDSIGNATURES)
"""
function NC4P(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},rl::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    r̄k = invR*(rk-origin_position)
    r̄l = invR*(rl-origin_position)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    w̄ = r̄l - r̄i
    X̄ = hcat(ū,v̄,w̄)
    np = 4
    nv = 0
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄)
end
