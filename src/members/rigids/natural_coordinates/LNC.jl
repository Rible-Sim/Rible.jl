
"""
local natural coords abstract type.
"""
abstract type LNC{T} end

"""
2D local natural coords abstract type.
"""
abstract type LNC2D{T} <: LNC{T} end

"""
3D local natural coords abstract type.
"""
abstract type LNC3D{T} <: LNC{T} end

"""
2 coords 2D local natural coords abstract type for point mass。
"""
abstract type LNC2D2C{T} <: LNC2D{T} end

"""
4 coords 2D local natural coords abstract type for rigid bars。
"""
abstract type LNC2D4C{T} <: LNC2D{T} end

"""
6 coords 2D local natural coords abstract type , for rigid bodies。
"""
abstract type LNC2D6C{T} <: LNC2D{T} end

"""
3 coords 3D local natural coords abstract type for point mass。
"""
abstract type LNC3D3C{T} <: LNC3D{T} end

"""
6 coords 3D local natural coords abstract type for rigid bars。
"""
abstract type LNC3D6C{T} <: LNC3D{T} end

"""
12 coords 3D local natural coords abstract type , for rigid bodies。
"""
abstract type LNC3D12C{T} <: LNC3D{T} end

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
Return the dimension of space.
$(TYPEDSIGNATURES)
"""
get_num_of_dims(nmcs::LNC) = get_num_of_dims(nmcs.data)
get_num_of_dims(::LNCData{N,M,T,L}) where {N,M,T,L} = N

"""
Return local dimension of natural coodinates.
$(TYPEDSIGNATURES)
"""
get_num_of_local_dims(nmcs::LNC) = get_num_of_local_dims(nmcs.data)
get_num_of_local_dims(::LNCData{N,M,T,L}) where {N,M,T,L} = M

"""
Return the number of coords。
$(TYPEDSIGNATURES)
"""
get_num_of_coords(nmcs::LNC) = get_num_of_coords(nmcs.data)
get_num_of_coords(::LNCData{N,M,T,L}) where {N,M,T,L} = N+L

"""
Return 约束方程个数。
$(TYPEDSIGNATURES)
"""
get_num_of_cstr(::LNC2D2C) = 0
get_num_of_cstr(::LNC3D3C) = 0
get_num_of_cstr(::LNC2D4C) = 1
get_num_of_cstr(::LNC2D6C) = 3
get_num_of_cstr(::LNC3D6C) = 1
get_num_of_cstr(::LNC3D12C) = 6

"""
Return 自由度数。
$(TYPEDSIGNATURES)
"""
get_num_of_dof(nmcs::LNC) =  get_num_of_coords(nmcs) - get_num_of_cstr(nmcs)

function Base.getproperty(nmcs::LNC,p::Symbol)
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


"""
3D local natural coords of a  point mass.
$(TYPEDEF)
"""
struct LNC3D1P{T} <: LNC3D3C{T}
    data::LNCData{3,0,T,0}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
2D local natural coords of a  point mass.
$(TYPEDEF)
"""
struct LNC2D1P{T} <: LNC2D2C{T}
    data::LNCData{2,0,T,0}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
4 coords 2D local natural coords, using 2 basic points.
$(TYPEDEF)
"""
struct LNC2D2P{T} <: LNC2D4C{T}
    data::LNCData{2,1,T,2}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
4 coords 2D local natural coords, using 1 basic point and 1 base vector.
$(TYPEDEF)
"""
struct LNC2D1P1V{T} <: LNC2D4C{T}
    data::LNCData{2,1,T,2}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coords 3D local natural coords, using 2 basic points.
$(TYPEDEF)
"""
struct LNC3D2P{T} <: LNC3D6C{T}
    data::LNCData{3,1,T,3}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coords 3D local natural coords, using 1 basic point and 1 base vector.
$(TYPEDEF)
"""
struct LNC3D1P1V{T} <: LNC3D6C{T}
    data::LNCData{3,1,T,3}
    conversion::SparseMatrixCSC{Int64,Int64}
end


"""
6 coords 2D local natural coords, using 1 basic point and 2 base vectors.
$(TYPEDEF)
"""
struct LNC1P2V{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coords 2D local natural coords, using 2 basic points and 1 base vector.
$(TYPEDEF)
"""
struct LNC2P1V{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coords 2D local natural coords, using 3 basic points.
$(TYPEDEF)
"""
struct LNC3P{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coords 3D local natural coords, using 1 basic point and 3 base vectors.
$(TYPEDEF)
"""
struct LNC1P3V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coords 3D local natural coords, using 2 basic points and 2 base vectors.
$(TYPEDEF)
"""
struct LNC2P2V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coords 3D local natural coords, using 3 basic points and 1 base vector.
$(TYPEDEF)
"""
struct LNC3P1V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coords 3D local natural coords, using 4 basic points.
$(TYPEDEF)
"""
struct LNC4P{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
！！！！！！！！！
$(TYPEDSIGNATURES)
"""
function get_conversion(ndim,np,nv)
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

"""
Return 3D point mass natural coodinates.
$(TYPEDSIGNATURES)
"""
function NC3D1P(ri::AbstractVector{T}) where T
    r̄i = @SVector zeros(T,3)
    nmcs = LNC3D1P(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,0,T}()
        ),
        get_conversion(3,1,0)
    )
    nmcs,ri
end

"""
Return 2D point mass natural coodinates.
$(TYPEDSIGNATURES)
"""
function NC2D1P(ri::AbstractVector{T}) where T
    r̄i = @SVector zeros(T,2)
    nmcs = LNC2D1P(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,0,T}()
        ),
        get_conversion(2,1,0)
    )
    nmcs,ri
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
    nmcs = LNC2D1P1V(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,1}(X̄)
        ),
        get_conversion(2,1,1)
    )
    q = vcat(ri,u)
    nmcs
end

function NC2D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
    ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    r̄j = invR*(rj-origin_position) #
    ū = r̄j-r̄i
    X̄ = ū
    nmcs = LNC2D2P(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,1}(X̄)
        ),
        get_conversion(3,2,2)
    )
    q = vcat(ri,rj)
    nmcs
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
    nmcs = LNC3D1P1V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,1}(X̄)
        ),
        get_conversion(3,1,1)
    )
    q = vcat(ri,u)
    nmcs
end

function NC3D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
    ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    r̄j = invR*(rj-origin_position)
    ū = r̄j-r̄i
    X̄ = ū
    nmcs = LNC3D2P(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,1}(X̄)
        ),
        get_conversion(3,2,0)
    )
    q = vcat(ri,rj)
    nmcs
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
    nmcs = LNC1P2V(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,2}(X̄)
        ),
        get_conversion(2,1,2)
    )
    q = vcat(ri,u,v)
    nmcs
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
    nmcs = LNC2P1V(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,2}(X̄)
        ),
        get_conversion(2,2,1)
    )
    q = vcat(ri,rj,v)
    nmcs
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
    nmcs = LNC3P(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,2}(X̄)
        ),
        get_conversion(2,3,0)
    )
    q = vcat(ri,rj,rk)
    nmcs
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
    nmcs = LNC1P3V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        ),
        get_conversion(3,1,3)
    )
    q = vcat(ri,u,v,w)
    nmcs
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
    nmcs = LNC1P3V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        ),
        get_conversion(3,1,3)
    )
    q = vcat(ri,u,v,w)
    nmcs
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
    nmcs = LNC2P2V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        ),
        get_conversion(3,2,2)
    )
    q = vcat(ri,rj,v,w)
    nmcs
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
    nmcs = LNC2P2V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        ),
        get_conversion(3,2,2)
    )
    q = vcat(ri,rj,v,w)
    nmcs
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
    nmcs = LNC3P1V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        ),
        get_conversion(3,3,1)
    )
    q = vcat(ri,rj,rk,w)
    nmcs
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
    nmcs = LNC4P(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        ),
        get_conversion(3,4,0)
    )
    q = vcat(ri,rj,rk,rl)
    nmcs
end


function get_uvw(nmcs::LNC3D,q)
    qstd = nmcs.conversion*q
    if     nmcs isa LNC3D12C
        u = @view qstd[4:6]
        v = @view qstd[7:9]
        w = @view qstd[10:12]
    elseif nmcs isa LNC3D6C
        u = @view qstd[4:6]
        u /= norm(u)
        v,w = HouseholderOrthogonalization(u)
    end
    SVector{3}(u),SVector{3}(v),SVector{3}(w)
end

function get_uv(nmcs::LNC2D,q)
    qstd = nmcs.conversion*q
    if     nmcs isa LNC2D6C
        u = @view qstd[3:4]
        v = @view qstd[5:6]
    elseif nmcs isa LNC2D4C
        u = @view qstd[3:4]
        v = rotation_matrix(π/2)*u
    end
    SVector{2}(u),SVector{2}(v)
end