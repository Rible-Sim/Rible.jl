module NCF
using LinearAlgebra
using StaticArrays
using SparseArrays
using LuxurySparse
using ForwardDiff
using DocStringExtensions
using SymmetricFormats

import ..Rible: HouseholderOrthogonalization

export get_num_of_constraints, get_num_of_coordinates, get_num_of_dof, get_num_of_local_dims

"""
local natural coordinates abstract type.
"""
abstract type LNC{T} end

"""
2D local natural coordinates abstract type.
"""
abstract type LNC2D{T} <: LNC{T} end

"""
3D local natural coordinates abstract type.
"""
abstract type LNC3D{T} <: LNC{T} end

"""
2 coordinates 2D local natural coordinates abstract type for point mass。
"""
abstract type LNC2D2C{T} <: LNC2D{T} end

"""
4 coordinates 2D local natural coordinates abstract type for rigid bars。
"""
abstract type LNC2D4C{T} <: LNC2D{T} end

"""
6 coordinates 2D local natural coordinates abstract type , for rigid bodies。
"""
abstract type LNC2D6C{T} <: LNC2D{T} end

"""
3 coordinates 3D local natural coordinates abstract type for point mass。
"""
abstract type LNC3D3C{T} <: LNC3D{T} end

"""
6 coordinates 3D local natural coordinates abstract type for rigid bars。
"""
abstract type LNC3D6C{T} <: LNC3D{T} end

"""
12 coordinates 3D local natural coordinates abstract type , for rigid bodies。
"""
abstract type LNC3D12C{T} <: LNC3D{T} end

"""
Data local natural coordinates。
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
Return the number of coordinates。
$(TYPEDSIGNATURES)
"""
get_num_of_coordinates(nmcs::LNC) = get_num_of_coordinates(nmcs.data)
get_num_of_coordinates(::LNCData{N,M,T,L}) where {N,M,T,L} = N+L

"""
Return 约束方程个数。
$(TYPEDSIGNATURES)
"""
get_num_of_constraints(::LNC2D2C) = 0
get_num_of_constraints(::LNC3D3C) = 0
get_num_of_constraints(::LNC2D4C) = 1
get_num_of_constraints(::LNC2D6C) = 3
get_num_of_constraints(::LNC3D6C) = 1
get_num_of_constraints(::LNC3D12C) = 6

"""
Return 自由度数。
$(TYPEDSIGNATURES)
"""
get_num_of_dof(nmcs::LNC) =  get_num_of_coordinates(nmcs) - get_num_of_constraints(nmcs)

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
3D local natural coordinates of a  point mass.
$(TYPEDEF)
"""
struct LNC3D1P{T} <: LNC3D3C{T}
    data::LNCData{3,0,T,0}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
2D local natural coordinates of a  point mass.
$(TYPEDEF)
"""
struct LNC2D1P{T} <: LNC2D2C{T}
    data::LNCData{2,0,T,0}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
4 coordinates 2D local natural coordinates, using 2 basic points.
$(TYPEDEF)
"""
struct LNC2D2P{T} <: LNC2D4C{T}
    data::LNCData{2,1,T,2}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
4 coordinates 2D local natural coordinates, using 1 basic point and 1 base vector.
$(TYPEDEF)
"""
struct LNC2D1P1V{T} <: LNC2D4C{T}
    data::LNCData{2,1,T,2}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coordinates 3D local natural coordinates, using 2 basic points.
$(TYPEDEF)
"""
struct LNC3D2P{T} <: LNC3D6C{T}
    data::LNCData{3,1,T,3}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coordinates 3D local natural coordinates, using 1 basic point and 1 base vector.
$(TYPEDEF)
"""
struct LNC3D1P1V{T} <: LNC3D6C{T}
    data::LNCData{3,1,T,3}
    conversion::SparseMatrixCSC{Int64,Int64}
end


"""
6 coordinates 2D local natural coordinates, using 1 basic point and 2 base vectors.
$(TYPEDEF)
"""
struct LNC1P2V{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coordinates 2D local natural coordinates, using 2 basic points and 1 base vector.
$(TYPEDEF)
"""
struct LNC2P1V{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
6 coordinates 2D local natural coordinates, using 3 basic points.
$(TYPEDEF)
"""
struct LNC3P{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coordinates 3D local natural coordinates, using 1 basic point and 3 base vectors.
$(TYPEDEF)
"""
struct LNC1P3V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coordinates 3D local natural coordinates, using 2 basic points and 2 base vectors.
$(TYPEDEF)
"""
struct LNC2P2V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coordinates 3D local natural coordinates, using 3 basic points and 1 base vector.
$(TYPEDEF)
"""
struct LNC3P1V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
    conversion::SparseMatrixCSC{Int64,Int64}
end

"""
12 coordinates 3D local natural coordinates, using 4 basic points.
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
    nmcs,q
end

function NC2D1P1V(ri,u,origin_position,θ::Number,origin_velocity,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC2D1P1V(ri,u,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    u̇ = origin_velocity + ω×u
    q̇ = vcat(ṙi,u̇)
    nmcs,SVector{4}(q),SVector{4}(q̇)
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
    nmcs,q
end

function NC2D2P(ri,rj,origin_position,θ::Number,origin_velocity,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC2D2P(ri,rj,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    q̇ = vcat(ṙi,ṙj)
    nmcs,SVector{4}(q),SVector{4}(q̇)
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
    nmcs,q
end

function NC3D1P1V(ri,u,origin_position,R::AbstractMatrix,origin_velocity,ω)
    nmcs,q = NC3D1P1V(ri,u,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    u̇ = origin_velocity + ω×u
    q̇ = vcat(ṙi,u̇)
    nmcs,SVector{6}(q),SVector{6}(q̇)
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
    nmcs,q
end

function NC3D2P(ri,rj,origin_position,R::AbstractMatrix,origin_velocity,ω)
    nmcs,q = NC3D2P(ri,rj,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    q̇ = vcat(ṙi,ṙj)
    nmcs,SVector{6}(q),SVector{6}(q̇)
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
    nmcs,q
end

function NC1P2V(ri,origin_position,θ,origin_velocity,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC1P2V(ri,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    u̇ = ω×q[3:4]
    v̇ = ω×q[5:6]
    q̇ = vcat(ṙi,u̇,v̇)
    nmcs,SVector{6}(q),SVector{6}(q̇)
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
    nmcs,q
end

function NC2P1V(ri,rj,origin_position,θ,origin_velocity,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC2P1V(ri,rj,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    v̇ = ω×q[5:6]
    q̇ = vcat(ṙi,ṙj,v̇)
    nmcs,SVector{6}(q),SVector{6}(q̇)
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
    nmcs,q
end

function NC3P(ri,rj,rk,origin_position,θ,origin_velocity,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC3P(ri,rj,rk,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    ṙk = origin_velocity + ω×(rk-origin_position)
    q̇ = vcat(ṙi,ṙj,ṙk)
    nmcs,SVector{6}(q),SVector{6}(q̇)
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
    nmcs,q
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
    nmcs,q
end

function NC1P3V(ri,origin_position,R,origin_velocity,ω)
    nmcs,q = NC1P3V(ri,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    u̇ = ω×q[4:6]
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

function NC1P3V(ri,u,v,w,origin_position,R,origin_velocity,ω)
    nmcs,q = NC1P3V(ri,u,v,w,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    u̇ = ω×u
    v̇ = ω×v
    ẇ = ω×w
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
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
    nmcs,q
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
    nmcs,q
end

function NC2P2V(ri,rj,origin_position,R,origin_velocity,ω)
    nmcs,q = NC2P2V(ri,rj,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

function NC2P2V(ri,rj,v,w,origin_position,R,origin_velocity,ω)
    nmcs,q = NC2P2V(ri,rj,v,w,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    v̇ = ω×v
    ẇ = ω×w
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
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
    nmcs,q
end

function NC3P1V(ri,rj,rk,origin_position,R,origin_velocity,ω)
    nmcs,q = NC3P1V(ri,rj,rk,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    ṙk = origin_velocity + ω×(rk-origin_position)
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,ṙk,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
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
    nmcs,q
end

function NC4P(ri,rj,rk,rl,origin_position,R,origin_velocity,ω)
    nmcs,q = NC4P(ri,rj,rk,rl,origin_position,R)
    ṙi = origin_velocity + ω×(ri-origin_position)
    ṙj = origin_velocity + ω×(rj-origin_position)
    ṙk = origin_velocity + ω×(rk-origin_position)
    ṙl = origin_velocity + ω×(rl-origin_position)
    q̇ = vcat(ṙi,ṙj,ṙk,ṙl)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

# """
# Return 2D or 3D local natural coordinates deformations for rigid bars。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds get_deform(ū::AbstractVector) = sqrt(ū⋅ū)

@inline @inbounds function get_deform(nmcs::Union{LNC2D1P,LNC3D1P})
    (;r̄i) = nmcs
    get_deform(r̄i)
end

@inline @inbounds function get_deform(nmcs::Union{LNC2D1P1V,LNC3D1P1V,LNC2D2P,LNC3D2P})
    (;ū) = nmcs
    get_deform(ū)
end

# """
# Return 2D natural coodinates deformations , for rigid bodies。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(ū,v̄)
    sqrt(ū⋅ū),sqrt(v̄⋅v̄),ū⋅v̄
end

@inline @inbounds function get_deform(nmcs::LNC2D6C)
    (;ū,v̄) = nmcs
    get_deform(ū,v̄)
end

# """
# Return 3D natural coodinates deformations , for rigid bodies。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(ū,v̄,w̄)
    u_sqrt = sqrt(ū⋅ū)
    v_sqrt = sqrt(v̄⋅v̄)
    w_sqrt = sqrt(w̄⋅w̄)
    vw = v̄⋅w̄
    uw = ū⋅w̄
    uv = ū⋅v̄
    u_sqrt,v_sqrt,w_sqrt,vw,uw,uv
end

@inline @inbounds function get_deform(nmcs::LNC3D12C)
    (;ū,v̄,w̄) = nmcs
    get_deform(ū,v̄,w̄)
end

# Intrinsic Constraints
## Intrinsic Constraints: Dispatch
"""
Return 2D or 3D intrinsic constraint(s) ，用于Dispatch。
$(TYPEDSIGNATURES)
"""
function make_constraints_function(nmcs::LNC,constraints_indices)
    deforms = get_deform(nmcs::LNC)
    make_constraints_function(nmcs,constraints_indices,deforms)
end

function make_inner_Φ(func,deforms)
    @inline @inbounds function ret_func(q)
        func(q,deforms)
    end
    @inline @inbounds function ret_func(q,d)
        func(q,d)
    end
    ret_func
end

"""
Return 2D or 3D intrinsic constraint(s) for point mass。
$(TYPEDSIGNATURES)
"""
function make_constraints_function(nmcs::Union{LNC2D2C,LNC3D3C},constraints_indices,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        nothing
    end
    make_inner_Φ(_inner_Φ,deforms)
end

"""
Return 2D or 3D intrinsic constraint(s) for rigid bars。
$(TYPEDSIGNATURES)
"""
function make_constraints_function(nmcs::Union{LNC2D4C,LNC3D6C},constraints_indices,deforms)
    ndim = get_num_of_dims(nmcs)
    cv = nmcs.conversion
    @inline @inbounds function _inner_Φ(q,d)
        qstd = cv*q
        u = @view qstd[ndim+1:2ndim]
        all = [u⋅u - d^2]
        all[constraints_indices]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

"""
Return 2D intrinsic constraint(s) , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_constraints_function(nmcs::LNC2D6C,constraints_indices,deforms)
    cv = nmcs.conversion
    @inline @inbounds function _inner_Φ(q,d)
        qstd = cv*q
        u = @view qstd[3:4]
        v = @view qstd[5:6]
        all = [
            (u⋅u - d[1]^2), 
            (v⋅v - d[2]^2), 
            (u⋅v - d[3])
        ]
        all[constraints_indices]
    end
    make_inner_Φ(_inner_Φ,deforms)
end


"""
Return 3D intrinsic constraint(s) , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_constraints_function(nmcs::LNC3D12C,constraints_indices,deforms)
    cv = nmcs.conversion
    @inline @inbounds function _inner_Φ(q,d)
        qstd = cv*q
        u = @view qstd[4:6]
        v = @view qstd[7:9]
        w = @view qstd[10:12]
        all = [
            u⋅u - d[1]^2, 
            v⋅v - d[2]^2, 
            w⋅w - d[3]^2, 
            v⋅w - d[4], 
            u⋅w - d[5], 
            u⋅v - d[6]
        ]
        @view all[constraints_indices]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

## Jacobians


"""
Return 2D or 3D Jacobian matrix for rigid bars。
$(TYPEDSIGNATURES)
"""
function make_constraints_jacobian(nmcs::Union{LNC2D2C,LNC3D3C},unconstrained_indices,constraints_indices)
    @inline @inbounds function inner_Φq(q)
        nothing
    end
end

"""
Return 2D or 3D Jacobian matrix for rigid bars。
$(TYPEDSIGNATURES)
"""
function make_constraints_jacobian(nmcs::Union{LNC2D4C,LNC3D6C},unconstrained_indices,constraints_indices)
    ndim = get_num_of_dims(nmcs)
    cv = nmcs.conversion
    @inline @inbounds function inner_Φq(q)
        qstd = cv*q
        u = @view qstd[ndim+1:2ndim]
        ret = zeros(eltype(q),1,2ndim)
        ret[ndim+1:2ndim] = 2u
        @view (ret*cv)[constraints_indices,unconstrained_indices]
    end
end

"""
Return 2D Jacobian matrix , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_constraints_jacobian(nmcs::LNC2D6C,unconstrained_indices,constraints_indices)
    cv = nmcs.conversion
    @inline @inbounds function inner_Φq(q)
        qstd = cv*q
        u = @view qstd[3:4]
        v = @view qstd[5:6]
        ret = zeros(eltype(q),3,6)
        ret[1,3:4] = 2u
        ret[2,5:6] = 2v
        ret[3,3:4] =  v
        ret[3,5:6] =  u
        @view (ret*cv)[constraints_indices,unconstrained_indices]
    end
end
"""
Return 3D Jacobian matrix , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_constraints_jacobian(nmcs::LNC3D12C,unconstrained_indices,constraints_indices)
    cv = nmcs.conversion
    @inline @inbounds function inner_Φq(q)
        qstd = cv*q
        u  = @view qstd[4:6]
        v  = @view qstd[7:9]
        w  = @view qstd[10:12]
        ret = zeros(eltype(q), 6, 12)
        ret[1,4:6]   = 2u
        ret[2,7:9]   = 2v
        ret[3,10:12] = 2w

        ret[4 ,7:9]  = w
        ret[4,10:12] = v

        ret[5, 4:6]  = w
        ret[5,10:12] = u

        ret[6,4:6] =  v
        ret[6,7:9] =  u

        @view (ret*cv)[constraints_indices,unconstrained_indices]
    end
end

"""
Return nullspace matrix
$(TYPEDSIGNATURES)
"""
function make_N(nmcs::LNC2D4C)
    cv = nmcs.conversion
    O2 = zero(I2_Bool)
    @inline @inbounds function inner_N(q)
        u,v = get_uv(nmcs,q)
        o2 = zero(u)
        ret = [
            I2_Bool  o2;
            O2       v;
        ]
        inv(cv)*ret
    end
end

function make_N(nmcs::LNC2D6C)
    cv = nmcs.conversion
    O2 = zero(I2_Bool)
    @inline @inbounds function inner_N(q)
        u,v = get_uv(nmcs,q)
        o2 = zero(u)
        ret = [
            I2_Bool    o2;
            O2         v;
            O2        -u;
        ]
        inv(cv)*ret
    end
end

function make_N(nmcs::LNC3D6C)
    cv = nmcs.conversion
    O3 = zero(I3_Bool)
    @inline @inbounds function inner_N(q)
        u,v,w = get_uvw(nmcs,q)
        o3 = zero(u)
        ret = [
            I3_Bool  o3 o3;
            O3      -w v;
        ]
        inv(cv)*ret
    end
end

function make_N(nmcs::LNC3D12C)
    cv = nmcs.conversion
    O3 = zero(I3_Bool)
    @inline @inbounds function inner_N(q)
        u,v,w = get_uvw(nmcs,q)
        # ret = [
        #     I3_Bool   O3;
        #     O3 -skew(u);
        #     O3 -skew(v);
        #     O3 -skew(w);
        # ]
        o3 = zero(u)
        ret = [
            I3_Bool    O3;
            O3  o3  -w  v;
            O3   w  o3 -u;
            O3  -v   u o3;
        ]
        inv(cv)*ret
    end
end
# Transformations
"""
Return 转换矩阵C。
$(TYPEDSIGNATURES)
"""
function to_local_coordinates(nmcs)
    (;r̄i,invX̄) = nmcs
    function c(r̄)
        invX̄*(r̄-r̄i)
    end
end

function to_transformation(nmcs::LNC)
    function C(c)
        ndim = get_num_of_dims(nmcs)
        ncoords = get_num_of_coordinates(nmcs)
        conversion = nmcs.conversion
        C_raw = hcat(1,transpose(c))
        SMatrix{ndim,ncoords}(kron(C_raw,IMatrix(ndim))*conversion)
    end
end

# CoordinateFunctions
"""
封装有函数的natural coodinates.
$(TYPEDEF)
"""
struct CoordinateFunctions{nmcsType,
                    cT<:Function,CT<:Function,
                    ΦT<:Function,ΦqT<:Function,
                    ∂Aᵀλ∂qT,∂Aq̇∂qT}
    nmcs::nmcsType
    c::cT
    C::CT
    Φ::ΦT
    Φq::ΦqT
    ∂Aᵀλ∂q::∂Aᵀλ∂qT
    ∂Aq̇∂q::∂Aq̇∂qT
end

get_idx(nmcs::Union{LNC2D4C,LNC3D6C}) = [
    [CartesianIndex(2,2),CartesianIndex(2,2)],
]

get_idx(nmcs::LNC2D6C) = [
    [CartesianIndex(2,2),CartesianIndex(2,2)],
    [CartesianIndex(3,3),CartesianIndex(3,3)],
    [CartesianIndex(2,3),CartesianIndex(3,2)]
]

get_idx(nmcs::LNC3D12C) = [
    [CartesianIndex(2,2),CartesianIndex(2,2)],
    [CartesianIndex(3,3),CartesianIndex(3,3)],
    [CartesianIndex(4,4),CartesianIndex(4,4)],
    [CartesianIndex(3,4),CartesianIndex(4,3)],
    [CartesianIndex(2,4),CartesianIndex(4,2)],
    [CartesianIndex(2,3),CartesianIndex(3,2)]
]

#todo cache Φqᵀq
function make_Φqᵀq(nmcs::LNC)
    cv = nmcs.conversion
    nld = get_num_of_local_dims(nmcs)
    ndim = get_num_of_dims(nmcs)
    I_Bool = IMatrix(ndim)
    idx = get_idx(nmcs)
    Φqᵀq = [
        begin
            ret_raw = zeros(Int,nld+1,nld+1)
            for ij in id
                ret_raw[ij] += 1
            end
            SymmetricPacked(transpose(cv)*kron(ret_raw,I_Bool)*cv)
        end
        for id in idx
    ]
end

#todo use SymmetricPacked to the end
function make_∂Aᵀλ∂q(nmcs::Union{LNC2D2C,LNC3D3C},unconstrained_indices,constraints_indices)
    function ∂Aᵀλ∂q(λ)
        nothing
    end
end

function make_∂Aᵀλ∂q(nmcs::LNC,unconstrained_indices,constraints_indices)
    Φqᵀq = make_Φqᵀq(nmcs)
    function ∂Aᵀλ∂q(λ)
        ret = [
            begin
                a = Φqᵀq[j][unconstrained_indices,unconstrained_indices] .* λ[i]
                # display(a)
                a 
            end
            for (i,j) in enumerate(constraints_indices)
        ]
        sum(ret)
    end
end

function make_∂Aq̇∂q(nmcs::Union{LNC2D2C,LNC3D3C},unconstrained_indices,constraints_indices)
    function ∂Aq̇∂q(q̇)
        nothing
    end
end

function make_∂Aq̇∂q(nmcs::LNC,unconstrained_indices,constraints_indices)
    Φqᵀq = make_Φqᵀq(nmcs)
    function ∂Aq̇∂q(q̇)
        q̇uc = @view q̇[unconstrained_indices]
        ret = [
            begin
                a = transpose(q̇uc)*Φqᵀq[j][unconstrained_indices,unconstrained_indices]
                # display(a)
                a 
            end
            for j in constraints_indices
        ]
        sum(ret)
    end
end
"""
Return ∂Aᵀλ∂q的前向自动微分结果。
$(TYPEDSIGNATURES)
"""
function make_∂Aᵀλ∂q_forwarddiff(Φq,nq,nuc)
    function ∂Aᵀλ∂q(λ)
        function ATλ(q)
            transpose(Φq(q))*λ
        end
        λT = eltype(λ)
        out = zeros(λT,nuc,nq)
        ForwardDiff.jacobian!(out,ATλ,ones(λT,nq))
    end
end


"""
Return ∂Aq̇∂q的前向自动微分结果。
$(TYPEDSIGNATURES)
"""
function make_∂Aq̇∂q_forwarddiff(Φq,nq,nλ)
    function ∂Aq̇∂q(q̇)
        function Aq̇(q)
            Φq(q)*q̇
        end
        q̇T = eltype(q̇)
        out = zeros(q̇T,nλ,nq)
        ForwardDiff.jacobian!(out,Aq̇,ones(q̇T,nq))
    end
end

"""
Return 未约束的natural coodinates 编号。
$(TYPEDSIGNATURES)
"""
function get_free_idx(nmcs::LNC,pres_idx)
    deleteat!(collect(1:get_num_of_coordinates(nmcs)),pres_idx)
end

"""
封装有函数的natural coodinates 构造子。
$(TYPEDSIGNATURES)
"""
function CoordinateFunctions(nmcs,unconstrained_indices,constraints_indices)
    c = to_local_coordinates(nmcs)
    C = to_transformation(nmcs)
    Φ = make_constraints_function(nmcs,constraints_indices)
    Φq = make_constraints_jacobian(nmcs,unconstrained_indices,constraints_indices)
    nq = get_num_of_coordinates(nmcs)
    nuc = length(unconstrained_indices)
    nΦ = length(constraints_indices)
    # ∂Aᵀλ∂q = make_∂Aᵀλ∂q_forwarddiff(Φq,nq,nuc)
    ∂Aᵀλ∂q = make_∂Aᵀλ∂q(nmcs,unconstrained_indices,constraints_indices)
    # ∂Aq̇∂q = make_∂Aq̇∂q_forwarddiff(Φq,nuc,nΦ)
    ∂Aq̇∂q = make_∂Aq̇∂q(nmcs,unconstrained_indices,constraints_indices)
    CoordinateFunctions(nmcs,c,C,Φ,Φq,∂Aᵀλ∂q,∂Aq̇∂q)
end


# Mass matrices
Ī2J̄(::Union{LNC2D2C,LNC3D3C,LNC2D6C,LNC3D6C},Ī)  = Ī
Ī2J̄(::LNC3D12C,Ī) = 1/2*tr(Ī)*I-Ī

function Īg2az(nmcs,m,Īg,mass_locus)
    (;r̄i,invX̄) = nmcs
    a = invX̄*(mass_locus-r̄i)
    J̄g = Ī2J̄(nmcs,Īg)
    J̄o = J̄g + m*mass_locus*transpose(mass_locus)
    J̄i = J̄o - m*r̄i*transpose(mass_locus) - m*mass_locus*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*J̄i*transpose(invX̄)
    #@assert issymmetric(z)
    a,Symmetric(z)
end

## Mass matrices: standard
function make_M(cf::CoordinateFunctions,m::T,Īg,mass_locus) where {T} # ami (area moment of inertia tensor)
    (;nmcs) = cf
    ndim = get_num_of_dims(nmcs)
    nld = get_num_of_local_dims(nmcs)
    ncoords = get_num_of_coordinates(nmcs)
    a,z = Īg2az(nmcs,m,Īg,mass_locus)
    M_raw = zeros(T,1+nld,1+nld)
    M_raw[1,1] = m
    M_raw[2:1+nld,1] = m*a
    M_raw[1,2:1+nld] = M_raw[2:1+nld,1]
    M_raw[2:1+nld,2:1+nld] .= z
    M_std = kron(M_raw,IMatrix(ndim))
    Y_nonstd = nmcs.conversion
    M = SMatrix{ncoords,ncoords}(transpose(Y_nonstd)*M_std*Y_nonstd)
end

make_X(q::AbstractVector,nmcs::LNC) = make_X(nmcs,q)

function make_X(nmcs::LNC,q::AbstractVector)
    Y = nmcs.conversion
    qstd = Y*q
    ndim = get_num_of_dims(nmcs)
    X = reshape(qstd[ndim+1:end],ndim,:)
    if (nmcs isa LNC2D2C) || (nmcs isa LNC3D3C)
        return X
    elseif  (nmcs isa LNC2D4C) || (nmcs isa LNC3D6C)
        # return find_rotation(nmcs,q)
        return SMatrix{ndim,1}(X)
    else
        return SMatrix{ndim,ndim}(X)
    end
end

find_rotation(q::AbstractVector, nmcs::LNC) = find_rotation(nmcs,q)
function find_rotation(nmcs::LNC,q::AbstractVector)
    (;invX̄) = nmcs
    ndim = get_num_of_dims(nmcs)
    if nmcs isa LNC2D4C
        (;r̄i,X̄) = nmcs
        ū,v̄ = get_uv(nmcs,vcat(r̄i,vec(X̄)))
        u,v = get_uv(nmcs,q)
        R = SMatrix{ndim,ndim}([u;;v]*inv([ū;;v̄]))
    elseif nmcs isa LNC3D6C
        (;r̄i,X̄) = nmcs
        ū,v̄,w̄ = get_uvw(nmcs,vcat(r̄i,vec(X̄)))
        u,v,w = get_uvw(nmcs,q)
        R = SMatrix{ndim,ndim}([u;;v;;w]*inv([ū;;v̄;;w̄]))
    else
        X = make_X(nmcs,q)
        R = SMatrix{ndim,ndim}(X*invX̄)
    end
    return R
end

find_angular_velocity(q::AbstractVector,q̇::AbstractVector,nmcs::LNC3D) = find_angular_velocity(nmcs,q,q̇)

function find_angular_velocity(nmcs::LNC,q::AbstractVector,q̇::AbstractVector)
    Ẋ = make_X(nmcs,q̇)
    X = make_X(nmcs,q)
    ndim = get_num_of_dims(nmcs)
    if ndim == 2
        u = X[:,1]
        u̇ = Ẋ[:,1]
        ω = SVector{1}([-u[2],u[1]]\u̇)
    else
        Ω = Ẋ*pinv(X)
        ω = SVector{3}(Ω[3,2],Ω[1,3],Ω[2,1])
    end
end

function get_uvw(nmcs::LNC3D,q)
    if     nmcs isa LNC1P3V
        u = @view q[4:6]
        v = @view q[7:9]
        w = @view q[10:12]
    elseif nmcs isa LNC2P2V
        ri = @view q[1:3]
        rj = @view q[4:6]
        v = @view q[7:9]
        w = @view q[10:12]
        u = rj - ri
    elseif nmcs isa LNC3P1V
        ri = @view q[1:3]
        rj = @view q[4:6]
        rk = @view q[7:9]
        w = @view q[10:12]
        u = rj - ri
        v = rk - ri
    elseif nmcs isa LNC4P
        ri = @view q[1:3]
        rj = @view q[4:6]
        rk = @view q[7:9]
        rl = @view q[10:12]
        u = rj - ri
        v = rk - ri
        w = rl - ri
    elseif nmcs isa LNC3D1P1V
        u = @view q[4:6]
        v,w = HouseholderOrthogonalization(u)
    elseif nmcs isa LNC3D2P
        ri = @view q[1:3]
        rj = @view q[4:6]
        u = rj -ri
        u ./= norm(u)
        v,w = HouseholderOrthogonalization(u)
    end
    SVector{3}(u),SVector{3}(v),SVector{3}(w)
end

function get_uv(nmcs::LNC2D,q)
    if     nmcs isa LNC1P2V
        u = @view q[3:4]
        v = @view q[5:6]
    elseif nmcs isa LNC2P1V
        ri = @view q[1:2]
        rj = @view q[3:4]
        v = @view q[5:6]
        u = rj - ri
    elseif nmcs isa LNC3P
        ri = @view q[1:2]
        rj = @view q[3:4]
        rk = @view q[5:6]
        u = rj - ri
        v = rk - ri
    elseif nmcs isa LNC2D1P1V
        u = @view q[3:4]
        v = rotation_matrix(π/2)*u
    elseif nmcs isa LNC2D2P
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj - ri
        v = rotation_matrix(π/2)*u
    end
    SVector{2}(u),SVector{2}(v)
end

"""
Return rigid body natural coodinates 
$(TYPEDSIGNATURES)
"""
function rigidstate2naturalcoords(nmcs::Union{LNC2D2C,LNC3D3C},origin_position,R)
    origin_position
end

function rigidstate2naturalcoords(nmcs::Union{LNC2D2C,LNC3D3C},origin_position,R,origin_velocity,ω)
    origin_position,origin_velocity
end

function rigidstate2naturalcoords(nmcs,origin_position,R)
    (;r̄i,X̄) = nmcs
    ri = origin_position + R*r̄i
    X = R*X̄
    qstd = vcat(ri,vec(X))
    Y = nmcs.conversion
    q = Y\qstd
    q
end

function rigidstate2naturalcoords(nmcs,origin_position,R,origin_velocity,ω)
    (;r̄i,X̄) = nmcs
    ri = origin_position + R*r̄i
    ṙi = origin_velocity + ω×(ri-origin_position)
    X = R*X̄
    Ẋ = reduce(hcat,Ref(ω) .× eachcol(X))
    qstd = vcat(ri,vec(X))
    q̇std = vcat(ṙi,vec(Ẋ))
    Y = nmcs.conversion
    q = Y\qstd
    q̇ = Y\q̇std
    q,q̇
end


end
