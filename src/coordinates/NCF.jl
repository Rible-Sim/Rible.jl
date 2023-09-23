module NCF
using LinearAlgebra
using StaticArrays
using Parameters
using ForwardDiff
using DocStringExtensions
using SymmetricFormats

export get_nconstraints, get_ncoords, get_ndof, get_nlocaldim
"""
所有局部自然坐标抽象超类。
"""
abstract type LNC{T} end
"""
所有二维局部自然坐标抽象超类。
"""
abstract type LNC2D{T} <: LNC{T} end
"""
所有三维局部自然坐标抽象超类。
"""
abstract type LNC3D{T} <: LNC{T} end
"""
所有坐标数目为4的二维局部自然坐标抽象超类，用于刚性直杆。
"""
abstract type LNC2D4C{T} <: LNC2D{T} end
"""
所有坐标数目为6的二维局部自然坐标抽象超类，用于任意形状刚体。
"""
abstract type LNC2D6C{T} <: LNC2D{T} end
"""
所有坐标数目为6的三维局部自然坐标抽象超类，用于刚性直杆。
"""
abstract type LNC3D6C{T} <: LNC3D{T} end
"""
所有坐标数目为12的三维局部自然坐标抽象超类，用于任意形状刚体。
"""
abstract type LNC3D12C{T} <: LNC3D{T} end

"""
二维或三维局部自然坐标类，用于任意一列点。
$(TYPEDEF)
"""
struct LNC1P{N,T} <: LNC{T}
    r̄i::SArray{Tuple{N},T,1,N}
end

struct LNCData{N,M,T,L}
    r̄i::SArray{Tuple{N},T,1,N}
    X̄::SArray{Tuple{N,M},T,2,L}
    invX̄::SArray{Tuple{M,N},T,2,L}
end

function LNCData(r̄i::SVector,X̄::SMatrix)
    LNCData(r̄i,X̄,(pinv(X̄)))
end

make_I(T,N) = SMatrix{N,N}(one(T)*I)
const I2_Int = make_I(Int,2)
const I3_Int = make_I(Int,3)

"""
返回空间维数。
$(TYPEDSIGNATURES)
"""
get_ndim(::LNC1P{N}) where N = N
get_ndim(nmcs::LNC) = get_ndim(nmcs.data)
get_ndim(::LNCData{N,M,T,L}) where {N,M,T,L} = N

"""
返回自然坐标所构成的坐标系的维数。
$(TYPEDSIGNATURES)
"""
get_nlocaldim(::LNC1P) = 0
get_nlocaldim(nmcs::LNC) = get_nlocaldim(nmcs.data)
get_nlocaldim(::LNCData{N,M,T,L}) where {N,M,T,L} = M
"""
返回坐标个数。
$(TYPEDSIGNATURES)
"""
get_ncoords(::LNC1P{N}) where N = N
get_ncoords(nmcs::LNC) = get_ncoords(nmcs.data)
get_ncoords(::LNCData{N,M,T,L}) where {N,M,T,L} = N+L

"""
返回约束方程个数。
$(TYPEDSIGNATURES)
"""
get_nconstraints(::LNC2D4C) = 1
get_nconstraints(::LNC2D6C) = 3
get_nconstraints(::LNC3D6C) = 1
get_nconstraints(::LNC3D12C) = 6
get_nconstraints(::LNC1P) = 0
"""
返回自由度数。
$(TYPEDSIGNATURES)
"""
get_ndof(nmcs::LNC) =  get_ncoords(nmcs) - get_nconstraints(nmcs)

# """
# 返回外积结果。
# $(TYPEDSIGNATURES)
# """
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

"""
返回旋转矩阵。
$(TYPEDSIGNATURES)
"""
rotation_matrix(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

# """
# ！！！！！！！！！
# $(TYPEDSIGNATURES)
# """
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

"""
返回斜对称矩阵。
$(TYPEDSIGNATURES)
"""
function skew(w)
    w1,w2,w3 = w
    o = zero(w1)
    @SMatrix [o -w3 w2;
              w3 o -w1;
             -w2 w1 o]
end

function skew(w::SVector{3})
    w1,w2,w3 = w
    o = zero(w1)
    @SMatrix [o -w3 w2;
              w3 o -w1;
             -w2 w1 o]
end

function skew(a::SVector{2})
	[-a[2],a[1]]
end

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
坐标数目为4的二维局部自然坐标类，使用2个基本点。
$(TYPEDEF)
"""
struct LNC2D2P{T} <: LNC2D4C{T}
    data::LNCData{2,1,T,2}
end

struct LNC2D1P1V{T} <: LNC2D4C{T}
    data::LNCData{2,1,T,2}
end

"""
坐标数目为6的三维局部自然坐标类，使用2个基本点。
$(TYPEDEF)
"""
struct LNC3D2P{T} <: LNC3D6C{T}
    data::LNCData{3,1,T,3}
end

struct LNC3D1P1V{T} <: LNC3D6C{T}
    data::LNCData{3,1,T,3}
end


"""
坐标数目为6的二维局部自然坐标类，使用1个基本点、2个基本向量。
$(TYPEDEF)
"""
struct LNC1P2V{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
end

"""
坐标数目为6的二维局部自然坐标类，使用2个基本点、1个基本向量。
$(TYPEDEF)
"""
struct LNC2P1V{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
end

"""
坐标数目为6的二维局部自然坐标类，使用3个基本点。
$(TYPEDEF)
"""
struct LNC3P{T} <: LNC2D6C{T}
    data::LNCData{2,2,T,4}
end

"""
坐标数目为12的三维局部自然坐标类，使用1个基本点、3个基本向量。
$(TYPEDEF)
"""
struct LNC1P3V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
end

"""
坐标数目为12的三维局部自然坐标类，使用2个基本点、2个基本向量。
$(TYPEDEF)
"""
struct LNC2P2V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
end

"""
坐标数目为12的三维局部自然坐标类，使用3个基本点、1个基本向量。
$(TYPEDEF)
"""
struct LNC3P1V{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
end

"""
坐标数目为12的三维局部自然坐标类，使用4个基本点。
$(TYPEDEF)
"""
struct LNC4P{T} <: LNC3D12C{T}
    data::LNCData{3,3,T,9}
end

"""
！！！！！！！！！
$(TYPEDSIGNATURES)
"""

function get_conversion(nmcs::Union{LNC2D1P1V,LNC3D1P1V})
    ndim = get_ndim(nmcs)
    kron(
        [ 1 0;
          0 1],
        make_I(Int,ndim)
    )
end

function get_conversion(nmcs::Union{LNC2D2P,LNC3D2P})
    ndim = get_ndim(nmcs)
    kron(
        [ 1 0;
         -1 1],
        make_I(Int,ndim)
    )
end

get_conversion(::LNC1P2V) = kron(
    [ 1 0 0;
      0 1 0;
      0 0 1],
    I2_Int
)

get_conversion(::LNC2P1V) = kron(
    [ 1 0 0;
     -1 1 0;
      0 0 1],
    I2_Int
)

get_conversion(::LNC3P) = kron(
    [ 1 0 0;
     -1 1 0;
     -1 0 1],
    I2_Int
)

get_conversion(::LNC1P3V) = kron(
    [ 1 0 0 0;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1],
    I3_Int
)

get_conversion(::LNC2P2V) = kron(
    [ 1 0 0 0;
     -1 1 0 0;
      0 0 1 0;
      0 0 0 1],
    I3_Int
)

get_conversion(::LNC3P1V) = kron(
    [ 1 0 0 0;
     -1 1 0 0;
     -1 0 1 0;
      0 0 0 1],
    I3_Int
)

get_conversion(::LNC4P) = kron(
    [ 1 0 0 0;
     -1 1 0 0;
     -1 0 1 0;
     -1 0 0 1],
    I3_Int
)

"""
返回刚体自然坐标
$(TYPEDSIGNATURES)
"""
function rigidstate2naturalcoords(nmcs,ro,R)
    (;r̄i,X̄) = nmcs
    ri = ro + R*r̄i
    X = R*X̄
    qstd = vcat(ri,vec(X))
    Y = get_conversion(nmcs)
    q = Y\qstd
    q
end

function rigidstate2naturalcoords(nmcs,ro,R,ṙo,ω)
    (;r̄i,X̄) = nmcs
    ri = ro + R*r̄i
    ṙi = ṙo + ω×(ri-ro)
    X = R*X̄
    Ẋ = reduce(hcat,Ref(ω) .× eachcol(X))
    qstd = vcat(ri,vec(X))
    q̇std = vcat(ṙi,vec(Ẋ))
    Y = get_conversion(nmcs)
    q = Y\qstd
    q̇ = Y\q̇std
    q,q̇
end


"""
返回二维刚性直杆自然坐标类。
$(TYPEDSIGNATURES)
"""
function NC2D1P1V(ri::AbstractVector{T},u::AbstractVector{T},
               ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
               ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    ū = invR*u
    X̄ = ū
    nmcs = LNC2D1P1V(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,1}(X̄)
        )
    )
    q = vcat(ri,u)
    nmcs,q
end

function NC2D1P1V(ri,u,ro,θ::Number,ṙo,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC2D1P1V(ri,u,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ṙo + ω×u
    q̇ = vcat(ṙi,u̇)
    nmcs,SVector{4}(q),SVector{4}(q̇)
end

function NC2D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
    ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    r̄j = invR*(rj-ro) #
    ū = r̄j-r̄i
    X̄ = ū
    nmcs = LNC2D2P(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,1}(X̄)
        )
    )
    q = vcat(ri,rj)
    nmcs,q
end

function NC2D2P(ri,rj,ro,θ::Number,ṙo,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC2D2P(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    q̇ = vcat(ṙi,ṙj)
    nmcs,SVector{4}(q),SVector{4}(q̇)
end

"""
返回三维刚性直杆自然坐标类。
$(TYPEDSIGNATURES)
"""
function NC3D1P1V(ri::AbstractVector{T},u::AbstractVector{T},
               ro=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
               ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) 
    ū = invR*u
    X̄ = ū
    nmcs = LNC3D1P1V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,1}(X̄)
        )
    )
    q = vcat(ri,u)
    nmcs,q
end

function NC3D1P1V(ri,u,ro,R::AbstractMatrix,ṙo,ω)
    nmcs,q = NC3D1P1V(ri,u,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ṙo + ω×u
    q̇ = vcat(ṙi,u̇)
    nmcs,SVector{6}(q),SVector{6}(q̇)
end

function NC3D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    ro=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
    ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    r̄j = invR*(rj-ro)
    ū = r̄j-r̄i
    X̄ = ū
    nmcs = LNC3D2P(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,1}(X̄)
        )
    )
    q = vcat(ri,rj)
    nmcs,q
end

function NC3D2P(ri,rj,ro,R::AbstractMatrix,ṙo,ω)
    nmcs,q = NC3D2P(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    q̇ = vcat(ṙi,ṙj)
    nmcs,SVector{6}(q),SVector{6}(q̇)
end

"""
返回二维任意形状刚体自然坐标类，使用1个基本点、2个基本向量。
$(TYPEDSIGNATURES)
"""
function NC1P2V(ri::AbstractVector{T},
               ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
               ) where T
    o = zero(T)
    i = one(T)
    ū = SVector(i,o)
    v̄ = SVector(o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    u = R*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    nmcs = LNC1P2V(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,2}(X̄)
        )
    )
    q = vcat(ri,u,v)
    nmcs,q
end

function NC1P2V(ri,ro,θ,ṙo,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC1P2V(ri,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×q[3:4]
    v̇ = ω×q[5:6]
    q̇ = vcat(ṙi,u̇,v̇)
    nmcs,SVector{6}(q),SVector{6}(q̇)
end

"""
返回二维任意形状刚体自然坐标类，使用2个基本点、1个基本向量。
$(TYPEDSIGNATURES)
"""
function NC2P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    ū = r̄j-r̄i
    v̄ = rotation_matrix(π/2)*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    nmcs = LNC2P1V(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,2}(X̄)
        )
    )
    q = vcat(ri,rj,v)
    nmcs,q
end

function NC2P1V(ri,rj,ro,θ,ṙo,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC2P1V(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×q[5:6]
    q̇ = vcat(ṙi,ṙj,v̇)
    nmcs,SVector{6}(q),SVector{6}(q̇)
end

"""
返回二维任意形状刚体自然坐标类，使用3个基本点。
$(TYPEDSIGNATURES)
"""
function NC3P(ri::AbstractVector{T},rj::AbstractVector{T},rk::AbstractVector{T},
              ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    X̄ = hcat(ū,v̄)
    nmcs = LNC3P(
        LNCData(
            SVector{2}(r̄i),
            SMatrix{2,2}(X̄)
        )
    )
    q = vcat(ri,rj,rk)
    nmcs,q
end
function NC3P(ri,rj,rk,ro,θ,ṙo,ω)
    R = rotation_matrix(θ)
    nmcs,q = NC3P(ri,rj,rk,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    q̇ = vcat(ṙi,ṙj,ṙk)
    nmcs,SVector{6}(q),SVector{6}(q̇)
end

## 3D Rigid body
"""
返回三维任意形状刚体自然坐标类，使用1个基本点、3个基本向量。
$(TYPEDSIGNATURES)
"""
function NC1P3V(ri::AbstractVector{T},
                ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                ) where T
    o = zero(T)
    i = one(T)
    ū = SVector(i,o,o)
    v̄ = SVector(o,i,o)
    w̄ = SVector(o,o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    u = R*ū
    v = R*v̄
    w = R*w̄
    X̄ = hcat(ū,v̄,w̄)
    nmcs = LNC1P3V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        )
    )
    q = vcat(ri,u,v,w)
    nmcs,q
end
function NC1P3V(ri::AbstractVector{T},u::AbstractVector{T},v::AbstractVector{T},w::AbstractVector{T},
                ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                ) where T
    o = zero(T)
    i = one(T)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    nmcs = LNC1P3V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        )
    )
    q = vcat(ri,u,v,w)
    nmcs,q
end

function NC1P3V(ri,ro,R,ṙo,ω)
    nmcs,q = NC1P3V(ri,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×q[4:6]
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

function NC1P3V(ri,u,v,w,ro,R,ṙo,ω)
    nmcs,q = NC1P3V(ri,u,v,w,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×u
    v̇ = ω×v
    ẇ = ω×w
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

"""
返回三维任意形状刚体自然坐标类，使用2个基本点、2个基本向量。
$(TYPEDSIGNATURES)
"""
function NC2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
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
        )
    )
    q = vcat(ri,rj,v,w)
    nmcs,q
end

function NC2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                v::AbstractVector{T},w::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    u = rj-ri
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    nmcs = LNC2P2V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        )
    )
    q = vcat(ri,rj,v,w)
    nmcs,q
end

function NC2P2V(ri,rj,ro,R,ṙo,ω)
    nmcs,q = NC2P2V(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

function NC2P2V(ri,rj,v,w,ro,R,ṙo,ω)
    nmcs,q = NC2P2V(ri,rj,v,w,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×v
    ẇ = ω×w
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

"""
返回三维任意形状刚体自然坐标类，使用3个基本点、1个基本向量。
$(TYPEDSIGNATURES)
"""
function NC3P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    u = rj-ri
    v = rk-ri
    w = u×v
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    nmcs = LNC3P1V(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        )
    )
    q = vcat(ri,rj,rk,w)
    nmcs,q
end

function NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
    nmcs,q = NC3P1V(ri,rj,rk,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,ṙk,ẇ)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

"""
返回三维任意形状刚体自然坐标类，使用4个基本点。
$(TYPEDSIGNATURES)
"""
function NC4P(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},rl::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    r̄l = invR*(rl-ro)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    w̄ = r̄l - r̄i
    X̄ = hcat(ū,v̄,w̄)
    nmcs = LNC4P(
        LNCData(
            SVector{3}(r̄i),
            SMatrix{3,3}(X̄)
        )
    )
    q = vcat(ri,rj,rk,rl)
    nmcs,q
end

function NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
    nmcs,q = NC4P(ri,rj,rk,rl,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    ṙl = ṙo + ω×(rl-ro)
    q̇ = vcat(ṙi,ṙj,ṙk,ṙl)
    nmcs,SVector{12}(q),SVector{12}(q̇)
end

# """
# 返回二维或三维局部自然坐标类的变形量，用于刚性直杆。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds get_deform(ū) = sqrt(ū⋅ū)
@inline @inbounds function get_deform(nmcs::Union{LNC2D1P1V,LNC3D1P1V})
    (;ū) = nmcs
    get_deform(ū)
end
@inline @inbounds function get_deform(nmcs::Union{LNC2D2P,LNC3D2P})
    (;r̄i,r̄j) = nmcs
    ū = r̄j-r̄i
    get_deform(ū)
end

# """
# 返回二维自然坐标类的变形量，用于任意形状刚体。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(ū,v̄)
    sqrt(ū⋅ū),sqrt(v̄⋅v̄),ū⋅v̄
end
@inline @inbounds function get_deform(nmcs::LNC1P2V)
    @unpack ū,v̄ = nmcs
    get_deform(ū,v̄)
end
@inline @inbounds function get_deform(nmcs::LNC2P1V)
    @unpack r̄i,r̄j,v̄ = nmcs
    ū = r̄j-r̄i
    get_deform(ū,v̄)
end
@inline @inbounds function get_deform(nmcs::LNC3P)
    @unpack r̄i,r̄j,r̄k = nmcs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    get_deform(ū,v̄)
end

# """
# 返回三维自然坐标类的变形量，用于任意形状刚体。
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
@inline @inbounds function get_deform(nmcs::LNC1P3V)
    @unpack ū,v̄,w̄ = nmcs
    get_deform(ū,v̄,w̄)
end
@inline @inbounds function get_deform(nmcs::LNC2P2V)
    @unpack r̄i,r̄j,v̄,w̄ = nmcs
    ū = r̄j-r̄i
    get_deform(ū,v̄,w̄)
end
@inline @inbounds function get_deform(nmcs::LNC3P1V)
    @unpack r̄i,r̄j,r̄k,w̄ = nmcs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    get_deform(ū,v̄,w̄)
end
@inline @inbounds function get_deform(nmcs::LNC4P)
    @unpack r̄i,r̄j,r̄k,r̄l = nmcs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = r̄l-r̄i
    get_deform(ū,v̄,w̄)
end

# Intrinsic Constraints
## Intrinsic Constraints: Dispatch
"""
返回二维或三维内部约束，用于Dispatch。
$(TYPEDSIGNATURES)
"""
function make_Φ(nmcs::LNC,Φi)
    deforms = get_deform(nmcs::LNC)
    make_Φ(nmcs,Φi,deforms)
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
返回二维或三维内部约束，用于刚性直杆。
$(TYPEDSIGNATURES)
"""
function make_Φ(nmcs::Union{LNC2D4C,LNC3D6C},Φi,deforms)
    ndim = get_ndim(nmcs)
    cv = get_conversion(nmcs)
    @inline @inbounds function _inner_Φ(q,d)
        qstd = cv*q
        u = @view qstd[ndim+1:2ndim]
        all = [u⋅u - d^2]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
"""
返回二维内部约束，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φ(nmcs::LNC2D6C,Φi,deforms)
    cv = get_conversion(nmcs)
    @inline @inbounds function _inner_Φ(q,d)
        qstd = cv*q
        u = @view qstd[3:4]
        v = @view qstd[5:6]
        all = [
            (u⋅u - d[1]^2), 
            (v⋅v - d[2]^2), 
            (u⋅v - d[3])
        ]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end


"""
返回三维内部约束，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φ(nmcs::LNC3D12C,Φi,deforms)
    cv = get_conversion(nmcs)
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
        @view all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

## Jacobians

"""
返回二维或三维雅可比矩阵，用于刚性直杆。
$(TYPEDSIGNATURES)
"""
function make_Φq(nmcs::Union{LNC2D4C,LNC3D6C},uci,Φi)
    ndim = get_ndim(nmcs)
    cv = get_conversion(nmcs)
    @inline @inbounds function inner_Φq(q)
        qstd = cv*q
        u = @view qstd[ndim+1:2ndim]
        ret = zeros(eltype(q),1,2ndim)
        ret[ndim+1:2ndim] = 2u
        @view (ret*cv)[Φi,uci]
    end
end

"""
返回二维雅可比矩阵，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φq(nmcs::LNC2D6C,uci,Φi)
    cv = get_conversion(nmcs)
    @inline @inbounds function inner_Φq(q)
        qstd = cv*q
        u = @view qstd[3:4]
        v = @view qstd[5:6]
        ret = zeros(eltype(q),3,6)
        ret[1,3:4] = 2u
        ret[2,5:6] = 2v
        ret[3,3:4] =  v
        ret[3,5:6] =  u
        @view (ret*cv)[Φi,uci]
    end
end
"""
返回三维雅可比矩阵，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φq(nmcs::LNC3D12C,uci,Φi)
    cv = get_conversion(nmcs)
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

        @view (ret*cv)[Φi,uci]
    end
end

"""
Return nullspace matrix
$(TYPEDSIGNATURES)
"""
function make_N(nmcs::LNC2D4C)
    cv = get_conversion(nmcs)
    O2 = zero(I2_Int)
    @inline @inbounds function inner_N(q)
        u,v = get_uv(nmcs,q)
        o2 = zero(u)
        ret = [
            I2_Int  o2;
            O2       v;
        ]
        inv(cv)*ret
    end
end

function make_N(nmcs::LNC2D6C)
    cv = get_conversion(nmcs)
    O2 = zero(I2_Int)
    @inline @inbounds function inner_N(q)
        u,v = get_uv(nmcs,q)
        o2 = zero(u)
        ret = [
            I2_Int    o2;
            O2         v;
            O2        -u;
        ]
        inv(cv)*ret
    end
end

function make_N(nmcs::LNC3D6C)
    cv = get_conversion(nmcs)
    O3 = zero(I3_Int)
    @inline @inbounds function inner_N(q)
        u,v,w = get_uvw(nmcs,q)
        o3 = zero(u)
        ret = [
            I3_Int  o3 o3;
            O3      -w v;
        ]
        inv(cv)*ret
    end
end

function make_N(nmcs::LNC3D12C)
    cv = get_conversion(nmcs)
    O3 = zero(I3_Int)
    @inline @inbounds function inner_N(q)
        u,v,w = get_uvw(nmcs,q)
        # ret = [
        #     I3_Int   O3;
        #     O3 -skew(u);
        #     O3 -skew(v);
        #     O3 -skew(w);
        # ]
        o3 = zero(u)
        ret = [
            I3_Int    O3;
            O3  o3  -w  v;
            O3   w  o3 -u;
            O3  -v   u o3;
        ]
        inv(cv)*ret
    end
end
# Transformations
"""
返回转换矩阵C。
$(TYPEDSIGNATURES)
"""
function make_c(nmcs)
    (;r̄i,invX̄) = nmcs
    function c(r̄)
        invX̄*(r̄-r̄i)
    end
end
function make_C(nmcs::LNC)
    function C(c)
        ndim = get_ndim(nmcs)
        ncoords = get_ncoords(nmcs)
        conversion = get_conversion(nmcs)
        C_raw = hcat(1,transpose(c))
        SMatrix{ndim,ncoords}(kron(C_raw,make_I(Int,ndim))*conversion)
    end
end

# CoordinateFunctions
"""
封装有函数的自然坐标类。
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
    cv = get_conversion(nmcs)
    nld = get_nlocaldim(nmcs)
    ndim = get_ndim(nmcs)
    I_int = make_I(Int,ndim)
    idx = get_idx(nmcs)
    Φqᵀq = [
        begin
            ret_raw = zeros(Int,nld+1,nld+1)
            for ij in id
                ret_raw[ij] += 1
            end
            SymmetricPacked(transpose(cv)*kron(ret_raw,I_int)*cv)
        end
        for id in idx
    ]
end

#todo use SymmetricPacked to the end
function make_∂Aᵀλ∂q(nmcs::LNC,uci,Φi)
    Φqᵀq = make_Φqᵀq(nmcs)
    function ∂Aᵀλ∂q(λ)
        ret = [
            begin
                a = Φqᵀq[j][uci,uci] .* λ[i]
                # display(a)
                a 
            end
            for (i,j) in enumerate(Φi)
        ]
        sum(ret)
    end
end

function make_∂Aq̇∂q(nmcs::LNC,uci,Φi)
    Φqᵀq = make_Φqᵀq(nmcs)
    function ∂Aq̇∂q(q̇)
        q̇uc = @view q̇[uci]
        ret = [
            begin
                a = transpose(q̇uc)*Φqᵀq[j][uci,uci]
                # display(a)
                a 
            end
            for j in Φi
        ]
        sum(ret)
    end
end
"""
返回∂Aᵀλ∂q的前向自动微分结果。
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
返回∂Aq̇∂q的前向自动微分结果。
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
返回未约束的自然坐标编号。
$(TYPEDSIGNATURES)
"""
function get_free_idx(nmcs::LNC,pres_idx)
    deleteat!(collect(1:get_ncoords(nmcs)),pres_idx)
end

"""
封装有函数的自然坐标类构造子。
$(TYPEDSIGNATURES)
"""
function CoordinateFunctions(nmcs,uci,Φi)
    c = make_c(nmcs)
    C = make_C(nmcs)
    Φ = make_Φ(nmcs,Φi)
    Φq = make_Φq(nmcs,uci,Φi)
    nq = get_ncoords(nmcs)
    nuc = length(uci)
    nΦ = length(Φi)
    # ∂Aᵀλ∂q = make_∂Aᵀλ∂q_forwarddiff(Φq,nq,nuc)
    ∂Aᵀλ∂q = make_∂Aᵀλ∂q(nmcs,uci,Φi)
    # ∂Aq̇∂q = make_∂Aq̇∂q_forwarddiff(Φq,nuc,nΦ)
    ∂Aq̇∂q = make_∂Aq̇∂q(nmcs,uci,Φi)
    CoordinateFunctions(nmcs,c,C,Φ,Φq,∂Aᵀλ∂q,∂Aq̇∂q)
end


# Mass matrices
Ī2J̄(::LNC2D4C,Ī)  = Ī
Ī2J̄(::LNC3D6C,Ī)  = Ī
Ī2J̄(::LNC2D6C,Ī)  = Ī
Ī2J̄(::LNC3D12C,Ī) = 1/2*tr(Ī)*I-Ī

function Īg2az(nmcs,m,Īg,r̄g)
    (;r̄i,invX̄) = nmcs
    a = invX̄*(r̄g-r̄i)
    J̄g = Ī2J̄(nmcs,Īg)
    J̄o = J̄g + m*r̄g*transpose(r̄g)
    J̄i = J̄o - m*r̄i*transpose(r̄g) - m*r̄g*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*J̄i*transpose(invX̄)
    #@assert issymmetric(z)
    a,Symmetric(z)
end

## Mass matrices: standard
function make_M(cf::CoordinateFunctions,m::T,Īg,r̄g) where {T} # ami (area moment of inertia tensor)
    (;nmcs) = cf
    ndim = get_ndim(nmcs)
    nld = get_nlocaldim(nmcs)
    ncoords = get_ncoords(nmcs)
    a,z = Īg2az(nmcs,m,Īg,r̄g)
    M_raw = zeros(T,1+nld,1+nld)
    M_raw[1,1] = m
    M_raw[2:1+nld,1] = m*a
    M_raw[1,2:1+nld] = M_raw[2:1+nld,1]
    M_raw[2:1+nld,2:1+nld] .= z
    M_std = kron(M_raw,make_I(T,ndim))
    Y_nonstd = get_conversion(nmcs)
    M = SMatrix{ncoords,ncoords}(transpose(Y_nonstd)*M_std*Y_nonstd)
end

make_X(q::AbstractVector,nmcs::LNC) = make_X(nmcs,q)

function make_X(nmcs::LNC,q::AbstractVector)
    Y = get_conversion(nmcs)
    qstd = Y*q
    ndim = get_ndim(nmcs)
    X = reshape(qstd[ndim+1:end],ndim,:)
    if (nmcs isa LNC2D4C) || (nmcs isa LNC3D6C)
        return find_R(nmcs,q)
    else
        return SMatrix{ndim,ndim}(X)
    end
end

find_R(q::AbstractVector, nmcs::LNC) = find_R(nmcs,q)
function find_R(nmcs::LNC,q::AbstractVector)
    (;invX̄) = nmcs
    ndim = get_ndim(nmcs)
    if nmcs isa LNC3D6C
        (;r̄i,X̄) = nmcs
        ū,v̄,w̄ = get_uvw(nmcs,vcat(r̄i,vec(X̄)))
        u,v,w = get_uvw(nmcs,q)
        R = SMatrix{ndim,ndim}([u;;v;;w]*inv([ū;;v̄;;w̄]))
    elseif nmcs isa LNC2D4C
        (;r̄i,X̄) = nmcs
        ū,v̄ = get_uv(nmcs,vcat(r̄i,vec(X̄)))
        u,v = get_uv(nmcs,q)
        R = SMatrix{ndim,ndim}([u;;v]*inv([ū;;v̄]))
    else
        X = make_X(nmcs,q)
        R = SMatrix{ndim,ndim}(X*invX̄)
    end
    return R
end

find_ω(q::AbstractVector,q̇::AbstractVector,nmcs::LNC3D) = find_ω(nmcs,q,q̇)

function find_ω(nmcs::LNC,q::AbstractVector,q̇::AbstractVector)
    Ẋ = make_X(nmcs,q̇)
    X = make_X(nmcs,q)
    ndim = get_ndim(nmcs)
    # @show Ẋ,X
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
    if nmcs isa LNC1P3V
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
    if nmcs isa LNC1P2V
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
    elseif nmcs isa LNC2D2P
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

end
