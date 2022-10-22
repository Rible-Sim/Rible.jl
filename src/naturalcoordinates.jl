module NaturalCoordinates
using LinearAlgebra
using StaticArrays
using Parameters
using ForwardDiff
using DocStringExtensions

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
struct LNCMP{M,L,T,N} <: LNC{T}
    r̄ps::SArray{Tuple{M,L},T,2,N}
end

make_I(T,N) = SMatrix{N,N}(one(T)*I)
I2_Int = make_I(Int,2)
I3_Int = make_I(Int,3)

"""
返回空间维数。
$(TYPEDSIGNATURES)
"""
get_ndim(::LNC2D) = 2
get_ndim(::LNC3D) = 3
get_ndim(::LNCMP{M}) where M = M
"""
返回自然坐标所构成的坐标系的维数。
$(TYPEDSIGNATURES)
"""
get_nlocaldim(::LNC2D4C) = 1
get_nlocaldim(::LNC2D6C) = 2
get_nlocaldim(::LNC3D6C) = 1
get_nlocaldim(::LNC3D12C) = 3
get_nlocaldim(::LNCMP{M,L}) where {M,L} = L
"""
返回坐标个数。
$(TYPEDSIGNATURES)
"""
get_ncoords(::LNC2D4C) = 4
get_ncoords(::LNC2D6C) = 6
get_ncoords(::LNC3D6C) = 6
get_ncoords(::LNC3D12C) = 12
get_ncoords(::LNCMP{M,L}) where {M,L}  = L*M
"""
返回约束方程个数。
$(TYPEDSIGNATURES)
"""
get_nconstraints(::LNC2D4C) = 1
get_nconstraints(::LNC2D6C) = 3
get_nconstraints(::LNC3D6C) = 1
get_nconstraints(::LNC3D12C) = 6
get_nconstraints(::LNCMP{M,L}) where {M,L}  = L*M
"""
返回自由度数。
$(TYPEDSIGNATURES)
"""
get_ndof(lncs::LNC) =  get_ncoords(lncs) - get_nconstraints(lncs)

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

"""
坐标数目为4的二维局部自然坐标类，使用2个基本点。
$(TYPEDEF)
"""
struct LNC2D2P{T} <: LNC2D4C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,1},T,2,2}
    invX̄::SArray{Tuple{1,2},T,2,2}
end

"""
坐标数目为6的三维局部自然坐标类，使用2个基本点。
$(TYPEDEF)
"""
struct LNC3D2P{T} <: LNC3D6C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,1},T,2,3}
    invX̄::SArray{Tuple{1,3},T,2,3}
end

"""
坐标数目为6的二维局部自然坐标类，使用1个基本点、2个基本向量。
$(TYPEDEF)
"""
struct LNC1P2V{T} <: LNC2D6C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    ū::SArray{Tuple{2},T,1,2}
    v̄::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

"""
坐标数目为6的二维局部自然坐标类，使用2个基本点、1个基本向量。
$(TYPEDEF)
"""
struct LNC2P1V{T} <: LNC2D6C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
    v̄::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

"""
坐标数目为6的二维局部自然坐标类，使用3个基本点。
$(TYPEDEF)
"""
struct LNC3P{T} <: LNC2D6C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
    r̄k::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

"""
坐标数目为12的三维局部自然坐标类，使用1个基本点、3个基本向量。
$(TYPEDEF)
"""
struct LNC1P3V{T} <: LNC3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    ū::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

"""
坐标数目为12的三维局部自然坐标类，使用2个基本点、2个基本向量。
$(TYPEDEF)
"""
struct LNC2P2V{T} <: LNC3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

"""
坐标数目为12的三维局部自然坐标类，使用3个基本点、1个基本向量。
$(TYPEDEF)
"""
struct LNC3P1V{T} <: LNC3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

"""
坐标数目为12的三维局部自然坐标类，使用4个基本点。
$(TYPEDEF)
"""
struct LNC4P{T} <: LNC3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    r̄l::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

"""
！！！！！！！！！
$(TYPEDSIGNATURES)
"""
get_conversion(lncs::LNCMP) = I
get_conversion(lncs::LNC1P2V) = I
get_conversion(lncs::LNC1P3V) = I
function get_conversion(lncs::Union{LNC2D2P,LNC3D2P})
    ndim = get_ndim(lncs)
    kron(
        [ 1 0;
         -1 1],
        make_I(Int,ndim)
    )
end
get_conversion(lncs::LNC2P1V) = kron(
    [ 1 0 0;
     -1 1 0;
      0 0 1],
    I2_Int
)
get_conversion(lncs::LNC3P) = kron(
    [ 1 0 0;
     -1 1 0;
     -1 0 1],
    I2_Int
)
get_conversion(lncs::LNC2P2V) = kron(
    [ 1 0 0 0;
     -1 1 0 0;
      0 0 1 0;
      0 0 0 1],
    I3_Int
)
get_conversion(lncs::LNC3P1V) = kron(
    [ 1 0 0 0;
     -1 1 0 0;
     -1 0 1 0;
      0 0 0 1],
    I3_Int
)
get_conversion(lncs::LNC4P) = kron(
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
function rigidstate2naturalcoords(lncs::LNCMP,ro,R,ṙo,ω)
    (;r̄ps) = lncs
    rps = Ref(ro) .+ Ref(R) .* eachcol(r̄ps)
    ṙps = Ref(ṙo) .+ Ref(ω) .× (rps .- Ref(ro))
    q = reduce(vcat,rps)
    q̇ = reduce(vcat,ṙps)
    q,q̇
end
function rigidstate2naturalcoords(lncs,ro,R,ṙo,ω)
    (;r̄i,X̄) = lncs
    ri = ro + R*r̄i
    ṙi = ṙo + ω×(ri-ro)
    X = R*X̄
    Ẋ = reduce(hcat,Ref(ω) .× eachcol(X))
    qstd = vcat(ri,vec(X))
    q̇std = vcat(ṙi,vec(Ẋ))
    Y = get_conversion(lncs)
    q = Y\qstd
    q̇ = Y\q̇std
    q,q̇
end


"""
返回任意列点自然坐标类、自然坐标。
$(TYPEDSIGNATURES)
"""
function NCMP(rps::AbstractVector{<:AbstractVector},ro,R,ṙo,ω)
    M = size(eltype(rps))[1]
    L = length(rps)
    invR = transpose(R)
    r̄ps = SMatrix{M,L}(reduce(hcat,Ref(invR).*(rps.-Ref(ro))))
    q = reduce(vcat,rps)
    LNCMP(r̄ps),q
end


"""
返回二维刚性直杆自然坐标类。
$(TYPEDSIGNATURES)
"""
function NC2D2P(ri::AbstractVector{T},rj::AbstractVector{T},
               ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)
               ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    r̄j = invR*(rj-ro) #
    ū = r̄j-r̄i
    X̄ = ū
    invX̄ = pinv(X̄)
    lncs = LNC2D2P(SVector{2}(r̄i),SVector{2}(r̄j),
                                        SMatrix{2,1}(X̄),SMatrix{1,2}(invX̄))
    q = vcat(ri,rj)
    lncs,q
end
function NC2D2P(ri,rj,ro,θ::Number,ṙo,ω)
    R = rotation_matrix(θ)
    lncs,q = NC2D2P(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    q̇ = vcat(ṙi,ṙj)
    lncs,SVector{4}(q),SVector{4}(q̇)
end

"""
返回三维刚性直杆自然坐标类。
$(TYPEDSIGNATURES)
"""
function NC3D2P(ri::AbstractVector{T},rj::AbstractVector{T},
               ro=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
               ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    r̄j = invR*(rj-ro)
    ū = r̄j-r̄i
    X̄ = ū
    invX̄ = pinv(X̄)
    lncs = LNC3D2P(SVector{3}(r̄i),SVector{3}(r̄j),
                                        SMatrix{3,1}(X̄),SMatrix{1,3}(invX̄))
    q = vcat(ri,rj)
    lncs,q
end
function NC3D2P(ri,rj,ro,R::AbstractMatrix,ṙo,ω)
    lncs,q = NC3D2P(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    q̇ = vcat(ṙi,ṙj)
    lncs,SVector{6}(q),SVector{6}(q̇)
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
    lncs = LNC1P2V(SVector{2}(r̄i),SVector{2}(ū),SVector{2}(v̄),
                                        SMatrix{2,2}(X̄),SMatrix{2,2}(inv(X̄)))
    q = vcat(ri,u,v)
    lncs,q
end
function NC1P2V(ri,ro,θ,ṙo,ω)
    R = rotation_matrix(θ)
    lncs,q = NC1P2V(ri,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×q[3:4]
    v̇ = ω×q[5:6]
    q̇ = vcat(ṙi,u̇,v̇)
    lncs,SVector{6}(q),SVector{6}(q̇)
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
    lncs = LNC2P1V(SVector{2}(r̄i),SVector{2}(r̄j),SVector{2}(v̄),
                                        SMatrix{2,2}(X̄),SMatrix{2,2}(inv(X̄)))
    q = vcat(ri,rj,v)
    lncs,q
end
function NC2P1V(ri,rj,ro,θ,ṙo,ω)
    R = rotation_matrix(θ)
    lncs,q = NC2P1V(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×q[5:6]
    q̇ = vcat(ṙi,ṙj,v̇)
    lncs,SVector{6}(q),SVector{6}(q̇)
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
    lncs = LNC3P(SVector{2}(r̄i),SVector{2}(r̄j),SVector{2}(r̄k),
                                        SMatrix{2,2}(X̄),SMatrix{2,2}(inv(X̄)))
    q = vcat(ri,rj,rk)
    lncs,q
end
function NC3P(ri,rj,rk,ro,θ,ṙo,ω)
    R = rotation_matrix(θ)
    lncs,q = NC3P(ri,rj,rk,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    q̇ = vcat(ṙi,ṙj,ṙk)
    lncs,SVector{6}(q),SVector{6}(q̇)
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
    lncs = LNC1P3V(SVector{3}(r̄i),SVector{3}(ū),SVector{3}(v̄),
                                        SVector{3}(w̄),SMatrix{3,3}(X̄),SMatrix{3,3}(inv(X̄)))
    q = vcat(ri,u,v,w)
    lncs,q
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
    lncs = LNC1P3V(SVector{3}(r̄i),SVector{3}(ū),SVector{3}(v̄),
                                        SVector{3}(w̄),SMatrix{3,3}(X̄),SMatrix{3,3}(inv(X̄)))
    q = vcat(ri,u,v,w)
    lncs,q
end
function NC1P3V(ri,ro,R,ṙo,ω)
    lncs,q = NC1P3V(ri,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×q[4:6]
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    lncs,SVector{12}(q),SVector{12}(q̇)
end
function NC1P3V(ri,u,v,w,ro,R,ṙo,ω)
    lncs,q = NC1P3V(ri,u,v,w,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×u
    v̇ = ω×v
    ẇ = ω×w
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    lncs,SVector{12}(q),SVector{12}(q̇)
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
    lncs = LNC2P2V(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(v̄),
                                        SVector{3}(w̄),SMatrix{3,3}(X̄),SMatrix{3,3}(inv(X̄)))
    q = vcat(ri,rj,v,w)
    lncs,q
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
    lncs = LNC2P2V(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(v̄),
                                        SVector{3}(w̄),SMatrix{3,3}(X̄),SMatrix{3,3}(inv(X̄)))
    q = vcat(ri,rj,v,w)
    lncs,q
end
function NC2P2V(ri,rj,ro,R,ṙo,ω)
    lncs,q = NC2P2V(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    lncs,SVector{12}(q),SVector{12}(q̇)
end
function NC2P2V(ri,rj,v,w,ro,R,ṙo,ω)
    lncs,q = NC2P2V(ri,rj,v,w,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×v
    ẇ = ω×w
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    lncs,SVector{12}(q),SVector{12}(q̇)
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
    lncs = LNC3P1V(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(r̄k),
                                        SVector{3}(w̄),SMatrix{3,3}(X̄),SMatrix{3,3}(inv(X̄)))
    q = vcat(ri,rj,rk,w)
    lncs,q
end
function NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
    lncs,q = NC3P1V(ri,rj,rk,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,ṙk,ẇ)
    lncs,SVector{12}(q),SVector{12}(q̇)
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
    lncs = LNC4P(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(r̄k),
                                        SVector{3}(r̄l),SMatrix{3,3}(X̄),SMatrix{3,3}(inv(X̄)))
    q = vcat(ri,rj,rk,rl)
    lncs,q
end
function NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
    lncs,q = NC4P(ri,rj,rk,rl,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    ṙl = ṙo + ω×(rl-ro)
    q̇ = vcat(ṙi,ṙj,ṙk,ṙl)
    lncs,SVector{12}(q),SVector{12}(q̇)
end

"""
返回二维或三维局部自然坐标类的变形量，用于任意一列点。
$(TYPEDSIGNATURES)
"""
function get_deform(lncs::LNCMP)
    reduce(vcat,lncs.r̄ps)
end

# """
# 返回二维或三维局部自然坐标类的变形量，用于刚性直杆。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds get_deform(ū) = sqrt(ū⋅ū)
@inline @inbounds function get_deform(lncs::Union{LNC2D2P,LNC3D2P})
    @unpack r̄i,r̄j = lncs
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
@inline @inbounds function get_deform(lncs::LNC1P2V)
    @unpack ū,v̄ = lncs
    get_deform(ū,v̄)
end
@inline @inbounds function get_deform(lncs::LNC2P1V)
    @unpack r̄i,r̄j,v̄ = lncs
    ū = r̄j-r̄i
    get_deform(ū,v̄)
end
@inline @inbounds function get_deform(lncs::LNC3P)
    @unpack r̄i,r̄j,r̄k = lncs
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
@inline @inbounds function get_deform(lncs::LNC1P3V)
    @unpack ū,v̄,w̄ = lncs
    get_deform(ū,v̄,w̄)
end
@inline @inbounds function get_deform(lncs::LNC2P2V)
    @unpack r̄i,r̄j,v̄,w̄ = lncs
    ū = r̄j-r̄i
    get_deform(ū,v̄,w̄)
end
@inline @inbounds function get_deform(lncs::LNC3P1V)
    @unpack r̄i,r̄j,r̄k,w̄ = lncs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    get_deform(ū,v̄,w̄)
end
@inline @inbounds function get_deform(lncs::LNC4P)
    @unpack r̄i,r̄j,r̄k,r̄l = lncs
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
function make_Φ(lncs::LNC,Φi)
    deforms = get_deform(lncs::LNC)
    make_Φ(lncs,Φi,deforms)
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
返回二维或三维内部约束，用于任意一列点。
$(TYPEDSIGNATURES)
"""
function make_Φ(lncs::LNCMP,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        all = q - d
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

"""
返回二维或三维内部约束，用于刚性直杆。
$(TYPEDSIGNATURES)
"""
function make_Φ(lncs::LNC2D2P,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        xi,yi,xj,yj = q
        all = [(xj-xi)^2 + (yj-yi)^2 - d^2]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
function make_Φ(lncs::LNC3D2P,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        xi,yi,zi,xj,yj,zj = q
        all = [(xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2- d^2]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

"""
返回二维内部约束，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φ(lncs::LNC1P2V,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        u = @view q[3:4]
        v = @view q[5:6]
        all = [(u⋅u - d[1]^2)/2, (v⋅v - d[2]^2)/2, √2/2*(u⋅v - d[3])]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
function make_Φ(lncs::LNC2P1V,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj-ri
        v = @view q[5:6]
        all = [(u⋅u - d[1]^2)/2, (v⋅v - d[2]^2)/2, √2/2*(u⋅v - d[3])]
        # all = [√2/4/d[1]*(u⋅u - d[1]^2), √2/4/d[2]*(v⋅v - d[2]^2), sqrt(1/(2d[1]+d[2]))*(u⋅v - d[3]^2)]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
function make_Φ(lncs::LNC3P,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        ri = @view q[1:2]
        rj = @view q[3:4]
        rk = @view q[5:6]
        u = rj-ri
        v = rk-ri
        all = [(u⋅u - d[1]^2), (v⋅v - d[2]^2), (u⋅v - d[3])]
        all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

"""
返回三维内部约束，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φ(lncs::LNC1P3V,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        u = @view q[4:6]
        v = @view q[7:9]
        w = @view q[10:12]
        all = [u⋅u - d[1]^2, v⋅v - d[2]^2, w⋅w - d[3]^2, v⋅w - d[4], u⋅w - d[5], u⋅v - d[6]]
        @view all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
function make_Φ(lncs::LNC2P2V,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u = rj-ri
        v = @view q[7:9]
        w = @view q[10:12]
        all = [u⋅u - d[1]^2, v⋅v - d[2]^2, w⋅w - d[3]^2, v⋅w - d[4], u⋅w - d[5], u⋅v - d[6]]
        @view all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
function make_Φ(lncs::LNC3P1V,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        ri = @view q[1:3]
        rj = @view q[4:6]
        rk = @view q[7:9]
        u = rj-ri
        v = rk-ri
        w = @view q[10:12]
        all = [u⋅u - d[1]^2, v⋅v - d[2]^2, w⋅w - d[3]^2, v⋅w - d[4], u⋅w - d[5], u⋅v - d[6]]
        @view all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end
function make_Φ(lncs::LNC4P,Φi,deforms)
    @inline @inbounds function _inner_Φ(q,d)
        ri = @view q[1:3]
        rj = @view q[4:6]
        rk = @view q[7:9]
        rl = @view q[10:12]
        u = rj-ri
        v = rk-ri
        w = rl-ri
        all = [u⋅u - d[1]^2, v⋅v - d[2]^2, w⋅w - d[3]^2, v⋅w - d[4], u⋅w - d[5], u⋅v - d[6]]
        @view all[Φi]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

## Jacobians
"""
返回二维或三维雅可比矩阵，用于任意一列点。
$(TYPEDSIGNATURES)
"""
function make_Φq(lncs::LNCMP,uci,Φi)
    ncoords = get_ncoords(lncs)
    @inline @inbounds function inner_Φq(q)
        ret = Matrix(one(eltype(q))*I,ncoords,ncoords)
        @view ret[Φi,uci]
    end
end

"""
返回二维或三维雅可比矩阵，用于刚性直杆。
$(TYPEDSIGNATURES)
"""
function make_Φq(lncs::Union{LNC2D2P,LNC3D2P},uci,Φi)
    ndim = get_ndim(lncs)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:ndim]
        rj = @view q[ndim+1:2ndim]
        ret = zeros(eltype(q),1,2ndim)
        ret[1,     1: ndim] = -2 .*(rj.-ri)
        ret[1,ndim+1:2ndim] = -ret[1,1:ndim]
        @view ret[Φi,uci]
    end
end

"""
返回二维雅可比矩阵，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φq(lncs::LNC1P2V,uci,Φi)
    @inline @inbounds function inner_Φq(q)
        u = @view q[3:4]
        v = @view q[5:6]
        ret = zeros(eltype(q),3,6)
        ret[1,3:4] = u
        ret[2,5:6] = v
        ret[3,3:4] = √2/2*v
        ret[3,5:6] = √2/2*u
        @view ret[Φi,uci]
    end
end
function make_Φq(lncs::LNC2P1V,uci,Φi)
    # d = get_deform(lncs)
    # a = sqrt(1/(2d[2]^2+d[1]^2))
    # weights = zeros(Rational{Int64},3,6)
    # weights[1,1:4] .= 1//4
    # weights[2,3:4] .= 1//2
    # weights[3,1:6] .= 1//6
    # wsum = sqrt.(inv.(sum(weights[:,uci],dims=2)))
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj-ri
        v = @view q[5:6]
        ret = zeros(eltype(q),3,6)
        ret[1,1:2] =-u
        ret[1,3:4] = u
        ret[2,5:6] = v
        ret[3,1:2] =-√2/2*v
        ret[3,3:4] = √2/2*v
        ret[3,5:6] = √2/2*u
        @view ret[Φi,uci]
    end
end
function make_Φq(lncs::LNC3P,uci,Φi)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:2]
        rj = @view q[3:4]
        rk = @view q[5:6]
        u = rj-ri
        v = rk-ri
        ret = zeros(eltype(q),3,6)
        ret[1,1:2] = -2u
        ret[1,3:4] =  2u
        ret[2,1:2] = -2v
        ret[2,5:6] =  2v
        ret[3,1:2] =  -v-u
        ret[3,3:4] =   v
        ret[3,5:6] =   u
        @view ret[Φi,uci]
    end
end

"""
返回三维雅可比矩阵，用于任意形状刚体。
$(TYPEDSIGNATURES)
"""
function make_Φq(lncs::LNC1P3V,uci,Φi)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:3]
        u  = @view q[4:6]
        v  = @view q[7:9]
        w  = @view q[10:12]
        ret = zeros(eltype(q), 6, 12)
        ret[1,4:6]   = 2u
        ret[2,7:9]   = 2v
        ret[3,10:12] = 2w

        ret[4,7:9]   = w
        ret[4,10:12] = v

        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,4:6] =  v
        ret[6,7:9] =  u

        @view ret[Φi,uci]
    end
end
function make_Φq(lncs::LNC2P2V,uci,Φi)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u  = rj-ri
        v  = @view q[7:9]
        w  = @view q[10:12]
        ret = zeros(eltype(q), 6, 12)
        ret[1,1:3] = -2u
        ret[1,4:6] =  2u

        ret[2,7:9]   = 2v
        ret[3,10:12] = 2w

        ret[4,7:9]   = w
        ret[4,10:12] = v

        ret[5,1:3] = -w
        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,1:3] = -v
        ret[6,4:6] =  v
        ret[6,7:9] =  u

        @view ret[Φi,uci]
    end
end
function make_Φq(lncs::LNC3P1V,uci,Φi)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u  = rj-ri
        rk = @view q[7:9]
        v  = rk-ri
        w  = @view q[10:12]
        ret = zeros(eltype(q), 6, 12)
        ret[1,1:3] = -2u
        ret[1,4:6] =  2u

        ret[2,1:3] = -2v
        ret[2,7:9] =  2v

        ret[3,10:12] = 2w

        ret[4,1:3]   =-w
        ret[4,7:9]   = w
        ret[4,10:12] = v

        ret[5,1:3] = -w
        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,1:3] = -v-u
        ret[6,4:6] =  v
        ret[6,7:9] =  u

        @view ret[Φi,uci]
    end
end
function make_Φq(lncs::LNC4P,uci,Φi)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u  = rj-ri
        rk = @view q[7:9]
        v  = rk-ri
        rl = @view q[10:12]
        w  = rl-ri
        ret = zeros(eltype(q), 6, 12)
        ret[1,1:3] = -2u
        ret[1,4:6] =  2u

        ret[2,1:3] = -2v
        ret[2,7:9] =  2v

        ret[3,1:3]   = -2w
        ret[3,10:12] =  2w

        ret[4,1:3]   =-w-v
        ret[4,7:9]   = w
        ret[4,10:12] = v

        ret[5,1:3] = -w-u
        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,1:3] = -v-u
        ret[6,4:6] =  v
        ret[6,7:9] =  u

        @view ret[Φi,uci]
    end
end

# Transformations
"""
返回转换矩阵C。
$(TYPEDSIGNATURES)
"""
function make_c(lncs)
    (;r̄i,invX̄) = lncs
    function c(r̄)
        invX̄*(r̄-r̄i)
    end
end
function make_C(lncs::LNC)
    function C(c)
        ndim = get_ndim(lncs)
        ncoords = get_ncoords(lncs)
        conversion = get_conversion(lncs)
        C_raw = hcat(1,transpose(c))
        SMatrix{ndim,ncoords}(kron(C_raw,make_I(Int,ndim))*conversion)
    end
end

## Transformations: Many points
function make_c(lncs::LNCMP)
    (;r̄ps) = lncs
    function c(r̄)
        dts = [norm(r̄-r̄p) for r̄p in eachcol(r̄ps)]
        c = zeros(eltype(r̄),size(r̄ps,2))
        c[argmin(dts)] = 1
        c
     end
end

function make_C(lncs::LNCMP)
    function C(c)
        ndim = get_ndim(lncs)
        ncoords = get_ncoords(lncs)
        conversion = get_conversion(lncs)
        C_raw = transpose(c)
        SMatrix{ndim,ncoords}(kron(C_raw,make_I(Int,ndim))*conversion)
    end
end

# CoordinateFunctions
"""
封装有函数的自然坐标类。
$(TYPEDEF)
"""
struct CoordinateFunctions{lncsType,
                    cT<:Function,CT<:Function,
                    ΦT<:Function,ΦqT<:Function,
                    ∂Aᵀλ∂qT,∂Aq̇∂qT}
    lncs::lncsType
    c::cT
    C::CT
    Φ::ΦT
    Φq::ΦqT
    ∂Aᵀλ∂q::∂Aᵀλ∂qT
    ∂Aq̇∂q::∂Aq̇∂qT
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
function get_unconstrained_indices(lncs::LNC,constrained_index)
    deleteat!(collect(1:get_ncoords(lncs)),constrained_index)
end

"""
封装有函数的自然坐标类构造子。
$(TYPEDSIGNATURES)
"""
function CoordinateFunctions(lncs,uci,Φi)
    c = make_c(lncs)
    C = make_C(lncs)
    Φ = make_Φ(lncs,Φi)
    Φq = make_Φq(lncs,uci,Φi)
    nq = get_ncoords(lncs)
    nuc = length(uci)
    nΦ = length(Φi)
    ∂Aᵀλ∂q = make_∂Aᵀλ∂q_forwarddiff(Φq,nq,nuc)
    ∂Aq̇∂q = make_∂Aq̇∂q_forwarddiff(Φq,nuc,nΦ)
    CoordinateFunctions(lncs,c,C,Φ,Φq,∂Aᵀλ∂q,∂Aq̇∂q)
end


# Mass matrices
Ī2J̄(::LNC2D2P,Ī)  = Ī
Ī2J̄(::LNC3D2P,Ī)  = Ī
Ī2J̄(::LNC2D6C,Ī)  = Ī
Ī2J̄(::LNC3D12C,Ī) = 1/2*tr(Ī)*I-Ī

function Īg2az(lncs,m,Īg,r̄g)
    (;r̄i,invX̄) = lncs
    a = invX̄*(r̄g-r̄i)
    J̄g = Ī2J̄(lncs,Īg)
    J̄o = J̄g + m*r̄g*transpose(r̄g)
    J̄i = J̄o - m*r̄i*transpose(r̄g) - m*r̄g*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*J̄i*transpose(invX̄)
    #@assert issymmetric(z)
    a,Symmetric(z)
end

## Mass matrices: standard
function make_M(cf::CoordinateFunctions,m::T,Īg,r̄g) where {T} # ami (area moment of inertia tensor)
    (;lncs) = cf
    ndim = get_ndim(lncs)
    nld = get_nlocaldim(lncs)
    ncoords = get_ncoords(lncs)
    a,z = Īg2az(lncs,m,Īg,r̄g)
    M_raw = zeros(T,1+nld,1+nld)
    M_raw[1,1] = m
    M_raw[2:1+nld,1] = m*a
    M_raw[1,2:1+nld] = M_raw[2:1+nld,1]
    M_raw[2:1+nld,2:1+nld] .= z
    M_std = kron(M_raw,make_I(T,ndim))
    Y_nonstd = get_conversion(lncs)
    M = SMatrix{ncoords,ncoords}(transpose(Y_nonstd)*M_std*Y_nonstd)
end

function make_M(cf::CoordinateFunctions{<:LNCMP},m::T,Īg,r̄g) where {T} # ami (area moment of inertia tensor)
    ncoords = get_ncoords(cf.lncs)
    M = SMatrix{ncoords,ncoords}(Matrix(m*I,ncoords,ncoords))
end

make_X(q::AbstractVector,lncs::LNC) = make_X(lncs,q)

function make_X(lncs::LNC,q::AbstractVector)
    Y = get_conversion(lncs)
    qstd = Y*q
    ndim = get_ndim(lncs)
    reshape(qstd[ndim+1:end],ndim,:)
end

find_R(q::AbstractVector, lncs::LNC) = find_R(lncs,q)
find_R(lncs::LNCMP,q::AbstractVector) = Matrix(one(eltype(q))*I,get_ndim(lncs),get_ndim(lncs))
function find_R(lncs::LNC,q::AbstractVector)
    X = make_X(lncs,q)
    (;invX̄) = lncs
    ndim = get_ndim(lncs)
    if lncs isa LNC3D2P
        (;r̄i,r̄j) = lncs
        ū = r̄j - r̄i
        v̄,w̄ = HouseholderOrthogonalization(ū)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u = rj -ri
        v,w = HouseholderOrthogonalization(u)
        R = SMatrix{ndim,ndim}([u;;v;;w]*inv([ū;;v̄;;w̄]))
    elseif lncs isa LNC2D2P
        (;r̄i,r̄j) = lncs
        ū = r̄j - r̄i
        v̄ = rotation_matrix(π/2)*ū
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj - ri
        v = rotation_matrix(π/2)*u
        R = SMatrix{ndim,ndim}([u;;v]*inv([ū;;v̄]))
    else
        R = SMatrix{ndim,ndim}(X*invX̄)
    end
    return R
end

find_ω(q::AbstractVector,q̇::AbstractVector,lncs::LNC3D) = find_ω(lncs,q,q̇)
function find_ω(lncs::LNCMP,q::AbstractVector,q̇::AbstractVector)
    ndim = get_ndim(lncs)
    if ndim == 3
        return zeros(eltype(q),ndim)
    else
        return zeros(eltype(q),1)
    end
end

function find_ω(lncs::LNC,q::AbstractVector,q̇::AbstractVector)
    Ẋ = make_X(lncs,q̇)
    X = make_X(lncs,q)
    ndim = get_ndim(lncs)
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

function get_uvw(lncs,q)
    if lncs isa LNC1P3V
        u = q[4:6]
        v = q[7:9]
        w = q[10:12]
    elseif lncs isa LNC2P2V
        ri = q[1:3]
        rj = q[4:6]
        v = q[7:9]
        w = q[10:12]
        u = rj - ri
    elseif lncs isa LNC3P1V
        ri = q[1:3]
        rj = q[4:6]
        rk = q[7:9]
        w = q[10:12]
        u = rj - ri
        v = rk - ri
    elseif lncs isa LNC4P
        ri = q[1:3]
        rj = q[4:6]
        rk = q[7:9]
        rl = q[10:12]
        u = rj - ri
        v = rk - ri
        w = rl - ri
    end
    SVector{3}(u),SVector{3}(v),SVector{3}(w)
end

function get_uv(lncs,q)
    if lncs isa LNC1P2V
        u = q[3:4]
        v = q[5:6]
    elseif lncs isa LNC2P1V
        ri = q[1:2]
        rj = q[3:4]
        v = q[5:6]
        u = rj - ri
    elseif lncs isa LNC3P
        ri = q[1:2]
        rj = q[3:4]
        rk = q[5:6]
        u = rj - ri
        v = rk - ri
    end
    SVector{2}(u),SVector{2}(v)
end

end
