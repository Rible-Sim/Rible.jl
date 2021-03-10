module NaturalCoordinates
using LinearAlgebra
using StaticArrays
using Parameters
using ForwardDiff

abstract type LocalNaturalCoordinates{T} end
abstract type LocalNaturalCoordinates2D{T} <: LocalNaturalCoordinates{T} end
abstract type LocalNaturalCoordinates3D{T} <: LocalNaturalCoordinates{T} end
abstract type LocalNaturalCoordinates2D4C{T} <: LocalNaturalCoordinates2D{T} end
abstract type LocalNaturalCoordinates2D6C{T} <: LocalNaturalCoordinates2D{T} end
abstract type LocalNaturalCoordinates3D12C{T} <: LocalNaturalCoordinates3D{T} end

make_I(T,N) = SMatrix{N,N}(one(T)*I)
get_ncoords(::LocalNaturalCoordinates2D4C) = 4
get_ncoords(::LocalNaturalCoordinates2D6C) = 6
get_ncoords(::LocalNaturalCoordinates3D12C) = 12

function LinearAlgebra.cross(a::Number,b::AbstractVector)
    ret = similar(b)
    ret[1] = -a*b[2]
    ret[2] =  a*b[1]
    ret
end

rotation_matrix(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

@inline @inbounds function HouseholderOrthogonalization(n)
    n̄ = norm(n)
    h1 = max(n[1]-n̄,n[1]+n̄)
    h2 = n[2]
    h3 = n[3]
    h = [h1,h2,h3]
    H = I - 2h*transpose(h)/(transpose(h)*h)
    t = H[1:3,2]
    b = H[1:3,3]
    t,b
end

function skew(w)
    w1,w2,w3 = w
    o = zero(w1)
    @SMatrix [o -w3 w2;
              w3 o -w1;
             -w2 w1 o]
end

struct LocalNaturalCoordinates1P1V{T} <: LocalNaturalCoordinates2D4C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    ū::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

struct LocalNaturalCoordinates2P{T} <: LocalNaturalCoordinates2D4C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

struct LocalNaturalCoordinates1P2V{T} <: LocalNaturalCoordinates2D6C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    ū::SArray{Tuple{2},T,1,2}
    v̄::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

struct LocalNaturalCoordinates2P1V{T} <: LocalNaturalCoordinates2D6C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
    v̄::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

struct LocalNaturalCoordinates3P{T} <: LocalNaturalCoordinates2D6C{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
    r̄k::SArray{Tuple{2},T,1,2}
    X̄::SArray{Tuple{2,2},T,2,4}
    invX̄::SArray{Tuple{2,2},T,2,4}
end

struct LocalNaturalCoordinates1P3V{T} <: LocalNaturalCoordinates3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    ū::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

struct LocalNaturalCoordinates2P2V{T} <: LocalNaturalCoordinates3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

struct LocalNaturalCoordinates3P1V{T} <: LocalNaturalCoordinates3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

struct LocalNaturalCoordinates4P{T} <: LocalNaturalCoordinates3D12C{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    r̄l::SArray{Tuple{3},T,1,3}
    X̄::SArray{Tuple{3,3},T,2,9}
    invX̄::SArray{Tuple{3,3},T,2,9}
end

# Constructors

# 2D
function LocalNaturalCoordinates1P1V{T}(r̄i = SVector{2}(zeros(T,2))) where T
    ū = SVector(one(T),zero(T))
    X̄ = hcat(ū,[-ū[2],ū[1]])
    LocalNaturalCoordinates1P1V(SVector{2}(r̄i),SVector{2}(ū),SMatrix{2,2}(X̄),SMatrix{2,2}(inv(X̄)))
end

function LocalNaturalCoordinates2P(Lij::Real)
    o = zero(Lij)
    r̄i = SVector(o,o)
    r̄j = SVector(Lij,o)
    ū = r̄j-r̄i
    X̄ = hcat(ū,[-ū[2],ū[1]])
    LocalNaturalCoordinates2P(SVector{2}(r̄i),SVector{2}(r̄j),SMatrix{2,2}(X̄),SMatrix{2,2}(inv(X̄)))
end

function transform_to_1P1V(lncs::LocalNaturalCoordinates2P{T}) where T
    I2 = make_I(T,2)
    Lij = norm(lncs.r̄i-lncs.r̄j)
    V = kron(
        [1     0;
        -1/Lij 1/Lij],
        I2
    )
end

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
    lncs = LocalNaturalCoordinates1P2V(SVector{2}(r̄i),SVector{2}(ū),SVector{2}(v̄),
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

function NC2P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    u = rj-ri
    ū = invR*u
    v̄ = rotation_matrix(π/2)*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    lncs = LocalNaturalCoordinates2P1V(SVector{2}(r̄i),SVector{2}(r̄j),SVector{2}(v̄),
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

function NC3P(ri::AbstractVector{T},rj::AbstractVector{T},rk::AbstractVector{T},
              ro=SVector{2}(zeros(T,2)),R=SMatrix{2,2}(one(T)*I)) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    X̄ = hcat(ū,v̄)
    lncs = LocalNaturalCoordinates3P(SVector{2}(r̄i),SVector{2}(r̄j),SVector{2}(r̄k),
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

# 3D

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
    lncs = LocalNaturalCoordinates1P3V(SVector{3}(r̄i),SVector{3}(ū),SVector{3}(v̄),
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
    lncs = LocalNaturalCoordinates2P2V(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(v̄),
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
    lncs = LocalNaturalCoordinates2P2V(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(v̄),
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

function NC3P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    u = rj-ri
    v = rk-ri
    # w_raw = cross(u,v)
    # w = w_raw/norm(w_raw)
    w = u×v
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    lncs = LocalNaturalCoordinates3P1V(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(r̄k),
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
    lncs = LocalNaturalCoordinates4P(SVector{3}(r̄i),SVector{3}(r̄j),SVector{3}(r̄k),
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

function transform_to_2P2V(lncs::LocalNaturalCoordinates3P1V{T}) where T
    I3 = make_I(T,3)
    Lij = norm(lncs.r̄i-lncs.r̄j)
    V_raw = [1      0       0        0;
             0      1       0        0;
             0      0       0        1;
             0     -1/Lij   1/Lij    0]
    V = kron(V_raw,I3)
end

function transform_to_2P2V(lncs::LocalNaturalCoordinates4P{T}) where T
    I3 = make_I(T,3)
    Lij = norm(lncs.r̄i-lncs.r̄j)
    V_raw = [1      0       0        0;
             0      1       0        0;
            -1/Lij  0       1/Lij   0;
             0     -1/Lij   0        1/Lij]
    V = kron(V_raw,I3)
end

function make_Φ(lncs::LocalNaturalCoordinates1P1V)
    ū = lncs.ū
    u_square = ū⋅ū
    @inline @inbounds function inner_Φ(q)
        u = q[3:4]
        u⋅u - u_square
    end
end

function make_Φ(lncs::LocalNaturalCoordinates2P)
    @unpack r̄i,r̄j = lncs
    ū = r̄j-r̄i
    u_square = ū⋅ū
    @inline @inbounds function inner_Φ(q)
        xi,yi,xj,yj = q
        (xj-xi)^2 + (yj-yi)^2 - u_square
    end
end

@inline @inbounds function get_deform(lncs::LocalNaturalCoordinates1P2V)
    @unpack ū,v̄ = lncs
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    uv_dotprod = ū⋅v̄
    u_square,v_square,uv_dotprod
end

@inline @inbounds function get_deform(lncs::LocalNaturalCoordinates2P1V)
    @unpack r̄i,r̄j,v̄ = lncs
    ū = r̄j-r̄i
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    uv_dotprod = ū⋅v̄
    u_square,v_square,uv_dotprod
end

@inline @inbounds function get_deform(lncs::LocalNaturalCoordinates3P)
    @unpack r̄i,r̄j,r̄k = lncs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    uv_dotprod = ū⋅v̄
    u_square,v_square,uv_dotprod
end

function make_Φ(lncs::LocalNaturalCoordinates2D6C)
    deforms = get_deform(lncs)
    make_Φ(lncs,deforms)
end

function make_inner_Φ(func,deforms)
    @inline @inbounds function ret_func(q)
        u_square,v_square,uv_dotprod = deforms
        func(q,u_square,v_square,uv_dotprod)
    end
    @inline @inbounds function ret_func(q,u)
        u_square,v_square,uv_dotprod = u
        func(q,u_square,v_square,uv_dotprod)
    end
    ret_func
end

function make_Φ(lncs::LocalNaturalCoordinates1P2V,deforms)
    @inline @inbounds function _inner_Φ(q,u_square,v_square,uv_dotprod)
        u = @view q[3:4]
        v = @view q[5:6]
        [u⋅u - u_square, v⋅v - v_square, u⋅v - uv_dotprod]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

function make_Φ(lncs::LocalNaturalCoordinates2P1V,deforms)
    @inline @inbounds function _inner_Φ(q,u_square,v_square,uv_dotprod)
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj-ri
        v = @view q[5:6]
        [u⋅u - u_square, v⋅v - v_square, u⋅v - uv_dotprod]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

function make_Φ(lncs::LocalNaturalCoordinates3P,deforms)
    @inline @inbounds function _inner_Φ(q,u_square,v_square,uv_dotprod)
        ri = @view q[1:2]
        rj = @view q[3:4]
        rk = @view q[5:6]
        u = rj-ri
        v = rk-ri
        [u⋅u - u_square, v⋅v - v_square, u⋅v - uv_dotprod]
    end
    make_inner_Φ(_inner_Φ,deforms)
end

function make_Φ(lncs::LocalNaturalCoordinates1P3V)
    @unpack r̄i,ū,v̄,w̄ = lncs
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    w_square = w̄⋅w̄
    uv_dotprod = ū⋅v̄
    uw_dotprod = ū⋅w̄
    vw_dotprod = v̄⋅w̄
    @inline @inbounds function inner_Φ(q)
        ri = @view q[1:3]
        u  = @view q[4:6]
        v  = @view q[7:9]
        w  = @view q[10:12]
        [
        u⋅u - u_square,
        v⋅v - v_square,
        w⋅w - w_square,
        u⋅v - uv_dotprod,
        u⋅w - uw_dotprod,
        v⋅w - vw_dotprod
        ]
    end
end

function make_Φ(lncs::LocalNaturalCoordinates2P2V)
    @unpack r̄i,r̄j,v̄,w̄ = lncs
    ū = r̄j-r̄i
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    w_square = w̄⋅w̄
    uv_dotprod = ū⋅v̄
    uw_dotprod = ū⋅w̄
    vw_dotprod = v̄⋅w̄
    @inline @inbounds function inner_Φ(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        v  = @view q[7:9]
        w  = @view q[10:12]
        u = rj-ri
        [
        u⋅u - u_square,
        v⋅v - v_square,
        w⋅w - w_square,
        u⋅v - uv_dotprod,
        u⋅w - uw_dotprod,
        v⋅w - vw_dotprod
        ]
    end
end

function make_Φ(lncs::LocalNaturalCoordinates3P1V)
    @unpack r̄i,r̄j,r̄k,w̄ = lncs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    w_square = w̄⋅w̄
    uv_dotprod = ū⋅v̄
    uw_dotprod = ū⋅w̄
    vw_dotprod = v̄⋅w̄
    @inline @inbounds function inner_Φ(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        rk = @view q[7:9]
        w  = @view q[10:12]
        u = rj-ri
        v = rk-ri
        [
        u⋅u - u_square,
        v⋅v - v_square,
        w⋅w - w_square,
        u⋅v - uv_dotprod,
        u⋅w - uw_dotprod,
        v⋅w - vw_dotprod
        ]
    end
end

function make_Φ(lncs::LocalNaturalCoordinates4P)
    @unpack r̄i,r̄j,r̄k,r̄l = lncs
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = r̄l-r̄i
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    w_square = w̄⋅w̄
    uv_dotprod = ū⋅v̄
    uw_dotprod = ū⋅w̄
    vw_dotprod = v̄⋅w̄
    @inline @inbounds function inner_Φ(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        rk  = @view q[7:9]
        rl  = @view q[10:12]
        u = rj-ri
        v = rk-ri
        w = rl-ri
        [
        u⋅u - u_square,
        v⋅v - v_square,
        w⋅w - w_square,
        u⋅v - uv_dotprod,
        u⋅w - uw_dotprod,
        v⋅w - vw_dotprod
        ]
    end
end

function make_Φq(lncs::LocalNaturalCoordinates1P1V)
    @inline @inbounds function inner_Φq(q)
        ret = zeros(eltype(q),1,4)
        u = q[3:4]
        ret[1,3:4] =  2u
        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates2P)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj - ri
        ret = zeros(eltype(q),1,4)
        ret[1,1:2] = -2u
        ret[1,3:4] =  2u
        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates1P2V)
    @inline @inbounds function inner_Φq(q)
        u = @view q[3:4]
        v = @view q[5:6]
        ret = zeros(eltype(q),3,6)
        ret[1,3:4] = 2u
        ret[2,5:6] = 2v
        ret[3,3:4] = v
        ret[3,5:6] = u
        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates2P1V)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj-ri
        v = @view q[5:6]
        ret = zeros(eltype(q),3,6)
        ret[1,1:2] = -2u
        ret[1,3:4] = 2u
        ret[2,5:6] = 2v
        ret[3,1:2] = -v
        ret[3,3:4] = v
        ret[3,5:6] = u
        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates3P)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:2]
        rj = @view q[3:4]
        u = rj-ri
        rk = @view q[5:6]
        v = rk-ri
        ret = zeros(eltype(q),3,6)
        ret[1,1:2] = -2u
        ret[1,3:4] = 2u
        ret[2,1:2] = -2v
        ret[2,5:6] = 2v
        ret[3,1:2] = -v-u
        ret[3,3:4] = v
        ret[3,5:6] = u
        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates1P3V)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:3]
        u  = @view q[4:6]
        v  = @view q[7:9]
        w  = @view q[10:12]
        ret = zeros(eltype(q), 6, 12)
        ret[1,4:6] =  2u
        ret[2,7:9]   = 2v
        ret[3,10:12] = 2w

        ret[4,4:6] =  v
        ret[4,7:9] =  u

        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,7:9]   = w
        ret[6,10:12] = v

        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates2P2V)
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

        ret[4,1:3] = -v
        ret[4,4:6] =  v
        ret[4,7:9] =  u

        ret[5,1:3] = -w
        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,7:9]   = w
        ret[6,10:12] = v

        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates3P1V)
    @inline @inbounds function inner_Φq(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u  = rj-ri
        rk  = @view q[7:9]
        v  = rk-ri
        w  = @view q[10:12]
        ret = zeros(eltype(q), 6, 12)
        ret[1,1:3] = -2u
        ret[1,4:6] =  2u

        ret[2,1:3] = -2v
        ret[2,7:9] =  2v

        ret[3,10:12] = 2w

        ret[4,1:3] = -v-u
        ret[4,4:6] =  v
        ret[4,7:9] =  u

        ret[5,1:3] = -w
        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,1:3]   =-w
        ret[6,7:9]   = w
        ret[6,10:12] = v

        ret
    end
end

function make_Φq(lncs::LocalNaturalCoordinates4P)
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

        ret[4,1:3] = -v-u
        ret[4,4:6] =  v
        ret[4,7:9] =  u

        ret[5,1:3] = -w-u
        ret[5,4:6] =  w
        ret[5,10:12] = u

        ret[6,1:3]   =-w-v
        ret[6,7:9]   = w
        ret[6,10:12] = v

        ret
    end
end

function make_c(lncs)
    @unpack r̄i,invX̄ = lncs
    function c(r̄)
        invX̄*(r̄-r̄i)
    end
end

function make_C(lncs::LocalNaturalCoordinates1P1V)
    function C(c)
        C_raw = Matrix{eltype(c)}(undef,2,4)
        C_raw[1,1] =    1
        C_raw[1,2] =    0
        C_raw[1,3] =  c[1]
        C_raw[1,4] = -c[2]
        C_raw[2,1] =    0
        C_raw[2,2] =    1
        C_raw[2,3] =  c[2]
        C_raw[2,4] =  c[1]
        SMatrix{2,4}(C_raw)
    end
end

function make_C(lncs::LocalNaturalCoordinates2P)
    function C(c)
        C_raw = Matrix{eltype(c)}(undef,2,4)
        C_raw[1,1] = 1-c[1]
        C_raw[1,2] =   c[2]
        C_raw[1,3] =   c[1]
        C_raw[1,4] =  -c[2]
        C_raw[2,1] =  -c[2]
        C_raw[2,2] = 1-c[1]
        C_raw[2,3] =   c[2]
        C_raw[2,4] =   c[1]
        SMatrix{2,4}(C_raw)
    end
end

function make_C(lncs::LocalNaturalCoordinates1P2V{T}) where T
    I2 = make_I(T,2)
    function C(c)
        C_raw = [1  c[1]  c[2]]
        SMatrix{2,6}(kron(C_raw,I2))
    end
end

function make_C(lncs::LocalNaturalCoordinates2P1V{T}) where T
    I2 = make_I(T,2)
    function C(c)
        C_raw = [1-c[1]  c[1]  c[2]]
        SMatrix{2,6}(kron(C_raw,I2))
    end
end

function make_C(lncs::LocalNaturalCoordinates3P{T}) where T
    I2 = make_I(T,2)
    function C(c)
        C_raw = [1-c[1]-c[2]  c[1]  c[2]]
        SMatrix{2,6}(kron(C_raw,I2))
    end
end

function make_C(lncs::LocalNaturalCoordinates1P3V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

function make_C(lncs::LocalNaturalCoordinates2P2V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1]  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

function make_C(lncs::LocalNaturalCoordinates3P1V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1]-c[2]  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

function make_C(lncs::LocalNaturalCoordinates4P{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1]-c[2]-c[3]  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

struct CoordinateFunctions{lncsType,cT,CT,ΦT,ΦqT}
    lncs::lncsType
    c::cT
    C::CT
    Φ::ΦT
    Φq::ΦqT
end

function CoordinateFunctions(lncs)
    c = make_c(lncs)
    C = make_C(lncs)
    Φ = make_Φ(lncs)
    Φq = make_Φq(lncs)
    CoordinateFunctions(lncs,c,C,Φ,Φq)
end


function ∂Φqᵀ∂q_forwarddiff(cf::CoordinateFunctions{<:LocalNaturalCoordinates})
    function ∂Aᵀλ∂q(λ)
        function ATλ(q)
            transpose(cf.Φq(q))*λ
        end
        ForwardDiff.jacobian(ATλ,ones(eltype(λ),get_ncoords(cf.lncs)))
    end
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates1P1V{T}},
                m::T,polar::T,r̄g) where {T}
    @unpack lncs = cf
    @unpack r̄i,ū,invX̄ = lncs
    u_square = ū⋅ū
    a = invX̄*(r̄g-r̄i)
    z = polar/u_square
    M = zeros(T,4,4)
    M[1,1] = m
    M[1,2] = 0.0
    M[1,3] =  m*a[1]
    M[1,4] = -m*a[2]
    M[2,2] = m
    M[2,3] = m*a[2]
    M[2,4] = m*a[1]
    M[3,3] = z
    M[3,4] = 0.0
    M[4,4] = z
    M_ret = SMatrix{4,4}(Symmetric(M))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates2P{T}},
                m::T,polar::T,r̄g) where {T}
    @unpack lncs = cf
    @unpack r̄i,r̄j,invX̄ = lncs
    ū = r̄j-r̄i
    u_square = ū⋅ū
    a = invX̄*(r̄g-r̄i)
    polar_o = polar + m*(r̄g⋅r̄g)
    z = polar_o/u_square
    M = zeros(T,4,4)
    M[1,1] = m - 2m*a[1] + z
    M[1,2] = 0.0
    M[1,3] =  m*a[1] - z
    M[1,4] = -m*a[2]
    M[2,2] = m - 2m*a[1] + z
    M[2,3] = m*a[2]
    M[2,4] = m*a[1] - z
    M[3,3] = z
    M[3,4] = 0.0
    M[4,4] = z
    M_ret = SMatrix{4,4}(Symmetric(M))
end

inertia2J(inertia) = _inertia2J(Size(inertia),inertia)
_inertia2J(::Size{(2,2)},inertia) =     tr(inertia)*I-inertia
_inertia2J(::Size{(3,3)},inertia) = 1/2*tr(inertia)*I-inertia

@inline @inbounds function inertia2z(m,inertia_o::AbstractMatrix{T},r̄g,r̄i,invX̄) where T
    Jo = inertia2J(inertia_o)
    Ji = Jo - m*r̄i*transpose(r̄g) - m*r̄g*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*Ji*transpose(invX̄)
    Symmetric(z)
end

@inline @inbounds function Jg2az(m,Jg::AbstractMatrix{T},r̄g,cf) where T
    @unpack r̄i,invX̄ = cf.lncs
    a = invX̄*(r̄g-r̄i)
    Jo = Jg + m*r̄g*transpose(r̄g)
    Ji = Jo - m*r̄i*transpose(r̄g) - m*r̄g*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*Ji*transpose(invX̄)
    a,Symmetric(z)
end

function compute_a_z(mass,inertia_o,r̄g,cf)
    @unpack r̄i,invX̄ = cf.lncs
    a = invX̄*(r̄g-r̄i)
    z = inertia2z(mass,inertia_o,r̄g,r̄i,invX̄)
    a,z
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates1P2V{T}},
                m::T,ami::AbstractMatrix{T},r̄g) where {T}
                # ami (area moment of inertia tensor)
    I2 = make_I(T,2)
    Jg = inertia2J(ami)
    a,z = Jg2az(m,Jg,r̄g,cf)
    M_raw = zeros(T,3,3)
    M_raw[1,1] = m
    M_raw[2:3,1] = m*a
    M_raw[1,2:3] = M_raw[2:3,1]
    M_raw[2:3,2:3] .= z
    M = SMatrix{6,6}(kron(M_raw,I2))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates2P1V{T}},
                m::T,ami::AbstractMatrix{T},r̄g) where {T}
    I2 = make_I(T,2)
    # ami (area moment of inertia tensor)
    Jg = inertia2J(ami)
    a,z = Jg2az(m,Jg,r̄g,cf)
    M_raw = zeros(T,3,3)
    M_raw[1,1] = m-2m*a[1]+z[1,1]
    M_raw[2:3,1] = m*a-z[1:2,1]
    M_raw[1,2:3] = M_raw[2:3,1]
    M_raw[2:3,2:3] .= z
    M = SMatrix{6,6}(kron(M_raw,I2))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates3P{T}},
                m::T,ami::AbstractMatrix{T},r̄g) where {T}
    I2 = make_I(T,2)
    # ami (area moment of inertia tensor)
    Jg = inertia2J(ami)
    a,z = Jg2az(m,Jg,r̄g,cf)
    M_raw = zeros(T,3,3)
    M_raw[1,1] = m-2m*a[1]-2m*a[2]+2z[1,2]+z[1,1]+z[2,2]
    M_raw[2:3,1] = m*a-z[1:2,1]-z[1:2,2]
    M_raw[1,2:3] = M_raw[2:3,1]
    M_raw[2:3,2:3] .= z
    M = SMatrix{6,6}(kron(M_raw,I2))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates1P3V{T}},
                m::T,inertia::AbstractMatrix{T},r̄g) where {T}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄g)^2
    a,z = compute_a_z(m,inertia_o,r̄g,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m
    M_raw[2:4,1] = m*a
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates2P2V{T}},
                m::T,inertia::AbstractMatrix{T},r̄g) where {T}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄g)^2
    a,z = compute_a_z(m,inertia_o,r̄g,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m-2m*a[1]+z[1,1]
    M_raw[2:4,1] = m*a-z[1:3,1]
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates3P1V{T}},
                m::T,inertia::AbstractMatrix{T},r̄g) where {T}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄g)^2
    a,z = compute_a_z(m,inertia_o,r̄g,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m-2m*a[1]-2m*a[2]+
                     2z[1,2]+
                      z[1,1]+z[2,2]
    M_raw[2:4,1] = m*a-z[1:3,1]-z[1:3,2]
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function make_M(cf::CoordinateFunctions{LocalNaturalCoordinates4P{T}},
                m::T,inertia::AbstractMatrix{T},r̄g) where {T}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄g)^2
    a,z = compute_a_z(m,inertia_o,r̄g,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m-2m*a[1]-2m*a[2]-2m*a[3]+
                     2z[1,2]+2z[1,3]+2z[2,3]+
                      z[1,1]+ z[2,2]+ z[3,3]
    M_raw[2:4,1] = m*a-z[1:3,1]-z[1:3,2]-z[1:3,3]
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

end
