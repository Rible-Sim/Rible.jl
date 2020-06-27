module NaturalCoordinates
using LinearAlgebra
using StaticArrays
using Parameters

abstract type BasicPoints end
make_I(T,N) = Matrix(one(T)*I,N,N)

struct BasicPoints2P{T} <: BasicPoints
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
end

struct BasicPoints1P1V{T} <: BasicPoints
    r̄i::SArray{Tuple{2},T,1,2}
    ū::SArray{Tuple{2},T,1,2}
end

function BasicPoints2P(Lij::Real)
    o = zero(Lij)
    r̄i = SVector(o,o)
    r̄j = SVector(Lij,o)
    BasicPoints2P(r̄i,r̄j)
end


function transform_to_1P1V(bps::BasicPoints2P{T}) where T
    I2 = make_I(T,2)
    Lij = norm(bps.r̄i-bps.r̄j)
    V = kron(
        [1     0;
        -1/Lij 1/Lij],
        I2
    )
end

struct BasicPoints1P3V{T} <: BasicPoints
    r̄i::SArray{Tuple{3},T,1,3}
    ū::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
end

struct BasicPoints2P2V{T} <: BasicPoints
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
end

struct BasicPoints3P1V{T} <: BasicPoints
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    w̄::SArray{Tuple{3},T,1,3}
end

struct BasicPoints4P{T} <: BasicPoints
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    r̄l::SArray{Tuple{3},T,1,3}
end

# Constructors
function BP1P3V(ri::AbstractVector{T},
               ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
               ) where T
    o = zero(T)
    i = one(T)
    u = SVector(i,o,o)
    v = SVector(o,i,o)
    w = SVector(o,o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro) #
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    bps = BasicPoints1P3V(SVector{3}(r̄i),SVector{3}(ū),
                          SVector{3}(v̄), SVector{3}(w̄))
    q = vcat(ri,u,v,w)
    bps,q
end

function BP1P3V(ri,ro,R,ṙo,ω)
    bps,q = BP1P3V(ri,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    u̇ = ω×q[4:6]
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,u̇,v̇,ẇ)
    bps,q,q̇
end

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

function BP2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    u = rj-ri
    v,w = HouseholderOrthogonalization(u)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    v̄ = invR*v
    w̄ = invR*w
    bps = BasicPoints2P2V(SVector{3}(r̄i),SVector{3}(r̄j),
                          SVector{3}(v̄), SVector{3}(w̄))
    q = vcat(ri,rj,v,w)
    bps,q
end

function BP2P2V(ri,rj,ro,R,ṙo,ω)
    bps,q = BP2P2V(ri,rj,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    v̇ = ω×q[7:9]
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,v̇,ẇ)
    bps,q,q̇
end

function BP3P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    u = rj-ri
    v = rk-ri
    w_raw = cross(u,v)
    w = w_raw/norm(w_raw)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    w̄ = invR*w
    bps = BasicPoints3P1V(SVector{3}(r̄i),SVector{3}(r̄j),
                          SVector{3}(r̄k),SVector{3}(w̄))
    q = vcat(ri,rj,rk,w)
    bps,q
end

function BP3P1V(ri,rj,rk,ro,R,ṙo,ω)
    bps,q = BP3P1V(ri,rj,rk,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    ẇ = ω×q[10:12]
    q̇ = vcat(ṙi,ṙj,ṙk,ẇ)
    bps,q,q̇
end

function BP4P(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},rl::AbstractVector{T},
                       ro=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-ro)
    r̄j = invR*(rj-ro)
    r̄k = invR*(rk-ro)
    r̄l = invR*(rl-ro)
    bps = BasicPoints4P(SVector{3}(r̄i),SVector{3}(r̄j),
                        SVector{3}(r̄k),SVector{3}(r̄l))
    q = vcat(ri,rj,rk,rl)
    bps,q
end

function BP4P(ri,rj,rk,rl,ro,R,ṙo,ω)
    bps,q = BP4P(ri,rj,rk,rl,ro,R)
    ṙi = ṙo + ω×(ri-ro)
    ṙj = ṙo + ω×(rj-ro)
    ṙk = ṙo + ω×(rk-ro)
    ṙl = ṙo + ω×(rl-ro)
    q̇ = vcat(ṙi,ṙj,ṙk,ṙl)
    bps,q,q̇
end

function transform_to_2P2V(bps::BasicPoints3P1V{T}) where T
    I3 = make_I(T,3)
    Lij = norm(bps.r̄i-bps.r̄j)
    V_raw = [1      0       0        0;
             0      1       0        0;
             0      0       0        1;
             0     -1/Lij   1/Lij    0]
    V = kron(V_raw,I3)
end

function transform_to_2P2V(bps::BasicPoints4P{T}) where T
    I3 = make_I(T,3)
    Lij = norm(bps.r̄i-bps.r̄j)
    V_raw = [1      0       0        0;
             0      1       0        0;
            -1/Lij  0       1/Lij   0;
             0     -1/Lij   0        1/Lij]
    V = kron(V_raw,I3)
end


function make_X̄(bps::BasicPoints2P{T}) where T
    @unpack r̄i,r̄j = bps
    ū = r̄j-r̄i
    X̄_raw = hcat(ū,[-ū[2],ū[1]])
    X̄ = SMatrix{2,2}(X̄_raw)
end

function make_X̄(bps::BasicPoints1P1V{T}) where T
    @unpack ū = bps
    X̄_raw = hcat(ū,[-ū[2],ū[1]])
    X̄ = SMatrix{2,2}(X̄_raw)
end

function make_X̄(bps::BasicPoints1P3V{T}) where T
    @unpack r̄i,ū,v̄,w̄ = bps
    X̄_raw = hcat(ū,v̄,w̄)
    X̄ = SMatrix{3,3}(X̄_raw)
end

function make_X̄(bps::BasicPoints2P2V{T}) where T
    @unpack r̄i,r̄j,v̄,w̄ = bps
    ū = r̄j-r̄i
    X̄_raw = hcat(ū,v̄,w̄)
    X̄ = SMatrix{3,3}(X̄_raw)
end

function make_X̄(bps::BasicPoints3P1V{T}) where T
    @unpack r̄i,r̄j,r̄k,w̄ = bps
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    X̄_raw = hcat(ū,v̄,w̄)
    X̄ = SMatrix{3,3}(X̄_raw)
end

function make_X̄(bps::BasicPoints4P{T}) where T
    @unpack r̄i,r̄j,r̄k,r̄l = bps
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = r̄l-r̄i
    X̄_raw = hcat(ū,v̄,w̄)
    X̄ = SMatrix{3,3}(X̄_raw)
end

function make_Φ(bps::BasicPoints2P)
    ū = bps.r̄j-bps.r̄i
    u_square = ū⋅ū
    @inline @inbounds function inner_Φ(q)
        xi,yi,xj,yj = q
        (xj-xi)^2 + (yj-yi)^2 - u_square
    end
end

function make_Φ(bps::BasicPoints1P1V)
    ū = bps.ū
    u_square = ū⋅ū
    @inline @inbounds function inner_Φ(q)
        u1,u2 = q[3],q[4]
        u1^2 + u2^2 - u_square
    end
end

function make_Φ(bps::BasicPoints1P3V)
    @unpack r̄i,ū,v̄,w̄ = bps
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

function make_Φ(bps::BasicPoints2P2V)
    @unpack r̄i,r̄j,v̄,w̄ = bps
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

function make_Φ(bps::BasicPoints3P1V)
    @unpack r̄i,r̄j,r̄k,w̄ = bps
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

function make_Φ(bps::BasicPoints4P)
    @unpack r̄i,r̄j,r̄k,r̄l = bps
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

function make_Φq(bps::BasicPoints2P)
    @inline @inbounds function inner_Φq(q)
        xi,yi,xj,yj = q
        ret = zeros(eltype(q),1,4)
        ret[1,1] = 2(xi-xj)
        ret[1,2] = 2(yi-yj)
        ret[1,3] = 2(xj-xi)
        ret[1,4] = 2(yj-yi)
        ret
    end
end

function make_Φq(bps::BasicPoints1P3V)
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

function make_Φq(bps::BasicPoints2P2V)
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

function make_Φq(bps::BasicPoints3P1V)
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

function make_Φq(bps::BasicPoints4P)
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

function make_c(bps,invX̄)
    @unpack r̄i = bps
    function c(r̄)
        invX̄*(r̄-r̄i)
    end
end

function make_C(bps::BasicPoints2P)
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

function make_C(bps::BasicPoints1P3V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

function make_C(bps::BasicPoints2P2V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1]  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

function make_C(bps::BasicPoints3P1V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1]-c[2]  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

function make_C(bps::BasicPoints4P{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1]-c[2]-c[3]  c[1]  c[2]  c[3]]
        SMatrix{3,12}(kron(C_raw,I3))
    end
end

struct CoordinateFunctions{bpsType,XT,cT,CT,ΦT,ΦqT}
    bps::bpsType
    invX̄::XT
    c::cT
    C::CT
    Φ::ΦT
    Φq::ΦqT
end

function CoordinateFunctions(bps)
    X̄ = make_X̄(bps)
    invX̄ = inv(X̄)
    c = make_c(bps,invX̄)
    C = make_C(bps)
    Φ = make_Φ(bps)
    Φq = make_Φq(bps)
    CoordinateFunctions(bps,invX̄,c,C,Φ,Φq)
end

function make_M(cf::CoordinateFunctions{BasicPoints1P1V{T},XT,cT,CT,ΦT,ΦqT},
                m::T,polar::T,r̄G) where {T,XT,cT,CT,ΦT,ΦqT}
    @unpack invX̄,bps = cf
    @unpack r̄i,ū = bps
    u_square = ū⋅ū
    a = invX̄*(r̄G-r̄i)
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

function make_M(cf::CoordinateFunctions{BasicPoints2P{T},XT,cT,CT,ΦT,ΦqT},
                m::T,polar::T,r̄G) where {T,XT,cT,CT,ΦT,ΦqT}
    @unpack invX̄,bps = cf
    @unpack r̄i,r̄j = bps
    ū = r̄j-r̄i
    u_square = ū⋅ū
    a = invX̄*(r̄G-r̄i)
    z = polar/u_square
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

@inline @inbounds function inertia2z(m,inertia_o::AbstractMatrix{T},
                                        r̄G,r̄i,invX̄) where T
    Jo = -MMatrix(inertia_o)
    Jo[1,1] = (inertia_o[2,2] + inertia_o[3,3] - inertia_o[1,1])/2
    Jo[2,2] = (inertia_o[1,1] + inertia_o[3,3] - inertia_o[2,2])/2
    Jo[3,3] = (inertia_o[1,1] + inertia_o[2,2] - inertia_o[3,3])/2
    Ji = Jo -
         m*r̄i*transpose(r̄G) -
         m*r̄G*transpose(r̄i) +
         m*r̄i*transpose(r̄i)
    z = invX̄*Ji*transpose(invX̄)
    Symmetric(z)
end

function compute_a_z(mass,inertia_o,r̄G,cf)
    @unpack bps,invX̄ = cf
    @unpack r̄i = bps
    a = invX̄*(r̄G-r̄i)
    z = inertia2z(mass,inertia_o,r̄G,r̄i,invX̄)
    a,z
end

function make_M(cf::CoordinateFunctions{BasicPoints1P3V{T},XT,cT,CT,ΦT,ΦqT},
                m::T,inertia::AbstractMatrix{T},r̄G) where {T,XT,cT,CT,ΦT,ΦqT}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄G)
    a,z = compute_a_z(m,inertia_o,r̄G,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m
    M_raw[2:4,1] = m*a
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function make_M(cf::CoordinateFunctions{BasicPoints2P2V{T},XT,cT,CT,ΦT,ΦqT},
                m::T,inertia::AbstractMatrix{T},r̄G) where {T,XT,cT,CT,ΦT,ΦqT}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄G)^2
    a,z = compute_a_z(m,inertia_o,r̄G,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m-2m*a[1]+z[1,1]
    M_raw[2:4,1] = m*a-z[1:3,1]
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function make_M(cf::CoordinateFunctions{BasicPoints3P1V{T},XT,cT,CT,ΦT,ΦqT},
                m::T,inertia::AbstractMatrix{T},r̄G) where {T,XT,cT,CT,ΦT,ΦqT}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄G)
    a,z = compute_a_z(m,inertia_o,r̄G,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m-2m*a[1]-2m*a[2]+
                     2z[1,2]+
                      z[1,1]+z[2,2]
    M_raw[2:4,1] = m*a-z[1:3,1]-z[1:3,2]
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function make_M(cf::CoordinateFunctions{BasicPoints4P{T},XT,cT,CT,ΦT,ΦqT},
                m::T,inertia::AbstractMatrix{T},r̄G) where {T,XT,cT,CT,ΦT,ΦqT}
    I3 = make_I(T,3)
    inertia_o = inertia - m*skew(r̄G)
    a,z = compute_a_z(m,inertia_o,r̄G,cf)
    M_raw = zeros(T,4,4)
    M_raw[1,1] = m-2m*a[1]-2m*a[2]-2m*a[3]+
                     2z[1,2]+2z[1,3]+2z[2,3]+
                      z[1,1]+ z[2,2]+ z[3,3]
    M_raw[2:4,1] = m*a-z[1:3,1]-z[1:3,2]-z[1:3,3]
    M_raw[1,2:4] = M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = SMatrix{12,12}(kron(M_raw,I3))
end

function skew(w)
    w1,w2,w3 = w
    o = zero(w1)
    [o -w3 w2;
     w3 o -w1;
    -w2 w1 o]
end
end
