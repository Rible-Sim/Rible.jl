module NaturalCoordinates
using LinearAlgebra
using StaticArrays
using Parameters

make_I(T,N) = Matrix(one(T)*I,N,N)

struct BasicPoints2P{T}
    r̄i::SArray{Tuple{2},T,1,2}
    r̄j::SArray{Tuple{2},T,1,2}
end

struct BasicPoints1P1V{T}
    r̄i::SArray{Tuple{2},T,1,2}
    ū::SArray{Tuple{2},T,1,2}
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

struct BasicPoints2P2V{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    ū::SArray{Tuple{3},T,1,3}
    v̄::SArray{Tuple{3},T,1,3}
end

struct BasicPoints4P{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    r̄l::SArray{Tuple{3},T,1,3}
end

struct BasicPoints3P1V{T}
    r̄i::SArray{Tuple{3},T,1,3}
    r̄j::SArray{Tuple{3},T,1,3}
    r̄k::SArray{Tuple{3},T,1,3}
    u::SArray{Tuple{3},T,1,3}
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
    ret = Matrix{T}(undef,2,2)
    ret[:,1] .= r̄j - r̄i
    ret[1,2] = -ret[2,1]
    ret[2,2] =  ret[1,1]
    ret
end

function make_X̄(bps::BasicPoints1P1V{T}) where T
    @unpack ū = bps
    ret = Matrix{T}(undef,2,2)
    ret[:,1] .= ū
    ret[1,2] = -ret[2,1]
    ret[2,2] =  ret[1,1]
    ret
end

function make_X̄(bps::BasicPoints2P2V{T}) where T
    @unpack r̄i,r̄j,ū,v̄ = bps
    ret = Matrix{T}(undef,3,3)
    ret[:,1] .= r̄j .- r̄i
    ret[:,2] .= ū
    ret[:,3] .= v̄
    ret
end

function make_Φ(bps::BasicPoints2P)
    Lij = norm(bps.r̄i-bps.r̄j)
    function inner_Φ(q)
        xi,yi,xj,yj = q
        (xj-xi)^2 + (yj-yi)^2 - Lij^2
    end
end

function make_Φ(bps::BasicPoints1P1V)
    function inner_Φ(q)
        u1,u2 = q[3],q[4]
        u1^2 + u2 - 1
    end
end

function make_Φ(bps::BasicPoints2P2V)
    Lij2 = norm(bps.r̄i-bps.r̄j)^2
    function inner_Φ(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u  = @view q[7:9]
        v  = @view q[10:12]
        [
        (rj-ri)⋅(rj-ri) - Lij2,
        (rj-ri)⋅u,
        (rj-ri)⋅v,
        u⋅v,
        u⋅u - 1,
        v⋅v - 1
        ]
    end
end

function make_Φq(bps::BasicPoints2P)
    function inner_Φq(q)
        xi,yi,xj,yj = q
        ret = similar(q)
        ret[1] = 2(xi-xj)
        ret[2] = 2(yi-yj)
        ret[3] = 2(xj-xi)
        ret[4] = 2(yj-yi)
        ret
    end
end

function make_Φq(bps::BasicPoints2P2V)
    function inner_Φq(q)
        ri = @view q[1:3]
        rj = @view q[4:6]
        u  = @view q[7:9]
        v  = @view q[10:12]
        ret = zeros(eltype(q), 6, 12)
        rirj = rj-ri
        ret[1,1:3] = -2rirj
        ret[1,4:6] =  2rirj

        ret[2,1:3] = -u
        ret[2,4:6] =  u
        ret[2,7:9] = rirj

        ret[3,1:3] = -v
        ret[3,4:6] =  v
        ret[3,10:12] = rirj

        ret[4,7:9]   = v
        ret[4,10:12] = u

        ret[5,7:9]   = 2u
        ret[6,10:12] = 2v
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
        ret = Matrix{eltype(c)}(undef,2,4)
        ret[1,1] = 1-c[1]
        ret[1,2] =   c[2]
        ret[1,3] =   c[1]
        ret[1,4] =  -c[2]
        ret[2,1] =  -c[2]
        ret[2,2] = 1-c[1]
        ret[2,3] =   c[2]
        ret[2,4] =   c[1]
        ret
    end
end

function make_C(bps::BasicPoints2P2V{T}) where T
    I3 = make_I(T,3)
    function C(c)
        C_raw = [1-c[1], c[1], c[2], c[3]]
        C_ret = kron(C_raw,I3)
    end
end

struct CoordinateFunctions{bpsType,XT,CT,cT,ΦT,ΦqT}
    bps::bpsType
    invX̄::XT
    c::cT
    C::CT
    ΦT::ΦT
    ΦqT::ΦqT
end

function CoordinateFunctions(bps)
    X̄ = make_X̄(bps)
    invX̄ = inv(X̄)
    c = make_c(bps,invX̄)
    C = make_C(bps)
    Φ = make_Φ(bps)
    Φq = make_Φq(bps)
    NaturalCoordinates(bps,invX̄,c,C,Φ,Φq)
end

function make_M(cf::CoordinateFunctions{T,BasicPoints2P{T},CT,cT,ΦT,ΦqT},
                m::T,polar::T,r̄G) where {T,CT,cT,ΦT,ΦqT}
    @unpack invX̄,bps = cf
    @unpack r̄i,r̄j = bps
    Lij = norm(r̄i-r̄j)
    a = invX̄*r̄G
    z = polar/Lij^2
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
    M_ret = Matrix(Symmetric(M))
end


function inertia2z(inertia::AbstractVector{T},bps) where T
    invX̄ = bps.invX̄
    Ji = -inertia
    Ji[1,1] = (inertia[2,2] + inertia[3,3] - inertia[1,1])/2
    Ji[2,2] = (inertia[1,1] + inertia[3,3] - inertia[2,2])/2
    Ji[3,3] = (inertia[1,1] + inertia[2,2] - inertia[3,3])/2
    invX̄*Ji*transpose(invX̄)
end

function make_M(cf::CoordinateFunctions{T,BasicPoints2P2V{T},CT,cT,ΦT,ΦqT},
                m::T,inertia::AbstractVector{T},r̄G) where {T,CT,cT,ΦT,ΦqT}
    I3 = make_I(T,3)
    @unpack bps,invX̄ = cf
    @unpack r̄i = bps
    a = invX̄*(r̄G-r̄i)
    z = inertia2z(inertia,bps)
    M_raw = zeros(T,4,4)
    M_raw[1:4,1] .= [m-2m*a[1]+z[1,1],
                        m*a[1]-z[1,1],
                        m*a[2]-z[2,1],
                        m*a[3]-z[3,1]]
    M_raw[1,2:4] .= M_raw[2:4,1]
    M_raw[2:4,2:4] .= z
    M = kron(M_raw,I3)
end


const DEFAULT_2P2V = BasicPoints2P2V(
    SVector(0.0,0.0,0.0),
    SVector(1.0,0.0,0.0),
    SVector(0.0,1.0,0.0),
    SVector(0.0,0.0,1.0)
)

end
