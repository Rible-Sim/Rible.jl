struct ContactFrame{T}
    n::SArray{Tuple{3},T,1,3}
    u::SArray{Tuple{3},T,1,3}
    w::SArray{Tuple{3},T,1,3}
end

function ContactFrame(n)
    n /= norm(n)
    u,w = NaturalCoordinates.HouseholderOrthogonalization(n)
    ContactFrame(SVector{3}(n),SVector{3}(u),SVector{3}(w))
end

struct GeneralizedDirections{N,T}
    Dn::SArray{Tuple{N},T,1,N}
    Du::SArray{Tuple{N},T,1,N}
    Dw::SArray{Tuple{N},T,1,N}
end

function GeneralizedDirections{N}(cf,C) where N
    @unpack n,u,w = cf
    Dn = SVector{N}(transpose(n)*C)
    Du = SVector{N}(transpose(u)*C)
    Dw = SVector{N}(transpose(w)*C)
    GeneralizedDirections(Dn,Du,Dw)
end

struct FrictionCone{N,T}
    μ::T
    e::T
    cf::ContactFrame{T}
    gd::GeneralizedDirections{N,T}
end

mutable struct ApproxFrictionalImpulse{T}
    active::Bool
    Λn::T
    β::Vector{T}
end

ApproxFrictionalImpulse{T}(m) where T = ApproxFrictionalImpulse(false,zero(T),zeros(T,m))

struct ApproxFrictionalContact{T}
    ϵ::T
    μ::T
    frame::ContactFrame{T}
    m::Int
    e::Vector{T}
    d::Vector{SArray{Tuple{3},T,1,3}}
    D::Array{T,2}
    impulse::ApproxFrictionalImpulse{T}
end

function ApproxFrictionalContact(ϵ::T,μ::T,m::Int) where T
    frame = ContactFrame{T}()
    e = ones(m)
    d = [SVector{3,T}(cos(2π*i/m),sin(2π*i/m),0.0) for i = 0:m-1]
    D = Matrix{T}(undef,3,m)
    for i = 1:m
        D[:,i] .= d[i]
    end
    impulse = ApproxFrictionalImpulse{T}(m)
    ApproxFrictionalContact(
    ϵ,μ,frame,m,e,d,D,impulse
    )
end
