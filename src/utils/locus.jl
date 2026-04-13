mutable struct ContactState{N,T,L}
    active::Bool
    persistent::Bool
    gap::T
    axes::Axes{N,T,L}
    relative_velocity::SVector{N,T}
    force::SVector{N,T}
end

function ContactState(normal::StaticArray{Tuple{N},T}) where {N,T}
    active = false
    persistent = true
    gap = typemax(T)
    axes = Axes(normal)
    relative_velocity = zeros(SVector{N,T})
    force = zeros(SVector{N,T})
    ContactState(
        active,persistent,
        gap,axes,
        relative_velocity,force
    )
end

function reset!(contact_state::ContactState{N,T}) where {N,T}
    contact_state.active = false
    contact_state.persistent = true
    contact_state.gap = typemax(T)
    contact_state.relative_velocity = zeros(SVector{N,T})
    contact_state.force = zeros(SVector{N,T})
end

function activate!(contact_state::ContactState,gap)
    contact_state.gap = gap
    pre_active = contact_state.active
    contact_state.active = (gap<=0)
    contact_state.persistent = (pre_active == contact_state.active)
end

struct Locus{N,T,L}
    position::SVector{N,T}
    axes::Axes{N,T,L}
    friction_coefficient::T
    restitution_coefficient::T
end

function Locus(position::SVector{N,T},
        friction_coefficient::T = zero(T),restitution_coefficient::T = zero(T)
    ) where {N,T}
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = Axes(normal)
    Locus(
        position,axes,
        friction_coefficient,
        restitution_coefficient
    )
end

function Locus(position::SVector{N,T},normal::SVector{N,T},
        friction_coefficient::T = zero(T),restitution_coefficient::T = zero(T)
    ) where {N,T}
    axes = Axes(normal)
    Locus(
        position,axes,
        friction_coefficient,
        restitution_coefficient
    )
end

mutable struct LocusState{N,M,T,L}
    frame::CartesianFrame{N,M,T,L}
    force::MVector{N,T}
    torque::MVector{M,T}
    contact_state::ContactState{N,T}
    curvature::MVector{N,T}
end

function LocusState(position::StaticArray{Tuple{N},T},velocity::StaticArray{Tuple{N},T},
        force = zeros(MVector{N,T}),torque = zeros(MVector{2N-3,T})
    ) where {N,T}
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = Axes(normal)
    frame = CartesianFrame(position,axes,velocity)
    contact_state = ContactState(normal)
    curvature = MVector{N,T}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    LocusState(
        frame,
        force,torque,
        contact_state,
        curvature
    )
end

function LocusState(locus::Locus{N,T,L},origin_frame::CartesianFrame{N,M,T,L},
        force = zeros(MVector{N,T}),torque = zeros(MVector{M,T})
    ) where {N,M,T,L}
    relative_position = origin_frame.axes*locus.position
    position = origin_frame.position+relative_position
    velocity = origin_frame.velocity+origin_frame.angular_velocity×relative_position
    axes = origin_frame.axes*locus.axes
    frame = CartesianFrame(position,axes,velocity,origin_frame.angular_velocity)
    contact_state = ContactState(axes.normal)
    curvature = MVector{N,T}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    LocusState(
        frame,
        force,torque,
        contact_state,
        curvature
    )
end

activate!(locus_state::LocusState,gap) = activate!(locus_state.contact_state,gap)

function update_contacts!(loci, loci_previous, v, Λ)
    rel_vel_split = split_by_lengths(v,3)
    forces_split = split_by_lengths(Λ,3)
    for (ac,acp,va,Λa) in zip(loci, loci_previous, rel_vel_split,forces_split)
        if acp.state.active # == ac.state.active
            ac.state.persistent = true
        end
        ac.state.active = true
        ac.state.relative_velocity = SVector{3}(va)
        Λa[1] /= ac.μ #denormalized
        ac.state.force             = SVector{3}(Λa)
    end
end

struct Contact{T}
    id::Int
    μ::T
    e::T
    state::ContactState{3,T}
end

Base.length(::Contact) = 1
Base.iterate(x::Contact) = (x, nothing)
Base.iterate(x::Contact, ::Nothing) = nothing

function Contact(id,μ,e)
    active = false
    persistent = true
    i = one(μ)
    o = zero(μ)
    gap = i
    n = SVector(i,o,o)
    frame = Axes(n)
    v = SVector(o,o,o)
    Λ = SVector(o,o,o)
    state = ContactState(active,persistent,gap,frame,v,Λ)
    Contact(id,μ,e,state)
end

mutable struct ApproxFrictionalImpulse{T}
    active::Bool
    Λn::T
    β::Vector{T}
end

ApproxFrictionalImpulse{T}(m) where T = ApproxFrictionalImpulse(false,zero(T),zeros(T,m))

struct ApproxFrictionalContact{T,L}
    ϵ::T
    μ::T
    frame::Axes{3,T,L}
    m::Int
    e::Vector{T}
    d::Vector{SArray{Tuple{3},T,1,3}}
    D::Array{T,2}
    impulse::ApproxFrictionalImpulse{T}
end

function ApproxFrictionalContact(ϵ::T,μ::T,m::Int) where T
    frame = spatial_axes{T}()
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

