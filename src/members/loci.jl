struct Axes{N,T}
    X::SMatrix{N,N,T}
end

function Base.getproperty(a::Axes{2}, p) 
    if (p == :n) || (p == :normal)
        return a.X[:,1]
    elseif (p == :t) || (p == :tangent)
        return  a.X[:,2]
    else # fallback to getfield
        return getfield(a, p)
    end
end

function Base.getproperty(a::Axes{3}, p::Symbol) 
    if (p == :n) || (p == :normal)
        return a.X[:,1]
    elseif (p == :t) || (p == :tangent)
        return  a.X[:,2]
    elseif (p == :b) || (p == :bitangent)
        return  a.X[:,3]
    else # fallback to getfield
        return getfield(a, p)
    end
end

function planar_frame(n)
    n /= norm(n)
    Axes(
        SMatrix{2,2,eltype{n}}(
            n[1], n[2], -n[2], n[1]
        )
    )
end

function spatial_frame(n)
    n |> SVector{3} |> get_orthonormal_axes |> Axes
end

function orthonormal_frame(n::SVector{N}) where {N}
    if N == 2
        axes = planar_frame(n)
    elseif N == 3
        axes = spatial_frame(n)
    end
    axes
end

function (Base.:*)(R::StaticArray,axes::Axes{N}) where N
    Axes(SMatrix{N,N}(R*axes.X))
end

struct Locus{N,T}
    position::SVector{N,T}
    axes::Axes{N,T}
    friction_coefficient::T
    restitution_coefficient::T
end

function Locus(
        position::SVector{N,T},axes::Axes,
        friction_coefficient=zero(T),
        restitution_coefficient=zero(T)
    ) where {N,T}
    Locus(
        position,
        axes,
        friction_coefficient,
        restitution_coefficient
    )
end

function Locus(
        position::SVector{N,T},
        friction_coefficient = zero(T),
        restitution_coefficient = zero(T)
    ) where {N,T}
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = orthonormal_frame(normal)
    Locus(
        position,axes,
        friction_coefficient,
        restitution_coefficient
    )
end


function Locus(
        position::SVector{N,T},
        normal::SVector{N,T},
        friction_coefficient = zero(T),
        restitution_coefficient = zero(T)
    ) where {N,T}
    axes = orthonormal_frame(normal)
    Locus(
        position,axes,
        friction_coefficient,
        restitution_coefficient
    )
end

mutable struct FrictionalContactState{N,T}
    active::Bool
    persistent::Bool
    gap::T
    frame::Axes{N,T}
    relative_velocity::SVector{N,T}
    force::SVector{N,T}
end

function FrictionalContactState(normal::StaticArray{Tuple{N},T}) where {N,T}
    active = false
    persistent = true
    gap = typemax(T)
    frame = orthonormal_frame(normal)
    relative_velocity = @SVector zeros(T,N)
    force = @SVector zeros(T,N)
    FrictionalContactState(
        active,persistent,
        gap,frame,
        relative_velocity,force
    )
end

function activate!(contact_state::FrictionalContactState,gap)
    contact_state.gap = gap
    pre_active = contact_state.active
    contact_state.active = (gap<=0)
    contact_state.persistent = (pre_active == contact_state.active)
end

mutable struct LocusState{N,M,T}
    position::MVector{N,T}
    velocity::MVector{N,T}
    axes::Axes{N,T}
    force::MVector{N,T}
    torque::MVector{M,T}
    contact_state::FrictionalContactState{N,T}
end

function LocusState(
        position::StaticArray{Tuple{N},T},
        velocity::StaticArray{Tuple{N},T}
    ) where {N,T}
    M = 2N-3
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = orthonormal_frame(normal)
    force = @MVector zeros(T,N)
    torque = @MVector zeros(T,M)
    contact_state = FrictionalContactState(normal)
    LocusState(
        position,velocity,
        axes,
        force,torque,
        contact_state
    )
end

function LocusState(
        locus::Locus{N,T},
        origin_position,R,
        origin_velocity,ω,
        force,torque
    ) where {N,T}
    relative_position = R*locus.position
    position = origin_position+relative_position
    velocity = origin_velocity+ω×relative_position
    axes = R*locus.axes
    contact_state = FrictionalContactState(axes.normal)
    LocusState(
        position,velocity,
        axes,
        force,torque,
        contact_state
    )
end

function LocusState(
        locus::Locus{N,T},
        origin_position,R,
        origin_velocity,ω::StaticArray{Tuple{M}}
    ) where {N,M,T}
    force = @MVector zeros(T,N)
    torque = @MVector zeros(T,M)
    LocusState(
        locus,
        origin_position,R,
        origin_velocity,ω,
        force,torque
    )
end

activate!(locus_state::LocusState,gap) = activate!(locus_state.contact_state,gap)

function update_contacts!(loci, loci_previous, v, Λ)
    rel_vel_split = split_by_lengths(v,3)
    forces_split = split_by_lengths(Λ,3)
    for (ac,acp,va,Λa) in zip(loci, loci_previous, rel_vel_split,forces_split)
        if acp.state.active
            ac.state.persistent = true
        end
        ac.state.active = true
        ac.state.relative_velocity = SVector{3}(va)
        ac.state.force             = SVector{3}(Λa)
    end
end