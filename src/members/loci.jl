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

mutable struct LocusState{N,M,T}
    position::MVector{N,T}
    velocity::MVector{N,T}
    axes::Axes{N,T}
    frame::Axes{N,T}
    force::MVector{N,T}
    torque::MVector{M,T}
end

function LocusState(
        position::StaticArray{Tuple{N},T},
        velocity::StaticArray{Tuple{N},T}
    ) where {N,T}
    M = 2N-3
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = orthonormal_frame(normal)
    frame = deepcopy(axes)
    force = @MVector zeros(T,N)
    torque = @MVector zeros(T,M)
    LocusState(
        position,velocity,
        axes,frame,
        force,torque
    )
end

function LocusState(
        lo::Locus{N,T},
        origin_position,R,
        origin_velocity,ω,
        force,torque
    ) where {N,T}
    relative_position = R*lo.position
    position = origin_position+relative_position
    velocity = origin_velocity+ω×relative_position
    axes = R*lo.axes
    frame = deepcopy(axes)
    LocusState(
        position,velocity,
        axes,frame,
        force,torque
    )
end

function LocusState(
        lo::Locus{N,T},
        origin_position,R,
        origin_velocity,ω::StaticArray{Tuple{M}}
    ) where {N,M,T}
    force = @MVector zeros(T,N)
    torque = @MVector zeros(T,M)
    LocusState(
        lo,
        origin_position,R,
        origin_velocity,ω,
        force,torque
    )
end
