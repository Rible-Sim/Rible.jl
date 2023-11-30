struct Axes{N,T}
    X::SMatrix{N,N,T}
end


function Base.getproperty(a::Axes{2}, p::Symbol) 
    if (p === :n) || (p === :normal)
        return a.X[:,2]
    elseif (p === :t) || (p === :tangent)
        return  a.X[:,1]
    else # fallback to getfield
        return getfield(a, p)
    end
end

function Base.getproperty(a::Axes{3}, p::Symbol) 
    if (p === :n) || (p === :normal)
        return a.X[:,3]
    elseif (p === :t) || (p === :tangent)
        return  a.X[:,1]
    elseif (p === :b) || (p === :bitangent)
        return  a.X[:,2]
    else # fallback to getfield
        return getfield(a, p)
    end
end

function planar_axes(n)
    n /= norm(n)
    Axes(
        SMatrix{2,2}(
            -n[2], n[1], n[1], n[2]
        )
    )
end

function spatial_axes(n)
    normal = SVector{3}(n) 
    normal /= norm(normal)
    tangent, bitangent = HouseholderOrthogonalization(normal)
    SMatrix{3,3}(
        tangent[1], tangent[2], tangent[3],
        bitangent[1], bitangent[2], bitangent[3],
        normal[1], normal[2], normal[3],
    ) |> Axes
end

orthonormal_axes(n::SVector{2}) = planar_axes(n)
orthonormal_axes(n::SVector{3}) = spatial_axes(n)

Axes(X::Axes) = X
Axes(n::SVector) = orthonormal_axes(n)

function (Base.:*)(R::StaticArray{Tuple{N,N}},axes::Axes{N}) where N
    Axes(SMatrix{N,N}(R*axes.X))
end

function (Base.:*)(axes1::Axes{N},axes2::Axes{N}) where N
    Axes(axes1.X*axes2.X)
end

# work around
function (Base.:*)(n::StaticArray{Tuple{1,N}},axes::Axes{N}) where N
    orthonormal_axes(SVector{N}(n)).X*axes
end

function (Base.:*)(axes::Axes{N},n::SVector{N}) where N
    axes.X*n
end

mutable struct CartesianFrame{N,M,T}
    position::SVector{N,T}
    axes::Axes{N,T}
    velocity::SVector{N,T}
    angular_velocity::SVector{M,T}
    local_angular_velocity::SVector{M,T}
end

function CartesianFrame(
        position::SVector{N,T},R,
        velocity::SVector{N,T}=zero(position),
        ω=zeros(SVector{2N-3,T})
    ) where {N,T}
    axes = Axes(R)
    M = 2N-3
    angular_velocity = SVector{M}(ω)
    if N == 2
        local_angular_velocity = SVector{M}(ω)
    else
        local_angular_velocity = transpose(R.X)*SVector{M}(ω) # cause 3 allocations
    end
    CartesianFrame(
        position,
        axes,
        velocity,
        angular_velocity,
        local_angular_velocity
    )
end

function CartesianFrame(r::SVector{N},v::SVector{N},R,ω) where {N}
    CartesianFrame(r,R,v,ω)
end

to_3D(frame::CartesianFrame{3}) = frame

function to_3D(frame::CartesianFrame{2,M,T}) where {M,T}
    (;position,velocity,axes,angular_velocity) = frame
    o = zero(T)
    i = one(T)
    position_3D = SVector{3}(position[1],position[2],o)
    velocity_3D = SVector{3}(velocity[1],velocity[2],o)
    axes_3D = Axes(
        SMatrix{3,3}(
            [
                 axes.X[1,1] -axes.X[2,1] o;
                -axes.X[1,2]  axes.X[2,2] o;
                        o      o     i;
            ]
        )
    )
    angular_velocity_3D = SVector{3}(o,o,angular_velocity[1])
    CartesianFrame(
        position_3D, 
        axes_3D, 
        velocity_3D, 
        angular_velocity_3D
    )
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
    axes = orthonormal_axes(normal)
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
    axes = orthonormal_axes(normal)
    Locus(
        position,axes,
        friction_coefficient,
        restitution_coefficient
    )
end

mutable struct ContactState{N,T}
    active::Bool
    persistent::Bool
    gap::T
    axes::Axes{N,T}
    relative_velocity::SVector{N,T}
    force::SVector{N,T}
end

function ContactState(normal::StaticArray{Tuple{N},T}) where {N,T}
    active = false
    persistent = true
    gap = typemax(T)
    axes = orthonormal_axes(normal)
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

mutable struct LocusState{N,M,T}
    frame::CartesianFrame{N,M,T}
    force::MVector{N,T}
    torque::MVector{M,T}
    contact_state::ContactState{N,T}
end

function LocusState(
        position::StaticArray{Tuple{N},T},
        velocity::StaticArray{Tuple{N},T},
        force = zeros(MVector{N,T}),
        torque = zeros(MVector{2N-3,T})
    ) where {N,T}
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = orthonormal_axes(normal)
    frame = CartesianFrame(position,axes,velocity)
    contact_state = ContactState(normal)
    LocusState(
        frame,
        force,torque,
        contact_state
    )
end

function LocusState(
        locus::Locus{N,T},
        origin_frame::CartesianFrame{N,M,T},
        force = zeros(MVector{N,T}),
        torque = zeros(MVector{M,T})
    ) where {N,M,T}
    relative_position = origin_frame.axes*locus.position
    position = origin_frame.position+relative_position
    velocity = origin_frame.velocity+origin_frame.angular_velocity×relative_position
    axes = origin_frame.axes*locus.axes
    frame = CartesianFrame(position,axes,velocity,origin_frame.angular_velocity)
    contact_state = ContactState(axes.normal)
    LocusState(
        frame,
        force,torque,
        contact_state
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
