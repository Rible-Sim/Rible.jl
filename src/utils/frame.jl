
mutable struct CartesianFrame{N,M,T,L}
    position::MVector{N,T}
    axes::Axes{N,T,L}
    velocity::MVector{N,T}
    angular_velocity::MVector{M,T}
    local_angular_velocity::MVector{M,T}
end

function CartesianFrame(
        position::MVector{N,T},R,
        velocity::MVector{N,T}=zero(position),
        ω=zeros(MVector{2N-3,T})
    ) where {N,T}
    axes = Axes(R)
    M = 2N-3
    angular_velocity = MVector{M}(ω)
    if N == 2
        local_angular_velocity = MVector{M}(ω)
    else
        local_angular_velocity = MVector{M}(transpose(R.X)*MVector{M}(ω)) # cause 3 allocations
    end
    CartesianFrame(
        MVector{N}(position),
        axes,
        MVector{N}(velocity),
        angular_velocity,
        local_angular_velocity
    )
end

function CartesianFrame(r::MVector{N},v::MVector{N},R,ω) where {N}
    CartesianFrame(r,R,v,ω)
end

to_3D(frame::CartesianFrame{3}) = frame

function CartesianFrame{3}(frame::CartesianFrame{2})
    to_3D(frame)
end

function to_3D(frame::CartesianFrame{2,M,T}) where {M,T}
    (;position,velocity,axes,angular_velocity) = frame
    o = zero(T)
    i = one(T)
    position_3D = MVector{3}(position[1],position[2],o)
    velocity_3D = MVector{3}(velocity[1],velocity[2],o)
    axes_3D = Axes(
        SMatrix{3,3}(
            [
                 axes.X[1,1] -axes.X[2,1] o;
                -axes.X[1,2]  axes.X[2,2] o;
                      o            o      i;
            ]
        )
    )
    angular_velocity_3D = MVector{3}(o,o,angular_velocity[1])
    CartesianFrame(
        position_3D, 
        axes_3D, 
        velocity_3D, 
        angular_velocity_3D
    )
end
