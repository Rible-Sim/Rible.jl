function get_uv(nmcs::NC2D4C{T}, q) where {T}
    C = nmcs.conversion_to_X
    X = C * SVector{4,T}(ntuple(i -> q[i], Val(4)))
    u = SVector{2,T}(X[3], X[4])
    v = SVector{2,T}(-u[2],u[1])
    u, v
end

function get_uv(nmcs::NC2D6C{T}, q) where {T}
    C = nmcs.conversion_to_X
    X = C * SVector{6,T}(ntuple(i -> q[i], Val(6)))
    u = SVector{2,T}(X[3], X[4])
    v = SVector{2,T}(X[5], X[6])
    u, v
end

function get_uv(nmcs::NC2D, q)
    X = nmcs.conversion_to_X * q
    if nmcs isa NC2D6C
        u = @view X[3:4]
        v = @view X[5:6]
    elseif nmcs isa NC2D4C
        u = @view X[3:4]
        v = rotation_matrix(π/2) * u
    end
    SVector{2}(u), SVector{2}(v)
end

function get_uvw(nmcs::NC3D6C{T}, q) where {T}
    C = nmcs.conversion_to_X
    X = C * SVector{6,T}(ntuple(i -> q[i], Val(6)))
    u = SVector{3,T}(X[4], X[5], X[6])
    u = u / norm(u)
    v, w = HouseholderOrthogonalization(u)
    u, v, w
end

function get_uvw(nmcs::NC3D12C{T}, q) where {T}
    C = nmcs.conversion_to_X
    X = C * SVector{12,T}(ntuple(i -> q[i], Val(12)))
    u = SVector{3,T}(X[4], X[5], X[6])
    v = SVector{3,T}(X[7], X[8], X[9])
    w = SVector{3,T}(X[10], X[11], X[12])
    u, v, w
end

function get_uvw(nmcs::NC3D, q)
    X = nmcs.conversion_to_X * q
    if nmcs isa NC3D12C
        u = @view X[4:6]
        v = @view X[7:9]
        w = @view X[10:12]
    elseif nmcs isa NC3D6C
        u = @view X[4:6]
        u /= norm(u)
        v, w = HouseholderOrthogonalization(u)
    end
    SVector{3}(u), SVector{3}(v), SVector{3}(w)
end

function get_X(nmcs::NC,q::AbstractVector)
    X = nmcs.conversion_to_X*q
    num_of_dim = get_num_of_dims(nmcs)
    nld = get_num_of_local_dims(nmcs)
    SMatrix{num_of_dim,nld}(X[num_of_dim+1:end])
end

get_X(q::AbstractVector,nmcs::NC) = get_X(nmcs,q)

function find_rotation(nmcs::NC,q::AbstractVector)
    num_of_dim = get_num_of_dims(nmcs)
    if (nmcs isa NC2D2C) || (nmcs isa NC3D3C)
        R = SMatrix{num_of_dim,num_of_dim,eltype(q)}(I(num_of_dim))
    elseif nmcs isa NC2D4C
        (;r̄i,X̄) = nmcs.data
        ū,v̄ = get_uv(nmcs,vcat(r̄i,vec(X̄)))
        u,v = get_uv(nmcs,q)
        R = SMatrix{num_of_dim,num_of_dim}([u;;v]*inv([ū;;v̄]))
    elseif nmcs isa NC3D6C
        (;r̄i,X̄) = nmcs.data
        ū,v̄,w̄ = get_uvw(nmcs,vcat(r̄i,vec(X̄)))
        u,v,w = get_uvw(nmcs,q)
        R = SMatrix{num_of_dim,num_of_dim}([u;;v;;w]*inv([ū;;v̄;;w̄]))
    else
        X = get_X(nmcs,q)
        (;invX̄) = nmcs.data
        R = SMatrix{num_of_dim,num_of_dim}(X*invX̄)
    end
    return R
end

find_rotation(q::AbstractVector, nmcs::NC) = find_rotation(nmcs,q)


function find_angular_velocity(nmcs::NC{N,M,T,L,NCOORDS,NCOORDS2},q::AbstractVector,q̇::AbstractVector) where {N,M,T,L,NCOORDS,NCOORDS2}
    Ẋ = get_X(nmcs,q̇)
    X = get_X(nmcs,q)
    num_of_dim = get_num_of_dims(nmcs)
    o = zero(T)
    if num_of_dim == 2
        if nmcs isa NC2D2C
            ω = SVector{1}(o)
        else
            u = X[:,1]
            u̇ = Ẋ[:,1]
            ω = SVector{1}([-u[2],u[1]]\u̇)
        end
    else
        if nmcs isa NC3D3C
            ω = SVector{3}(o,o,o)
        else
            nld = get_num_of_local_dims(nmcs)
            if nld == 1
                ω = pinv(-skew(X[:,1]))*Ẋ[:,1]
                Ω = skew(ω)
            else
                Ω = Ẋ*pinv(X)
                ω = SVector{3}(Ω[3,2],Ω[1,3],Ω[2,1])
            end
        end
    end
    return ω
end

find_angular_velocity(q::AbstractVector,q̇::AbstractVector,nmcs::NC3D) = find_angular_velocity(nmcs,q,q̇)

# Transformations
"""
Return transformation matrix。
$(TYPEDSIGNATURES)
"""
function to_local_coords(nmcs::NC,r̄)
    (;r̄i,invX̄) = nmcs.data
    invX̄*(r̄-r̄i)
end

function to_position(nmcs::NC,q,c)
    to_position_jacobian(nmcs,c)*q
end

function to_position_jacobian(nmcs::NC,c)
    num_of_dim = get_num_of_dims(nmcs)
    nlds = get_num_of_local_dims(nmcs)
    ncoords = get_num_of_coords(nmcs)
    conversion = nmcs.conversion_to_std
    T = eltype(c)
    C_raw = SMatrix{1, nlds+1, T}(i == 1 ? one(T) : T(c[i-1]) for i in 1:(nlds+1))
    kron_block = SMatrix{num_of_dim, num_of_dim*(nlds+1), T}(
        begin
            j = ((col-1) ÷ num_of_dim) + 1
            r = ((col-1) % num_of_dim) + 1
            (r == row ? one(T) : zero(T)) * C_raw[1, j]
        end for row in 1:num_of_dim, col in 1:num_of_dim*(nlds+1)
    )
    return kron_block * conversion
end

function to_position_jacobian(nmcs::NC,q,c)
    to_position_jacobian(nmcs,c)
end

function to_velocity_jacobian(nmcs::NC,q,q̇,c)
    num_of_dim = get_num_of_dims(nmcs)
    nlds = get_num_of_local_dims(nmcs)
    ncoords = get_num_of_coords(nmcs)
    T = promote_type(eltype(nmcs.conversion_to_std), eltype(c))
    return zeros(SMatrix{num_of_dim, ncoords, T})
end

"""
Return rigid body natural coordinates 
$(TYPEDSIGNATURES)
"""
function cartesian_frame2coords(nmcs::Union{NC2D2C,NC3D3C},origin_frame)
    (;position,velocity) = origin_frame
    position, velocity
end

function cartesian_frame2coords(nmcs::NC,origin_frame)
    (;position,velocity,axes,angular_velocity) = origin_frame
    (;r̄i,X̄) = nmcs.data
    ri = position + axes*r̄i
    ṙi = velocity + angular_velocity×(ri-position)
    X = (axes.X*X̄)
    Ẋ = reduce(hcat,Ref(angular_velocity) .× eachcol(X))
    qstd = vcat(ri,vec(X))
    q̇std = vcat(ṙi,vec(Ẋ))
    Y = nmcs.conversion_to_std
    q = Y\qstd
    q̇ = Y\q̇std
    q,q̇
end
