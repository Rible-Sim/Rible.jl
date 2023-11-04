function get_uv(nmcs::NC2D,q)
    X = nmcs.conversion_to_X*q
    if     nmcs isa NC2D6C
        u = @view X[3:4]
        v = @view X[5:6]
    elseif nmcs isa NC2D4C
        u = @view X[3:4]
        v = rotation_matrix(π/2)*u
    end
    SVector{2}(u),SVector{2}(v)
end

function get_uvw(nmcs::NC3D,q)
    X = nmcs.conversion_to_X*q
    if     nmcs isa NC3D12C
        u = @view X[4:6]
        v = @view X[7:9]
        w = @view X[10:12]
    elseif nmcs isa NC3D6C
        u = @view X[4:6]
        u /= norm(u)
        v,w = HouseholderOrthogonalization(u)
    end
    SVector{3}(u),SVector{3}(v),SVector{3}(w)
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
    if nmcs isa NC2D4C
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


function find_angular_velocity(nmcs::NC,q::AbstractVector,q̇::AbstractVector)
    Ẋ = get_X(nmcs,q̇)
    X = get_X(nmcs,q)
    num_of_dim = get_num_of_dims(nmcs)
    if num_of_dim == 2
        u = X[:,1]
        u̇ = Ẋ[:,1]
        ω = SVector{1}([-u[2],u[1]]\u̇)
    else
        Ω = Ẋ*pinv(X)
        ω = SVector{3}(Ω[3,2],Ω[1,3],Ω[2,1])
    end
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

function to_transformation(nmcs::NC,c)
    num_of_dim = get_num_of_dims(nmcs)
    nlds = get_num_of_local_dims(nmcs)
    ncoords = get_num_of_coords(nmcs)
    conversion = nmcs.conversion_to_std
    C_raw = hcat(1,transpose(SVector{nlds}(c)))
    SMatrix{num_of_dim,ncoords}(kron(C_raw,IMatrix(num_of_dim))*conversion)
end

"""
Return rigid body natural coodinates 
$(TYPEDSIGNATURES)
"""
function cartesian_frame2coords(nmcs::Union{NC2D2C,NC3D3C},origin_position,R)
    origin_position
end

function cartesian_frame2coords(nmcs::Union{NC2D2C,NC3D3C},origin_position,R,origin_velocity,ω)
    origin_position,origin_velocity
end

function cartesian_frame2coords(nmcs,origin_position,R)
    (;r̄i,X̄) = nmcs.data
    ri = origin_position + R*r̄i
    X = R*X̄
    qstd = vcat(ri,vec(X))
    Y = nmcs.conversion_to_std
    q = Y\qstd
    q
end

function cartesian_frame2coords(nmcs,origin_position,R,origin_velocity,ω)
    (;r̄i,X̄) = nmcs.data
    ri = origin_position + R*r̄i
    ṙi = origin_velocity + ω×(ri-origin_position)
    X = R*X̄
    Ẋ = reduce(hcat,Ref(ω) .× eachcol(X))
    qstd = vcat(ri,vec(X))
    q̇std = vcat(ṙi,vec(Ẋ))
    Y = nmcs.conversion_to_std
    q = Y\qstd
    q̇ = Y\q̇std
    q,q̇
end