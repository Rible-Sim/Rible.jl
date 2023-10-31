
"""
Return 未约束的natural coodinates 编号。
$(TYPEDSIGNATURES)
"""
function get_free_idx(nmcs::LNC,pres_idx)
    deleteat!(collect(1:get_num_of_coords(nmcs)),pres_idx)
end


make_X(q::AbstractVector,nmcs::LNC) = make_X(nmcs,q)

function make_X(nmcs::LNC,q::AbstractVector)
    qstd = nmcs.conversion*q
    ndim = get_num_of_dims(nmcs)
    X = reshape(qstd[ndim+1:end],ndim,:)
    if (nmcs isa LNC2D2C) || (nmcs isa LNC3D3C)
        return X
    elseif  (nmcs isa LNC2D4C) || (nmcs isa LNC3D6C)
        # return find_rotation(nmcs,q)
        return SMatrix{ndim,1}(X)
    else
        return SMatrix{ndim,ndim}(X)
    end
end

find_rotation(q::AbstractVector, nmcs::LNC) = find_rotation(nmcs,q)
function find_rotation(nmcs::LNC,q::AbstractVector)
    ndim = get_num_of_dims(nmcs)
    if nmcs isa LNC2D4C
        (;r̄i,X̄) = nmcs
        ū,v̄ = get_uv(nmcs,vcat(r̄i,vec(X̄)))
        u,v = get_uv(nmcs,q)
        R = SMatrix{ndim,ndim}([u;;v]*inv([ū;;v̄]))
    elseif nmcs isa LNC3D6C
        (;r̄i,X̄) = nmcs
        ū,v̄,w̄ = get_uvw(nmcs,vcat(r̄i,vec(X̄)))
        u,v,w = get_uvw(nmcs,q)
        R = SMatrix{ndim,ndim}([u;;v;;w]*inv([ū;;v̄;;w̄]))
    else
        X = make_X(nmcs,q)
        (;invX̄) = nmcs
        R = SMatrix{ndim,ndim}(X*invX̄)
    end
    return R
end

find_angular_velocity(q::AbstractVector,q̇::AbstractVector,nmcs::LNC3D) = find_angular_velocity(nmcs,q,q̇)

function find_angular_velocity(nmcs::LNC,q::AbstractVector,q̇::AbstractVector)
    Ẋ = make_X(nmcs,q̇)
    X = make_X(nmcs,q)
    ndim = get_num_of_dims(nmcs)
    if ndim == 2
        u = X[:,1]
        u̇ = Ẋ[:,1]
        ω = SVector{1}([-u[2],u[1]]\u̇)
    else
        Ω = Ẋ*pinv(X)
        ω = SVector{3}(Ω[3,2],Ω[1,3],Ω[2,1])
    end
end


# Transformations
"""
Return transformation matrix。
$(TYPEDSIGNATURES)
"""
function to_local_coords(nmcs::LNC,r̄)
    (;r̄i,invX̄) = nmcs
    invX̄*(r̄-r̄i)
end

function to_transformation(nmcs::LNC,c)
    ndim = get_num_of_dims(nmcs)
    nlds = get_num_of_local_dims(nmcs)
    ncoords = get_num_of_coords(nmcs)
    conversion = nmcs.conversion
    C_raw = hcat(1,transpose(SVector{nlds}(c)))
    SMatrix{ndim,ncoords}(kron(C_raw,IMatrix(ndim))*conversion)
end


"""
Return rigid body natural coodinates 
$(TYPEDSIGNATURES)
"""
function cartesian_frame2coords(nmcs::Union{LNC2D2C,LNC3D3C},origin_position,R)
    origin_position
end

function cartesian_frame2coords(nmcs::Union{LNC2D2C,LNC3D3C},origin_position,R,origin_velocity,ω)
    origin_position,origin_velocity
end

function cartesian_frame2coords(nmcs,origin_position,R)
    (;r̄i,X̄) = nmcs
    ri = origin_position + R*r̄i
    X = R*X̄
    qstd = vcat(ri,vec(X))
    Y = nmcs.conversion
    q = Y\qstd
    q
end

function cartesian_frame2coords(nmcs,origin_position,R,origin_velocity,ω)
    (;r̄i,X̄) = nmcs
    ri = origin_position + R*r̄i
    ṙi = origin_velocity + ω×(ri-origin_position)
    X = R*X̄
    Ẋ = reduce(hcat,Ref(ω) .× eachcol(X))
    qstd = vcat(ri,vec(X))
    q̇std = vcat(ṙi,vec(Ẋ))
    Y = nmcs.conversion
    q = Y\qstd
    q̇ = Y\q̇std
    q,q̇
end