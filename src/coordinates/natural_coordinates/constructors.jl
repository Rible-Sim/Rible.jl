"""
Natural Coordinates Constructors

This module provides constructors for various natural coordinate systems.
"""

# Helper functions for constructors
"""
！！！！！！！！！
$(TYPEDSIGNATURES)
"""
function get_conversion_core(np, nv)
    nb = np + nv
    B = Matrix(1I, nb, nb)
    for i = 2:np
        B[i, 1] = -1
    end
    B
end

function get_conversion_core(nmcs::NC)
    (;np, nv) = nmcs
    get_conversion_core(np, nv)
end

function get_conversion_to_X(ndim,np,nv)
    nb = np+nv
    ncoords = ndim * nb
    B = Matrix(1I,nb,nb)
    B[1] = 0
    for i = 2:np
        B[i,1] = -1
    end
    C = kron(
        sparse(B),
        IMatrix(ndim)
    )
    SMatrix{ncoords,ncoords,Int64,ncoords*ncoords}(C)
end

"""
！！！！！！！！！
$(TYPEDSIGNATURES)
"""
function get_conversion_to_std(ndim,np,nv)
    nb = np + nv
    ncoords = ndim * nb
    B = get_conversion_core(np, nv)
    C = kron(
        sparse(B),
        IMatrix(ndim)
    )
    SMatrix{ncoords,ncoords,Int64,ncoords*ncoords}(C)
end


function get_hessians_idx(::LNCData{N,0,T,L}) where {N,T,L}
    [
        CartesianIndex{2}[]
    ]
end

function get_hessians_idx(::LNCData{N,1,T,L}) where {N,T,L}
    [
        [CartesianIndex(2,2),CartesianIndex(2,2)],
    ]
end

function get_hessians_idx(::LNCData{2,2,T,4}) where {T}
    [
        [CartesianIndex(2,2),CartesianIndex(2,2)],
        [CartesianIndex(3,3),CartesianIndex(3,3)],
        [CartesianIndex(2,3),CartesianIndex(3,2)]
    ]
end

function get_hessians_idx(::LNCData{3,3,T,9}) where {T}
    [
        [CartesianIndex(2,2),CartesianIndex(2,2)],
        [CartesianIndex(3,3),CartesianIndex(3,3)],
        [CartesianIndex(4,4),CartesianIndex(4,4)],
        [CartesianIndex(3,4),CartesianIndex(4,3)],
        [CartesianIndex(2,4),CartesianIndex(4,2)],
        [CartesianIndex(2,3),CartesianIndex(3,2)]
    ]
end

#todo cache cstr_hessians
function get_cstr_hessians(data::LNCData,cv)
    ndim = get_num_of_dims(data)
    nld = get_num_of_local_dims(data)
    I_Bool = IMatrix(ndim)
    idx = get_hessians_idx(data)
    [
        begin
            ret_raw = zeros(Int,nld+1,nld+1)
            for ij in id
                ret_raw[ij] += 1
            end
            Symmetric(sparse(transpose(cv)*kron(ret_raw,I_Bool)*cv))
        end
        for id in idx
    ]
end

function _normalize_cstr_idx(num_of_intrinsic_cstr::Int, cstr_idx)
    idx = if cstr_idx isa Colon || cstr_idx === nothing
        collect(1:num_of_intrinsic_cstr)
    elseif cstr_idx isa AbstractVector{<:Integer}
        collect(Int, cstr_idx)
    else
        throw(ArgumentError("cstr_idx must be `:`, `nothing`, or an integer vector, got $(typeof(cstr_idx))."))
    end
    if any(i -> i < 1 || i > num_of_intrinsic_cstr, idx)
        throw(ArgumentError("Invalid cstr_idx=$idx. Valid intrinsic constraint indices are in 1:$num_of_intrinsic_cstr."))
    end
    if !allunique(idx)
        throw(ArgumentError("Invalid cstr_idx=$idx. Constraint indices must be unique."))
    end
    idx
end

function NC(np,nv,ndim,nld,r̄i::AbstractVector{T},X̄; cstr_idx=:) where T
    data = LNCData(
        SVector{ndim}(r̄i),
        SMatrix{ndim,nld,T}(X̄)
    )
    cv = get_conversion_to_std(ndim,np,nv)
    nb = np + nv
    NCOORDS = ndim * nb
    NCOORDS2 = NCOORDS * NCOORDS
    L = ndim*nld

    ndof = if nld == 0
        ndim
    elseif ndim == 2
        3
    elseif (ndim == 3) && (nld == 1)
        # 3D rigid bar: translation (3) + direction on S^2 (2)
        5
    else
        6
    end

    num_of_intrinsic_cstr = NCOORDS - ndof
    idx = _normalize_cstr_idx(num_of_intrinsic_cstr, cstr_idx)
    num_of_cstr = length(idx)

    NC{ndim,nld,T,L,NCOORDS,NCOORDS2}(
        np,
        nv,
        data,
        cv,
        get_conversion_to_X(ndim,np,nv),
        get_cstr_hessians(data,cv),
        num_of_cstr,
        idx
    )
end

function NC(nmcs::NC; cstr_idx=:)
    (; np, nv, data, conversion_to_std, conversion_to_X, hessians) = nmcs
    num_of_intrinsic_cstr = get_num_of_intrinsic_cstr(nmcs)
    idx = _normalize_cstr_idx(num_of_intrinsic_cstr, cstr_idx)
    num_of_cstr = length(idx)
    NC(
        np,
        nv,
        data,
        conversion_to_std,
        conversion_to_X,
        hessians,
        num_of_cstr,
        idx
    )
end

"""
Return 2D point mass natural coordinates.
$(TYPEDSIGNATURES)
"""
function NC2D1P(ri::AbstractVector{T}; cstr_idx=:) where T
    r̄i = @SVector zeros(T,2)
    np = 1
    nv = 0
    ndim = 2
    nld = 0
    X̄ = @SMatrix zeros(T,2,0)
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 2D rigid bar natural coordinates, using 1 basic point and 1 base vector.
$(TYPEDSIGNATURES)
"""
function NC2D1P1V(ri::AbstractVector{T},u::AbstractVector{T},
    origin_position=SVector{2}(zeros(T,2)),α=SMatrix{2,2}(one(T)*I);
    cstr_idx=:) where T
    R = rotation_matrix(α)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    ū = invR*u
    X̄ = ū
    np = 1
    nv = 1
    ndim = 2
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 2D rigid bar natural coordinates, starting from two points.
$(TYPEDSIGNATURES)
"""
function NC2D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{2}(zeros(T,2)),α=SMatrix{2,2}(one(T)*I);
    cstr_idx=:) where T
    R = rotation_matrix(α)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    r̄j = invR*(rj-origin_position) #
    ū = r̄j-r̄i
    X̄ = ū
    np = 2
    nv = 0
    ndim = 2
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 2D rigid bodies natural coordinates, using 1 basic point and 2 base vectors.
$(TYPEDSIGNATURES)
"""
function NC1P2V(ri::AbstractVector{T},
    origin_position=SVector{2}(zeros(T,2)),α=SMatrix{2,2}(one(T)*I)
    ; cstr_idx=:) where T
    R = rotation_matrix(α)
    o = zero(T)
    i = one(T)
    ū = SVector(i,o)
    v̄ = SVector(o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    u = R*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    np = 1
    nv = 2
    ndim = 2
    nld = 2
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 2D rigid bodies natural coordinates, using 2 basic points and 1 base vector.
$(TYPEDSIGNATURES)
"""
function NC2P1V(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{2}(zeros(T,2)),α=SMatrix{2,2}(one(T)*I);
    cstr_idx=:) where T
    R = rotation_matrix(α)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    ū = r̄j-r̄i
    v̄ = rotation_matrix(π/2)*ū
    v = R*v̄
    X̄ = hcat(ū,v̄)
    np = 2
    nv = 1
    ndim = 2
    nld = 2
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 2D rigid bodies natural coordinates, using 3 basic points
$(TYPEDSIGNATURES)
"""
function NC3P(ri::AbstractVector{T},rj::AbstractVector{T},rk::AbstractVector{T},
              origin_position=SVector{2}(zeros(T,2)),α=SMatrix{2,2}(one(T)*I);
              cstr_idx=:) where T
    R = rotation_matrix(α)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    r̄k = invR*(rk-origin_position)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    X̄ = hcat(ū,v̄)
    np = 3
    nv = 0
    ndim = 2
    nld = 2
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D point mass natural coordinates.
$(TYPEDSIGNATURES)
"""
function NC3D1P(ri::AbstractVector{T}; cstr_idx=:) where T
    r̄i = @SVector zeros(T,3)
    np = 1
    nv = 0
    ndim = 3
    nld = 0
    X̄ = @SMatrix zeros(T,3,0)
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bar natural coordinates.
$(TYPEDSIGNATURES)
"""
function NC3D1P1V(ri::AbstractVector{T},u::AbstractVector{T},
               origin_position=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
               ; cstr_idx=:) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) 
    ū = invR*u
    X̄ = ū
    np = 1
    nv = 1
    ndim = 3
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bar natural coordinates, starting from two points.
$(TYPEDSIGNATURES)
"""
function NC3D2P(ri::AbstractVector{T},rj::AbstractVector{T},
    origin_position=SVector{3}(zeros(T,3)),R=SMatrix{3,3}(one(T)*I)
    ; cstr_idx=:) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    r̄j = invR*(rj-origin_position)
    ū = r̄j-r̄i
    X̄ = ū
    np = 2
    nv = 0
    ndim = 3
    nld = 1
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

## 3D Rigid body
"""
Return 3D rigid bodies natural coordinates, using 1 basic point and 3 base vectors
$(TYPEDSIGNATURES)
"""
function NC1P3V(ri::AbstractVector{T},
                origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                ; cstr_idx=:) where T
    o = zero(T)
    i = one(T)
    ū = SVector(i,o,o)
    v̄ = SVector(o,i,o)
    w̄ = SVector(o,o,i)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    u = R*ū
    v = R*v̄
    w = R*w̄
    X̄ = hcat(ū,v̄,w̄)
    np = 1
    nv = 3
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bodies natural coordinates, using 1 basic point and 3 base vectors.
$(TYPEDSIGNATURES)
"""
function NC1P3V(ri::AbstractVector{T},u::AbstractVector{T},v::AbstractVector{T},w::AbstractVector{T},
                origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                ; cstr_idx=:) where T
    o = zero(T)
    i = one(T)
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position) #
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    np = 1
    nv = 3
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bodies natural coordinates, using 2 basic points and 2 base vectors.
$(TYPEDSIGNATURES)
"""
function NC2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ; cstr_idx=:) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    u = rj-ri
    ū = invR*u
    v̄,w̄ = HouseholderOrthogonalization(ū./norm(ū)).*norm(ū)
    v = R*v̄
    w = R*w̄
    X̄ = hcat(ū,v̄,w̄)
    np = 2
    nv = 2
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bodies natural coordinates, using 2 basic points and 2 base vectors.
$(TYPEDSIGNATURES)
"""
function NC2P2V(ri::AbstractVector{T},rj::AbstractVector{T},
                v::AbstractVector{T},w::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ; cstr_idx=:) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    u = rj-ri
    ū = invR*u
    v̄ = invR*v
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    np = 2
    nv = 2
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bodies natural coordinates, using 3 basic points and 1 base vector.
$(TYPEDSIGNATURES)
"""
function NC3P1V(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ; cstr_idx=:) where T
    u = rj-ri
    v = rk-ri
    w = u×v
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    r̄k = invR*(rk-origin_position)
    ū = r̄j-r̄i
    v̄ = r̄k-r̄i
    w̄ = invR*w
    X̄ = hcat(ū,v̄,w̄)
    np = 3
    nv = 1
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end

"""
Return 3D rigid bodies natural coordinates, using 4 basic points
$(TYPEDSIGNATURES)
"""
function NC4P(ri::AbstractVector{T},rj::AbstractVector{T},
                       rk::AbstractVector{T},rl::AbstractVector{T},
                       origin_position=zeros(T,3),R=Matrix(one(T)*I,3,3)
                       ; cstr_idx=:) where T
    invR = transpose(R) # because this is a rotation matrix
    r̄i = invR*(ri-origin_position)
    r̄j = invR*(rj-origin_position)
    r̄k = invR*(rk-origin_position)
    r̄l = invR*(rl-origin_position)
    ū = r̄j - r̄i
    v̄ = r̄k - r̄i
    w̄ = r̄l - r̄i
    X̄ = hcat(ū,v̄,w̄)
    np = 4
    nv = 0
    ndim = 3
    nld = 3
    NC(np,nv,ndim,nld,r̄i,X̄; cstr_idx)
end
