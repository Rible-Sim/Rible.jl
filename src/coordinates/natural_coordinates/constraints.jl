# """
# Return 2D or 3D local natural coords deformations for rigid bars.
# $(TYPEDSIGNATURES)
# """
@inline @inbounds get_deform(ū::AbstractVector) = SVector(sqrt(ū⋅ū))

@inline @inbounds function get_deform(nmcs::Union{NC2D2C,NC3D3C})
    (;r̄i) = nmcs.data
    get_deform(r̄i)
end

@inline @inbounds function get_deform(nmcs::Union{NC2D4C,NC3D6C})
    (;X̄) = nmcs.data
    ū = X̄[:,1]
    get_deform(ū)
end

# """
# Return 2D natural coordinates deformations , for rigid bodies.
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(nmcs::NC2D6C)
    (;X̄) = nmcs.data
    ū = X̄[:,1]
    v̄ = X̄[:,2]
    SVector(sqrt(ū⋅ū),sqrt(v̄⋅v̄),ū⋅v̄)
end

# """
# Return 3D natural coordinates deformations , for rigid bodies.
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(nmcs::NC3D12C)
    (;X̄) = nmcs.data
    ū = X̄[:,1]
    v̄ = X̄[:,2]
    w̄ = X̄[:,3]
    u_sqrt = sqrt(ū⋅ū)
    v_sqrt = sqrt(v̄⋅v̄)
    w_sqrt = sqrt(w̄⋅w̄)
    vw = v̄⋅w̄
    uw = ū⋅w̄
    uv = ū⋅v̄
    SVector(u_sqrt,v_sqrt,w_sqrt,vw,uw,uv)
end

# Intrinsic Constraints

"""
Return 2D or 3D intrinsic cstr(s) for point mass.
$(TYPEDSIGNATURES)
"""
function cstr_function!(ret, nmcs::Union{NC2D2C,NC3D3C}, q::AbstractVector, deforms = get_deform(nmcs))
    #Not constrained for a point mass
end

"""
Return 2D or 3D intrinsic cstr(s) for rigid bars.
$(TYPEDSIGNATURES)
"""
function cstr_function!(ret::AbstractVector{T}, nmcs::NC2D4C{T}, q::AbstractVector{T}, d::SVector{1,T} = get_deform(nmcs)) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{4,T}(ntuple(i -> q[i], Val(4)))
    u = SVector{2,T}(qstd[3], qstd[4])
    ret[1] = (u ⋅ u) - d[1]^2
    @view ret[nmcs.cstr_idx]
end

function cstr_function!(ret::AbstractVector{T}, nmcs::NC3D6C{T}, q::AbstractVector{T}, d::SVector{1,T} = get_deform(nmcs)) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{6,T}(ntuple(i -> q[i], Val(6)))
    u = SVector{3,T}(qstd[4], qstd[5], qstd[6])
    ret[1] = (u ⋅ u) - d[1]^2
    @view ret[nmcs.cstr_idx]
end

"""
Return 2D intrinsic cstr(s) , for rigid bodies.
$(TYPEDSIGNATURES)
"""
function cstr_function!(ret::AbstractVector{T}, nmcs::NC2D6C, q::AbstractVector{T}, d::SVector{3,T} = get_deform(nmcs)) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{6,T}(ntuple(i -> q[i], Val(6)))
    u = SVector{2,T}(qstd[3], qstd[4])
    v = SVector{2,T}(qstd[5], qstd[6])
    ret[1] = (u⋅u - d[1]^2)
    ret[2] = (v⋅v - d[2]^2)
    ret[3] = (u⋅v - d[3])
    @view ret[nmcs.cstr_idx]
end


"""
Return 3D intrinsic cstr(s) , for rigid bodies.
$(TYPEDSIGNATURES)
"""
function cstr_function!(ret::AbstractVector{T}, nmcs::NC3D12C, q::AbstractVector{T}, d::SVector{6,T} = get_deform(nmcs)) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{12,T}(ntuple(i -> q[i], Val(12)))
    u = SVector{3,T}(qstd[4], qstd[5], qstd[6])
    v = SVector{3,T}(qstd[7], qstd[8], qstd[9])
    w = SVector{3,T}(qstd[10], qstd[11], qstd[12])
    ret[1] = u⋅u - d[1]^2
    ret[2] = v⋅v - d[2]^2
    ret[3] = w⋅w - d[3]^2
    ret[4] = v⋅w - d[4]
    ret[5] = u⋅w - d[5]
    ret[6] = u⋅v - d[6]
    @view ret[nmcs.cstr_idx]
end

## Jacobians
"""
Return 2D or 3D Jacobian matrix for point mass.
$(TYPEDSIGNATURES)
"""
function cstr_jacobian!(ret, nmcs::Union{NC2D2C,NC3D3C}, jac, q::AbstractVector)
    # not constrained for a point mass
end

"""
Return 2D or 3D Jacobian matrix for rigid bars.
$(TYPEDSIGNATURES)
"""
function cstr_jacobian!(ret::AbstractMatrix{T}, nmcs::NC2D4C{T}, jac::AbstractMatrix{T}, q::AbstractVector{T}) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{4,T}(ntuple(i -> q[i], Val(4)))
    fill!(jac, zero(T))
    @inbounds begin
        jac[1,3] = 2qstd[3]
        jac[1,4] = 2qstd[4]
    end
    jac_out = SMatrix{1,4,T,4}(jac) * C
    @inbounds for j in 1:4
        ret[1, j] = jac_out[1, j]
    end
    ret
end

"""
Return 2D or 3D Jacobian matrix for rigid bars.
$(TYPEDSIGNATURES)
"""
function cstr_jacobian!(ret::AbstractMatrix{T}, nmcs::NC3D6C{T}, jac::AbstractMatrix{T}, q::AbstractVector{T}) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{6,T}(ntuple(i -> q[i], Val(6)))
    fill!(jac, zero(T))
    @inbounds begin
        jac[1,4] = 2qstd[4]
        jac[1,5] = 2qstd[5]
        jac[1,6] = 2qstd[6]
    end
    jac_out = SMatrix{1,6,T,6}(jac) * C
    @inbounds for j in 1:6
        ret[1, j] = jac_out[1, j]
    end
    ret
end

"""
Return 2D Jacobian matrix , for rigid bodies.
$(TYPEDSIGNATURES)
"""
function cstr_jacobian!(ret::AbstractMatrix{T}, nmcs::NC2D6C, jac::AbstractMatrix{T}, q::AbstractVector{T}) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{6,T}(ntuple(i -> q[i], Val(6)))
    fill!(jac, zero(T))
    @inbounds begin
        jac[1,3] = 2qstd[3]; jac[1,4] = 2qstd[4]
        jac[2,5] = 2qstd[5]; jac[2,6] = 2qstd[6]
        jac[3,3] =  qstd[5]; jac[3,4] =  qstd[6]
        jac[3,5] =  qstd[3]; jac[3,6] =  qstd[4]
    end
    jac_out = SMatrix{3,6,T,18}(jac) * C
    @inbounds for j in 1:6, i in 1:3
        ret[i, j] = jac_out[i, j]
    end
    ret
end
"""
Return 3D Jacobian matrix , for rigid bodies.
$(TYPEDSIGNATURES)
"""
function cstr_jacobian!(ret::AbstractMatrix{T}, nmcs::NC3D12C, jac::AbstractMatrix{T}, q::AbstractVector{T}) where {T}
    C = nmcs.conversion_to_std
    qstd = C * SVector{12,T}(ntuple(i -> q[i], Val(12)))
    fill!(jac, zero(T))
    @inbounds begin
        # u
        jac[1,4] = 2qstd[4]; jac[1,5] = 2qstd[5]; jac[1,6] = 2qstd[6]
        # v
        jac[2,7] = 2qstd[7]; jac[2,8] = 2qstd[8]; jac[2,9] = 2qstd[9]
        # w
        jac[3,10] = 2qstd[10]; jac[3,11] = 2qstd[11]; jac[3,12] = 2qstd[12]

        jac[4,7] = qstd[10]; jac[4,8] = qstd[11]; jac[4,9]  = qstd[12]
        jac[4,10] = qstd[7]; jac[4,11] = qstd[8]; jac[4,12] = qstd[9]

        jac[5,4] = qstd[10]; jac[5,5] = qstd[11]; jac[5,6]  = qstd[12]
        jac[5,10] = qstd[4]; jac[5,11] = qstd[5]; jac[5,12] = qstd[6]

        jac[6,4] = qstd[7];  jac[6,5] = qstd[8];  jac[6,6]  = qstd[9]
        jac[6,7] = qstd[4];  jac[6,8] = qstd[5];  jac[6,9]  = qstd[6]
    end
    jac_out = SMatrix{6,12,T,72}(jac) * C
    @inbounds for j in 1:12, i in 1:6
        ret[i, j] = jac_out[i, j]
    end
    ret
end

"""
Return nullspace matrix
$(TYPEDSIGNATURES)
"""
function nullspace_mat(nmcs::NC2D4C,q)
    cv = nmcs.conversion_to_std
    u,v = get_uv(nmcs,q)
    o2 = zero(u)
    ret = [
        I2_Bool  o2;
        o2 o2     v;
    ]
    cv\ret
end

function nullspace_mat(nmcs::NC2D6C,q)
    cv = nmcs.conversion_to_std
    u,v = get_uv(nmcs,q)
    o2 = zero(u)
    ret = [
        I2_Bool    o2;
        o2 o2       v;
        o2 o2      -u;
    ]
    cv\ret
end

function nullspace_mat(nmcs::NC3D3C,q)
    ret = SMatrix{3,3}(IMatrix(3))
end

function nullspace_mat(nmcs::NC3D6C,q)
    cv = nmcs.conversion_to_std
    u,v,w = get_uvw(nmcs,q)
    o3 = zero(u)
    O3 = [o3 o3 o3;]
    ret = [
        I3_Bool  o3 o3;
        O3       -w  v;
    ]
    cv\ret
end

function nullspace_mat(nmcs::NC3D12C,q)
    cv = nmcs.conversion_to_std
    u,v,w = get_uvw(nmcs,q)
    o3 = zero(u)
    O3 = [o3 o3 o3;]
    ret = [
        I3_Bool   O3;
        O3 -skew(u);
        O3 -skew(v);
        O3 -skew(w);
    ]
    # ret = [
    #     I3_Bool    O3;
    #     O3  o3  -w  v;
    #     O3   w  o3 -u;
    #     O3  -v   u o3;
    # ]
    cv\ret
end


#todo use SymmetricPacked to the end
function add_cstr_forces_jacobian!(ret, nmcs::Union{NC2D2C,NC3D3C}, λ)
    # Not constrained for a point mass
end

function add_cstr_forces_jacobian!(ret, nmcs::NC, λ)
    hessians = nmcs.hessians
    for (λj,hessian_j) in zip(λ,hessians)
        @. ret += λj * hessian_j
    end
end

#todo use SymmetricPacked to the end
function cstr_velocity_jacobian!(ret::AbstractMatrix, nmcs::Union{NC2D2C,NC3D3C},q̇)
    # Not constrained for a point mass
end

function cstr_velocity_jacobian!(ret::AbstractMatrix, nmcs::Union{NC2D2C,NC3D3C}, jac, q̇)
    cstr_velocity_jacobian!(ret, nmcs, q̇)
end

function cstr_velocity_jacobian!(ret::AbstractMatrix, nmcs::NC,q̇)
    num_of_cstr = get_num_of_cstr(nmcs)
    @inbounds for j in 1:num_of_cstr
        # q̇' * H = (H' * q̇)'
        mul!((@view ret[j,:]), transpose(nmcs.hessians[j]), q̇)
    end
    ret
end


"""
Return the ForwardDiff results for ∂Aq̇∂q
$(TYPEDSIGNATURES)
"""
function make_∂Aq̇∂q_forwarddiff(Φq,nq,nλ)
    function ∂Aq̇∂q(q̇)
        function Aq̇(q)
            Φq(q)*q̇
        end
        q̇T = eltype(q̇)
        out = zeros(q̇T,nλ,nq)
        ForwardDiff.jacobian!(out,Aq̇,ones(q̇T,nq))
    end
end

function find_independent_free_idx(nmcs::NC,q)
    ncoords = get_num_of_coords(nmcs)
    ncstr = get_num_of_cstr(nmcs)
    A = MMatrix{ncstr,ncoords, eltype(q)}(undef,)
    jac = MMatrix{ncstr,ncoords, eltype(q)}(undef,)
    cstr_jacobian!(A, nmcs, jac, q)
    col_index = GECP(A)
    col_index[size(A,1)+1:end] |> sort
end
