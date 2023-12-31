# """
# Return 2D or 3D local natural coords deformations for rigid bars。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds get_deform(ū::AbstractVector) = SVector(sqrt(ū⋅ū))

@inline @inbounds function get_deform(nmcs::Union{NC2D2C,NC3D3C})
    (;r̄i) = nmcs.data
    get_deform(r̄i)
end

@inline @inbounds function get_deform(nmcs::Union{NC2D4C,NC3D6C})
    (;ū) = nmcs.data
    get_deform(ū)
end

# """
# Return 2D natural coodinates deformations , for rigid bodies。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(nmcs::NC2D6C)
    (;ū,v̄) = nmcs.data
    SVector(sqrt(ū⋅ū),sqrt(v̄⋅v̄),ū⋅v̄)
end

# """
# Return 3D natural coodinates deformations , for rigid bodies。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(nmcs::NC3D12C)
    (;ū,v̄,w̄) = nmcs.data
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
Return 2D or 3D intrinsic cstr(s) for point mass。
$(TYPEDSIGNATURES)
"""
function cstr_function(nmcs::Union{NC2D2C,NC3D3C},cstr_idx,q,deforms = get_deform(nmcs))
    eltype(q)[]
end

"""
Return 2D or 3D intrinsic cstr(s) for rigid bars。
$(TYPEDSIGNATURES)
"""
function cstr_function(nmcs::Union{NC2D4C,NC3D6C},cstr_idx,q,d = get_deform(nmcs))
    ndim = get_num_of_dims(nmcs)
    cv = nmcs.conversion_to_std
    qstd = cv*q
    u = @view qstd[ndim+1:2ndim]
    all = [u⋅u - d[1]^2]
    @view all[cstr_idx]
end

"""
Return 2D intrinsic cstr(s) , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function cstr_function(nmcs::NC2D6C,cstr_idx,q,d = get_deform(nmcs))
    cv = nmcs.conversion_to_std
    qstd = cv*q
    u = @view qstd[3:4]
    v = @view qstd[5:6]
    all = [
        (u⋅u - d[1]^2), 
        (v⋅v - d[2]^2), 
        (u⋅v - d[3])
    ]
    @view all[cstr_idx]
end


"""
Return 3D intrinsic cstr(s) , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function cstr_function(nmcs::NC3D12C,cstr_idx,q,d = get_deform(nmcs))
    cv = nmcs.conversion_to_std
    qstd = cv*q
    u = @view qstd[4:6]
    v = @view qstd[7:9]
    w = @view qstd[10:12]
    all = [
        u⋅u - d[1]^2, 
        v⋅v - d[2]^2, 
        w⋅w - d[3]^2, 
        v⋅w - d[4], 
        u⋅w - d[5], 
        u⋅v - d[6]
    ]
    @view all[cstr_idx]
end

## Jacobians
"""
Return 2D or 3D Jacobian matrix for rigid bars。
$(TYPEDSIGNATURES)
"""
function cstr_jacobian(nmcs::Union{NC2D2C,NC3D3C},free_idx,cstr_idx,q)
    ndim = get_num_of_dims(nmcs)
    zeros(eltype(q),0,ndim)
end

"""
Return 2D or 3D Jacobian matrix for rigid bars。
$(TYPEDSIGNATURES)
"""
function cstr_jacobian(nmcs::Union{NC2D4C,NC3D6C},free_idx,cstr_idx,q)
    ndim = get_num_of_dims(nmcs)
    cv = nmcs.conversion_to_std
    qstd = cv*q
    u = @view qstd[ndim+1:2ndim]
    ret = zeros(eltype(q),1,2ndim)
    ret[ndim+1:2ndim] = 2u
    @view (ret*cv)[cstr_idx,free_idx]
end

"""
Return 2D Jacobian matrix , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function cstr_jacobian(nmcs::NC2D6C,free_idx,cstr_idx,q)
    cv = nmcs.conversion_to_std
    u,v = get_uv(nmcs,q)
    ret = zeros(eltype(q),3,6)
    ret[1,3:4] = 2u
    ret[2,5:6] = 2v
    ret[3,3:4] =  v
    ret[3,5:6] =  u
    @view (ret*cv)[cstr_idx,free_idx]
end
"""
Return 3D Jacobian matrix , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function cstr_jacobian(nmcs::NC3D12C,free_idx,cstr_idx,q)
    cv = nmcs.conversion_to_std
    u,v,w = get_uvw(nmcs,q)
    ret = zeros(eltype(q), 6, 12)
    ret[1,4:6]   = 2u
    ret[2,7:9]   = 2v
    ret[3,10:12] = 2w

    ret[4 ,7:9]  = w
    ret[4,10:12] = v

    ret[5, 4:6]  = w
    ret[5,10:12] = u

    ret[6,4:6] =  v
    ret[6,7:9] =  u

    @view (ret*cv)[cstr_idx,free_idx]
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
function cstr_forces_jacobian(nmcs::Union{NC2D2C,NC3D3C},free_idx,cstr_idx,λ)
    nothing
end

function cstr_forces_jacobian(nmcs::NC,free_idx,cstr_idx,λ)
    ret = [
        begin
            a = -λ[i] .* nmcs.hessians[j][free_idx,free_idx]
            # display(a)
            a 
        end
        for (i,j) in enumerate(cstr_idx)
    ]
    sum(ret)
end

#todo use SymmetricPacked to the end
function cstr_velocity_jacobian(nmcs::Union{NC2D2C,NC3D3C},free_idx,cstr_idx,q̇)
    nothing
end

function cstr_velocity_jacobian(nmcs::NC,free_idx,cstr_idx,q̇)
    reduce(vcat,
        [
            transpose(q̇)*nmcs.hessians[j][:,free_idx]
            for j in cstr_idx
        ]
    )
end


"""
Return ∂Aq̇∂q的前向自动微分结果。
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
    free_idx = collect(1:get_num_of_coords(nmcs))
    cstr_idx = collect(1:get_num_of_cstr(nmcs))
    A = cstr_jacobian(nmcs,free_idx,cstr_idx,q)
    col_index = GECP(A)
    col_index[size(A,1)+1:end] |> sort
end