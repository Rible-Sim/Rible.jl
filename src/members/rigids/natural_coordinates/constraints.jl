
# """
# Return 2D or 3D local natural coords deformations for rigid bars。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds get_deform(ū::AbstractVector) = sqrt(ū⋅ū)

@inline @inbounds function get_deform(nmcs::Union{LNC2D1P,LNC3D1P})
    (;r̄i) = nmcs
    get_deform(r̄i)
end

@inline @inbounds function get_deform(nmcs::Union{LNC2D1P1V,LNC3D1P1V,LNC2D2P,LNC3D2P})
    (;ū) = nmcs
    get_deform(ū)
end

# """
# Return 2D natural coodinates deformations , for rigid bodies。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(ū,v̄)
    sqrt(ū⋅ū),sqrt(v̄⋅v̄),ū⋅v̄
end

@inline @inbounds function get_deform(nmcs::LNC2D6C)
    (;ū,v̄) = nmcs
    get_deform(ū,v̄)
end

# """
# Return 3D natural coodinates deformations , for rigid bodies。
# $(TYPEDSIGNATURES)
# """
@inline @inbounds function get_deform(ū,v̄,w̄)
    u_sqrt = sqrt(ū⋅ū)
    v_sqrt = sqrt(v̄⋅v̄)
    w_sqrt = sqrt(w̄⋅w̄)
    vw = v̄⋅w̄
    uw = ū⋅w̄
    uv = ū⋅v̄
    u_sqrt,v_sqrt,w_sqrt,vw,uw,uv
end

@inline @inbounds function get_deform(nmcs::LNC3D12C)
    (;ū,v̄,w̄) = nmcs
    get_deform(ū,v̄,w̄)
end

# Intrinsic Constraints
## Intrinsic Constraints: Dispatch
"""
Return 2D or 3D intrinsic cstr(s) ，用于Dispatch。
$(TYPEDSIGNATURES)
"""
function make_cstr_function(nmcs::LNC,cstr_idx)
    deforms = get_deform(nmcs::LNC)
    make_cstr_function(nmcs,cstr_idx,deforms)
end

function make_inner_cstr_function(func,deforms)
    @inline @inbounds function ret_func(q)
        func(q,deforms)
    end
    @inline @inbounds function ret_func(q,d)
        func(q,d)
    end
    ret_func
end

"""
Return 2D or 3D intrinsic cstr(s) for point mass。
$(TYPEDSIGNATURES)
"""
function make_cstr_function(nmcs::Union{LNC2D2C,LNC3D3C},cstr_idx,deforms)
    @inline @inbounds function _inner_cstr_function(q,d)
        nothing
    end
    make_inner_cstr_function(_inner_cstr_function,deforms)
end

"""
Return 2D or 3D intrinsic cstr(s) for rigid bars。
$(TYPEDSIGNATURES)
"""
function make_cstr_function(nmcs::Union{LNC2D4C,LNC3D6C},cstr_idx,deforms)
    ndim = get_num_of_dims(nmcs)
    cv = nmcs.conversion
    @inline @inbounds function _inner_cstr_function(q,d)
        qstd = cv*q
        u = @view qstd[ndim+1:2ndim]
        all = [u⋅u - d^2]
        all[cstr_idx]
    end
    make_inner_cstr_function(_inner_cstr_function,deforms)
end

"""
Return 2D intrinsic cstr(s) , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_cstr_function(nmcs::LNC2D6C,cstr_idx,deforms)
    cv = nmcs.conversion
    @inline @inbounds function _inner_cstr_function(q,d)
        qstd = cv*q
        u = @view qstd[3:4]
        v = @view qstd[5:6]
        all = [
            (u⋅u - d[1]^2), 
            (v⋅v - d[2]^2), 
            (u⋅v - d[3])
        ]
        all[cstr_idx]
    end
    make_inner_cstr_function(_inner_cstr_function,deforms)
end


"""
Return 3D intrinsic cstr(s) , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_cstr_function(nmcs::LNC3D12C,cstr_idx,deforms)
    cv = nmcs.conversion
    @inline @inbounds function _inner_cstr_function(q,d)
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
    make_inner_cstr_function(_inner_cstr_function,deforms)
end

## Jacobians


"""
Return 2D or 3D Jacobian matrix for rigid bars。
$(TYPEDSIGNATURES)
"""
function make_cstr_jacobian(nmcs::Union{LNC2D2C,LNC3D3C},free_coords_idx,cstr_idx)
    @inline @inbounds function inner_cstr_jacobian(q)
        nothing
    end
end

"""
Return 2D or 3D Jacobian matrix for rigid bars。
$(TYPEDSIGNATURES)
"""
function make_cstr_jacobian(nmcs::Union{LNC2D4C,LNC3D6C},free_coords_idx,cstr_idx)
    ndim = get_num_of_dims(nmcs)
    cv = nmcs.conversion
    @inline @inbounds function inner_cstr_jacobian(q)
        qstd = cv*q
        u = @view qstd[ndim+1:2ndim]
        ret = zeros(eltype(q),1,2ndim)
        ret[ndim+1:2ndim] = 2u
        @view (ret*cv)[cstr_idx,free_coords_idx]
    end
end

"""
Return 2D Jacobian matrix , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_cstr_jacobian(nmcs::LNC2D6C,free_coords_idx,cstr_idx)
    @inline @inbounds function inner_cstr_jacobian(q)
        u,v = get_uv(nmcs,q)
        ret = zeros(eltype(q),3,6)
        ret[1,3:4] = 2u
        ret[2,5:6] = 2v
        ret[3,3:4] =  v
        ret[3,5:6] =  u
        @view (ret*cv)[cstr_idx,free_coords_idx]
    end
end
"""
Return 3D Jacobian matrix , for rigid bodies。
$(TYPEDSIGNATURES)
"""
function make_cstr_jacobian(nmcs::LNC3D12C,free_coords_idx,cstr_idx)
    @inline @inbounds function inner_cstr_jacobian(q)
        u,v,w = get_uvw(nmcs,q)
        cv = nmcs.conversion
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

        @view (ret*cv)[cstr_idx,free_coords_idx]
    end
end

"""
Return nullspace matrix
$(TYPEDSIGNATURES)
"""
function make_nullspace(nmcs::LNC2D4C)
    cv = nmcs.conversion
    @inline @inbounds function inner_nullspace(q)
        u,v = get_uv(nmcs,q)
        o2 = zero(u)
        ret = [
            I2_Bool  o2;
            o2 o2     v;
        ]
        cv\ret
    end
end

function make_nullspace(nmcs::LNC2D6C)
    cv = nmcs.conversion
    @inline @inbounds function inner_nullspace(q)
        u,v = get_uv(nmcs,q)
        o2 = zero(u)
        ret = [
            I2_Bool    o2;
            o2 o2       v;
            o2 o2      -u;
        ]
        cv\ret
    end
end

function make_nullspace(nmcs::LNC3D6C)
    cv = nmcs.conversion
    @inline @inbounds function inner_nullspace(q)
        u,v,w = get_uvw(nmcs,q)
        o3 = zero(u)
        O3 = [o3 o3 o3;]
        ret = [
            I3_Bool  o3 o3;
            O3       -w  v;
        ]
        cv\ret
    end
end

function make_nullspace(nmcs::LNC3D12C)
    cv = nmcs.conversion
    @inline @inbounds function inner_nullspace(q)
        u,v,w = get_uvw(nmcs,q)
        o3 = zero(u)
        O3 = [o3 o3 o3;]
        # ret = [
        #     I3_Bool   O3;
        #     O3 -skew(u);
        #     O3 -skew(v);
        #     O3 -skew(w);
        # ]
        ret = [
            I3_Bool    O3;
            O3  o3  -w  v;
            O3   w  o3 -u;
            O3  -v   u o3;
        ]
        cv\ret
    end
end


get_idx(nmcs::Union{LNC2D4C,LNC3D6C}) = [
    [CartesianIndex(2,2),CartesianIndex(2,2)],
]

get_idx(nmcs::LNC2D6C) = [
    [CartesianIndex(2,2),CartesianIndex(2,2)],
    [CartesianIndex(3,3),CartesianIndex(3,3)],
    [CartesianIndex(2,3),CartesianIndex(3,2)]
]

get_idx(nmcs::LNC3D12C) = [
    [CartesianIndex(2,2),CartesianIndex(2,2)],
    [CartesianIndex(3,3),CartesianIndex(3,3)],
    [CartesianIndex(4,4),CartesianIndex(4,4)],
    [CartesianIndex(3,4),CartesianIndex(4,3)],
    [CartesianIndex(2,4),CartesianIndex(4,2)],
    [CartesianIndex(2,3),CartesianIndex(3,2)]
]

#todo cache cstr_hessians
function make_cstr_hessians(nmcs::LNC)
    cv = nmcs.conversion
    nld = get_num_of_local_dims(nmcs)
    ndim = get_num_of_dims(nmcs)
    I_Bool = IMatrix(ndim)
    idx = get_idx(nmcs)
    cstr_hessians = [
        begin
            ret_raw = zeros(Int,nld+1,nld+1)
            for ij in id
                ret_raw[ij] += 1
            end
            SymmetricPacked(transpose(cv)*kron(ret_raw,I_Bool)*cv)
        end
        for id in idx
    ]
end

#todo use SymmetricPacked to the end
function make_cstr_forces_jacobian(nmcs::Union{LNC2D2C,LNC3D3C},free_coords_idx,cstr_idx)
    function cstr_forces_jacobian(λ)
        nothing
    end
end


# CoordinateFunctions
"""
封装有函数的natural coodinates.
$(TYPEDEF)
"""
struct CoordinateFunctions{nmcsType,IndicesType}
    nmcs::nmcsType
    free_coords_idx::IndicesType
    cstr_idx::IndicesType
end

"""
封装有函数的natural coodinates 构造子。
$(TYPEDSIGNATURES)
"""
function CoordinateFunctions(nmcs,free_coords_idx,cstr_idx)
    # # cstr_forces_jacobian = make_cstr_forces_jacobian_forwarddiff(Φq,nq,nuc)
    # cstr_forces_jacobian = make_cstr_forces_jacobian(nmcs,free_coords_idx,cstr_idx)
    # # ∂Aq̇∂q = make_∂Aq̇∂q_forwarddiff(Φq,nuc,num_of_cstr)
    # ∂Aq̇∂q = make_∂Aq̇∂q(nmcs,free_coords_idx,cstr_idx)
    CoordinateFunctions(
        nmcs,
        free_coords_idx,
        cstr_idx
    )
end

function make_cstr_function(cf::CoordinateFunctions)
    make_cstr_function(cf.nmcs,cf.cstr_idx)
end

function make_cstr_jacobian(cf::CoordinateFunctions)
    make_cstr_jacobian(
        cf.nmcs,
        cf.free_coords_idx,
        cf.cstr_idx,
    )
end

function make_cstr_forces_jacobian(
        cf::CoordinateFunctions,
        cstr_hessians
    )
    (;nmcs::LNC,free_coords_idx,cstr_idx) = cf
    function cstr_forces_jacobian(λ)
        ret = [
            begin
                a = cstr_hessians[j][free_coords_idx,free_coords_idx] .* λ[i]
                # display(a)
                a 
            end
            for (i,j) in enumerate(cstr_idx)
        ]
        sum(ret)
    end
end

function make_∂Aq̇∂q(nmcs::Union{LNC2D2C,LNC3D3C},free_coords_idx,cstr_idx)
    function ∂Aq̇∂q(q̇)
        nothing
    end
end

# this is wrong
function make_∂Aq̇∂q(nmcs::LNC,free_coords_idx,cstr_idx)
    cstr_hessians = make_cstr_hessians(nmcs)
    function ∂Aq̇∂q(q̇)
        q̇uc = @view q̇[free_coords_idx]
        ret = [
            begin
                a = transpose(q̇uc)*cstr_hessians[j][free_coords_idx,free_coords_idx]
                # display(a)
                a 
            end
            for j in cstr_idx
        ]
        sum(ret)
    end
end

"""
Return cstr_forces_jacobian的前向自动微分结果。
$(TYPEDSIGNATURES)
"""
function make_cstr_forces_jacobian_forwarddiff(Φq,nq,nuc)
    function cstr_forces_jacobian(λ)
        function ATλ(q)
            transpose(Φq(q))*λ
        end
        λT = eltype(λ)
        out = zeros(λT,nuc,nq)
        ForwardDiff.jacobian!(out,ATλ,ones(λT,nq))
    end
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