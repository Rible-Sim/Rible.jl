"""
Nonminimal Coordinates State Type.
$(TYPEDEF)
"""
mutable struct PresFreeCoordinatesState{T,qT,qviewT,pT,sT} <: AbstractCoordinatesState
    t::T
    q::qT
    q̇::qT
    q̈::qT
    F::qT
    λ::qT
    s::sT
    q̌::qviewT
    q̌̇::qviewT
    q̌̈::qviewT
    q̃::qviewT
    q̃̇::qviewT
    q̃̈::qviewT
    F̌::qviewT
    p::pT
    p̌::pT
    c::qT
    """
    Nonminimal Coordinates State Constructor.
    $(TYPEDSIGNATURES)
    """
    function PresFreeCoordinatesState(t,q,q̇,q̈,F,p,p̌,λ,s,free_idx,pres_idx,c)
        q̌ = @view q[free_idx]
        q̌̇ = @view q̇[free_idx]
        q̌̈ = @view q̈[free_idx]
        q̃ = @view q[pres_idx]
        q̃̇ = @view q̇[pres_idx]
        q̃̈ = @view q̈[pres_idx]
        F̌ = @view F[free_idx]
        new{typeof(t),typeof(q),typeof(q̌),typeof(p),typeof(s)}(t,q,q̇,q̈,F,λ,s,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌,p,p̌,c)
    end
end


"""
PresFreeCoordinates wrapper.
$(TYPEDEF)
"""
struct PresFreeCoordinates{N,T} <: AbstractCoordinates{N,T}
    nmcs::AbstractNonminimalCoordinates{N,T}
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
end


function PresFreeCoordinates(nmcs,pres_idx=Int[],)
    free_idx = get_free_idx(nmcs,pres_idx)
    PresFreeCoordinates(nmcs,pres_idx,free_idx)
end

function CoordinatesState(inst_state::PresFreeCoordinatesState)
    (;t,q,q̇,q̈,F,λ,s,p,c) = inst_state
    CoordinatesState(t,q,q̇,q̈,F,λ,s,p,c)
end

function RigidBodyCache(prop::RigidBodyProperty{N,T},coords::PresFreeCoordinates{N,T}) where {N,T}
    RigidBodyCache(prop,coords.nmcs)
end

function cartesian_frame2coords(coords::PresFreeCoordinates,origin_frame)
    cartesian_frame2coords(coords.nmcs,origin_frame)
end

function find_rotation(coords::PresFreeCoordinates,q)
    find_rotation(coords.nmcs,q)
end

function find_angular_velocity(coords::PresFreeCoordinates,q,v)
    find_angular_velocity(coords.nmcs,q,v)
end

# Default implementations for Coordinates
function cstr_function!(ret, coords::PresFreeCoordinates, q::AbstractVector, d = get_deform(coords.nmcs))
    ret_full = zeros(eltype(ret),NCF.get_num_of_intrinsic_cstr(coords.nmcs))
    cstr_function!(ret_full, coords.nmcs, q ,d)
    ret .= ret_full[get_cstr_idx(coords.nmcs)]
end

function cstr_function!(ret, coords::PresFreeCoordinates, cache::AbstractBodyCache, q::AbstractVector, d = get_deform(coords.nmcs))
    cstr_function!(ret, coords, q, d)
end

function cstr_jacobian!(ret, coords::PresFreeCoordinates, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    #todo make use cache.jac_cache as ret_full
    ret_full = zeros(eltype(ret),get_num_of_intrinsic_cstr(coords.nmcs),get_num_of_coords(coords.nmcs))
    cstr_jacobian!(ret_full, coords.nmcs, cache.jac_cache, inst_state.q)
    ret .= ret_full[get_cstr_idx(coords.nmcs),:]
end

function add_cstr_forces_jacobian!(ret, coords::PresFreeCoordinates, λ)
    hessians = coords.nmcs.hessians
    for (i,j) in enumerate(get_cstr_idx(coords.nmcs))
        ret .+= λ[i] * hessians[j]
    end
end

function cstr_velocity_jacobian(coords::PresFreeCoordinates,q̇)
    cstr_velocity_jacobian(coords.nmcs,q̇)
end


"""
Return free indices of natural coordinates.
$(TYPEDSIGNATURES)
"""
function get_free_idx(nmcs,pres_idx)
    deleteat!(collect(1:get_num_of_coords(nmcs)),pres_idx)
end

get_num_of_cstr(coords::PresFreeCoordinates) = get_num_of_cstr(coords.nmcs)
get_num_of_intrinsic_cstr(coords::PresFreeCoordinates) = get_num_of_intrinsic_cstr(coords.nmcs)
get_num_of_local_dims(coords::PresFreeCoordinates) = get_num_of_local_dims(coords.nmcs)
get_num_of_coords(coords::PresFreeCoordinates) = get_num_of_coords(coords.nmcs)
get_num_of_dof(coords::PresFreeCoordinates) = get_num_of_dof(coords.nmcs)

to_local_coords(coords::PresFreeCoordinates,r̄p) = to_local_coords(coords.nmcs,r̄p)
has_constant_mass_matrix(coords::PresFreeCoordinates) =  has_constant_mass_matrix(coords.nmcs)



get_coords(coords::PresFreeCoordinatesState) = coords.q
get_params(coords::PresFreeCoordinatesState) = coords.c
get_auxiliary(coords::PresFreeCoordinatesState) = coords.s
get_velocs(coords::PresFreeCoordinatesState) = coords.q̇
get_accels(coords::PresFreeCoordinatesState) = coords.q̈
get_free_coords(coords::PresFreeCoordinatesState) = coords.q̌
get_free_velocs(coords::PresFreeCoordinatesState) = coords.q̌̇
get_free_accels(coords::PresFreeCoordinatesState) = coords.q̌̈
get_pres_coords(coords::PresFreeCoordinatesState) = coords.q̃
get_pres_velocs(coords::PresFreeCoordinatesState) = coords.q̃̇
get_pres_accels(coords::PresFreeCoordinatesState) = coords.q̃̈
get_multipliers(coords::PresFreeCoordinatesState) = coords.λ

function build_joint_cache(
        coords_hen::PresFreeCoordinates,
        coords_egg::PresFreeCoordinates,
        args...
    )
    build_joint_cache(
        coords_hen.nmcs,
        coords_egg.nmcs,
        args...
    )
end

function get_joint_violations!(
        ret,
        coords_hen::PresFreeCoordinates, 
        coords_egg::PresFreeCoordinates,
        args...
    )
    get_joint_violations!(
        ret,
        coords_hen.nmcs, 
        coords_egg.nmcs,
        args...
    )
end

function get_joint_jacobian!(ret,coords_hen::PresFreeCoordinates, coords_egg::PresFreeCoordinates,
        args...
    )
    get_joint_jacobian!(ret,coords_hen.nmcs, coords_egg.nmcs,
        args...
    )
end

function add_joint_forces_jacobian!(ret,num_of_cstr,
        coords_hen::PresFreeCoordinates, coords_egg::PresFreeCoordinates,
        args...
    )
    add_joint_forces_jacobian!(ret,num_of_cstr,
        coords_hen.nmcs, coords_egg.nmcs,
        args...
    )
end

function add_joint_velocity_jacobian!(
        ret,
        num_of_cstr,
        coords_hen::PresFreeCoordinates, coords_egg::PresFreeCoordinates,
        args...
    )
    add_joint_velocity_jacobian!(
        ret,
        num_of_cstr,
        coords_hen.nmcs, coords_egg.nmcs,
        args...
    )
end

function to_position_jacobian(coords::PresFreeCoordinates,q,c)
    to_position_jacobian(coords.nmcs,q,c)
end

function nullspace_mat(coords::PresFreeCoordinates, q)
    nullspace_mat(coords.nmcs, q)
end
