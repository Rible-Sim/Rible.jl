abstract type AbstractBody{N,T} end
abstract type AbstractBodyProperty{N,T} end
abstract type AbstractBodyState{N,T} end
abstract type AbstractBodyCache end

update_cache!(body::AbstractBody,q,q̇) = update_cache!(body.cache,body.coords,body.prop,q,q̇)
lazy_update_state!(body::AbstractBody,q,q̇) = lazy_update_state!(body.state,body.coords,body.cache,body.prop,q,q̇)
update_state!(body::AbstractBody,q,q̇) = update_state!(body.state,body.coords,body.cache,body.prop,q,q̇)
update_loci_states!(body::AbstractBody,q,q̇) = update_loci_states!(body.state,body.coords,body.cache,body.prop,q,q̇)
stretch_loci!(body::AbstractBody,c) = stretch_loci!(body.coords,body.cache,body.prop,c)
update_transformations!(body::AbstractBody,q) = update_transformations!(body.coords,body.cache,body.state,body.prop,q)
function body_state2coords_state(body::AbstractBody)
    body_state2coords_state(body.state,body.coords)
end

"""
Nonminimal Coordinates State Type.
$(TYPEDEF)
"""
mutable struct NonminimalCoordinatesState{T,qT,qviewT,pT}
    t::T
    q::qT
    q̇::qT
    q̈::qT
    F::qT
    λ::qT
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
    function NonminimalCoordinatesState(t,q,q̇,q̈,F,p,p̌,λ,free_idx,pres_idx,c)
        q̌ = @view q[free_idx]
        q̌̇ = @view q̇[free_idx]
        q̌̈ = @view q̈[free_idx]
        q̃ = @view q[pres_idx]
        q̃̇ = @view q̇[pres_idx]
        q̃̈ = @view q̈[pres_idx]
        F̌ = @view F[free_idx]
        new{typeof(t),typeof(q),typeof(q̌),typeof(p)}(t,q,q̇,q̈,F,λ,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌,p,p̌,c)
    end
end

struct NonminimalCoordinates{nmcsType}
    nmcs::nmcsType
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
    cstr_idx::Vector{Int}
end

function NonminimalCoordinates(nmcs::NCF.NC,
        pres_idx=Int[],
        cstr_idx=get_cstr_idx(nmcs)
    )
    free_idx = NCF.get_free_idx(nmcs,pres_idx)
    NonminimalCoordinates(nmcs,pres_idx,free_idx,cstr_idx)
end

function make_cstr_function(coords::NonminimalCoordinates)
    make_cstr_function(coords.nmcs,coords.cstr_idx)
end

function make_cstr_jacobian(coords::NonminimalCoordinates)
    make_cstr_jacobian(coords.nmcs,coords.free_idx,coords.cstr_idx)
end


struct NonminimalCoordinatesCache{MType,JType,HType,GType}
    cstr_hessians::Vector{HType}
    M::MType
    M⁻¹::MType
    ∂Mq̇∂q::JType
    ∂M⁻¹p∂q::JType
    "Coriolis force"
    Ṁq̇::GType
    "Centripetal force"
    ∂T∂qᵀ::GType
end