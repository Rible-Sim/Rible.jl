
@inline function _zero_matrix!(M::SparseMatrixCSC)
    fill!(nonzeros(M), zero(eltype(M)))
    M
end
@inline function _zero_matrix!(M)
    M .= zero(eltype(M))
end

@inline function _add_block!(dest::SparseMatrixCSC{T,Int}, idx::AbstractVector{<:Integer}, blk::SparseMatrixCSC{T,Int}) where {T}
    @assert length(idx) == size(blk,1) == size(blk,2)
    for lcol in 1:size(blk,2)
        gcol = idx[lcol]
        for k in nzrange(blk, lcol)
            grow = idx[blk.rowval[k]]
            dest[grow, gcol] += blk.nzval[k]
        end
    end
    dest
end
@inline function _add_block!(dest::AbstractMatrix, idx::AbstractVector{<:Integer}, blk::AbstractMatrix)
    @assert length(idx) == size(blk,1) == size(blk,2)
    dest[idx, idx] .+= blk
    dest
end

function assemble_M!(M,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    _zero_matrix!(M)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        _add_block!(M, memfull, body.cache.inertia_cache.M)
    end
    # @assert issymmetric(M)
    M
end

function assemble_M(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    T = get_numbertype(st)
    M = spzeros(T,num_of_full_coords,num_of_full_coords)
    assemble_M!(M,st)
    M
end

function assemble_M̌(st::AbstractStructure)
    (;sys_free_coords_idx) = st.connectivity
    M = assemble_M(st)
    M̌ = Symmetric(M[sys_free_coords_idx,sys_free_coords_idx])
end

function assemble_∂Mq̇∂q!(∂Mq̇∂q,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    _zero_matrix!(∂Mq̇∂q)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        _add_block!(∂Mq̇∂q, memfull, body.cache.inertia_cache.∂Mq̇∂q)
    end
end

function assemble_∂Mq̇∂q(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    T = get_numbertype(st)
    ∂Mq̇∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    assemble_∂Mq̇∂q!(∂Mq̇∂q,st::AbstractStructure)
    ∂Mq̇∂q
    # symsparsecsr(M;symmetrize=true)
end

function assemble_∂T∂qᵀ!(∂T∂qᵀ,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    ∂T∂qᵀ .= 0
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂T∂qᵀ[memfull] .+= body.state.cache.∂T∂qᵀ
    end
end

function update_inertia_cache!(st::AbstractStructure)
    update_inertia_cache!(st, has_constant_mass_matrix(st))
end

function update_inertia_cache!(st::AbstractStructure, ::Val{true})
    # constant mass matrix, do nothing
    nothing 
end

function update_inertia_cache!(st::AbstractStructure, ::Val{false})
    cache = st.cache.system
    if !cache.dirty
        return cache
    end

    assemble_M!(cache.M, st)
    
    # Invert M globally 
    M⁻¹_dense = inv(Matrix(cache.M))
    cache.M⁻¹ = sparse(M⁻¹_dense)

    assemble_∂Mq̇∂q!(cache.∂Mq̇∂q, st)

    # Correct global calculation of ∂M⁻¹p∂q
    # ∂(M⁻¹p)/∂q = -M⁻¹ * (∂M/∂q * q̇) = -M⁻¹ * ∂Mq̇∂q
    cache.∂M⁻¹p∂q = sparse(-M⁻¹_dense * cache.∂Mq̇∂q)
    
    cache.dirty = false
    return cache
end

function assemble_M⁻¹!(M⁻¹, st::AbstractStructure)
    update_inertia_cache!(st)
    M⁻¹ .= st.cache.system.M⁻¹
    M⁻¹
end

function assemble_M⁻¹(st::AbstractStructure)
    update_inertia_cache!(st)
    copy(st.cache.system.M⁻¹)
end

function assemble_∂M⁻¹p∂q!(∂M⁻¹p∂q, st::AbstractStructure)
    update_inertia_cache!(st)
    ∂M⁻¹p∂q .= st.cache.system.∂M⁻¹p∂q
    ∂M⁻¹p∂q
end

function assemble_∂M⁻¹p∂q(st::AbstractStructure)
    update_inertia_cache!(st)
    copy(st.cache.system.∂M⁻¹p∂q)
end

function check_stale_inertia_cache(st::AbstractStructure)
    if st.cache.system.dirty
        @warn "Inertia cache is dirty! The computed M⁻¹, ∂Mq̇∂q, ∂M⁻¹p∂q might be stale. Recomputing now... (Consider calling update_inertia_cache! earlier)"
        update_inertia_cache!(st)
    end
end

function assemble_∂T∂qᵀ(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    T = get_numbertype(st)
    ∂T∂qᵀ = zeros(T,num_of_full_coords)
    assemble_∂T∂qᵀ!(∂T∂qᵀ,st)
    ∂T∂qᵀ
end

"""
Return System mass matrices
$(TYPEDSIGNATURES)
"""
function build_mass_matrices(structure::AbstractStructure,)
    sys_free_coords_idx = get_free_coords_idx(structure.connectivity)
    sys_pres_coords_idx = get_pres_coords_idx(structure.connectivity)
    M = assemble_M(structure,) |> sparse |> Symmetric
    M⁻¹ = M |> Matrix |> inv |> sparse |> Symmetric
    M̌   = M[sys_free_coords_idx,sys_free_coords_idx] |> sparse |> Symmetric
    M̌⁻¹ = M̌ |> Matrix |> inv |> sparse |> Symmetric
    Ḿ = M[sys_free_coords_idx,:]            |> sparse
    M̄ = M[sys_free_coords_idx,sys_pres_coords_idx] |> sparse
    @eponymtuple(M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)
end

function mass_center(structure)
    N = get_num_of_dims(structure)
    T = get_numbertype(structure)
    rg = zeros(T,N)
    M = Ref(zero(T))
    foreach(structure.bodies) do body
        M[] += body.prop.mass
        rg .+= body.prop.mass * body.state.mass_locus_state.frame.position
    end
    rg ./= M[]

end
