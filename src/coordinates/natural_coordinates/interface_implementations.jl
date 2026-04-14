"""
Natural Coordinates Interface Implementations

This module implements the coordinate interface methods for natural coordinates.
"""

# Core interface implementations
"""
Return the dimension of space.
$(TYPEDSIGNATURES)
"""
get_num_of_dims(::LNCData{N,M,T,L}) where {N,M,T,L} = N
get_num_of_dims(nmcs::NC) = get_num_of_dims(nmcs.data)

"""
Return local dimension of natural coordinates.
$(TYPEDSIGNATURES)
"""
get_num_of_local_dims(::LNCData{N,M,T,L}) where {N,M,T,L} = M
get_num_of_local_dims(nmcs::NC) = get_num_of_local_dims(nmcs.data)

"""
Return the number of coordinates.
$(TYPEDSIGNATURES)
"""
get_num_of_coords(::LNCData{N,M,T,L}) where {N,M,T,L} = N+L
get_num_of_coords(nmcs::NC) = get_num_of_coords(nmcs.data)

"""
Return the number of degrees of freedom
$(TYPEDSIGNATURES)
"""
get_num_of_dof(nmcs::NC) =  get_num_of_coords(nmcs) - get_num_of_cstr(nmcs)
get_num_of_cstr(nmcs::NC) = nmcs.num_of_cstr
"""
Return the number of constraints.
$(TYPEDSIGNATURES)
"""
get_num_of_intrinsic_cstr(::NC2D2C) = 0
get_num_of_intrinsic_cstr(::NC3D3C) = 0
get_num_of_intrinsic_cstr(::NC2D4C) = 1
get_num_of_intrinsic_cstr(::NC2D6C) = 3
get_num_of_intrinsic_cstr(::NC3D6C) = 1
get_num_of_intrinsic_cstr(::NC3D12C) = 6

# Property access for LNCData
function Base.getproperty(data::LNCData,p::Symbol)
    if     (p == :ū) 
        return data.X̄[:,1]
    elseif (p == :v̄)
        return data.X̄[:,2]
    elseif (p == :w̄)
        return data.X̄[:,3]
    elseif (p === :r̄j)
        return data.r̄i + data.ū
    elseif  (p === :r̄k)
        return data.r̄i + data.v̄
    elseif  (p === :r̄l)
        return data.r̄i + data.w̄
    else # fallback to getfield
        return getfield(data, p)
    end
end

# Mass matrix interface
has_constant_mass_matrix(::NC) = Val{true}()


function get_cstr_idx(nmcs::NC)
    nmcs.cstr_idx
end
