struct Apparatus{jointType,forceType}
    id::Int
    joint::jointType
    force::forceType
    num_of_add_var::Int
    full_coords_idx::Vector{Int}
    free_coords_idx::Vector{Int}
end

function Base.isless(a::Apparatus,b::Apparatus)
    isless(a.id,b.id)
end

get_id(appar::Apparatus) = appar.id
