


"""
Return system parameters (e.g., loci coordinates).
$(TYPEDSIGNATURES)
"""
get_params(st::AbstractStructure) = get_params(st.state.system)

"""
Return system auxiliary variables.
$(TYPEDSIGNATURES)
"""
get_auxiliary(st::AbstractStructure) = get_auxiliary(st.state.system)

"""
Return system generalized coordinates.
$(TYPEDSIGNATURES)
"""
get_coords(st::AbstractStructure) = get_coords(st.state.system)

"""
Return system generalized velocities.
$(TYPEDSIGNATURES)
"""
get_velocs(st::AbstractStructure) = get_velocs(st.state.system)

"""
Return system generalized accelerations.
$(TYPEDSIGNATURES)
"""
get_accels(st::AbstractStructure) = get_accels(st.state.system)

"""
Return system free coordinates.
$(TYPEDSIGNATURES)
"""
get_free_coords(st::AbstractStructure) = get_free_coords(st.state.system)

"""
Return system free velocities.
$(TYPEDSIGNATURES)
"""
get_free_velocs(st::AbstractStructure) = get_free_velocs(st.state.system)

"""
Return system free accelerations.
$(TYPEDSIGNATURES)
"""
get_free_accels(st::AbstractStructure) = get_free_accels(st.state.system)

"""
Return system prescribed coordinates.
$(TYPEDSIGNATURES)
"""
get_pres_coords(st::AbstractStructure) = get_pres_coords(st.state.system)

"""
Return system prescribed velocities.
$(TYPEDSIGNATURES)
"""
get_pres_velocs(st::AbstractStructure) = get_pres_velocs(st.state.system)

"""
Return system prescribed accelerations.
$(TYPEDSIGNATURES)
"""
get_pres_accels(st::AbstractStructure) = get_pres_accels(st.state.system)

"""
Return system Lagrange multipliers.
$(TYPEDSIGNATURES)
"""
get_multipliers(st::AbstractStructure) = get_multipliers(st.state.system)

"""
Return the initial state of the structure, including coordinates, velocities, and parameters.
$(TYPEDSIGNATURES)
"""
function get_initial(st::AbstractStructure)
    _,λ = check_static_equilibrium_output_multipliers(st::AbstractStructure)
    q̌ = get_free_coords(st::AbstractStructure)
    q = get_coords(st::AbstractStructure)
    c = get_params(st::AbstractStructure)
    @eponymtuple(q̌,q,s,λ,c)
end


"""
Get the parameters of a specific locus (identified by body ID and locus ID) from the structure.
$(TYPEDSIGNATURES)
"""
function get_params(st::AbstractStructure,bodyid,pid)
    (;bodyid2sys_locus_id,sys_locus_id2coords_idx) = st.connectivity
    (;c) = st.state.system
    cidx = bodyid2sys_locus_id[bodyid][pid]
    c[sys_locus_id2coords_idx[cidx]]
end

"""
Build the system parameters vector from bodies and apparatuses.
$(TYPEDSIGNATURES)
"""
function get_params(bodies,appars,cnt::AbstractConnectivity)
    num_of_dim = get_num_of_dims(bodies)
    T = get_numbertype(bodies)
    (;
        bodyid2sys_locus_id,
        sys_locus_id2coords_idx,
        num_of_params,
        num_of_sys_loci_coords,
        apparid2params_idx
    ) = cnt
    ret = zeros(T,num_of_params)
    foreach(bodies) do body
        bodyid = body.prop.id
        for i in eachindex(body.prop.loci)
            ip = bodyid2sys_locus_id[bodyid][i]
            get_params!(
                (@view ret[sys_locus_id2coords_idx[ip]]),
                body,body.prop.loci[i].position
            )
        end
    end
    foreach(appars) do appar
        get_params!(
            (@view ret[apparid2params_idx[appar.id]]),
            appar,
        )
    end
    ret
end



get_num_of_params(structure::AbstractStructure) = structure.connectivity.num_of_params


function set_C!(st::AbstractStructure,c)
    T = get_numbertype(st::AbstractStructure)
    cnt = st.connectivity
    (;bodyid2sys_locus_id,sys_locus_id2coords_idx,nc) = cnt
    foreach(st.bodies) do body
        bodyid = body.prop.id
        for i in eachindex(body.prop.loci)
            ip = bodyid2sys_locus_id[bodyid][i]
            body.state.cache.Cps[i] = body.state.cache.funcs.C(c[sys_locus_id2coords_idx[ip]])
        end
    end
end

"""
Return the violation vector (deformations or joint violations) for the entire structure.
$(TYPEDSIGNATURES)
"""
function get_d(st::AbstractStructure)
    (;num_of_cstr,bodies) = st
    cnt = st.connectivity
    (;bodyid2sys_intrinsic_cstr_idx,num_of_intrinsic_cstr) = cnt
    T = get_numbertype(st::AbstractStructure)
    d = Vector{T}(undef,num_of_cstr)
    is = Ref(num_of_intrinsic_cstr)
    foreach(bodies) do body
        bodyid = body.prop.id
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        if !isempty(memincst)
            d[memincst] .= NCF.get_deform(body.coords)
        end
    end
    foreach(st.apparatuses) do appar
        nc = appar.num_of_cstr
        d[is[]+1:is[]+nc] .= appar.violations
        is[] += nc
    end
    d
end

"""
Return System dimensions
$(TYPEDSIGNATURES)
"""
get_num_of_dims(st::AbstractStructure) = get_num_of_dims(st.bodies)
get_num_of_dims(bodies::AbstractVector) = get_num_of_dims(bodies[begin])
get_num_of_dims(bodies::TypeSortedCollection) = get_num_of_dims(bodies.data[1])


get_numbertype(st::AbstractStructure) = get_numbertype(st.bodies)
get_numbertype(bodies::AbstractVector) = get_numbertype(bodies[begin])
get_numbertype(bodies::TypeSortedCollection) = get_numbertype(bodies.data[1])


"""
Return the system's number of constraints.
$(TYPEDSIGNATURES)
"""
get_num_of_cstr(st::AbstractStructure) = st.connectivity.num_of_cstr


get_num_of_free_coords(st::AbstractStructure) = get_num_of_free_coords(st.connectivity)
get_num_of_coords(st::AbstractStructure) = get_num_of_coords(st.connectivity)


function get_num_of_cstr(rbs::TypeSortedCollection)
    num_of_intrinsic_cstr = mapreduce(get_num_of_intrinsic_cstr,+,rbs,init=0)
end

function get_num_of_aux_var(st::AbstractStructure)
    st.connectivity.num_of_aux_var
end

has_constant_mass_matrix(st::AbstractStructure) = has_constant_mass_matrix(st.bodies)

function has_constant_mass_matrix(rbs::Union{AbstractVector{<:AbstractBody},TypeSortedCollection})
    mapreduce(
        has_constant_mass_matrix,
        both_true,
        rbs;
        init = Val{true}()
    )
end

both_true(::Val{true},::Val{true}) = Val{true}()
both_true(::Val{true},::Val{false}) = Val{false}()
both_true(::Val{false},::Val{false}) = Val{false}()


"""
Return gravity
$(TYPEDSIGNATURES)
"""
get_gravity(st::AbstractStructure) = get_gravity(st.bodies)
get_gravity(rbs::TypeSortedCollection) = get_gravity(eltype(rbs.data[1]))

function _ordered_by_id(xs)
    ys = Any[]
    foreach(xs) do x
        push!(ys, x)
    end
    return _ordered_by_id(ys)
end

function _ordered_by_id(xs::AbstractVector)
    n = length(xs)
    n == 0 && return collect(xs)

    out = eltype(xs) === Any ? Vector{Any}(undef, n) : Vector{eltype(xs)}(undef, n)
    assigned = falses(n)

    for x in xs
        id = get_id(x)
        if !(1 <= id <= n) || assigned[id]
            # Fallback without sorting to avoid Julia 1.12 sort inference crash on Vector{Any}.
            return collect(xs)
        end
        out[id] = x
        assigned[id] = true
    end

    all(assigned) || return collect(xs)
    return out
end

get_bodies(st::AbstractStructure) = _ordered_by_id(st.bodies)
get_apparatuses(st::AbstractStructure) = _ordered_by_id(st.apparatuses)
