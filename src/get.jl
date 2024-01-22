
get_local_coords(bot::Robot) = get_local_coords(bot.structure)
get_s(bot::Robot) = get_s(bot.structure)

get_coords(st::Structure) = copy(st.state.system.q)
get_velocs(st::Structure) = copy(st.state.system.q̇)
get_free_coords(st::Structure) = copy(st.state.system.q̌)
get_free_velocs(st::Structure) = copy(st.state.system.q̌̇)

function get_λ(st::Structure)
    st.state.system.λ
end

"""
Return System Initial State.
$(TYPEDSIGNATURES)
"""
function get_initial(st::Structure)
    _,λ = check_static_equilibrium_output_multipliers(st::Structure)
    q̌ = get_free_coords(st::Structure)
    q = get_coords(st::Structure)
    ℓ = get_cables_len(st::Structure)
    s = 1 ./ℓ
    d = get_d(st::Structure)
    c = get_local_coords(st::Structure)
    k = get_cables_stiffness(st::Structure)
    μ = get_cables_restlen(st::Structure)
    @eponymtuple(q̌,q,s,λ,d,c,k,μ,)
end

function get_polyvar(st::Structure)
    (;num_of_cstr,connectivity,apparatuses) = st
    (;cables) = apparatuses
    (;indexed,numbered,) = connectivity
    (;nc) = numbered
    (;num_of_intrinsic_cstr) = indexed
    (;num_of_free_coords) = indexed
    ncables = length(cables)
    # state variables
    @polyvar q̌[1:num_of_free_coords]
    @polyvar s[1:ncables]
    @polyvar λ[1:num_of_cstr]
    # parameters
    @polyvar d[1:num_of_cstr]
    @polyvar c[1:nc]
    @polyvar k[1:ncables]
    @polyvar μ[1:ncables]
    @eponymtuple(q̌,s,λ,d,c,k,μ,)
end

function get_s(st::Structure)
    1 ./get_cables_len(st::Structure)
end

function get_local_coords(st::Structure,bodyid,pid)
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = st.connectivity.numbered
    (;c) = st.state.system
    cidx = bodyid2sys_loci_idx[bodyid][pid]
    c[sys_loci2coords_idx[cidx]]
end

get_local_coords(st::Structure) = copy(st.state.system.c)

function get_local_coords(bodies,numbered::Numbered)
    num_of_dim = get_num_of_dims(bodies)
    T = get_numbertype(bodies)
    (;bodyid2sys_loci_idx,sys_loci2coords_idx,nc) = numbered
    ret = zeros(T,nc)
    foreach(bodies) do body
        bodyid = body.prop.id
        for i in eachindex(body.prop.loci)
            ip = bodyid2sys_loci_idx[bodyid][i]
            ret[sys_loci2coords_idx[ip]] .= get_local_coords(body,body.prop.loci[i].position)
        end
    end
    ret
end

function get_local_coords(body::RigidBody,r̄p)
    to_local_coords(body.coords.nmcs,r̄p)
end

function get_local_coords(body::FlexibleBody,r̄p)
    body.cache.funcs.x(r̄p)
end

function set_C!(st::Structure,c)
    T = get_numbertype(st::Structure)
    (;numbered,indexed) = st.connectivity
    (;bodyid2sys_loci_idx,sys_loci2coords_idx,nc) = numbered
    foreach(st.bodies) do body
        bodyid = body.prop.id
        for i in eachindex(body.prop.loci)
            ip = bodyid2sys_loci_idx[bodyid][i]
            body.state.cache.Cps[i] = body.state.cache.funcs.C(c[sys_loci2coords_idx[ip]])
        end
    end
end

function get_d(st::Structure)
    (;num_of_cstr,bodies) = st
    (;indexed,jointed) = st.connectivity
    (;bodyid2sys_intrinsic_cstr_idx,num_of_intrinsic_cstr) = indexed
    T = get_numbertype(st::Structure)
    d = Vector{T}(undef,num_of_cstr)
    is = Ref(num_of_intrinsic_cstr)
    foreach(bodies) do body
        bodyid = body.prop.id
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        if !isempty(memincst)
            d[memincst] .= NCF.get_deform(body.coords.nmcs)
        end
    end
    foreach(jointed.joints) do joint
        nc = joint.num_of_cstr
        d[is[]+1:is[]+nc] .= joint.violations
        is[] += nc
    end
    d
end

"""
Return System dimensions
$(TYPEDSIGNATURES)
"""
get_num_of_dims(bot::Robot) = get_num_of_dims(bot.structure)
get_num_of_dims(st::AbstractStructure) = get_num_of_dims(st.bodies)
get_num_of_dims(bodies::AbstractVector) = get_num_of_dims(bodies[begin])
get_num_of_dims(bodies::TypeSortedCollection) = get_num_of_dims(bodies.data[1])
get_num_of_dims(body::BodyType) where {BodyType<:AbstractBody{N,T}} where {N,T} = N

get_numbertype(bot::Robot) = get_numbertype(bot.structure)
get_numbertype(st::AbstractStructure) = get_numbertype(st.bodies)
get_numbertype(bodies::AbstractVector) = get_numbertype(bodies[begin])
get_numbertype(bodies::TypeSortedCollection) = get_numbertype(bodies.data[1])
get_numbertype(body::BodyType) where {BodyType<:AbstractBody{N,T}} where {N,T} = T

"""
Return the system's number of constraints.
$(TYPEDSIGNATURES)
"""
get_num_of_cstr(st::Structure) = st.connectivity.indexed.num_of_cstr

get_num_of_free_coords(st::Structure) = st.connectivity.indexed.num_of_free_coords


function get_num_of_cstr(rbs::TypeSortedCollection)
    num_of_intrinsic_cstr = mapreduce(get_num_of_intrinsic_cstr,+,rbs,init=0)
end
get_num_of_intrinsic_cstr(body::AbstractRigidBody) = get_num_of_cstr(body.coords.nmcs)
get_num_of_coords(body::AbstractRigidBody) = get_num_of_coords(body.coords.nmcs)
get_num_of_dof(body::AbstractRigidBody) = get_num_of_dof(body.coords.nmcs)
get_num_of_local_dims(body::AbstractRigidBody) = get_num_of_local_dims(body.coords.nmcs)

get_num_of_intrinsic_cstr(body::AbstractFlexibleBody) = get_num_of_cstr(body.coords.nmcs)
get_num_of_coords(body::AbstractFlexibleBody) = get_num_of_coords(body.coords.nmcs)
get_num_of_dof(body::AbstractFlexibleBody) = get_num_of_dof(body.coords.nmcs)
get_num_of_local_dims(body::AbstractFlexibleBody) = get_num_of_local_dims(body.coords.nmcs)


has_constant_mass_matrix(bot::Robot) = has_constant_mass_matrix(bot.structure)
has_constant_mass_matrix(st::Structure) = has_constant_mass_matrix(st.bodies)
has_constant_mass_matrix(rbs::AbstractVector{<:AbstractRigidBody}) = has_constant_mass_matrix(eltype(rbs))


has_constant_mass_matrix(body::AbstractBody) = has_constant_mass_matrix(body.coords.nmcs)
has_constant_mass_matrix(::NCF.NC) = Val{true}()
has_constant_mass_matrix(::QCF.QC) = Val{false}()
has_constant_mass_matrix(::ANCF.ANC) = Val{true}()

function has_constant_mass_matrix(rbs::TypeSortedCollection)
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
get_gravity(bot::Robot) = get_gravity(bot.structure)
get_gravity(st::Structure) = get_gravity(st.bodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rbs::TypeSortedCollection) = get_gravity(eltype(rbs.data[1]))
get_gravity(body::AbstractBody) = get_gravity(typeof(body))
get_gravity(::Type{<:AbstractBody{2,T}}) where {T} = SVector(zero(T),        -9.81*one(T))
get_gravity(::Type{<:AbstractBody{3,T}}) where {T} = SVector(zero(T),zero(T),-9.81*one(T))

get_cables_len(bot::Robot) = get_cables_len(bot.structure)
get_cables_deform(bot::Robot) = get_cables_deform(bot.structure)
get_cables_restlen(bot::Robot) = get_cables_restlen(bot.structure)
get_cables_len_dot(bot::Robot) = get_cables_len_dot(bot.structure)
get_cables_tension(bot::Robot) = get_cables_tension(bot.structure)
get_cables_stiffness(bot::Robot) = get_cables_stiffness(bot.structure)
get_cables_force_density(bot::Robot) = get_cables_force_density(bot.structure)

get_bodies(bot::Robot) = get_bodies(bot.structure)
get_bodies(st::AbstractStructure) = sort(st.bodies)

get_rigidbars(bot::Robot) = get_rigidbars(bot.structure)

function get_rigidbars(st::AbstractStructure)
    filter(get_bodies(st::Structure)) do body
        body.state.cache.funcs.nmcs isa Union{NCF.NC2D4C,NCF.NC3D6C}
    end
end

get_connected(bot::Robot) = get_connected(bot.structure)
get_connected(st::AbstractStructure) = sort(st.connectivity.tensioned.connected)

function get_cables_len!(st::Structure,q)
    update_bodies!(st::Structure,q,zero(q))
    update_cables_apply_forces!(st::Structure)
    get_cables_len(st::Structure)
end

"""
Return System DistanceSpringDamper stiffness coefficients
$(TYPEDSIGNATURES)
"""
function get_cables_stiffness(st::Structure)
    [s.k for s in st.apparatuses.cables]
end

"""
Return System DistanceSpringDamper current Length.
$(TYPEDSIGNATURES)
"""
function get_cables_len(st::Structure)
    [s.state.length for s in st.apparatuses.cables]
end

function get_cables_len_dot(st::Structure)
    [s.state.lengthdot for s in st.apparatuses.cables]
end

"""
Return System DistanceSpringDamper deformation。
$(TYPEDSIGNATURES)
"""
function get_cables_deform(st::Structure)
    [s.state.length - s.state.restlen for s in st.apparatuses.cables]
end

"""
Return System DistanceSpringDamper Restlength。
$(TYPEDSIGNATURES)
"""
function get_cables_restlen(st::Structure)
    [s.state.restlen for s in st.apparatuses.cables]
end

"""
Return System DistanceSpringDamper Tension.
$(TYPEDSIGNATURES)
"""
function get_cables_tension(st::Structure)
    [s.state.tension for s in st.apparatuses.cables]
end

"""
Set cables' tension
$(TYPEDSIGNATURES)
"""
function set_cables_tension!(st::Structure,fs)
    for (s,f) in zip(st.apparatuses.cables,fs)
        s.state.tension = f
    end
end


"""
Return System DistanceSpringDamper force density。
$(TYPEDSIGNATURES)
"""
function get_cables_force_density(st::Structure)
    [s.state.tension/s.state.length for s in st.apparatuses.cables]
end

"""
Return System DistanceSpringDamper Initial Length.
$(TYPEDSIGNATURES)
"""
function get_original_restlen(botinput::Robot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_cables_restlen(bot.structure)
end

function force_densities_to_restlen(st::Structure,γs)
    [
    begin
        l = s.state.length
        l̇ = s.state.lengthdot
        k = s.k
        c = s.c
        u = l-(γ*l-c*l̇)/k
    end
        for (γ,s) in zip(γs,st.apparatuses.cables)]
end
