
"""
StructureState Type.
$(TYPEDEF)
"""
struct StructureState{sysT, msT} <: AbstractStructureState
    system::sysT
    members::msT
end

"""
State Constructor.
$(TYPEDSIGNATURES)
"""
function StructureState(bodies,apparatuses,cnt::Connectivity)
    (;
        num_of_bodies,
        num_of_full_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        num_of_aux_var,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_full_coords,
        apparid2sys_aux_var_idx,
        bodyid2sys_loci_coords_idx,
    ) = cnt
    
    T = get_numbertype(bodies)
    t = zero(T)
    q = zeros(T,num_of_full_coords)
    q̇ = zero(q)
    q̈ = zero(q)
    F = zero(q)
    λ = zeros(T,num_of_cstr)
    s = zeros(T,num_of_aux_var)
    p = zero(q)
    c = get_params(bodies,apparatuses,cnt)
    system = CoordinatesState(t,q,q̇,q̈,F,λ,s,p,c)
    members = [
        begin
            qmem = @view q[bodyid2sys_full_coords[bodyid]]
            q̇mem = @view q̇[bodyid2sys_full_coords[bodyid]]
            q̈mem = @view q̈[bodyid2sys_full_coords[bodyid]]
            Fmem = @view F[bodyid2sys_full_coords[bodyid]]
            λmem = @view λ[bodyid2sys_intrinsic_cstr_idx[bodyid]]
            smem = @view s[Int[]]
            cmem = @view c[bodyid2sys_loci_coords_idx[bodyid]]
            pmem = zero(p[bodyid2sys_full_coords[bodyid]])
            CoordinatesState(
                t,
                qmem,q̇mem,q̈mem,Fmem,
                λmem,smem,pmem,cmem
            )
        end
        for bodyid = 1:num_of_bodies
    ]
    M = spzeros(T,num_of_full_coords,num_of_full_coords)
    foreach(bodies) do body
        q, q̇ = body_state2coords_state(body)
        members[body.prop.id].q .= q
        members[body.prop.id].q̇ .= q̇
        memfull = bodyid2sys_full_coords[body.prop.id]
        M[memfull,memfull] .+= body.cache.inertia_cache.M
    end
    system.p = M*system.q̇
    
    foreach(apparatuses) do appar
        prepare_cache!(appar,cnt)
        s[apparid2sys_aux_var_idx[appar.id]] = get_auxilary(appar)
    end
    StructureState(system,members)
end


"""
Structure Type.
$(TYPEDEF)
"""
struct Structure{BodyType,TenType,CntType,StateType,CRType,CacheType} <: AbstractStructure{BodyType,TenType,CntType}
    num_of_dim::Int
    bodies::BodyType
    apparatuses::TenType
    connectivity::CntType
    state::StateType
    contacts_related::CRType
    cache::CacheType
end

"""
Structure Constructor.
$(TYPEDSIGNATURES)
"""
function Structure(bodies,apparatuses,cnt::AbstractConnectivity)
    num_of_dims = get_num_of_dims(bodies)
    state = StructureState(bodies,apparatuses,cnt)
    cache = StructureCache(bodies,cnt)
    (;bodyid2sys_locus_id) = cnt
    num_of_sys_loci = length.(bodyid2sys_locus_id) |> sum
    activated_bits = BitVector(undef,num_of_sys_loci)
    persistent_bits = BitVector(undef,num_of_sys_loci)
    T = get_numbertype(bodies)
    friction_coefficients = ones(T,num_of_sys_loci)
    restitution_coefficients = zeros(T,num_of_sys_loci)
    gaps = fill(typemax(T),num_of_sys_loci)

    # initilize
    foreach(bodies) do body
        (;prop,) = body
        bid = prop.id
        (;loci) = prop
        friction_coefficients[bodyid2sys_locus_id[bid]] .= [locus.friction_coefficient for locus in loci]
        restitution_coefficients[bodyid2sys_locus_id[bid]] .= [locus.restitution_coefficient for locus in loci]
    end

    contacts_related = @eponymtuple(
        activated_bits,
        persistent_bits,
        friction_coefficients,
        restitution_coefficients,
        gaps
    )
    structure = Structure(
        num_of_dims,
        bodies,
        apparatuses,
        cnt,
        state,
        contacts_related,
        cache
    )
    check_jacobian_singularity(structure)
    check_constraints_consistency(structure)
    structure
end

"""
VizStructure Type.
$(TYPEDEF)
"""
struct VizStructure{BodyType,TenType,CntType,StateType,CacheType} <: AbstractStructure{BodyType,TenType,CntType}
    num_of_dim::Int
    bodies::BodyType
    apparatuses::TenType
    connectivity::CntType
    state::StateType
    cache::CacheType
end


function VizStructure(structure::AbstractStructure)
    bodies = get_bodies(structure) .|> Observable
    apparatuses = get_apparatuses(structure) .|> Observable
    (;num_of_dim, state, connectivity, cache) = structure
    viz_st = VizStructure(
        num_of_dim,
        bodies,
        apparatuses,
        connectivity,
        state,
        cache
    )
end



function update!(st::AbstractStructure, field=NoField(); )
    clear_forces!(st)
    stretch!(st)
    update_bodies!(st)
    update_apparatuses!(st)
    apply_field!(st, field)
    assemble_forces!(st)
end

function lazy_update!(st::AbstractStructure, field=NoField(); )
    clear_forces!(st)
    lazy_update_bodies!(st)
    update_apparatuses!(st)
    apply_field!(st, field)
    assemble_forces!(st)
end

function update!(st::AbstractStructure,inst_state::AbstractCoordinatesState)
    (;q,q̇,t,s) = inst_state
    st.state.system.t = t
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    st.state.system.s .= s
    update!(st)
end

function build_F(st,bodyid,pid,f)
    T = get_numbertype(st)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity
    F = zeros(T,num_of_full_coords)
    foreach(st.bodies) do body
        if body.prop.id == bodyid
            C = body.cache.Cps[pid]
            memfull = bodyid2sys_full_coords[bodyid]
            F[memfull] = transpose(C)*f
        end
    end
    reshape(F,:,1)
end

struct StructureCache{sysType}
    system::sysType
    function StructureCache(st::AbstractStructure)
        system_cache = InertiaCache(st,)
        new{typeof(system_cache)}(
            system_cache,
        )
    end
    function StructureCache(bodies,cnt::AbstractConnectivity)
        system_cache = InertiaCache(bodies, cnt)
        new{typeof(system_cache)}(
            system_cache,
        )
    end
end

function InertiaCache(st::AbstractStructure)
    (;M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄) = build_mass_matrices(st)
    T = get_numbertype(st)
    (;num_of_full_coords) = st.connectivity
    ∂Mq̇∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    ∂M⁻¹p∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    ∂T∂qᵀ = spzeros(T,num_of_full_coords)
    Ṁq̇ = spzeros(T,num_of_full_coords)
    InertiaCache(
        M,M⁻¹,
        # M̌,M̌⁻¹,
        # Ḿ,M̄,
        ∂Mq̇∂q,∂M⁻¹p∂q,
        Ṁq̇,∂T∂qᵀ,
        true
    )
end

"""
check constraints jacobian singularity
$(TYPEDSIGNATURES)
"""
function check_jacobian_singularity(st::AbstractStructure)
    (;bodies,state) = st
    A = cstr_jacobian(st,st.state.system)
    sys_rank = rank(A)
    if sys_rank < minimum(size(A))
        @show A
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(A)))"
        foreach(bodies) do body
            bodyid = body.prop.id
            q_body = state.members[bodyid].q
            A_body = zeros(eltype(q_body),get_num_of_cstr(body.coords),get_num_of_coords(body.coords))
            cstr_jacobian!(A_body, body, state.members[bodyid])
            body_rank = rank(A_body)
            if body_rank < minimum(size(A_body))
                @warn "The $(bodyid)th body's Jacobian is singular: rank(A(q))=$(body_rank)<$(minimum(size(A_body)))"
            end
        end
    end
end

function check_constraints_consistency(st::AbstractStructure;tol=1e-14)
    (;bodies,apparatuses,state,connectivity) = st
    cnt = connectivity
    inst_state = state.system
    T = get_numbertype(st)
    Φ = cstr_function(st,inst_state)
    norm_position = norm(Φ,Inf)
    if norm_position > tol
        @warn "System's position-level constraints are inconsistent: norm_position=$(norm_position)"
        foreach(bodies) do body
            bodyid = body.prop.id
            q_body = state.members[bodyid].q
            Φ_body = zeros(T,get_num_of_cstr(body.coords))
            cstr_function!(Φ_body, body.coords, q_body)
            norm_position_body = norm(Φ_body,Inf)
            if norm_position_body > tol
                @warn "The $(bodyid)th body's position-level constraints are inconsistent" norm_position_body
                @debug "The $(bodyid)th body's position-level constraints" Φ_body
            end
        end
        foreach(apparatuses) do appar
            Φ_joint = zeros(T,length(cnt.apparid2sys_extrinsic_cstr_idx[appar.id]))
            cstr_function!(Φ_joint, appar,cnt,inst_state)
            norm_position_joint = norm(Φ_joint,Inf)
            if norm_position_joint > tol
                @warn "The $(appar.id)th apparatus's position-level constraints are inconsistent" norm_position_joint
                @debug "The $(appar.id)th apparatus's position-level constraints" Φ_joint
            end
        end
    end
    Aq̇ = cstr_jacobian(st,inst_state)*get_free_velocs(inst_state)
    norm_velocity = norm(Aq̇,Inf)
    if norm_velocity > tol
        @warn "System's velocity-level constraints are inconsistent" norm_velocity
        foreach(bodies) do body
            bodyid = body.prop.id
            q_body = state.members[bodyid].q
            q̇_body = state.members[bodyid].q̇
            A_body = zeros(eltype(q_body),get_num_of_cstr(body.coords),get_num_of_coords(body.coords))
            cstr_jacobian!(A_body, body, state.members[bodyid])
            Aq̇_body = A_body*q̇_body
            norm_velocity_body = norm(Aq̇_body,Inf)
            if norm_velocity_body > tol
                @warn "The $(bodyid)th body's velcoity-level constraints are inconsistent" norm_velocity_body
                @debug "The $(bodyid)th body's velcoity-level constraints" Aq̇_body
            end
        end
        foreach(apparatuses) do appar
            appar_sys_full_idx = cnt.apparid2sys_full_coords_idx[appar.id]
            q̇_jointed = state.system.q̇[appar_sys_full_idx]
            A_jointed = zeros(T,length(cnt.apparid2sys_extrinsic_cstr_idx[appar.id]),length(appar_sys_full_idx) )
            cstr_jacobian!(A_jointed, appar,cnt,inst_state)
            Aq̇_joint = A_jointed*q̇_jointed
            norm_velocity_joint = norm(Aq̇_joint,Inf)
            if norm_velocity_joint > tol
                @warn "The $(appar.id)th apparatus's velcoity-level constraints are inconsistent" Aq̇_joint
                @debug "The $(appar.id)th apparatus's velcoity-level constraints" Aq̇_joint
            end
        end
    end
end


"""
Return System Kinetic energy
$(TYPEDSIGNATURES)
"""
function kinetic_energy(st::AbstractStructure)
    M = assemble_M(st)
    (;q̇) = st.state.system
    T = 1/2*transpose(q̇)*M*q̇
end

"""
Return System potential energy
$(TYPEDSIGNATURES)
"""
function potential_energy_field(st::AbstractStructure, field = NoField();)
    V = Ref(zero(get_numbertype(st)))
    foreach(st.bodies) do body
        V[] += potential_energy_field(body, field, st.state.members[body.prop.id].q)
    end
    V[]
end

function potential_energy_strain(st::AbstractStructure)
    V = Ref(zero(get_numbertype(st)))
    foreach(st.bodies) do body
        V[] += potential_energy_strain(body)
    end
    V[]
end

function potential_energy(st::AbstractStructure, field = NoField();)
    V = Ref(zero(get_numbertype(st)))
    V[] += potential_energy_strain(st)
    V[] += potential_energy_field(st, field)
    foreach(st.apparatuses) do appar
        if !(appar.force isa NoForce)
            V[] += potential_energy(appar.force)
        end
    end
    V[]
end

"""
Return System mechanical_energy
$(TYPEDSIGNATURES)
"""
function mechanical_energy(st::AbstractStructure, field=NoField();)
    T = kinetic_energy(st)
    V = potential_energy(st, field)
    E = T+V
    @eponymtuple(T,V,E)
end

function mechanical_energy!(st::AbstractStructure, field=NoField())
    update!(st, field)
    mechanical_energy(st)
end

function build_T(st,i)
    cnt = st.connectivity
    (;num_of_full_coords,bodyid2sys_full_coords) = cnt
    T = spzeros(Int,length(bodyid2sys_full_coords[i]),num_of_full_coords)
    # Ť = zeros(Int,length(bodyid2sys_full_coords[i]),num_of_free_coords)
    # T̃ = zeros(Int,length(bodyid2sys_full_coords[i]),num_of_pres_coords)
    for (i,j) in enumerate(bodyid2sys_full_coords[i])
        T[i,j] = 1
    end
    T
end

