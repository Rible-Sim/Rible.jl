abstract type AbstractStructure end

"""
StructureState Type.
$(TYPEDEF)
"""
struct StructureState{sysT, msT}
    system::sysT
    members::msT
end

"""
State Constructor.
$(TYPEDSIGNATURES)
"""
function StructureState(bodies,apparatuses,cnt::Connectivity)
    (;numbered,indexed) = cnt
    (;bodyid2sys_loci_coords_idx) = numbered
    (;
        num_of_bodies,
        num_of_full_coords,
        num_of_free_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        sys_free_coords_idx,
        sys_pres_coords_idx,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_full_coords,
        bodyid2sys_free_coords,
    ) = indexed
    pres_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    free_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    foreach(bodies) do body
        bodyid = body.prop.id
        pres_idx_by_body[bodyid] = body.coords.pres_idx
        free_idx_by_body[bodyid] = body.coords.free_idx
    end
    T = get_numbertype(bodies)
    t = zero(T)
    q = Vector{T}(undef,num_of_full_coords)
    q̇ = zero(q)
    q̈ = zero(q)
    F = zero(q)
    λ = Vector{T}(undef,num_of_cstr)
    c = get_local_coords(bodies,numbered)
    p = zero(q)
    p̌ = Vector{T}(undef,num_of_free_coords)
    system = NonminimalCoordinatesState(t,q,q̇,q̈,F,p,p̌,λ,sys_free_coords_idx,sys_pres_coords_idx,c)
    members = [
        begin
            qmem = @view q[bodyid2sys_full_coords[bodyid]]
            q̇mem = @view q̇[bodyid2sys_full_coords[bodyid]]
            q̈mem = @view q̈[bodyid2sys_full_coords[bodyid]]
            Fmem = @view F[bodyid2sys_full_coords[bodyid]]
            λmem = @view λ[bodyid2sys_intrinsic_cstr_idx[bodyid]]
            cmem = @view c[bodyid2sys_loci_coords_idx[bodyid]]
            pmem = zero(p[bodyid2sys_full_coords[bodyid]])
            p̌mem = zero(p[bodyid2sys_free_coords[bodyid]])
            NonminimalCoordinatesState(
                t,
                qmem,q̇mem,q̈mem,Fmem,pmem,p̌mem,λmem,
                free_idx_by_body[bodyid],
                pres_idx_by_body[bodyid],
                cmem
            )
        end
        for bodyid = 1:num_of_bodies
    ]
    foreach(bodies) do body
        q, q̇ = body_state2coords_state(body)
        members[body.prop.id].q .= q
        members[body.prop.id].q̇ .= q̇
    end
    StructureState(system,members)
end

"""
Structure Type.
$(TYPEDEF)
"""
struct Structure{BodyType,TenType,CntType,StateType} <: AbstractStructure
    num_of_dim::Int
    bodies::BodyType
    apparatuses::TenType
    connectivity::CntType
    state::StateType
end

"""
Structure Constructor.
$(TYPEDSIGNATURES)
"""
function Structure(bodies,apparatuses,cnt = Connectivity(index(bodies,apparatuses),number(bodies,apparatuses),))
    num_of_dim = get_num_of_dims(bodies)
    state = StructureState(bodies,apparatuses,cnt)
    st = Structure(
        num_of_dim,
        bodies,
        apparatuses,
        cnt,
        state
    )
    check_jacobian_singularity(st)
    check_constraints_consistency(st)
    st
end

function Structure(bodies,cnt::Connectivity)
    apparatuses = (cables = Int[],)
    Structure(bodies,apparatuses,cnt)
end


function update!(st::AbstractStructure; gravity=false)
    clear_forces!(st)
    stretch_rigids!(st)
    update_bodies!(st)
    update_apparatuses!(st)
    if gravity
        apply_gravity!(st)
    end
    assemble_forces!(st)
end

function lazy_update!(st::AbstractStructure; gravity=false)
    clear_forces!(st)
    lazy_update_bodies!(st)
    update_apparatuses!(st)
    if gravity
        apply_gravity!(st)
    end
    assemble_forces!(st)
end

function update!(st::Structure,q::AbstractVector,q̇::AbstractVector=zero(q))
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    update!(st)
end

function assemble_M!(M,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    M .= 0
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        M[memfull,memfull] .+= body.cache.inertia_cache.M
    end
    # @assert issymmetric(M)
    M
end

function assemble_M(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    M = spzeros(T,num_of_full_coords,num_of_full_coords)
    assemble_M!(M,st)
    M
end

function assemble_M⁻¹!(M⁻¹,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    M⁻¹ .= 0
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        M⁻¹[memfull,memfull] .+= body.cache.inertia_cache.M⁻¹
    end
    # @assert issymmetric(M⁻¹)
    M⁻¹
end

function assemble_M⁻¹(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    M⁻¹ = spzeros(T,num_of_full_coords,num_of_full_coords)
    assemble_M⁻¹!(M⁻¹,st)
    M⁻¹
end

function assemble_M̌(st::AbstractStructure)
    (;sys_free_coords_idx) = st.connectivity.indexed
    M = assemble_M(st)
    M̌ = Symmetric(M[sys_free_coords_idx,sys_free_coords_idx])
end

function assemble_∂Mq̇∂q!(∂Mq̇∂q,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    ∂Mq̇∂q .= 0
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂Mq̇∂q[memfull,memfull] .+= body.cache.inertia_cache.∂Mq̇∂q
    end
end

function assemble_∂Mq̇∂q(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂Mq̇∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    assemble_∂Mq̇∂q!(∂Mq̇∂q,st::AbstractStructure)
    ∂Mq̇∂q
    # symsparsecsr(M;symmetrize=true)
end

function assemble_∂M⁻¹p∂q!(∂M⁻¹p∂q,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    ∂M⁻¹p∂q .= 0
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂M⁻¹p∂q[memfull,memfull] .+= body.cache.inertia_cache.∂M⁻¹p∂q
    end
    # symsparsecsr(M;symmetrize=true)
end

function assemble_∂M⁻¹p∂q(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂M⁻¹p∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    assemble_∂M⁻¹p∂q!(∂M⁻¹p∂q,st)
    ∂M⁻¹p∂q
end

function assemble_∂T∂qᵀ!(∂T∂qᵀ,st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    ∂T∂qᵀ .= 0
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂T∂qᵀ[memfull] .+= body.state.cache.∂T∂qᵀ
    end
end

function assemble_∂T∂qᵀ(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂T∂qᵀ = zeros(T,num_of_full_coords)
    assemble_∂T∂qᵀ!(∂T∂qᵀ,st)
    ∂T∂qᵀ
end

function make_M!(st)
    function inner_M!(M,q)
        update_bodies!(st,q)
        assemble_M!(M,st)
    end
end


function make_M⁻¹!(st)
    function inner_M⁻¹!(M⁻¹,q)
        update_bodies!(st,q)
        assemble_M⁻¹!(M⁻¹,st)
    end
end

function make_M_and_Jac_M!(st)
    function inner_M_and_Jac_M!(M,∂Mq̇∂q,q,q̇)
        update_bodies!(st,q,q̇)
        assemble_M!(M,st)
        assemble_∂Mq̇∂q!(∂Mq̇∂q,st)
    end
end

function make_M⁻¹_and_Jac_M⁻¹!(st)
    function inner_M⁻¹_and_Jac_M⁻¹!(M⁻¹,∂M⁻¹p∂q,q,q̇)
        update_bodies!(st,q,q̇)
        assemble_M⁻¹!(M⁻¹,st)
        assemble_∂M⁻¹p∂q!(∂M⁻¹p∂q,st)
    end
end

"""
Return System mass matrices
$(TYPEDSIGNATURES)
"""
function build_mass_matrices(structure::Structure,)
    (;sys_free_coords_idx,sys_pres_coords_idx) = structure.connectivity.indexed
    M = assemble_M(structure,) |> sparse |> Symmetric
    M⁻¹ = M |> Matrix |> inv |> sparse |> Symmetric
    M̌   = M[sys_free_coords_idx,sys_free_coords_idx] |> sparse |> Symmetric
    M̌⁻¹ = M̌ |> Matrix |> inv |> sparse |> Symmetric
    Ḿ = M[sys_free_coords_idx,:]            |> sparse
    M̄ = M[sys_free_coords_idx,sys_pres_coords_idx] |> sparse
    @eponymtuple(M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)
end

function cstr_function(structure::AbstractStructure,q,c=get_local_coords(structure))
    (;bodies,apparatuses,connectivity) = structure
    (;numbered,indexed,) = connectivity
    (;  
        num_of_cstr,
        num_of_free_coords,
        num_of_full_coords,
        sys_pres_coords_idx,
        sys_free_coords_idx,
        bodyid2sys_full_coords,
        bodyid2sys_pres_coords,
        bodyid2sys_free_coords,
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        apparid2sys_extrinsic_cstr_idx
    ) = indexed
    
    ret = Vector{eltype(q)}(undef,num_of_cstr)
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memfree = bodyid2sys_free_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        if !isempty(memincst)
            ret[memincst] .= cstr_function(
                body.coords,q[memfull]#,d[memincst]
            )
        end
    end
    foreach(apparatuses) do appar
        jointexcst = num_of_intrinsic_cstr.+apparid2sys_extrinsic_cstr_idx[appar.id]
        ret[jointexcst] .= cstr_function(appar,structure,
            q,
            c
        )
    end
    ret
end

function cstr_jacobian(structure::AbstractStructure,q,c=get_local_coords(structure))
    (;bodies,apparatuses,connectivity) = structure
    (;indexed,numbered) = connectivity
    (;
        num_of_cstr,
        num_of_free_coords,
        num_of_full_coords,
        sys_pres_coords_idx,
        sys_free_coords_idx,
        bodyid2sys_full_coords,
        bodyid2sys_pres_coords,
        bodyid2sys_free_coords,
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        apparid2sys_free_coords_idx,
        apparid2sys_extrinsic_cstr_idx
    ) = indexed
    ret = zeros(eltype(q),num_of_cstr,num_of_free_coords)
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memfree = bodyid2sys_free_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        if !isempty(memincst)
            ret[memincst,memfree] .= cstr_jacobian(
                body.coords,q[memfull]
            )
        end
    end
    foreach(apparatuses) do appar
        (;id,) = appar
        ext = apparid2sys_extrinsic_cstr_idx[id]
        if !isempty(ext)
            jointexcst = num_of_intrinsic_cstr.+apparid2sys_extrinsic_cstr_idx[id]
            jointed_sys_free_idx = apparid2sys_free_coords_idx[id]
            ret[jointexcst,jointed_sys_free_idx] .= cstr_jacobian(appar,structure,q)
        end
    end
    ret
end

make_cstr_jacobian(structure::AbstractStructure) = (q) -> cstr_jacobian(structure,q)
make_cstr_function(structure::AbstractStructure) = (q) -> cstr_function(structure,q)

function make_cstr_function(structure::AbstractStructure,q0::AbstractVector)
    (;connectivity) = structure
    (;indexed,) = connectivity
    (;
        num_of_full_coords,
        sys_pres_coords_idx,
        sys_free_coords_idx,
    ) = indexed
    function inner_cstr_function(q̌,d=get_d(structure),c=get_local_coords(structure))
        q = Vector{eltype(q̌)}(undef,num_of_full_coords)
        q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
        q[sys_free_coords_idx] .= q̌
        cstr_function(q,d,c)
    end
end

function make_cstr_jacobian(structure::AbstractStructure,q0::AbstractVector)
    (;connectivity) = structure
    (;indexed,) = connectivity
    (;
        num_of_full_coords,
        sys_pres_coords_idx,
        sys_free_coords_idx,
    ) = indexed
    function inner_cstr_jacobian(q̌,c=get_local_coords(structure))
        q = Vector{eltype(q̌)}(undef,num_of_full_coords)
        q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
        q[sys_free_coords_idx] .= q̌
        cstr_jacobian(structure,q,c)
    end
    inner_cstr_jacobian
end

function build_F̌(st,bodyid,pid,f)
    T = get_numbertype(st)
    (;num_of_free_coords,bodyid2sys_free_coords) = st.connectivity.indexed
    F̌ = zeros(T,num_of_free_coords)
    foreach(st.bodies) do body
        if body.prop.id == bodyid
            C = body.state.cache.Cps[pid]
            memfree = bodyid2sys_free_coords[bodyid]
            free_idx = body.state.cache.free_idx
            F̌[memfree] = (transpose(C)*f)[free_idx,:]
        end
    end
    reshape(F̌,:,1)
end

struct StructureCache{sysType}
    system::sysType
    function StructureCache(st::Structure)
        system_cache = InertiaCache(st,)
        new{typeof(system_cache)}(
            system_cache,
        )
    end
end

function InertiaCache(st::Structure)
    (;M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄) = build_mass_matrices(st)
    T = get_numbertype(st)
    (;num_of_full_coords) = st.connectivity.indexed
    ∂Mq̇∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    ∂M⁻¹p∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    ∂T∂qᵀ = spzeros(T,num_of_full_coords)
    Ṁq̇ = spzeros(T,num_of_full_coords)
    InertiaCache(
        M,M⁻¹,
        # M̌,M̌⁻¹,
        # Ḿ,M̄,
        ∂Mq̇∂q,∂M⁻¹p∂q,
        Ṁq̇,∂T∂qᵀ
    )
end

"""
check constraints jacobian singularity
$(TYPEDSIGNATURES)
"""
function check_jacobian_singularity(st)
    (;bodies,state) = st
    q = get_coords(st)
    A = cstr_jacobian(st,q)
    sys_rank = rank(A)
    if sys_rank < minimum(size(A))
        @show A
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(A)))"
        foreach(bodies) do body
            bodyid = body.prop.id
            free_idx = body.coords.free_idx
            q_body = state.members[bodyid].q
            A_body = cstr_jacobian(body.coords,q_body)
            body_rank = rank(A_body)
            if body_rank < minimum(size(A_body))
                @warn "The $(bodyid)th body's Jacobian is singular: rank(A(q))=$(body_rank)<$(minimum(size(A_rb)))"
            end
        end
    end
end

function check_constraints_consistency(st;tol=1e-14)
    (;bodies,apparatuses,state,connectivity) = st
    (;indexed,numbered) = connectivity
    q = get_coords(st)
    q̌̇ = get_free_velocs(st)
    Φ = cstr_function(st,q)
    norm_position = norm(Φ)
    if norm_position > tol
        @warn "System's position-level constraints are inconsistent: norm_position=$(norm_position)"
        foreach(bodies) do body
            bodyid = body.prop.id
            free_idx = body.coords.free_idx
            q_body = state.members[bodyid].q
            Φ_body = cstr_function(body.coords,q_body)
            norm_position_body = norm(Φ_body)
            if norm_position_body > tol
                @warn "The $(bodyid)th body's position-level constraints are inconsistent: norm_position_body=$(norm_position_body)"
            end
        end
        foreach(apparatuses) do appar
            Φ_joint = cstr_jacobian(appar,st,q)
            norm_position_joint = norm(Φ_joint)
            if norm_position_joint > tol
                @warn "The $(appar.id)th apparatus's position-level constraints are inconsistent: norm_position_joint=$(norm_position_joint)"
            end
        end
    end
    Aq̇ = cstr_jacobian(st,q)*q̌̇
    norm_velocity = norm(Aq̇)
    if norm_velocity > tol
        @warn "System's velocity-level constraints are inconsistent: norm_velocity=$(norm_velocity)"
        foreach(bodies) do body
            bodyid = body.prop.id
            q_body = state.members[bodyid].q
            q̇_body = state.members[bodyid].q̇
            Aq̇_body = cstr_jacobian(body.coords,q_body)*q̇_body
            norm_velocity_body = norm(Aq̇_body)
            if norm_velocity_body > tol
                @warn "The $(bodyid)th body's velcoity-level constraints are inconsistent: norm_velocity_body=$(norm_velocity_body)"
            end
        end
        foreach(apparatuses) do appar
            sys_free_coords_idx = indexed.apparid2sys_free_coords_idx[appar.id]
            q̇_jointed = state.system.q̌̇[sys_free_coords_idx]
            Aq̇_joint = cstr_jacobian(appar,st,q)*q̇_jointed
            norm_velocity_joint = norm(Aq̇_joint)
            if norm_velocity_joint > tol
                @warn "The $(appar.id)th apparatus's velcoity-level constraints are inconsistent: Aq̇_joint=$(Aq̇_joint)"
            end
        end
    end
end

function analyse_slack(st::AbstractStructure,verbose=false)
    (;cables) = st.apparatuses
    slackcases = [cable.id for cable in cables if cable.slack && (cable.state.length <= cable.state.restlen)]
    if verbose && !isempty(slackcases)
        @show slackcases
    end
    slackcases
end

"""
Return System Kinetic energy
$(TYPEDSIGNATURES)
"""
function kinetic_energy(st::Structure)
    M = assemble_M(st)
    (;q̇) = st.state.system
    T = 1/2*transpose(q̇)*M*q̇
end

"""
Return System potential energy
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(st::Structure;factor=1)
    V = Ref(zero(get_numbertype(st)))
    foreach(st.bodies) do body
        V[] += factor*potential_energy_gravity(body)
    end
    V[]
end

function potential_energy_strain(st::Structure)
    V = Ref(zero(get_numbertype(st)))
    foreach(st.bodies) do body
        V[] += potential_energy_strain(body)
    end
    V[]
end

function potential_energy(st::Structure;gravity=false)
    V = Ref(zero(get_numbertype(st)))
    V[] += potential_energy_strain(st)
    if gravity
        V[] += potential_energy_gravity(st)
    end
    foreach(st.apparatuses) do appar
        if !(appar.force isa Nothing)
            V[] += potential_energy(appar.force)
        end
    end
    V[]
end

"""
Return System mechanical_energy
$(TYPEDSIGNATURES)
"""
function mechanical_energy(st::Structure;gravity=false)
    T = kinetic_energy(st)
    V = potential_energy(st;gravity)
    E = T+V
    @eponymtuple(T,V,E)
end

function mechanical_energy!(st::Structure)
    update!(st)
    mechanical_energy(st)
end
