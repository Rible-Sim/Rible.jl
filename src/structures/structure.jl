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
function StructureState(bodies,tensiles,cnt::Connectivity{<:Any,<:Any,<:NamedTuple{(:connected, )},<:Any})
    (;numbered,indexed,jointed) = cnt
    (;bodyid2sys_loci_coords_idx) = numbered
    (;num_of_full_coords,num_of_intrinsic_cstr,sys_free_idx,sys_pres_idx) = indexed
    (;bodyid2sys_intrinsic_cstr_idx,bodyid2sys_full_coords) = indexed
    (;num_of_extrinsic_cstr) = jointed
    num_of_cstr = num_of_intrinsic_cstr + num_of_extrinsic_cstr
    nb = length(bodies)
    pres_idx_by_mem = Vector{Vector{Int}}(undef,nb)
    free_idx_by_mem = Vector{Vector{Int}}(undef,nb)
    foreach(bodies) do body
        bodyid = body.prop.id
        pres_idx_by_mem[bodyid] = body.state.cache.pres_idx
        free_idx_by_mem[bodyid] = body.state.cache.free_idx
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
    p̌ = zero(q)
    system = NonminimalCoordinatesState(t,q,q̇,q̈,F,p,p̌,λ,sys_free_idx,sys_pres_idx,c)
    members = [
        begin
            qmem = @view q[bodyid2sys_full_coords[bodyid]]
            q̇mem = @view q̇[bodyid2sys_full_coords[bodyid]]
            q̈mem = @view q̈[bodyid2sys_full_coords[bodyid]]
            Fmem = @view F[bodyid2sys_full_coords[bodyid]]
            λmem = @view λ[bodyid2sys_intrinsic_cstr_idx[bodyid]]
            cmem = @view c[bodyid2sys_loci_coords_idx[bodyid]]
            pmem = zero(p[bodyid2sys_full_coords[bodyid]])
            p̌mem = zero(p[bodyid2sys_full_coords[bodyid]])
            NonminimalCoordinatesState(
                t,
                qmem,q̇mem,q̈mem,Fmem,pmem,p̌mem,λmem,
                free_idx_by_mem[bodyid],
                pres_idx_by_mem[bodyid],
                cmem
            )
        end
        for bodyid = 1:nb
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
    num_of_dof::Int
    num_of_cstr::Int
    num_of_bodies::Int
    num_of_tensiles::Int
    bodies::BodyType
    tensiles::TenType
    connectivity::CntType
    state::StateType
end

"""
Structure Constructor.
$(TYPEDSIGNATURES)
"""
function Structure(bodies,tensiles,cnt::Connectivity)
    num_of_dim = get_num_of_dims(bodies)
    num_of_bodies = length(bodies)
    num_of_tensiles = sum(map(length,tensiles))
    (;num_of_free_coords,num_of_intrinsic_cstr) = cnt.indexed
    (;num_of_extrinsic_cstr) = cnt.jointed
    num_of_cstr = num_of_intrinsic_cstr + num_of_extrinsic_cstr
    num_of_dof = num_of_free_coords - num_of_cstr
    if num_of_dof <= 0
        @warn "Non positive degree of freedom: $num_of_dof."
    end
    state = StructureState(bodies,tensiles,cnt)
    st = Structure(
        num_of_dim,num_of_dof,num_of_cstr,
        num_of_bodies,num_of_tensiles,
        bodies,tensiles,
        cnt,
        state
    )
    # check_jacobian_singularity(st)
    st
end

function Structure(bodies,cnt::Connectivity)
    tensiles = (cables = Int[],)
    Structure(bodies,tensiles,cnt)
end


function update!(st::AbstractStructure; gravity=false)
    clear_forces!(st)
    stretch_rigids!(st)
    update_rigids!(st)
    update_tensiles!(st)
    if gravity
        apply_gravity!(st)
    end
    generate_forces!(st)
end

function update!(st::Structure,q::AbstractVector,q̇::AbstractVector=zero(q))
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    update!(st)
end

function build_M(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    M = spzeros(T,num_of_full_coords,num_of_full_coords)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        M[memfull,memfull] .+= body.state.cache.M
    end
    # @assert issymmetric(M)
    M
end

function build_M⁻¹(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    M⁻¹ = spzeros(T,num_of_full_coords,num_of_full_coords)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        M⁻¹[memfull,memfull] .+= body.state.cache.M⁻¹
    end
    # @assert issymmetric(M⁻¹)
    M⁻¹
end

function build_M̌(st::AbstractStructure)
    (;sys_free_coords_idx) = st.connectivity.indexed
    M = build_M(st)
    M̌ = Symmetric(M[sys_free_coords_idx,sys_free_coords_idx])
end

function build_∂Mq̇∂q(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂Mq̇∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂Mq̇∂q[memfull,memfull] .+= body.state.cache.∂Mq̇∂q
    end
    ∂Mq̇∂q
    # symsparsecsr(M;symmetrize=true)
end

function build_∂M⁻¹p∂q(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂M⁻¹p∂q = spzeros(T,num_of_full_coords,num_of_full_coords)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂M⁻¹p∂q[memfull,memfull] .+= body.state.cache.∂M⁻¹p∂q
    end
    ∂M⁻¹p∂q
    # symsparsecsr(M;symmetrize=true)
end

function make_M!(st)
    function inner_M!(M,q)
        update_rigids!(st,q)
        M .= build_M(st)
    end
end

function make_M⁻¹!(st)
    function inner_M⁻¹!(M⁻¹,q)
        update_rigids!(st,q)
        M⁻¹ .= build_M⁻¹(st)
    end
end

function make_Jac_M!(st)
    function Jac_M!(∂Mq̇∂q,q,q̇)
        update_rigids!(st,q,q̇)
        ∂Mq̇∂q .= build_∂Mq̇∂q(st)
    end
end

function make_Jac_M⁻¹!(st)
    function Jac_M⁻¹!(∂M⁻¹p∂q,q,q̇)
        update_rigids!(st,q,q̇)
        ∂M⁻¹p∂q .= build_∂M⁻¹p∂q(st)
    end
end

function build_∂T∂qᵀ(st::AbstractStructure)
    (;num_of_full_coords,bodyid2sys_full_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂T∂qᵀ = zeros(T,num_of_full_coords)
    foreach(st.bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        ∂T∂qᵀ[memfull] .+= body.state.cache.∂T∂qᵀ
    end
    ∂T∂qᵀ
end

"""
Return System mass matrices
$(TYPEDSIGNATURES)
"""
function build_mass_matrices(structure::Structure,)
    (;sys_free_idx,sys_pres_idx) = structure.connectivity.indexed
    M = assemble_M(structure,) |> sparse |> Symmetric
    M⁻¹ = M |> Matrix |> inv |> sparse |> Symmetric
    M̌   = M[sys_free_idx,sys_free_idx] |> sparse |> Symmetric
    M̌⁻¹ = M̌ |> Matrix |> inv |> sparse |> Symmetric
    Ḿ = M[sys_free_idx,:]            |> sparse
    M̄ = M[sys_free_idx,sys_pres_idx] |> sparse
    @eponymtuple(M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)
end

function make_cstr_function(st::AbstractStructure,q0::AbstractVector)
    (;bodies,num_of_cstr) = st
    (;numbered,indexed,jointed) = st.connectivity
    (;num_of_free_coords,num_of_full_coords,sys_pres_idx,sys_free_idx,bodyid2sys_full_coords,bodyid2sys_pres_coords,bodyid2sys_free_coords,num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx) = indexed

    function _inner_cstr_function(q̌,d,c)
        q = Vector{eltype(q̌)}(undef,num_of_full_coords)
        q[sys_pres_idx] .= q0[sys_pres_idx]
        q[sys_free_idx] .= q̌
        ret = Vector{eltype(q̌)}(undef,num_of_cstr)
        foreach(bodies) do body
            bodyid = body.prop.id
            memfull = bodyid2sys_full_coords[bodyid]
            memfree = bodyid2sys_free_coords[bodyid]
            memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
            if !isempty(memincst)
                ret[memincst] .= make_cstr_function(
                    body.state.cache.funcs
                )(q[memfull],d[memincst])
            end
        end
        is = Ref(num_of_intrinsic_cstr)
        foreach(jointed.joints) do joint
            nc = joint.num_of_cstr
            ret[is[]+1:is[]+nc] .= make_cstr_function(joint,st)(q,d[is[]+1:is[]+nc],c)
            is[] += nc
        end
        ret
    end

    function inner_cstr_function(q̌)
        _inner_cstr_function(q̌,get_d(st),get_local_coords(st))
    end
    function inner_cstr_function(q̌,d)
        _inner_cstr_function(q̌,d,get_local_coords(st))
    end
    function inner_cstr_function(q̌,d,c)
        _inner_cstr_function(q̌,d,c)
    end

    inner_cstr_function
end

function make_cstr_function(st::AbstractStructure)
    (;bodies,num_of_cstr) = st
    (;indexed,jointed,numbered) = st.connectivity
    (;num_of_free_coords,bodyid2sys_full_coords,bodyid2sys_free_coords,num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx) = indexed
    @inline @inbounds function inner_cstr_function(q)
        ret = Vector{eltype(q)}(undef,num_of_cstr)
        is = Ref(num_of_intrinsic_cstr)
        foreach(bodies) do body
            bodyid = body.prop.id
            memfull = bodyid2sys_full_coords[bodyid]
            memfree = bodyid2sys_free_coords[bodyid]
            memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
            if !isempty(memincst)
                ret[memincst] .= make_cstr_function(
                    body.state.cache.funcs
                )(q[memfull])
            end
        end
        foreach(jointed.joints) do joint
            nc = joint.num_of_cstr
            ret[is[]+1:is[]+nc] .= make_cstr_function(joint,st)(q)
            is[] += nc
        end
        ret
    end
    inner_cstr_function
end

function make_cstr_jacobian(st::AbstractStructure,q0::AbstractVector)
    (;bodies,num_of_cstr) = st
    (;numbered,indexed,jointed) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_pres_idx,sys_free_idx,bodyid2sys_full_coords,bodyid2sys_free_coords,num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx) = indexed

    function _inner_cstr_jacobian(q̌,c)
        q = Vector{eltype(q̌)}(undef,num_of_full_coords)
        q[sys_pres_idx] .= q0[sys_pres_idx]
        q[sys_free_idx] .= q̌
        ret = zeros(eltype(q̌),num_of_cstr,num_of_free_coords)
        foreach(bodies) do body
            bodyid = body.prop.id
            memfull = bodyid2sys_full_coords[bodyid]
            memfree = bodyid2sys_free_coords[bodyid]
            memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
            if !isempty(memincst)
                ret[memincst,memfree] .= make_cstr_jacobian(
                    body.state.cache.funcs
                )(q[memfull])
            end
        end
        is = Ref(num_of_intrinsic_cstr)
        foreach(jointed.joints) do joint
            nc = joint.num_of_cstr
            ret[is[]+1:is[]+nc,:] .= make_cstr_jacobian(joint,st)(q,c)
            is[] += nc
        end
        ret
    end
    function inner_cstr_jacobian(q̌)
        _inner_cstr_jacobian(q̌,get_local_coords(st))
    end
    function inner_cstr_jacobian(q̌,c)
        _inner_cstr_jacobian(q̌,c)
    end
    inner_cstr_jacobian
end

function make_cstr_jacobian(st::AbstractStructure)
    (;bodies,num_of_cstr) = st
    (;numbered,indexed,jointed) = st.connectivity
    (;num_of_free_coords,bodyid2sys_full_coords,bodyid2sys_free_coords,num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx) = indexed
    @inline @inbounds function inner_cstr_jacobian(q)
        ret = zeros(eltype(q),num_of_cstr,num_of_free_coords)
        is = Ref(num_of_intrinsic_cstr)
        foreach(bodies) do body
            bodyid = body.prop.id
            memfull = bodyid2sys_full_coords[bodyid]
            memfree = bodyid2sys_free_coords[bodyid]
            memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
            if !isempty(memincst)
                ret[memincst,memfree] .= make_cstr_jacobian(
                    body.state.cache.funcs
                )(q[memfull])
            end
        end
        foreach(jointed.joints) do joint
            nc = joint.num_of_cstr
            ret[is[]+1:is[]+nc,:] .= make_cstr_jacobian(joint,st)(q)
            is[] += nc
        end
        ret
    end
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


"""
检查雅可比矩阵奇异性
$(TYPEDSIGNATURES)
"""
function check_jacobian_singularity(st)
    (;bodies,state) = st
    q = get_coords(st)
    A = make_cstr_jacobian(st)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    foreach(bodies) do body
        bodyid = body.prop.id
        free_idx = body.state.cache.free_idx
        q_rb = state.members[bodyid].q
        Aq_rb = body.state.cache.funcs.Φq(q_rb)
        rb_rank = rank(Aq_rb)
        if rb_rank < minimum(size(Aq_rb))
            @warn "The $(bodyid)th rigid body's Jacobian is singular: rank(A(q))=$(rb_rank)<$(minimum(size(Aq_rb)))"
        end
    end
end

function analyse_slack(st::AbstractStructure,verbose=false)
    (;cables) = st.tensiles
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
    M = build_M(st)
    (;q̇) = st.state.system
    T = 1/2*transpose(q̇)*M*q̇
end

"""
Return System potential energy
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(st::Structure)
    V = Ref(zero(get_numbertype(st)))
    foreach(st.bodies) do body
        V[] += potential_energy_gravity(body)
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
    V = zero(get_numbertype(st))
    V += potential_energy_strain(st)
    if gravity
        V += potential_energy_gravity(st)
    end
    if !isempty(st.tensiles.cables)
        V += sum(potential_energy.(st.tensiles.cables))
    end
    V
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
