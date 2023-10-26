abstract type AbstractStructure end


"""
Rigid Body Natural Coordinates State Type 。
$(TYPEDEF)
"""
mutable struct NonminimalCoordinatesState{T,qT,qviewT}
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
    c::qT
end

"""
Natural Coordinates State Constructor 。
$(TYPEDSIGNATURES)
"""
function NonminimalCoordinatesState(t,q,q̇,q̈,F,λ,freei,presi,c)
    q̌ = @view q[freei]
    q̌̇ = @view q̇[freei]
    q̌̈ = @view q̈[freei]
    q̃ = @view q[presi]
    q̃̇ = @view q̇[presi]
    q̃̈ = @view q̈[presi]
    F̌ = @view F[freei]
    NonminimalCoordinatesState(t,q,q̇,q̈,F,λ,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌,c)
end


"""
State Type 。
$(TYPEDEF)
"""
struct StructureState{sysT, msT}
    system::sysT
    parts::msT
end

"""
State Constructor 。
$(TYPEDSIGNATURES)
"""
function StructureState(bodies,tensiles,cnt::Connectivity{<:Any,<:Any,<:NamedTuple{(:connected, )},<:Any})
    (;numbered,indexed,jointed) = cnt
    (;mem2sys) = numbered
    (;nfull,ninconstraints,sysfree,syspres) = indexed
    (;mem2sysincst,mem2sysfull) = indexed
    (;nexconstraints) = jointed
    nconstraints = ninconstraints + nexconstraints
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
    q = Vector{T}(undef,nfull)
    q̇ = zero(q)
    q̈ = zero(q)
    F = zero(q)
    λ = Vector{T}(undef,nconstraints)
    c = get_c(bodies,numbered)
    system = NonminimalCoordinatesState(t,q,q̇,q̈,F,λ,sysfree,syspres,c)
    parts = [
        begin
            qmem = @view q[mem2sysfull[bodyid]]
            q̇mem = @view q̇[mem2sysfull[bodyid]]
            q̈mem = @view q̈[mem2sysfull[bodyid]]
            Fmem = @view F[mem2sysfull[bodyid]]
            λmem = @view λ[mem2sysincst[bodyid]]
            cmem = @view c[mem2sys[bodyid]]
            NonminimalCoordinatesState(
                t,qmem,q̇mem,q̈mem,Fmem,λmem,
                free_idx_by_mem[bodyid],
                pres_idx_by_mem[bodyid],
                cmem
            )
        end
        for bodyid = 1:nb
    ]
    foreach(bodies) do body
        q, q̇ = body2coordinates(body)
        parts[body.prop.id].q .= q
        parts[body.prop.id].q̇ .= q̇
    end
    StructureState(system,parts)
end


"""
Structure Type 。
$(TYPEDEF)
"""
struct Structure{BodyType,TenType,CntType,StateType} <: AbstractStructure
    ndim::Int
    ndof::Int
    nconstraints::Int
    nbodies::Int
    ntensiles::Int
    bodies::BodyType
    tensiles::TenType
    connectivity::CntType
    state::StateType
end

"""
Structure Constructor 。
$(TYPEDSIGNATURES)
"""
function Structure(bodies,tensiles,cnt::Connectivity)
    ndim = get_ndim(bodies)
    nbodies = length(bodies)
    ntensiles = sum(map(length,tensiles))
    (;nfree,ninconstraints) = cnt.indexed
    (;nexconstraints) = cnt.jointed
    nconstraints = ninconstraints + nexconstraints
    ndof = nfree - nconstraints
    if ndof <= 0
        @warn "Non positive degree of freedom: $ndof."
    end
    state = StructureState(bodies,tensiles,cnt)

    (;mem2num) = cnt.numbered
    contacts = reduce(
        vcat,
        map(sort(bodies)) do body
            (;prop,) = body
            (;μs,es) = prop
            rbid = prop.id
            [
                Contact(id,μ,e)
                for (id,μ,e) in  zip(mem2num[rbid],μs,es)
            ]
        end
    )
    st = Structure(
        ndim,ndof,nconstraints,
        nbodies,ntensiles,
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
    # update_clustercables_apply_forces!(st)
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
    (;nfull,mem2sysfull) = st.connectivity.indexed
    T = get_numbertype(st)
    M = spzeros(T,nfull,nfull)
    foreach(st.bodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        M[memfull,memfull] .+= rb.state.cache.M
    end
    # @assert issymmetric(M)
    M
end

function build_M⁻¹(st::AbstractStructure)
    (;nfull,mem2sysfull) = st.connectivity.indexed
    T = get_numbertype(st)
    M⁻¹ = spzeros(T,nfull,nfull)
    foreach(st.bodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        M⁻¹[memfull,memfull] .+= rb.state.cache.M⁻¹
    end
    # @assert issymmetric(M⁻¹)
    M⁻¹
end

function build_M̌(st::AbstractStructure)
    (;sysfree) = st.connectivity.indexed
    M = build_M(st)
    M̌ = Symmetric(M[sysfree,sysfree])
end

function build_∂Mq̇∂q(st::AbstractStructure)
    (;nfull,mem2sysfull) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂Mq̇∂q = spzeros(T,nfull,nfull)
    foreach(st.bodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        ∂Mq̇∂q[memfull,memfull] .+= rb.state.cache.∂Mq̇∂q
    end
    ∂Mq̇∂q
    # symsparsecsr(M;symmetrize=true)
end

function build_∂M⁻¹p∂q(st::AbstractStructure)
    (;nfull,mem2sysfull) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂M⁻¹p∂q = spzeros(T,nfull,nfull)
    foreach(st.bodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        ∂M⁻¹p∂q[memfull,memfull] .+= rb.state.cache.∂M⁻¹p∂q
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
    (;nfull,mem2sysfull) = st.connectivity.indexed
    T = get_numbertype(st)
    ∂T∂qᵀ = zeros(T,nfull)
    foreach(st.bodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        ∂T∂qᵀ[memfull] .+= rb.state.cache.∂T∂qᵀ
    end
    ∂T∂qᵀ
end


function make_Φ(st::AbstractStructure,q0::AbstractVector)
    (;bodies,nconstraints) = st
    (;numbered,indexed,jointed) = st.connectivity
    (;nfree,nfull,syspres,sysfree,mem2sysfull,mem2syspres,mem2sysfree,ninconstraints,mem2sysincst) = indexed

    function _inner_Φ(q̌,d,c)
        q = Vector{eltype(q̌)}(undef,nfull)
        q[syspres] .= q0[syspres]
        q[sysfree] .= q̌
        ret = Vector{eltype(q̌)}(undef,nconstraints)
        foreach(bodies) do rb
            rbid = rb.prop.id
            memfull = mem2sysfull[rbid]
            memfree = mem2sysfree[rbid]
            memincst = mem2sysincst[rbid]
            if !isempty(memincst)
                ret[memincst] .= rb.state.cache.funcs.Φ(q[memfull],d[memincst])
            end
        end
        is = Ref(ninconstraints)
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(joint,st)(q,d[is[]+1:is[]+nc],c)
            is[] += nc
        end
        ret
    end

    function inner_Φ(q̌)
        _inner_Φ(q̌,get_d(st),get_c(st))
    end
    function inner_Φ(q̌,d)
        _inner_Φ(q̌,d,get_c(st))
    end
    function inner_Φ(q̌,d,c)
        _inner_Φ(q̌,d,c)
    end

    inner_Φ
end

function make_Φ(st::AbstractStructure)
    (;bodies,nconstraints) = st
    (;indexed,jointed,numbered) = st.connectivity
    (;nfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed
    @inline @inbounds function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraints)
        is = Ref(ninconstraints)
        foreach(bodies) do rb
            rbid = rb.prop.id
            memfull = mem2sysfull[rbid]
            memfree = mem2sysfree[rbid]
            memincst = mem2sysincst[rbid]
            if !isempty(memincst)
                ret[memincst] .= rb.state.cache.funcs.Φ(q[memfull])
            end
        end
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(joint,st)(q)
            is[] += nc
        end
        ret
    end
    inner_Φ
end


function make_A(st::AbstractStructure,q0::AbstractVector)
    (;bodies,nconstraints) = st
    (;numbered,indexed,jointed) = st.connectivity
    (;nfull,nfree,syspres,sysfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed

    function _inner_A(q̌,c)
        q = Vector{eltype(q̌)}(undef,nfull)
        q[syspres] .= q0[syspres]
        q[sysfree] .= q̌
        ret = zeros(eltype(q̌),nconstraints,nfree)
        foreach(bodies) do rb
            rbid = rb.prop.id
            memfull = mem2sysfull[rbid]
            memfree = mem2sysfree[rbid]
            memincst = mem2sysincst[rbid]
            if !isempty(memincst)
                ret[memincst,memfree] .= rb.state.cache.funcs.Φq(q[memfull])
            end
        end
        is = Ref(ninconstraints)
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc,:] .= make_A(joint,st)(q,c)
            is[] += nc
        end
        ret
    end
    function inner_A(q̌)
        _inner_A(q̌,get_c(st))
    end
    function inner_A(q̌,c)
        _inner_A(q̌,c)
    end
    inner_A
end

function make_A(st::AbstractStructure)
    (;bodies,nconstraints) = st
    (;numbered,indexed,jointed) = st.connectivity
    (;nfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed
    @inline @inbounds function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        is = Ref(ninconstraints)
        foreach(bodies) do rb
            rbid = rb.prop.id
            memfull = mem2sysfull[rbid]
            memfree = mem2sysfree[rbid]
            memincst = mem2sysincst[rbid]
            if !isempty(memincst)
                ret[memincst,memfree] .= rb.state.cache.funcs.Φq(q[memfull])
            end
        end
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc,:] .= make_A(joint,st)(q)
            is[] += nc
        end
        ret
    end
end

function build_F̌(st,rbid,pid,f)
    T = get_numbertype(st)
    (;nfree,mem2sysfree) = st.connectivity.indexed
    F̌ = zeros(T,nfree)
    foreach(st.bodies) do rb
        if rb.prop.id == rbid
            C = rb.state.cache.Cps[pid]
            memfree = mem2sysfree[rbid]
            uci = rb.state.cache.free_idx
            F̌[memfree] = (transpose(C)*f)[uci,:]
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
    q = get_q(st)
    A = make_A(st)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    foreach(bodies) do rb
        rbid = rb.prop.id
        uci = rb.state.cache.free_idx
        q_rb = state.parts[rbid].q
        Aq_rb = rb.state.cache.funcs.Φq(q_rb)
        rb_rank = rank(Aq_rb)
        if rb_rank < minimum(size(Aq_rb))
            @warn "The $(rbid)th rigid body's Jacobian is singular: rank(A(q))=$(rb_rank)<$(minimum(size(Aq_rb)))"
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
    foreach(st.bodies) do rb
        V[] += potential_energy_gravity(rb)
    end
    V[]
end

function potential_energy_strain(st::Structure)
    V = Ref(zero(get_numbertype(st)))
    foreach(st.bodies) do rb
        V[] += potential_energy_strain(rb)
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
