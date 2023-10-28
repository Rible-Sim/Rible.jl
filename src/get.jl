
get_c(bot::Robot) = get_c(bot.st)
get_s(bot::Robot) = get_s(bot.st)

get_q(st::Structure) = copy(st.state.system.q)
get_q̇(st::Structure) = copy(st.state.system.q̇)
get_q̌(st::Structure) = copy(st.state.system.q̌)
get_q̌̇(st::Structure) = copy(st.state.system.q̌̇)

function get_λ(st::Structure)
    st.state.system.λ
end

"""
Return System Initial State 。
$(TYPEDSIGNATURES)
"""
function get_initial(st::Structure)
    _,λ = check_static_equilibrium_output_multipliers(st::Structure)
    q̌ = get_q̌(st::Structure)
    q = get_q(st::Structure)
    ℓ = get_cables_len(st::Structure)
    s = 1 ./ℓ
    d = get_d(st::Structure)
    c = get_c(st::Structure)
    k = get_cables_stiffness(st::Structure)
    μ = get_cables_restlen(st::Structure)
    @eponymtuple(q̌,q,s,λ,d,c,k,μ,)
end

function get_polyvar(st::Structure)
    (;nconstraints,connectivity,tensiles) = st
    (;cables) = tensiles
    (;indexed,numbered,) = connectivity
    (;nc) = numbered
    (;ninconstraints) = indexed
    (;nfree) = indexed
    ncables = length(cables)
    # state variables
    @polyvar q̌[1:nfree]
    @polyvar s[1:ncables]
    @polyvar λ[1:nconstraints]
    # parameters
    @polyvar d[1:nconstraints]
    @polyvar c[1:nc]
    @polyvar k[1:ncables]
    @polyvar μ[1:ncables]
    @eponymtuple(q̌,s,λ,d,c,k,μ,)
end

function get_s(st::Structure)
    1 ./get_cables_len(st::Structure)
end

function get_c(st::Structure,bodyid,pid)
    (;mem2num,num2sys) = st.connectivity.numbered
    (;c) = st.state.system
    cidx = mem2num[bodyid][pid]
    c[num2sys[cidx]]
end

get_c(st::Structure) = copy(st.state.system.c)

function get_c(bodies,numbered::NumberedPoints)
    ndim = get_ndim(bodies)
    T = get_numbertype(bodies)
    (;mem2num,num2sys,nc) = numbered
    ret = zeros(T,nc)
    foreach(bodies) do body
        bodyid = body.prop.id
        for i in eachindex(body.prop.loci)
            ip = mem2num[bodyid][i]
            ret[num2sys[ip]] .= get_c(body,body.prop.loci[i].position)
        end
    end
    ret
end

function get_c(body::RigidBody,r̄p)
    body.state.cache.funcs.c(r̄p)
end

function get_c(body::FlexibleBody,r̄p)
    body.state.cache.funcs.x(r̄p)
end

function set_C!(st::Structure,c)
    T = get_numbertype(st::Structure)
    (;numbered,indexed) = st.connectivity
    (;mem2num,num2sys,nc) = numbered
    foreach(st.bodies) do body
        bodyid = body.prop.id
        for i in eachindex(body.prop.loci)
            ip = mem2num[bodyid][i]
            body.state.cache.Cps[i] = body.state.cache.funcs.C(c[num2sys[ip]])
        end
    end
end

function get_d(st::Structure)
    (;nconstraints,bodies) = st
    (;indexed,jointed) = st.connectivity
    (;mem2sysincst,ninconstraints) = indexed
    T = get_numbertype(st::Structure)
    d = Vector{T}(undef,nconstraints)
    is = Ref(ninconstraints)
    foreach(bodies) do body
        bodyid = body.prop.id
        memincst = mem2sysincst[bodyid]
        if !isempty(memincst)
            d[memincst] .= NCF.get_deform(body.state.cache.funcs.nmcs)
        end
    end
    foreach(jointed.joints) do joint
        nc = joint.nconstraints
        d[is[]+1:is[]+nc] .= joint.values
        is[] += nc
    end
    d
end

"""
Return System 维度。
$(TYPEDSIGNATURES)
"""
get_ndim(bot::Robot) = get_ndim(bot.st)
get_ndim(st::AbstractStructure) = get_ndim(st.bodies)
get_ndim(bodies::AbstractVector{<:AbstractBody}) = get_ndim(eltype(bodies))
get_ndim(bodies::TypeSortedCollection) = get_ndim(eltype(bodies.data[1]))
get_ndim(body::AbstractBody) = get_ndim(typeof(body))
get_ndim(::Type{<:AbstractBody{N,T}}) where {N,T} = N
get_ndim(::AbstractBodyProperty{N}) where {N} = N

get_numbertype(bot::Robot) = get_numbertype(bot.st)
get_numbertype(st::AbstractStructure) = get_numbertype(st.bodies)
get_numbertype(bodies::AbstractVector{<:AbstractBody}) = get_numbertype(eltype(bodies))
get_numbertype(bodies::TypeSortedCollection) = get_numbertype(eltype(bodies.data[1]))
get_numbertype(body::AbstractBody) = get_numbertype(typeof(body))
get_numbertype(::Type{<:AbstractBody{N,T}}) where {N,T} = T
get_numbertype(::AbstractBodyProperty{N,T}) where {N,T} = T

"""
Return System 约束数量。
$(TYPEDSIGNATURES)
"""
get_nconstraints(st::Structure) = st.nconstraints

function get_nconstraints(rbs::TypeSortedCollection)
    ninconstraints = mapreduce(get_ninconstraints,+,rbs,init=0)
end
get_ninconstraints(body::AbstractRigidBody) = get_nconstraints(body.state.cache.funcs.nmcs)
get_nbodycoords(body::AbstractRigidBody) = get_nbodycoords(body.state.cache.funcs.nmcs)
get_ndof(body::AbstractRigidBody) = get_ndof(body.state.cache.funcs.nmcs)
get_nlocaldim(body::AbstractRigidBody) = get_nlocaldim(body.state.cache)
get_nlocaldim(cache::NonminimalCoordinatesCache) = get_nlocaldim(cache.funcs.nmcs)

get_ninconstraints(fb::AbstractFlexibleBody) = get_nconstraints(fb.state.cache.funcs.ancs)
get_nbodycoords(fb::AbstractFlexibleBody) = get_nbodycoords(fb.state.cache.funcs.ancs)
get_ndof(body::AbstractFlexibleBody) = get_ndof(body.state.cache.funcs.ancs)
get_nlocaldim(fb::AbstractFlexibleBody) = get_nlocaldim(fb.state.cache)
get_nlocaldim(cache::FlexibleBodyCoordinatesCache) = get_nlocaldim(cache.funcs.ancs)

get_nconstraints(nmcs::NCF.LNC) = NCF.get_nconstraints(nmcs)
get_nbodycoords(nmcs::NCF.LNC) = NCF.get_ncoords(nmcs)
get_ndof(nmcs::NCF.LNC) = NCF.get_ndof(nmcs)
get_nlocaldim(nmcs::NCF.LNC) = NCF.get_nlocaldim(nmcs)

get_nconstraints(nmcs::QBF.QC) = QBF.get_nconstraints(nmcs)
get_nbodycoords(nmcs::QBF.QC) = QBF.get_ncoords(nmcs)
get_ndof(nmcs::QBF.QC) = QBF.get_ndof(nmcs)
get_nlocaldim(nmcs::QBF.QC) = QBF.get_nlocaldim(nmcs)


get_nconstraints(ancs::ANCF.ANC) = ANCF.get_nconstraints(ancs)
get_nbodycoords(ancs::ANCF.ANC) = ANCF.get_ncoords(ancs)
get_ndof(ancs::ANCF.ANC) = ANCF.get_ndof(ancs)
get_nlocaldim(ancs::ANCF.ANC) = ANCF.get_nlocaldim(ancs)


"""
Return System 重力。
$(TYPEDSIGNATURES)
"""
get_gravity(bot::Robot) = get_gravity(bot.st)
get_gravity(st::Structure) = get_gravity(st.bodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rbs::TypeSortedCollection) = get_gravity(eltype(rbs.data[1]))
get_gravity(body::AbstractBody) = get_gravity(typeof(body))
get_gravity(::Type{<:AbstractBody{2,T}}) where {T} = SVector(zero(T),        -9.81*one(T))
get_gravity(::Type{<:AbstractBody{3,T}}) where {T} = SVector(zero(T),zero(T),-9.81*one(T))
# get_gravity(::Type{<:AbstractBody{3,T}}) where {T} = SVector(zero(T),-9.81*one(T),zero(T))
# get_gravity(::Type{<:AbstractBody{3,T}}) where {T} = SVector(-9.81*one(T),zero(T),zero(T))

get_cables_len(bot::Robot) = get_cables_len(bot.st)
get_cables_deform(bot::Robot) = get_cables_deform(bot.st)
get_cables_restlen(bot::Robot) = get_cables_restlen(bot.st)
get_cables_len_dot(bot::Robot) = get_cables_len_dot(bot.st)
get_cables_tension(bot::Robot) = get_cables_tension(bot.st)
get_cables_stiffness(bot::Robot) = get_cables_stiffness(bot.st)
get_cables_force_density(bot::Robot) = get_cables_force_density(bot.st)

get_bodies(bot::Robot) = get_bodies(bot.st)
get_bodies(st::AbstractStructure) = sort(st.bodies)

get_rigidbars(bot::Robot) = get_rigidbars(bot.st)

function get_rigidbars(st::AbstractStructure)
    filter(get_bodies(st::Structure)) do body
        body.state.cache.funcs.nmcs isa Union{NCF.LNC2D4C,NCF.LNC3D6C}
    end
end

get_connected(bot::Robot) = get_connected(bot.st)
get_connected(st::AbstractStructure) = sort(st.connectivity.tensioned.connected)

function get_cables_len!(st::Structure,q)
    update_rigids!(st::Structure,q,zero(q))
    update_cables_apply_forces!(st::Structure)
    get_cables_len(st::Structure)
end

"""
Return System Cable 刚度。
$(TYPEDSIGNATURES)
"""
function get_cables_stiffness(st::Structure)
    [s.k for s in st.tensiles.cables]
end

"""
Return System Cable 当前Length 。
$(TYPEDSIGNATURES)
"""
function get_cables_len(st::Structure)
    [s.state.length for s in st.tensiles.cables]
end

function get_cables_len_dot(st::Structure)
    [s.state.lengthdot for s in st.tensiles.cables]
end

"""
Return System Cable 变形量。
$(TYPEDSIGNATURES)
"""
function get_cables_deform(st::Structure)
    [s.state.length - s.state.restlen for s in st.tensiles.cables]
end

"""
Return System Cable Restlength。
$(TYPEDSIGNATURES)
"""
function get_cables_restlen(st::Structure)
    [s.state.restlen for s in st.tensiles.cables]
end

"""
Return System Cable Tension 。
$(TYPEDSIGNATURES)
"""
function get_cables_tension(st::Structure)
    [s.state.tension for s in st.tensiles.cables]
end

"""
Set cables' tension
$(TYPEDSIGNATURES)
"""
function set_cables_tension!(st::Structure,fs)
    for (s,f) in zip(st.tensiles.cables,fs)
        s.state.tension = f
    end
end


"""
Return System Cable 力密度。
$(TYPEDSIGNATURES)
"""
function get_cables_force_density(st::Structure)
    [s.state.tension/s.state.length for s in st.tensiles.cables]
end

"""
Return System Cable Initial Length 。
$(TYPEDSIGNATURES)
"""
function get_original_restlen(botinput::Robot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_cables_restlen(bot.st)
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
        for (γ,s) in zip(γs,st.tensiles.cables)]
end
