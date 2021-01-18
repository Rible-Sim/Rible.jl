
struct ID
    rbid::Int
    apid::Int
end

struct Connectivity{bType,sType,cType}
    body2q::bType
    string2ap::sType
    contacts::cType
end

Connectivity(b) = Connectivity(b,nothing,nothing)
Connectivity(b,s) = Connectivity(b,s,nothing)

struct TensegrityStructure{BodyType,StringType,ActuatorType,ConnectType}
    ndim::Int
    ncoords::Int
    nconstraint::Int
    ndof::Int
    nbodies::Int
    nmvbodies::Int
    mvbodyindex::Vector{Int}
    nfixbodies::Int
    fixbodyindex::Vector{Int}
    nstrings::Int
    npoints::Int
    nmvpoints::Int
    rigidbodies::Vector{BodyType}
    strings::Vector{StringType}
    actuators::Vector{ActuatorType}
    connectivity::ConnectType
end

function TensegrityStructure(rbs::Vector{rbT},ss::Vector{sT},cnt::Connectivity
                ) where {sT<:SString,rbT<:AbstractRigidBody{N,T,CType}} where {N,T,CType}
    acs = [Actuator(SVector{1}(s)) for s in ss]
    TensegrityStructure(rbs,ss,acs,cnt)
end

function TensegrityStructure(rbs::Vector{rbT},ss::Vector{sT},acs::Vector{aT},cnt::Connectivity
                ) where {sT<:SString,aT<:Actuator,rbT<:AbstractRigidBody{N,T,CType}} where {N,T,CType}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbodies = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbodies = length(fixbodyindex)
    nstrings = length(ss)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.naps
    end
    nmvpoints = 0
    for rbid in mvbodyindex
        nmvpoints += rbs[rbid].prop.naps
    end
    ncoords = maximum(maximum.(cnt.body2q))
    nconstraint = get_nconstraint(rbs,mvbodyindex,nmvbodies,nfixbodies)
    ndof = ncoords - nconstraint
    TensegrityStructure(ndim,
                    ncoords,nconstraint,ndof,
                    nbodies,
                    nmvbodies,mvbodyindex,
                    nfixbodies,fixbodyindex,
                    nstrings,
                    npoints,nmvpoints,
                    rbs,ss,acs,cnt)
end


function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function reset_forces!(tgst::TensegrityStructure)
    reset_forces!.(tgst.rigidbodies)
end


function update_strings_apply_forces!(tgstruct)
    rbs = tgstruct.rigidbodies
    ss = tgstruct.strings
    cnt = tgstruct.connectivity
    for (istr,sstring) in enumerate(ss)
        @unpack k,c = sstring
        sstate = sstring.state
        a,b = cnt.string2ap[istr]
        state1 = rbs[a.rbid].state
        p1 = state1.p[a.apid]
        ṗ1 = state1.cache.Cp[a.apid]*state1.coords.q̇
        f1 = state1.Faps[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.p[b.apid]
        ṗ2 = state2.cache.Cp[b.apid]*state2.coords.q̇
        f2 = state2.Faps[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        sstate.direction = τ
        l̇ = (transpose(Δr)*Δṙ)/l
        Δl = l - sstate.restlen
        f = k*Δl + c*l̇
        if Δl < 0
            sstate.tension = 0.0
        elseif f < 0
            sstate.tension = 0.0
        else
            sstate.tension = f
        end
        𝐟 = τ*sstate.tension
        f1 .+=  𝐟
        f2 .+= -𝐟
    end
end
distribute_q_to_rbs!(tgstruct,globalq) = distribute_q_to_rbs!(tgstruct,globalq,zero(globalq))
function distribute_q_to_rbs!(tgstruct,globalq,globalq̇)
    rbs = tgstruct.rigidbodies
    cnt = tgstruct.connectivity
    for rbid in tgstruct.mvbodyindex
        pindex = cnt.body2q[rbid]
        @unpack q, q̇ = rbs[rbid].state.coords
        q .= globalq[pindex]
        q̇ .= globalq̇[pindex]
        @unpack cache,p,r = rbs[rbid].state
        @unpack CG,Cp = cache
        r .= CG*q
        for (i,ap) in enumerate(p)
            ap .= Cp[i]*q
        end
    end
end

function generate_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        @unpack Faps = state
        @unpack Cp,CG = state.cache
        @unpack Q = state.coords
        Q .= 0.0
        for (pid,f) in enumerate(Faps)
            Q .+= transpose(Cp[pid])*f
        end
        Q .+= transpose(CG)*state.F
    end
end

function assemble_forces!(F,tgstruct;factor=1.0)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    generate_forces!(rbs)
    F .= 0.0
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        F[pindex] .+= factor*rbs[rbid].state.coords.Q
    end
end

function assemble_forces(tgstruct;factor=1.0)
    T = get_numbertype(tgstruct)
    @unpack body2q = tgstruct.connectivity
    F = zeros(T,tgstruct.ncoords)
    assemble_forces!(F,tgstruct,factor=factor)
    F
end

function apply_gravity!(tgstruct)
    rbs = tgstruct.rigidbodies
    gravity_acceleration = get_gravity(tgstruct)
    for (rbid,rb) in enumerate(rbs)
        @unpack prop, state = rb
        rb.state.F .+= gravity_acceleration*prop.mass
    end
end

function kinetic_energy_coords(rb::RigidBody)
    @unpack q̇ = rb.state.coords
    @unpack M = rb.state.cache
    ke = 1/2*transpose(q̇)*M*q̇
end

function gravity_potential_energy(rb)
    q = rb.state.coords.q
    gravity_potential_energy(rb,q)
end

function gravity_potential_energy(rb::RigidBody,q)
    @unpack CG = rb.state.cache
    r = CG*q
    gravity_acceleration = get_gravity(rb)
    -transpose(r)*gravity_acceleration*rb.prop.mass
end

function potential_energy(s::SString)
    pe = 0.0
    @unpack k,state = s
    Δlen = s.state.length-s.state.restlen
    if Δlen > 0.0
        pe += 1/2*k*Δlen^2
    end
    pe
end

function kinetic_energy_coords(tgstruct::TensegrityStructure,q,q̇)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    ke = sum(kinetic_energy_coords.(tgstruct.rigidbodies))
end

function gravity_potential_energy(tgstruct::TensegrityStructure,q)
    distribute_q_to_rbs!(tgstruct,q)
    sum(gravity_potential_energy.(tgstruct.rigidbodies))
end

function elastic_potential_energy(tgstruct::TensegrityStructure)
    update_strings_apply_forces!(tgstruct)
    pe = sum(potential_energy.(tgstruct.strings))
end

function elastic_potential_energy(tgstruct::TensegrityStructure,q)
    distribute_q_to_rbs!(tgstruct,q)
    elastic_potential_energy(tgstruct)
end

function elastic_potential_energy(tgstruct::TensegrityStructure,q,a)
    actuate!(tgstruct,a)
    elastic_potential_energy(tgstruct,q)
end

function energy(tgstruct,q,q̇;gravity=false)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    ke = sum(kinetic_energy_coords.(tgstruct.rigidbodies))
    update_strings_apply_forces!(tgstruct)
    epe = sum(potential_energy.(tgstruct.strings))
    if gravity
        gpe = gravity_potential_energy(tgstruct,q)
    else
        gpe = 0
    end
    ke + epe + gpe
end

function build_body2q(rbs::Vector{rbType}) where rbType<:AbstractRigidBody{N,T,CType} where {N,T,CType}
    bps = Vector{Vector{T}}()
    bp_number = Vector{Int}()
    push!(bp_number,0)
    body2q = Vector{Vector{Int}}()
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        xi,yi,xj,yj = state.coords.q
        bp1 = [xi,yi]
        bp2 = [xj,yj]
        bp1_find = findfirst(x->x==bp1,bps)
        if bp1_find === nothing
            push!(bps,bp1)
            push!(bp_number,bp_number[end]+1)
            bp1_number = bp_number[end]
        else
            bp1_number = bp1_find
        end
        bp2_find = findfirst(x->x==bp2,bps)
        if bp2_find === nothing
            push!(bps,bp2)
            push!(bp_number,bp_number[end]+1)
            bp2_number = bp_number[end]
        else
            bp2_number = bp2_find
        end
        push!(body2q,[2bp1_number-1,2bp1_number,
                      2bp2_number-1,2bp2_number])
    end
    body2q
end

function build_massmatrix(tgstruct::TensegrityStructure)
    body2q = tgstruct.connectivity.body2q
    ncoords = tgstruct.ncoords
    mass_matrix = zeros(ncoords,ncoords)
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= tgstruct.rigidbodies[rbid].state.cache.M
    end
    mass_matrix
end

function get_nconstraint(rbs,mvbodyindex,nmvbodies,nfixbodies)
    nbodyconstraint = get_nbodyconstraint(rbs)
    nbodydof = get_nbodydof(rbs)
    ninconstraint = nbodyconstraint*nmvbodies
    nexconstraint = 0  #nbodydof*nfixbodies
    for id in mvbodyindex
        nexconstraint += rbs[id].state.cache.nc
    end
    nconstraint = ninconstraint + nexconstraint
end

function build_Φ(tgstruct,q0)
    rbs = tgstruct.rigidbodies
    #q0,q̇0 = get_q(tgstruct)
    @unpack body2q = tgstruct.connectivity
    nfixbodies = tgstruct.nfixbodies
    nconstraint = tgstruct.nconstraint
    nbodyc = get_nbodyconstraint(tgstruct)
    nbodydof = get_nbodydof(tgstruct)
    @inline @inbounds function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = 0
        #is += nbodydof*nfixbodies
        for rbid in tgstruct.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            if nc > 0
                ret[is+1:is+nc] = rb.state.cache.cfuncs.Φ(q[pindex])
                is += nc
            end
            ret[is+1:is+nbodyc] .=rb.state.cache.funcs.Φ(q[pindex])
            is += nbodyc
        end
        ret
    end
end

function build_A(tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    nfixbodies = tgstruct.nfixbodies
    nconstraint = tgstruct.nconstraint
    nbodyc = get_nbodyconstraint(tgstruct)
    nbodydof = get_nbodydof(tgstruct)
    ncoords = tgstruct.ncoords
    @inline @inbounds function inner_A(q)
        ret = zeros(eltype(q),nconstraint,ncoords)
        is = 0
        for rbid in tgstruct.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            if nc > 0
                ret[is+1:is+nc,pindex] = rb.state.cache.cfuncs.Φq(q[pindex])
                is += nc
            end
            ret[is+1:is+nbodyc,pindex] .= rb.state.cache.funcs.Φq(q[pindex])
            is += nbodyc
        end
        ret
    end
end

function get_q(tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    ncoords = tgstruct.ncoords
    q = zeros(ncoords)
    q̇ = zeros(ncoords)
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        q[pindex] .= rbs[rbid].state.coords.q
        q̇[pindex] .= rbs[rbid].state.coords.q̇
    end
    return q,q̇
end

get_λ(tgstruct) = zeros(tgstruct.nconstraint)
function get_initial(tgstruct)
    q0,q̇0 = get_q(tgstruct)
    λ0 = get_λ(tgstruct)
    q0,q̇0,λ0
end

function get_u(tgstruct)
end

get_ndim(tg::TensegrityStructure) = get_ndim(tg.rigidbodies)
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = N

get_numbertype(tg::TensegrityStructure) = get_numbertype(tg.rigidbodies)
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = T

get_nbodyconstraint(tg::TensegrityStructure) = get_nbodyconstraint(tg.rigidbodies)
get_nbodyconstraint(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodyconstraint(eltype(rbs))
get_nbodyconstraint(rb::AbstractRigidBody) = get_nbodyconstraint(typeof(rb))
get_nbodyconstraint(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = get_nbodyconstraint(lncsType)
get_nbodyconstraint(::Type{<:NaturalCoordinates.LocalNaturalCoordinates2D4C}) = 1
get_nbodyconstraint(::Type{<:NaturalCoordinates.LocalNaturalCoordinates2D6C}) = 3
get_nbodyconstraint(::Type{<:NaturalCoordinates.LocalNaturalCoordinates3D12C}) = 6

get_nbodycoords(tg::TensegrityStructure) = get_nbodycoords(tg.rigidbodies)
get_nbodycoords(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodycoords(eltype(rbs))
get_nbodycoords(rb::AbstractRigidBody) = get_nbodycoords(typeof(rb))
get_nbodycoords(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = get_nbodycoords(lncsType)
get_nbodycoords(::Type{<:NaturalCoordinates.LocalNaturalCoordinates2D4C}) = 4
get_nbodycoords(::Type{<:NaturalCoordinates.LocalNaturalCoordinates2D6C}) = 6
get_nbodycoords(::Type{<:NaturalCoordinates.LocalNaturalCoordinates3D12C}) = 12

get_nbodydof(tg::TensegrityStructure) = get_nbodydof(tg.rigidbodies)
get_nbodydof(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodydof(eltype(rbs))
get_nbodydof(rb::AbstractRigidBody) = get_nbodydof(typeof(rb))
get_nbodydof(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = 3
get_nbodydof(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = 6

get_gravity(tg::TensegrityStructure) = get_gravity(tg.rigidbodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rb::AbstractRigidBody) = get_gravity(typeof(rb))
get_gravity(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),-9.81*one(T)]
get_gravity(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),zero(T),-9.81*one(T)]

function get_strings_len(tg::TensegrityStructure,q)
    distribute_q_to_rbs!(tg,q,zero(q))
    update_strings_apply_forces!(tg)
    get_strings_len(tg)
end
function get_strings_len(tg::TensegrityStructure)
    [s.state.length for s in tg.strings]
end

function get_original_restlen(tg::TensegrityStructure)
    [s.original_restlen for s in tg.strings]
end

function find_remaining_index(body2q,rbs)
    original_nq = maximum(maximum.(body2q))
    switch_index = zeros(Int,original_nq)
    for (rbid,rb) in enumerate(rbs)
        qindex = body2q[rbid]
        if rb.prop.movable
            for i in qindex
                switch_index[i] = i
            end
        end
    end
    remaining_index = findall((x)->x!=0,switch_index)
end

function filter_body2q(rbs)
    body2q_raw = build_body2q(rbs)
    body2q = filter_body2q(body2q_raw,rbs)
end

function filter_body2q(body2q,rbs)
    original_nq = maximum(maximum.(body2q))
    remaining_index = find_remaining_index(body2q,rbs)
    qpointer = collect(1:original_nq)[remaining_index]
    filtered_body2q = Vector{Vector{Int}}()
    for (rbid,rb) in enumerate(rbs)
        qindex = body2q[rbid]
        filtered_index = zero(qindex)
        if rb.prop.movable
            for (j,i) in enumerate(qindex)
                filtered_index[j] = findfirst((x)->x==i,qpointer)
            end
        end
        push!(filtered_body2q,filtered_index)
    end
    filtered_body2q
end
