
struct ID
    rbid::Int
    apid::Int
end

struct Connectivity{BPConnectType,StringConnectType}
    body2q::BPConnectType
    string2ap::StringConnectType
end


struct TensegrityStructure{BodyType,StringType,ActuatorType,ConnectType}
    ndim::Int
    nbody::Int
    nmovablebody::Int
    mvbodyindex::Vector{Int}
    nfixbody::Int
    fixbodyindex::Vector{Int}
    nstring::Int
    npoints::Int
    rigidbodies::Vector{BodyType}
    strings::Vector{StringType}
    actuators::Vector{ActuatorType}
    connectivity::ConnectType
end

function TensegrityStructure(rbs::Vector{rbType},ss,acs,cnt
                ) where rbType<:AbstractRigidBody{N,T,CType} where {N,T,CType}
    ndim = N
    nbody = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbody = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbody = length(fixbodyindex)
    nstring = length(ss)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.naps
    end
    TensegrityStructure(ndim,nbody,nmvbody,mvbodyindex,
                    nfixbody,fixbodyindex,
                    nstring,npoints,
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
        l̇ = 1/l*(Δr[1]*Δṙ[1] + Δr[2]*Δṙ[2])
        f_raw = k*(l - sstate.restlen) +
                c*l̇
        sstate.tension = ifelse(f_raw > 0.0, f_raw, 0.0)
        f = τ*sstate.tension
        f1 .+= f
        f2 .+= -f
    end
end

function distribute_q_to_rbs!(tgstruct,globalq,globalq̇)
    rbs = tgstruct.rigidbodies
    cnt = tgstruct.connectivity
    for rbid in tgstruct.mvbodyindex
        pindex = cnt.body2q[rbid]
        @unpack q, q̇ = rbs[rbid].state.coords
        q .= globalq[pindex]
        q̇ .= globalq̇[pindex]
        @unpack cache,p = rbs[rbid].state
        for (i,ap) in enumerate(p)
            ap .= cache.Cp[i]*q
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

function assemble_forces!(F,tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    generate_forces!(rbs)
    F .= 0.0
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        F[pindex] .+= rbs[rbid].state.coords.Q
    end
end

function apply_gravity!(tgstruct)
    rbs = tgstruct.rigidbodies
    gravity_force = get_gravity(tgstruct)
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        rb.state.F .= gravity_force
    end
end
function kinetic_energy_coords(rb::RigidBody)
    @unpack q̇ = rb.state.coords
    @unpack M = rb.state.cache
    ke = 1/2*transpose(q̇)*M*q̇
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

function potential_energy(tgstruct::TensegrityStructure,q,q̇)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    update_strings_apply_forces!(tgstruct)
    pe = sum(potential_energy.(tgstruct.strings))
end

function energy(tgstruct,q,q̇)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    ke = sum(kinetic_energy_coords.(tgstruct.rigidbodies))
    update_strings_apply_forces!(tgstruct)
    pe = sum(potential_energy.(tgstruct.strings))
    ke + pe
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
    nq = body2q[end][end]
    mass_matrix = zeros(nq,nq)
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= tgstruct.rigidbodies[rbid].state.cache.M
    end
    mass_matrix
end

function get_nconstraint(tgstruct)
    @unpack nmovablebody,nfixbody = tgstruct
    nbodyconstraint = get_nbodyconstraint(tgstruct)
    nbodydof = get_nbodydof(tgstruct)
    ninconstraint = nbodyconstraint*nmovablebody
    nexconstraint = 0  #nbodydof*nfixbody
    for rb in tgstruct.rigidbodies
        nexconstraint += rb.state.cache.nc
    end
    nconstraint = ninconstraint + nexconstraint
end

function build_Φ(tgstruct,q0)
    rbs = tgstruct.rigidbodies
    #q0,q̇0 = get_q(tgstruct)
    @unpack body2q = tgstruct.connectivity
    nfixbody = tgstruct.nfixbody
    nconstraint = get_nconstraint(tgstruct)
    nbodyc = get_nbodyconstraint(tgstruct)
    nbodydof = get_nbodydof(tgstruct)
    fixindex = get_fixindex(tgstruct)
    function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = 0
        #is += nbodydof*nfixbody
        for rbid in tgstruct.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            ret[is+1:nc] = rb.state.cache.cfuncs.Φ(q[pindex])
            is += nc
            ret[is+1:is+nbodyc] .=rb.state.cache.funcs.Φ(q[pindex])
            is += nbodyc
        end
        # for (fixid,rbid) in enumerate(tgstruct.fixbodyindex)
        #     pindex = body2q[rbid]
        #     ret[nbodydof*(fixid-1)+1:nbodydof*fixid] .= q[pindex[fixindex]] - q0[pindex[fixindex]]
        # end
        ret
    end
end

function build_A(tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    nfixbody = tgstruct.nfixbody
    nconstraint = get_nconstraint(tgstruct)
    nbodyc = get_nbodyconstraint(tgstruct)
    nbodydof = get_nbodydof(tgstruct)
    fixindex = get_fixindex(tgstruct)
    fixA = get_fixA(tgstruct)
    nq = body2q[end][end]
    function inner_A(q)
        ret = zeros(eltype(q),nconstraint,nq)
        is = 0
        #is += nbodydof*nfixbody
        for rbid in tgstruct.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            ret[is+1:nc,pindex] = rb.state.cache.cfuncs.Φq(q[pindex])
            is += nc
            ret[is+1:is+nbodyc,pindex] .= rb.state.cache.funcs.Φq(q[pindex])
            is += nbodyc
        end
        # for (fixid,rbid) in enumerate(tgstruct.fixbodyindex)
        #     pindex = body2q[rbid]
        #     ret[nbodydof*(fixid-1)+1:nbodydof*fixid,pindex] .= fixA
        # end
        ret
    end
end

function get_q(tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    nq = body2q[end][end]
    q = zeros(nq)
    q̇ = zeros(nq)
    for rbid in tgstruct.mvbodyindex
        pindex = body2q[rbid]
        q[pindex] .= rbs[rbid].state.coords.q
        q̇[pindex] .= rbs[rbid].state.coords.q̇
    end
    return q,q̇
end

function get_initial(tgstruct)
    q0,q̇0 = get_q(tgstruct)
    λ0 = zeros(get_nconstraint(tgstruct))
    q0,q̇0,λ0
end

function get_u(tgstruct)
end

get_nbodyconstraint(tg::TensegrityStructure) = get_nbodyconstraint(tg.rigidbodies[1])
get_nbodyconstraint(rb::AbstractRigidBody{3,T,CType}) where {T,CType} = 6
get_nbodyconstraint(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = 1

get_nbodydof(tg::TensegrityStructure) = get_nbodydof(tg.rigidbodies[1])
get_nbodydof(rb::AbstractRigidBody{3,T,CType}) where {T,CType} = 6
get_nbodydof(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = 3

get_nbodycoords(tg::TensegrityStructure) = get_nbodycoords(tg.rigidbodies[1])
get_nbodycoords(rb::AbstractRigidBody{3,T,CType}) where {T,CType} = 12
get_nbodycoords(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = 4

get_fixindex(tg::TensegrityStructure) = get_fixindex(tg.rigidbodies[1])
get_fixindex(rb::AbstractRigidBody{3,T,CType}) where {T,CType} = 2:2:12
get_fixindex(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = [2,3,4]

get_fixA(tg::TensegrityStructure) = get_fixA(tg.rigidbodies[1])
function get_fixA(rb::AbstractRigidBody{3,T,CType}) where {T,CType}
    A = zeros(6,12)
    for (i,j) in zip(1:6,2:2:12)
        A[i,j] = 1
    end
    A
end
get_fixA(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = [0.0 1.0 0.0 0.0;
                                                              0.0 0.0 1.0 0.0;
                                                              0.0 0.0 0.0 1.0]
get_gravity(tg::TensegrityStructure) = get_gravity(tg.rigidbodies[1])
get_gravity(::AbstractRigidBody{3,T,CType}) where {T,CType} = [zero(T),zero(T),-9.8*one(T)]
get_gravity(::AbstractRigidBody{2,T,CType}) where {T,CType} = [zero(T),-9.8*one(T)]


function get_strings_len(tg::TensegrityStructure,q)
    distribute_q_to_rbs!(tg,q,zero(q))
    update_strings_apply_forces!(tg)
    get_strings_len(tg)
end
function get_strings_len(tg::TensegrityStructure)
    [s.state.length for s in tg.strings]
end

function find_remaining_index(body2q,rbs)
    original_nq = body2q[end][end]
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

function filter_body2q(body2q,rbs)
    original_nq = body2q[end][end]
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
