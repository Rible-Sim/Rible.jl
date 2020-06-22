
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
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt.body2q[rbid]
        @unpack q, q̇ = rb.state.coords
        q .= globalq[pindex]
        q̇ .= globalq̇[pindex]
        @unpack cache,p = rb.state
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
    for (rbid,rb) in enumerate(rbs)
        pindex = body2q[rbid]
        F[pindex] .+= rb.state.coords.Q
    end
end

function apply_gravity!(tgstruct)
    rbs = tgstruct.rigidbodies
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        g = 9.8
        rb.state.F .= [0.0,-g]
    end
end
function kineticenergy(rbs)
    ke = 0.0
    for (rbid,rb) in enumerate(rbs)
        @unpack q̇ = rb.state.coords
        @unpack M = rb.state.cache
        ke += 1/2*transpose(q̇)*M*q̇
    end
    ke
end

function potentialenergy(ss)
    pe = 0.0
    for (ssid,sstring) in enumerate(ss)
        @unpack k,original_restlen = sstring
        sstate = sstring.state
        Δ1 = sstate.length-original_restlen
        if Δ1 > 0.0
            pe += 1/2*k*Δ1^2
        end
    end
    pe
end

function energy(q,q̇,tgstruct)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    ss = tgstruct.strings
    rbs = tgstruct.rigidbodies
    update_strings_apply_forces!(tgstruct)
    ke = kineticenergy(rbs)
    pe = potentialenergy(ss)
    ke + pe
end

function build_body2q(rbs::Vector{RigidBody{T,CT,AT}}) where {T,CT,AT}
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
    rbs = tgstruct.rigidbodies
    body2q = tgstruct.connectivity.body2q
    build_massmatrix(rbs,body2q)
end
function build_massmatrix(rbs,body2q)
    nq = body2q[end][end]
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= rbs[rbid].state.cache.M
    end
    mass_matrix
end

function get_nconstraint(tgstruct)
    @unpack nbody,nfixbody = tgstruct
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
end

function build_Φ(tgstruct)
    rbs = tgstruct.rigidbodies
    q0,q̇0 = get_q(tgstruct)
    @unpack body2q = tgstruct.connectivity
    nfixbody = tgstruct.nfixbody
    nconstraint = get_nconstraint(tgstruct)
    nbodyc = get_nbodyconstraint(tgstruct)
    nbodydof = get_nbodydof(tgstruct)
    fixindex = get_fixindex(tgstruct)
    function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        for (rbid,rb) in enumerate(rbs)
            pindex = body2q[rbid]
            is = nbodydof*nfixbody
            ret[is+nbodyc*(rbid-1)+1:is+nbodyc*rbid] .= rb.state.cache.funcs.Φ(q[pindex])
        end
        for (fixid,rbid) in enumerate(tgstruct.fixbodyindex)
            pindex = body2q[rbid]
            ret[nbodydof*(fixid-1)+1:nbodydof*fixid] .= q[pindex[fixindex]] - q0[pindex[fixindex]]
        end
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
        for (rbid,rb) in enumerate(rbs)
            pindex = body2q[rbid]
            is = nbodydof*nfixbody
            ret[is+nbodyc*(rbid-1)+1:is+nbodyc*rbid,pindex] .= rb.state.cache.funcs.Φq(q[pindex])
        end
        for (fixid,rbid) in enumerate(tgstruct.fixbodyindex)
            pindex = body2q[rbid]
            ret[nbodydof*(fixid-1)+1:nbodydof*fixid,pindex] .= fixA
        end
        ret
    end
end

function get_q(tgstruct)
    rbs = tgstruct.rigidbodies
    @unpack body2q = tgstruct.connectivity
    nq = body2q[end][end]
    q = zeros(nq)
    q̇ = zeros(nq)
    for (rbid,rb) in enumerate(rbs)
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
get_fixindex(rb::AbstractRigidBody{3,T,CType}) where {T,CType} = 1:6
get_fixindex(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = [1,2,4]

get_fixA(tg::TensegrityStructure) = get_fixA(tg.rigidbodies[1])
get_fixA(rb::AbstractRigidBody{3,T,CType}) where {T,CType} = hcat(Matrix(1.0I,6,6),zeros(6,6))
get_fixA(rb::AbstractRigidBody{2,T,CType}) where {T,CType} = [1.0 0.0 0.0 0.0;
                                                                0.0 1.0 0.0 0.0;
                                                                0.0 0.0 0.0 1.0]
