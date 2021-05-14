
struct TensegrityRobot{tgT,hubT,trajT}
    tg::tgT
    hub::hubT
    traj::trajT
end

abstract type TensegrityRobotTrajectory{T} end

struct NaturalCoordinatesTrajectory{T} <: TensegrityRobotTrajectory{T}
    ts::Vector{T}
    qs::Vector{Vector{T}}
    q̇s::Vector{Vector{T}}
    λs::Vector{Vector{T}}
end

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

struct TensegrityStructure{BodyType,ST,TT,ConnectType}
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
    strings::ST
    tensiles::TT
    connectivity::ConnectType
end

function TensegrityStructure(rbs::Vector{rbT},tensiles::TT,cnt::Connectivity) where {TT,rbT<:AbstractRigidBody{N,T}} where {N,T}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbodies = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbodies = length(fixbodyindex)
    nstrings = sum(map(length,tensiles))
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
    strings = tensiles.strings
    TensegrityStructure(ndim,
                    ncoords,nconstraint,ndof,
                    nbodies,
                    nmvbodies,mvbodyindex,
                    nfixbodies,fixbodyindex,
                    nstrings,
                    npoints,nmvpoints,
                    rbs,strings,tensiles,cnt)
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
    for sstring in ss
        @unpack id,k,c = sstring
        sstate = sstring.state
        a,b = cnt.string2ap[id]
        state1 = rbs[a.rbid].state
        p1 = state1.rps[a.apid]
        ṗ1 = state1.ṙps[a.apid]
        f1 = state1.Faps[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.rps[b.apid]
        ṗ2 = state2.ṙps[b.apid]
        f2 = state2.Faps[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        sstate.direction = τ
        sstate.lengthdot = (transpose(Δr)*Δṙ)/l
        Δl = sstate.length - sstate.restlen
        f = k*Δl + c*sstate.lengthdot
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

function update_SMA_strings_apply_forces!(tgstruct)
    rbs = tgstruct.rigidbodies
    SMA_strings = tgstruct.tensiles.SMA_strings
    cnt = tgstruct.connectivity
    for SMA_string in SMA_strings
        @unpack id,law = SMA_string
        sstate = SMA_string.state
        a,b = cnt.string2ap[id]
        state1 = rbs[a.rbid].state
        p1 = state1.rps[a.apid]
        ṗ1 = state1.ṙps[a.apid]
        f1 = state1.Faps[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.rps[b.apid]
        ṗ2 = state2.ṙps[b.apid]
        f2 = state2.Faps[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        sstate.direction = τ
        sstate.lengthdot = (transpose(Δr)*Δṙ)/l
        Δl = sstate.length - sstate.restlen
        f = law(Δl)
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
        @unpack cache,rps,ṙps,ro,ṙo,rg,ṙg = rbs[rbid].state
        @unpack Co,Cg,Cp = cache
        mul!(ro, Co, q)
        mul!(ṙo, Co, q̇)
        mul!(rg, Cg, q)
        mul!(ṙg, Cg, q̇)
        for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
            mul!(rp, Cp[i], q)
            mul!(ṙp, Cp[i], q̇)
        end
    end
end

function generate_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        @unpack Faps = state
        @unpack Cp,Cg = state.cache
        @unpack Q = state.coords
        Q .= 0.0
        for (pid,f) in enumerate(Faps)
            Q .+= transpose(Cp[pid])*f
        end
        Q .+= transpose(Cg)*state.F
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

function apply_gravity!(tgstruct;factor=1)
    rbs = tgstruct.rigidbodies
    gravity_acceleration = factor*get_gravity(tgstruct)
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
    @unpack Cg = rb.state.cache
    r = Cg*q
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

potential_energy(rb::AbstractRigidBody) = gravity_potential_energy(rb)

function kinetic_energy_coords(tgstruct::TensegrityStructure,q,q̇)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    ke = sum(kinetic_energy_coords.(tgstruct.rigidbodies))
end

function gravity_potential_energy(tgstruct::TensegrityStructure,q)
    distribute_q_to_rbs!(tgstruct,q)
    sum(gravity_potential_energy.(tgstruct.rigidbodies))
end

function elastic_potential_energy(tgstruct::TensegrityStructure)
    reset_forces!(tgstruct)
    update_strings_apply_forces!(tgstruct)
    pe = sum(potential_energy.(tgstruct.strings))
end

function elastic_potential_energy(tgstruct::TensegrityStructure,q)
    distribute_q_to_rbs!(tgstruct,q)
    elastic_potential_energy(tgstruct)
end

function elastic_potential_energy(tr::TensegrityRobot,q,a)
    actuate!(tr,a)
    elastic_potential_energy(tr.tg,q)
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
    lncs = Vector{Vector{T}}()
    bp_number = Vector{Int}()
    push!(bp_number,0)
    body2q = Vector{Vector{Int}}()
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        xi,yi,xj,yj = state.coords.q
        bp1 = [xi,yi]
        bp2 = [xj,yj]
        bp1_find = findfirst(x->x==bp1,lncs)
        if bp1_find === nothing
            push!(lncs,bp1)
            push!(bp_number,bp_number[end]+1)
            bp1_number = bp_number[end]
        else
            bp1_number = bp1_find
        end
        bp2_find = findfirst(x->x==bp2,lncs)
        if bp2_find === nothing
            push!(lncs,bp2)
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

get_nconstraint(tg) = tg.nconstraint

function build_Φ(tgstruct)
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
    @inline @inbounds function inner_Φ(q,d)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = 0
        #is += nbodydof*nfixbodies
        for rbid in tgstruct.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            if nc > 0
                ret[is+1:is+nc] = rb.state.cache.cfuncs.Φ(q[pindex],d[is+1:is+nc])
                is += nc
            end
            ret[is+1:is+nbodyc] .=rb.state.cache.funcs.Φ(q[pindex],d[is+1:is+nbodyc])
            is += nbodyc
        end
        ret
    end
    inner_Φ
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

function get_q(tg)
    rbs = tg.rigidbodies
    @unpack body2q = tg.connectivity
    ncoords = tg.ncoords
    T = get_numbertype(tg)
    q = zeros(T,ncoords)
    q̇ = zeros(T,ncoords)
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        q[pindex] .= rbs[rbid].state.coords.q
        q̇[pindex] .= rbs[rbid].state.coords.q̇
    end
    return q,q̇
end

get_λ(tg) = zeros(get_numbertype(tg),tg.nconstraint)

function get_initial(tgstruct)
    q0,q̇0 = get_q(tgstruct)
    λ0 = get_λ(tgstruct)
    q0,q̇0,λ0
end

function lucompletepiv!(A)
    n=size(A, 1)
    rowpiv=zeros(Int, n)
    colpiv=zeros(Int, n)
    for k=1:n
        Asub = abs.(A[k:n, k:n])#Search for next pivot
        _, index_max = findmax(Asub)
        μ,λ = index_max.I
        μ += k-1; λ += k-1
        rowpiv[k] = μ
        A[[k,μ], 1:n] = A[[μ,k], 1:n]
        colpiv[k] = λ
        A[1:n, [k,λ]] = A[1:n, [λ,k]]
        if A[k,k]≠0
            ρ = k+1:n
            A[ρ,k] = A[ρ,k]./A[k,k]
            A[ρ,ρ] = A[ρ,ρ] - A[ρ,k:k]*A[k:k,ρ]
        end
    end
    return (rowpiv, colpiv)
end

function jac_singularity_check(tg)
    q,_ = get_q(tg)
    A = build_A(tg)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    for (rbid,rb) in enumerate(tg.rigidbodies)
        if rb.prop.movable && rb.prop.constrained
            q_rb = rb.state.coords.q
            Aq_rb = vcat(rb.state.cache.cfuncs.Φq(q_rb),
                         rb.state.cache.funcs.Φq(q_rb))
            rb_rank = rank(Aq_rb)
            intrinsic_Aq = rb.state.cache.funcs.Φq(q_rb)
            # @show rbid,lucompletepiv!(copy(intrinsic_Aq))
            # col_index = GECP(intrinsic_Aq)
            # @show rbid,col_index
            # @show rank(intrinsic_Aq[:,col_index[1:6]])
            if rb_rank < minimum(size(Aq_rb))
                @warn "The $(rbid)th rigid body's Jacobian is singular: rank(A(q))=$(rb_rank)<$(minimum(size(Aq_rb)))"
            end
        end
    end
end


function get_d(tg)
    @unpack nconstraint = tg
    rbs = tg.rigidbodies
    T = get_numbertype(tg)
    d = Vector{T}(undef,nconstraint)
    nbodyc = get_nbodyconstraint(tg)
    is = 0
    for rbid in tg.mvbodyindex
        rb = rbs[rbid]
        nc = rb.state.cache.nc
        if nc > 0
            d[is+1:is+nc] .= get_index_Φ_q0(rb)
            is += nc
        end
        d[is+1:is+nbodyc] .= NaturalCoordinates.get_deform(rb.state.cache.funcs.lncs)
        is += nbodyc
    end
    d
end

get_ndim(tr::TensegrityRobot) = get_ndim(tr.tg)
get_ndim(tg::TensegrityStructure) = get_ndim(tg.rigidbodies)
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = N

get_numbertype(tr::TensegrityRobot) = get_numbertype(tr.tg)
get_numbertype(tg::TensegrityStructure) = get_numbertype(tg.rigidbodies)
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = T

get_nbodyconstraint(tr::TensegrityRobot) = get_nbodyconstraint(tr.tg)
get_nbodyconstraint(tg::TensegrityStructure) = get_nbodyconstraint(tg.rigidbodies)
get_nbodyconstraint(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodyconstraint(eltype(rbs))
get_nbodyconstraint(rb::AbstractRigidBody) = get_nbodyconstraint(typeof(rb))
get_nbodyconstraint(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_nconstraint(lncsType)

get_nbodycoords(tr::TensegrityRobot) = get_nbodycoords(tr.tg)
get_nbodycoords(tg::TensegrityStructure) = get_nbodycoords(tg.rigidbodies)
get_nbodycoords(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodycoords(eltype(rbs))
get_nbodycoords(rb::AbstractRigidBody) = get_nbodycoords(typeof(rb))
get_nbodycoords(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_ncoords(lncsType)

get_nbodydof(tr::TensegrityRobot) = get_nbodydof(tr.tg)
get_nbodydof(tg::TensegrityStructure) = get_nbodydof(tg.rigidbodies)
get_nbodydof(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodydof(eltype(rbs))
get_nbodydof(rb::AbstractRigidBody) = get_nbodydof(typeof(rb))
get_nbodydof(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = 3
get_nbodydof(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = 6

get_gravity(tr::TensegrityRobot) = get_gravity(tr.tg)
get_gravity(tg::TensegrityStructure) = get_gravity(tg.rigidbodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rb::AbstractRigidBody) = get_gravity(typeof(rb))
get_gravity(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),-9.81*one(T)]
get_gravity(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),zero(T),-9.81*one(T)]

get_strings_len(bot::TensegrityRobot) = get_strings_len(bot.tg)
get_strings_deform(bot::TensegrityRobot) = get_strings_deform(bot.tg)
get_strings_restlen(bot::TensegrityRobot) = get_strings_restlen(bot.tg)
get_strings_len_dot(bot::TensegrityRobot) = get_strings_len_dot(bot.tg)

function get_strings_len!(tg::TensegrityStructure,q)
    distribute_q_to_rbs!(tg,q,zero(q))
    update_strings_apply_forces!(tg)
    get_strings_len(tg)
end

function get_strings_len(tg::TensegrityStructure)
    [s.state.length for s in tg.strings]
end

function get_strings_len_dot(tg::TensegrityStructure)
    [s.state.lengthdot for s in tg.strings]
end

function get_strings_deform(tg::TensegrityStructure)
    [s.state.length - s.state.restlen for s in tg.strings]
end

function get_strings_restlen(tg::TensegrityStructure)
    [s.state.restlen for s in tg.strings]
end

function get_original_restlen(tr_input::TensegrityRobot)
    tr = deepcopy(tr_input)
    T = get_numbertype(tr)
    actuate!(tr,zeros(T,length(tr.hub.actuators)))
    u0 = get_strings_restlen(tr.tg)
end

function force_densities_to_restlen(tg::TensegrityStructure,γs)
    [
    begin
        l = s.state.length
        l̇ = s.state.lengthdot
        k = s.k
        c = s.c
        u = l-(γ*l-c*l̇)/k
    end
        for (γ,s) in zip(γs,tg.strings)]
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


function new_trajectory(tg::TensegrityStructure)
    t0 = zero(get_numbertype(tg))
    q0, q̇0  = get_initial(tg)
    λ0 = get_λ(tg)
    NaturalCoordinatesTrajectory([t0], [q0], [q̇0], [λ0])
end

TensegrityRobot(tg,hub) = TensegrityRobot(tg,hub,new_trajectory(tg))

function reset!(tr::TensegrityRobot)
    @unpack tg, traj = tr
    @unpack qs, q̇s = traj
    distribute_q_to_rbs!(tg,qs[begin],q̇s[begin])
    reset!(traj)
end

function reset!(traj::TensegrityRobotTrajectory)
    @unpack ts, qs, q̇s, λs = traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
end

function set_new_initial!(tr::TensegrityRobot,q,q̇=zero(q))
    @unpack tg, traj = tr
    traj.qs[begin] .= q
    traj.q̇s[begin] .= q̇
    reset!(tr)
end
