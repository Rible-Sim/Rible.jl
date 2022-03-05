
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

struct ID{RBType,APType}
    rbsig::RBType
    pid::APType
end

struct Connectivity{bType,iType,sType,aType,cType}
    body2fullq::bType
    body2q::bType
    mvindices::iType
    fixindices::iType
    string2ap::sType
    apnb::aType
    contacts::cType
end

Connectivity(bf,b,mi,fi) = Connectivity(bf,b,mi,fi,nothing,nothing,nothing)
Connectivity(bf,b,mi,fi,s) = Connectivity(bf,b,mi,fi,s,nothing,nothing)
Connectivity(bf,b,mi,fi,s,a) = Connectivity(bf,b,mi,fi,s,a,nothing)

struct TensegrityStructure{BodyType,StrType,TenType,CntType,CstType}
    ndim::Int
    ncoords::Int
    nfullcoords::Int
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
    rigidbodies::BodyType
    strings::StrType
    tensiles::TenType
    connectivity::CntType
    constraints::CstType
end

function TensegrityStructure{N,T}(rbs,tensiles,cnt::Connectivity,
                            constraints = [EmptyConstraint()]) where {N,T}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = Vector{Int}()
    fixbodyindex = Vector{Int}()
    foreach(rbs) do rb
        if rb.prop.movable
            push!(mvbodyindex,rb.prop.id)
        else
            push!(fixbodyindex,rb.prop.id)
        end
    end
    nmvbodies = length(mvbodyindex)
    nfixbodies = length(fixbodyindex)
    nstrings = sum(map(length,tensiles))
    npoints = 0
    foreach(rbs) do rb
        npoints += rb.prop.naps
    end
    nmvpoints = 0
    foreach(rbs) do rb
        if rb.prop.id in mvbodyindex
            nmvpoints += rb.prop.naps
        end
    end
    ncoords = maximum(maximum.(cnt.body2q[mvbodyindex]))
    nfullcoords = maximum(maximum.(cnt.body2fullq))
    nconstraint = get_nconstraint(rbs,mvbodyindex,nmvbodies,nbodies,constraints)
    ndof = ncoords - nconstraint
    strings = tensiles.strings
    tg = TensegrityStructure(ndim,
                    ncoords,nfullcoords,nconstraint,ndof,
                    nbodies,
                    nmvbodies,mvbodyindex,
                    nfixbodies,fixbodyindex,
                    nstrings,
                    npoints,nmvpoints,
                    rbs,strings,tensiles,cnt,constraints)
    check_jacobian_singularity(tg)
    tg
end

function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function reset_forces!(tg::TensegrityStructure)
    reset_forces!(tg.rigidbodies)
end

function reset_forces!(rigidbodies::AbstractVector)
    foreach(reset_forces!,rigidbodies)
end

function reset_forces!(rigidbodies::TypeSortedCollection)
    foreach(reset_forces!,rigidbodies)
end

function update!(tg::TensegrityStructure)
    q,_ = get_q(tg)
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,q)
    update_strings_apply_forces!(tg)
end

function update_strings_apply_forces!(tg)
    rbs = tg.rigidbodies
    ss = tg.strings
    cnt = tg.connectivity
    (;string2ap) = cnt
    foreach(string2ap) do scnt
        sstring = ss[scnt.id]
        @unpack id,k,c = sstring
        sstate = sstring.state
        state1 = scnt.end1.rbsig.state
        state2 = scnt.end2.rbsig.state
        apid1 = scnt.end1.pid
        apid2 = scnt.end2.pid
        p1 = state1.rps[apid1]
        ṗ1 = state1.ṙps[apid1]
        f1 = state1.Faps[apid1]
        p2 = state2.rps[apid2]
        ṗ2 = state2.ṙps[apid2]
        f2 = state2.Faps[apid2]
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

distribute_q_to_rbs!(tg,globalq) = distribute_q_to_rbs!(tg,globalq,zero(globalq))
function distribute_q_to_rbs!(tg,globalq,globalq̇)
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        @unpack q, q̇ = rb.state.coords
        if rbid in tg.mvbodyindex
            pindex = cnt.body2q[rbid]
            uci = rb.state.cache.unconstrained_index
            q[uci] .= globalq[pindex]
            q̇[uci] .= globalq̇[pindex]
        end
        @unpack cache,rps,ṙps,ro,ṙo,rg,ṙg = rb.state
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

function update_rbs_states!(tg,q,q̇=zero(q))
    distribute_q_to_rbs!(tg,q,q̇)
    rbs = tg.rigidbodies
    for rbid in tg.mvbodyindex
        rb = rbs[rbid]
        lncs = rb.state.cache.funcs.lncs
        @unpack q, q̇ = rb.state.coords
        R = NaturalCoordinates.find_R(lncs,q)
        Ω = NaturalCoordinates.find_ω(lncs,q,q̇)
        rb.state.R .= R
        # @show Ω
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

function assemble_forces!(F,tg;factor=1.0)
    rbs = tg.rigidbodies
    @unpack body2q = tg.connectivity
    generate_forces!(rbs)
    F .= 0.0
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        F[pindex] .+= factor*rbs[rbid].state.coords.Q
    end
end

function assemble_forces(tg;factor=1.0)
    T = get_numbertype(tg)
    @unpack body2q = tg.connectivity
    F = zeros(T,tg.ncoords)
    assemble_forces!(F,tg,factor=factor)
    F
end

function apply_gravity!(tg;factor=1)
    rbs = tg.rigidbodies
    gravity_acceleration = factor*get_gravity(tg)
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

function renumbering!(body2q_raw,sync=false,rbs=nothing)
    imax = 0
    mvindices = Vector{Int}()
    fixindices = Vector{Int}()
    for (rbid,pindex) in enumerate(body2q_raw)
        if !(rbs isa Nothing)
            if !rbs[rbid].prop.movable
                foreach(pindex) do index
                    if !(index in fixindices)
                        push!(fixindices,index)
                    end
                end
                resize!(pindex,0)
                continue
            else
                ci = Vector(rbs[rbid].state.cache.constrained_index)
                foreach(pindex[ci]) do index
                    if !(index in fixindices)
                        push!(fixindices,index)
                    end
                end
                deleteat!(pindex,ci)
            end
        end
        for (pid,index) in enumerate(pindex)
            # @show pindex
            if !(index in mvindices)
                push!(mvindices,index)
            end
            indexdiff = index - imax
            if  indexdiff == 1
                imax += 1
            elseif indexdiff ≥ 1
                imax += 1
                pindex[pid] = imax
                if sync
                    for restpindex in body2q_raw[rbid+1:end]
                        ifind = findfirst((x)->x==index,restpindex)
                        if !(ifind isa Nothing)
                            restpindex[ifind]=imax
                        end
                    end
                end
            end
        end
    end
    mvindices,fixindices
end

function make_body2q(rbs,shared_index;move_only=true)
    nbodycoords = get_nbodycoords(rbs[1])
    nb = length(rbs)
    body2q_raw = collect.([(i-1)*nbodycoords+1:i*nbodycoords for i = 1:nb])
    # @show body2q_raw
    for ishare in shared_index
        rbid1 = ishare.end1.rbsig
        rbid2 = ishare.end2.rbsig
        pid1 = ishare.end1.pid
        pid2 = ishare.end2.pid
        body2q_raw[rbid2][pid2] = body2q_raw[rbid1][pid1]
        # @show body2q_raw
        renumbering!(body2q_raw)
    end
    # @show body2q_raw
    # if move_only
    body2fullq = deepcopy(body2q_raw)
    mvindices,fixindices = renumbering!(deepcopy(body2q_raw),false,rbs)
    renumbering!(body2q_raw,true,rbs)
    # end
    # @show body2q_raw
    body2fullq, body2q_raw, mvindices, fixindices
end

function make_apnb(rigidbodies)
	ret = Vector{Vector{Int}}()
	is = 0
	foreach(rigidbodies) do rb
		push!(ret,collect(is+1:is+rb.prop.naps))
		is += rb.prop.naps
	end
	ret
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

function get_nconstraint(rbs,mvbodyindex,nmvbodies,nbodies,constraints)
    nbodyconstraint = get_nbodyconstraint(rbs)
    nbodydof = get_nbodydof(rbs)
    ninconstraint = nbodyconstraint*nmvbodies
    nexconstraint = 0  #nbodydof*nfixbodies
    foreach(constraints) do cst
        nexconstraint += cst.nconstraints
    end
    nconstraint = ninconstraint + nexconstraint
end

get_nconstraint(tg) = tg.nconstraint

function build_Φ(tg)
    rbs = tg.rigidbodies
    csts = tg.constraints
    #q0,q̇0 = get_q(tg)
    @unpack body2q = tg.connectivity
    nfixbodies = tg.nfixbodies
    nconstraint = tg.nconstraint
    nbodyc = get_nbodyconstraint(tg)
    nbodydof = get_nbodydof(tg)
    @inline @inbounds function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = Ref(0)
        #is[] += nbodydof*nfixbodies
        foreach(rbs) do rb
            rbid = rb.prop.id
            mvrbid = findfirst((x)->x==rbid,tg.mvbodyindex)
            if !(mvrbid isa Nothing)
                pindex = body2q[rbid]
                # nc = rb.state.cache.nc
                # if nc > 0
                #     ret[is[]+1:is[]+nc] = rb.state.cache.cfuncs.Φ(q[pindex])
                #     is[] += nc
                # end
                ret[nbodyc*(mvrbid-1)+1:nbodyc*mvrbid] .= rb.state.cache.funcs.Φ(q[pindex])
                is[] += nbodyc
            end
        end
        foreach(csts) do cst
            nc = cst.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(cst)(q)
            is[] += nc
        end
        ret
    end
    @inline @inbounds function inner_Φ(q,d)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = Ref(0)
        #is += nbodydof*nfixbodies
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            # nc = rb.state.cache.nc
            # if nc > 0
            #     ret[is[]+1:is[]+nc] = rb.state.cache.cfuncs.Φ(q[pindex],d[is[]+1:is[]+nc])
            #     is[] += nc
            # end
            ret[is[]+1:is[]+nbodyc] .=rb.state.cache.funcs.Φ(q[pindex],d[is[]+1:is[]+nbodyc])
            is[] += nbodyc
        end
        foreach(csts) do cst
            nc = cst.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(cst)(q,d[is[]+1:is[]+nc])
            is[] += nc
        end
        ret
    end
    inner_Φ
end

function build_A(tg)
    rbs = tg.rigidbodies
    csts = tg.constraints
    @unpack body2q = tg.connectivity
    nfixbodies = tg.nfixbodies
    nconstraint = tg.nconstraint
    nbodyc = get_nbodyconstraint(tg)
    nbodydof = get_nbodydof(tg)
    ncoords = tg.ncoords
    @inline @inbounds function inner_A(q)
        ret = zeros(eltype(q),nconstraint,ncoords)
        is = Ref(0)
        foreach(rbs) do rb
            rbid = rb.prop.id
            mvrbid = findfirst((x)->x==rbid,tg.mvbodyindex)
            if !(mvrbid isa Nothing)
                pindex = body2q[rbid]
                # nc = rb.state.cache.nc
                # if nc > 0
                #     ret[is[]+1:is[]+nc,pindex] = rb.state.cache.cfuncs.Φq(q[pindex])
                #     is[] += nc
                # end
                ret[nbodyc*(mvrbid-1)+1:nbodyc*mvrbid,pindex] .= rb.state.cache.funcs.Φq(q[pindex])
                is[] += nbodyc
            end
        end
        foreach(csts) do cst
            nc = cst.nconstraints
            ret[is[]+1:is[]+nc,:] .= make_A(cst)(q)
            is[] += nc
        end
        ret
    end
end

function build_F(tg,rbid,pid,f)
    rbs = tg.rigidbodies
    Ti = build_Ti(tg,rbid)
    C = rbs[rbid].state.cache.Cp[pid]
    F = transpose(C*Ti)*f
    reshape(F,:,1)
end

function get_fullq(tg)
    rbs = tg.rigidbodies
    @unpack body2fullq = tg.connectivity
    nfullcoords = maximum.(body2fullq) |> maximum
    T = get_numbertype(tg)
    q = zeros(T,nfullcoords)
    q̇ = zeros(T,nfullcoords)
    foreach(rbs) do rb
        pindex = body2fullq[rb.prop.id]
        q[pindex] .= rb.state.coords.q
        q̇[pindex] .= rb.state.coords.q̇
    end
    return q,q̇
end

function get_q(tg)
    rbs = tg.rigidbodies
    @unpack body2q = tg.connectivity
    ncoords = tg.ncoords
    T = get_numbertype(tg)
    q = zeros(T,ncoords)
    q̇ = zeros(T,ncoords)
    foreach(rbs) do rb
        if rb.prop.id in tg.mvbodyindex
            pindex = body2q[rb.prop.id]
            uci = rb.state.cache.unconstrained_index
            q[pindex] .= rb.state.coords.q[uci]
            q̇[pindex] .= rb.state.coords.q̇[uci]
        end
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

function check_jacobian_singularity(tg)
    q,_ = get_q(tg)
    A = build_A(tg)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        if rb.prop.movable && rb.prop.constrained
            uci = rb.state.cache.unconstrained_index
            q_rb = rb.state.coords.q[uci]
            Aq_rb = rb.state.cache.funcs.Φq(q_rb)
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

get_s(bot::TensegrityRobot) = get_s(bot.tg)

function get_s(tg::TensegrityStructure)
    1 ./get_strings_len(tg)
end

get_c(bot::TensegrityRobot) = get_c(bot.tg)
function get_c(tg)
    (;ndim,npoints) = tg
    cnt = tg.connectivity
    (;apnb) = cnt
    T = get_numbertype(tg)
    ret = Vector{T}(undef,npoints*ndim)
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        for i in 1:rb.prop.naps
            ip = apnb[rbid][i]
            ret[(ip-1)*ndim+1:ip*ndim] = rb.state.cache.funcs.c(rb.prop.aps[i])
        end
    end
    ret
end


function set_C!(tg,c)
    (;ndim,npoints) = tg
    cnt = tg.connectivity
    (;apnb) = cnt
    T = get_numbertype(tg)
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        for i in 1:rb.prop.naps
            ip = apnb[rbid][i]
            ci = c[(ip-1)*ndim+1:ip*ndim]
            (;r̄i,X̄) = rb.state.cache.funcs.lncs
            # @set rb.prop.aps[i] = r̄i + X̄*ci
            # @show i
            # @show r̄i + X̄*ci
            # @show rb.prop.aps[i]
            rb.state.cache.Cp[i] = rb.state.cache.funcs.C(ci)
            # rb.state.cache.funcs.c(rb.prop.aps[i])
        end
    end
end

function get_d(tg)
    @unpack nconstraint = tg
    rbs = tg.rigidbodies
    csts = tg.constraints
    T = get_numbertype(tg)
    d = Vector{T}(undef,nconstraint)
    nbodyc = get_nbodyconstraint(tg)
    is = Ref(0)
    for rbid in tg.mvbodyindex
        rb = rbs[rbid]
        d[is[]+1:is[]+nbodyc] .= NaturalCoordinates.get_deform(rb.state.cache.funcs.lncs)
        is[] += nbodyc
    end
    foreach(csts) do cst
        nc = cst.nconstraints
        d[is[]+1:is[]+nc] .= cst.values
        is[] += nc
    end
    d
end

get_ndim(bot::TensegrityRobot) = get_ndim(bot.tg)
get_ndim(tg::TensegrityStructure) = get_ndim(tg.rigidbodies)
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rbs::TypeSortedCollection) = get_ndim(eltype(rbs.data[1]))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = N

get_numbertype(bot::TensegrityRobot) = get_numbertype(bot.tg)
get_numbertype(tg::TensegrityStructure) = get_numbertype(tg.rigidbodies)
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rbs::TypeSortedCollection) = get_numbertype(eltype(rbs.data[1]))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = T

get_nbodyconstraint(bot::TensegrityRobot) = get_nbodyconstraint(bot.tg)
get_nbodyconstraint(tg::TensegrityStructure) = get_nbodyconstraint(tg.rigidbodies)
get_nbodyconstraint(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodyconstraint(eltype(rbs))
get_nbodyconstraint(rbs::TypeSortedCollection) = get_nbodyconstraint(eltype(rbs.data[1]))
get_nbodyconstraint(rb::AbstractRigidBody) = get_nbodyconstraint(typeof(rb))
get_nbodyconstraint(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_nconstraint(lncsType)

get_nbodycoords(bot::TensegrityRobot) = get_nbodycoords(bot.tg)
get_nbodycoords(tg::TensegrityStructure) = get_nbodycoords(tg.rigidbodies)
get_nbodycoords(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodycoords(eltype(rbs))
get_nbodycoords(rbs::TypeSortedCollection) = get_nbodycoords(eltype(rbs.data[1]))
get_nbodycoords(rb::AbstractRigidBody) = get_nbodycoords(typeof(rb))
get_nbodycoords(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_ncoords(lncsType)

get_nbodydof(bot::TensegrityRobot) = get_nbodydof(bot.tg)
get_nbodydof(tg::TensegrityStructure) = get_nbodydof(tg.rigidbodies)
get_nbodydof(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodydof(eltype(rbs))
get_nbodydof(rbs::TypeSortedCollection) = get_nbodydof(eltype(rbs.data[1]))
get_nbodydof(rb::AbstractRigidBody) = get_nbodydof(typeof(rb))
get_nbodydof(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = 3
get_nbodydof(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = 6

get_gravity(bot::TensegrityRobot) = get_gravity(bot.tg)
get_gravity(tg::TensegrityStructure) = get_gravity(tg.rigidbodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rb::AbstractRigidBody) = get_gravity(typeof(rb))
get_gravity(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),-9.81*one(T)]
get_gravity(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),zero(T),-9.81*one(T)]

get_strings_len(bot::TensegrityRobot) = get_strings_len(bot.tg)
get_strings_deform(bot::TensegrityRobot) = get_strings_deform(bot.tg)
get_strings_restlen(bot::TensegrityRobot) = get_strings_restlen(bot.tg)
get_strings_len_dot(bot::TensegrityRobot) = get_strings_len_dot(bot.tg)
get_strings_tension(bot::TensegrityRobot) = get_strings_tension(bot.tg)
get_strings_stiffness(bot::TensegrityRobot) = get_strings_stiffness(bot.tg)
get_strings_force_density(bot::TensegrityRobot) = get_strings_force_density(bot.tg)

function get_strings_len!(tg::TensegrityStructure,q)
    distribute_q_to_rbs!(tg,q,zero(q))
    update_strings_apply_forces!(tg)
    get_strings_len(tg)
end

function get_strings_stiffness(tg::TensegrityStructure)
    [s.k for s in tg.strings]
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

function get_strings_tension(tg::TensegrityStructure)
    [s.state.tension for s in tg.strings]
end

function get_strings_force_density(tg::TensegrityStructure)
    [s.state.tension/s.state.length for s in tg.strings]
end

function get_original_restlen(botinput::TensegrityRobot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_strings_restlen(bot.tg)
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

function get_indices(tg)
    (;ncoords) = tg
    (;body2q) = tg.connectivity
    mvindices_bool = ones(Bool,ncoords)
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        if rbid in tg.fixbodyindex
            mvindices_bool[body2q[rbid]] .= false
        end
    end
    # @show mvindices_bool
    findall(mvindices_bool),findall(.!mvindices_bool)
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

function build_Y(bot)
	@unpack tg, hub = bot
	@unpack actuators = hub
    @unpack nstrings,strings = tg
    nact = length(actuators)
    ret = spzeros(Int,nstrings,nact)
    for (i,iact) in enumerate(actuators)
		if typeof(iact)<:ManualActuator
			is1 = iact.reg.id_string
	        ret[is1,i] = 1
		elseif typeof(iact)<:ManualGangedActuators
	        is1, is2 = iact.regs.id_strings
	        ret[is1,i] = 1
	        ret[is2,i] = -1
		else
			error("Unknown actuator type")
		end
    end
    ret
end

function new_trajectory(tg::TensegrityStructure)
    t0 = zero(get_numbertype(tg))
    q0, q̇0  = get_initial(tg)
    λ0 = get_λ(tg)
    NaturalCoordinatesTrajectory([t0], [q0], [q̇0], [λ0])
end

function TensegrityRobot(tg,hub)
	reset_forces!(tg)
    update_strings_apply_forces!(tg)
	# check_jacobian_singularity(tg)
	# check_stability(tg)
    TensegrityRobot(tg,hub,new_trajectory(tg))
end

function reset!(bot::TensegrityRobot)
    @unpack tg, traj = bot
    @unpack qs, q̇s = traj
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,qs[begin],q̇s[begin])
    update_strings_apply_forces!(tg)
    reset!(traj)
end

function reset!(traj::TensegrityRobotTrajectory)
    @unpack ts, qs, q̇s, λs = traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
end

function set_new_initial!(bot::TensegrityRobot,q,q̇=zero(q))
    @unpack tg, traj = bot
    traj.qs[begin] .= q
    traj.q̇s[begin] .= q̇
    reset!(bot)
end
