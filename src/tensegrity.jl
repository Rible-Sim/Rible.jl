struct TensegrityRobot{tgT,hubT,trajT}
    tg::tgT
    hub::hubT
    traj::trajT
end
abstract type AbstractTensegrity end
abstract type TensegrityRobotTrajectory{T} end

struct ConstrainedCoordinatesTrajectory{T} <: TensegrityRobotTrajectory{T}
    ts::Vector{T}
    qs::Vector{Vector{T}}
    q̇s::Vector{Vector{T}}
    λs::Vector{Vector{T}}
end

struct SlidingConstrainedCoordinatesTrajectory{T} <: TensegrityRobotTrajectory{T}
    ts::Vector{T}
    qs::Vector{Vector{T}}
    q̇s::Vector{Vector{T}}
    λs::Vector{Vector{T}}
    s̄s::Vector{Vector{T}}
end

struct ID
    rbid::Int
    apid::Int
end

# struct Connectivity{bType,sType,cType}
#     body2q::bType
#     string2ap::sType
#     contacts::cType
# end
#
# struct Cluster_Connectivity{bType,sType,cType}
#     body2q::bType
#     string2ap::sType
#     clusterstring2ap::Vector{sType}
#     contacts::cType
# end

Connectivity(b) = (body2q=b,)
Connectivity(b,s) = (body2q=b,string2ap=s)
Connectivity(b,s,c) = (body2q=b,string2ap=s,clusterstring2ap=c)

struct TensegrityStructure{BodyType,StrType,TenType,CntType,CstType} <: AbstractTensegrity
    ndim::Int
    ncoords::Int
    nconstraint::Int
    ndof::Int
    nbodies::Int
    nmvbodies::Int
    mvbodyindex::Vector{Int}
    nfixbodies::Int
    fixbodyindex::Vector{Int}
    npoints::Int
    nmvpoints::Int
    nstrings::Int
    rigidbodies::Vector{BodyType}
    strings::StrType
    tensiles::TenType
    connectivity::CntType
    constraints::CstType
end

function TensegrityStructure(rbs::Vector{rbT},tensiles::TT,cnt,
                            constraints = [EmptyConstraint()]) where {TT,rbT<:AbstractRigidBody{N,T}} where {N,T}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbodies = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbodies = length(fixbodyindex)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.naps
    end
    nmvpoints = 0
    for rbid in mvbodyindex
        nmvpoints += rbs[rbid].prop.naps
    end
    ncoords = maximum(maximum.(cnt.body2q))
    strings = tensiles.strings
    nstrings = length(strings)
    nconstraint = get_nconstraint(rbs,mvbodyindex,nmvbodies,constraints)
    ndof = ncoords - nconstraint
    tg = TensegrityStructure(
                    ndim,
                    ncoords,nconstraint,ndof,
                    nbodies,nmvbodies,mvbodyindex,nfixbodies,fixbodyindex,
                    npoints,nmvpoints,
                    nstrings,
                    rbs,strings,tensiles,
                    cnt,constraints)
    check_jacobian_singularity(tg)
    tg
end

struct ClusterTensegrityStructure{BodyType,StrType,CStrType,TenType,CntType,CstType} <: AbstractTensegrity
    ndim::Int
    ncoords::Int
    nconstraint::Int
    ndof::Int
    nbodies::Int
    nmvbodies::Int
    mvbodyindex::Vector{Int}
    nfixbodies::Int
    fixbodyindex::Vector{Int}
    npoints::Int
    nmvpoints::Int
    nstrings::Int
    nclusterstrings::Int
    nslidings::Int
    rigidbodies::Vector{BodyType}
    strings::StrType
    clusterstrings::CStrType
    tensiles::TenType
    connectivity::CntType
    constraints::CstType
end

function ClusterTensegrityStructure(rbs::Vector{rbT},tensiles::TT,cnt,
                            constraints = [EmptyConstraint()]) where {TT,rbT<:AbstractRigidBody{N,T}} where {N,T}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbodies = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbodies = length(fixbodyindex)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.naps
    end
    nmvpoints = 0
    for rbid in mvbodyindex
        nmvpoints += rbs[rbid].prop.naps
    end
    ncoords = maximum(maximum.(cnt.body2q))
    nconstraint = get_nconstraint(rbs,mvbodyindex,nmvbodies,constraints)
    ndof = ncoords - nconstraint
    strings = tensiles.strings
    nstrings = length(strings)
    clusterstrings = tensiles.clusterstrings
    nclusterstrings = length(clusterstrings)
    nslidings = sum(length(cs.sps.s) for cs in clusterstrings)
    tg = ClusterTensegrityStructure(
                    ndim,
                    ncoords,nconstraint,ndof,
                    nbodies,nmvbodies,mvbodyindex,nfixbodies,fixbodyindex,
                    npoints,nmvpoints,
                    nstrings,nclusterstrings,nslidings,
                    rbs,strings,clusterstrings,tensiles,
                    cnt,constraints)
    check_jacobian_singularity(tg)
    tg
end

function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function reset_forces!(tg::AbstractTensegrity)
    reset_forces!.(tg.rigidbodies)
end

function update_strings!(tg)
    update_strings!(tg, tg.tensiles)
end

function update_strings!(tg, @eponymargs(strings,))
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    for sstring in strings
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

function update_strings!(tg,@eponymargs(SMA_strings,))
    rbs = tg.rigidbodies
    cnt = tg.connectivity
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

function update_strings!(tg, @eponymargs(clusterstrings,))
    rbs = tg.rigidbodies
    cnt = tg.connectivity

    for clusterstring in clusterstrings
        s = clusterstring.sps.s
        for (segid, seg) in enumerate(clusterstring.segs)
            @unpack k,c,prestress,original_restlen = seg
            u0 = original_restlen
            segstate = seg.state
            a,b = cnt.clusterstring2ap[clusterstring.ID][segid]
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
            segstate.length,segstate.direction = lengthdir(p2-p1)
            l = segstate.length
            τ = segstate.direction
            segstate.lengthdot = (transpose(Δr)*Δṙ)/l
            if segid == 1
                u = u0 + s[segid]
            elseif segid == length(clusterstring.segs)
                u = u0 - s[segid-1]
            else
                u = u0 + s[segid] - s[segid-1]
            end
            segstate.tension = k*u0/u*(l-u)
            𝐟 = τ*segstate.tension
            f1 .+=  𝐟
            f2 .+= -𝐟
        end
    end
end

function update_strings!(tg, @eponymargs(strings,clusterstrings))
    update_strings!(tg,@eponymtuple(strings))
    update_strings!(tg,@eponymtuple(clusterstrings))
end


distribute_q_to_rbs!(tg,globalq) = distribute_q_to_rbs!(tg,globalq,zero(globalq))
function distribute_q_to_rbs!(tg,globalq,globalq̇)
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    for rbid in tg.mvbodyindex
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

function apply_gravity_y!(tgstruct;factor=1)
    rbs = tgstruct.rigidbodies
    gravity_acceleration = factor*get_gravity_y(tgstruct)
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

function kinetic_energy_coords(tg::AbstractTensegrity,q,q̇)
    distribute_q_to_rbs!(tg,q,q̇)
    ke = sum(kinetic_energy_coords.(tg.rigidbodies))
end

function gravity_potential_energy(tg::AbstractTensegrity,q)
    distribute_q_to_rbs!(tg,q)
    sum(gravity_potential_energy.(tg.rigidbodies))
end

function elastic_potential_energy(tg::TensegrityStructure)
    reset_forces!(tg)
    update_strings_apply_forces!(tg)
    pe = sum(potential_energy.(tg.strings))
end

function elastic_potential_energy(tg::TensegrityStructure,q)
    distribute_q_to_rbs!(tg,q)
    elastic_potential_energy(tg)
end

function elastic_potential_energy(bot::TensegrityRobot,q,a)
    actuate!(bot,a)
    elastic_potential_energy(bot.tg,q)
end

function energy(tg,q,q̇;gravity=false)
    distribute_q_to_rbs!(tg,q,q̇)
    ke = sum(kinetic_energy_coords.(tg.rigidbodies))
    update_strings_apply_forces!(tg)
    epe = sum(potential_energy.(tg.strings))
    if gravity
        gpe = gravity_potential_energy(tg,q)
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

function build_massmatrix(tg::AbstractTensegrity)
    body2q = tg.connectivity.body2q
    ncoords = tg.ncoords
    T = get_numbertype(tg)
    mass_matrix = zeros(T,ncoords,ncoords)
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= tg.rigidbodies[rbid].state.cache.M
    end
    mass_matrix
end

function get_nconstraint(rbs,mvbodyindex,nmvbodies,constraints)
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
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            # nc = rb.state.cache.nc
            # if nc > 0
            #     ret[is[]+1:is[]+nc] = rb.state.cache.cfuncs.Φ(q[pindex])
            #     is[] += nc
            # end
            ret[is[]+1:is[]+nbodyc] .= rb.state.cache.funcs.Φ(q[pindex])
            is[] += nbodyc
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
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            # if nc > 0
            #     ret[is[]+1:is[]+nc,pindex] = rb.state.cache.cfuncs.Φq(q[pindex])
            #     is[] += nc
            # end
            ret[is[]+1:is[]+nbodyc,pindex] .= rb.state.cache.funcs.Φq(q[pindex])
            is[] += nbodyc
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

function check_jacobian_singularity(tg)
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
get_ndim(tg::AbstractTensegrity) = get_ndim(tg.rigidbodies)
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = N

get_numbertype(bot::TensegrityRobot) = get_numbertype(bot.tg)
get_numbertype(tg::AbstractTensegrity) = get_numbertype(tg.rigidbodies)
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = T

get_nbodyconstraint(bot::TensegrityRobot) = get_nbodyconstraint(bot.tg)
get_nbodyconstraint(tg::AbstractTensegrity) = get_nbodyconstraint(tg.rigidbodies)
get_nbodyconstraint(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodyconstraint(eltype(rbs))
get_nbodyconstraint(rb::AbstractRigidBody) = get_nbodyconstraint(typeof(rb))
get_nbodyconstraint(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_nconstraint(lncsType)

get_nbodycoords(bot::TensegrityRobot) = get_nbodycoords(bot.tg)
get_nbodycoords(tg::AbstractTensegrity) = get_nbodycoords(tg.rigidbodies)
get_nbodycoords(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodycoords(eltype(rbs))
get_nbodycoords(rb::AbstractRigidBody) = get_nbodycoords(typeof(rb))
get_nbodycoords(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_ncoords(lncsType)

get_nbodydof(bot::TensegrityRobot) = get_nbodydof(bot.tg)
get_nbodydof(tg::AbstractTensegrity) = get_nbodydof(tg.rigidbodies)
get_nbodydof(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodydof(eltype(rbs))
get_nbodydof(rb::AbstractRigidBody) = get_nbodydof(typeof(rb))
get_nbodydof(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = 3
get_nbodydof(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = 6

get_gravity(bot::TensegrityRobot) = get_gravity(bot.tg)
get_gravity(tg::AbstractTensegrity) = get_gravity(tg.rigidbodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rb::AbstractRigidBody) = get_gravity(typeof(rb))
get_gravity(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),-9.81*one(T)]
get_gravity(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),zero(T),-9.81*one(T)]


get_gravity_y(tr::TensegrityRobot) = get_gravity_y(tr.tg)
get_gravity_y(tg::AbstractTensegrity) = get_gravity_y(tg.rigidbodies)
get_gravity_y(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity_y(eltype(rbs))
get_gravity_y(rb::AbstractRigidBody) = get_gravity_y(typeof(rb))
get_gravity_y(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),9.81*one(T)]
get_gravity_y(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),9.81*one(T),zero(T)]

get_strings_len(bot::TensegrityRobot) = get_strings_len(bot.tg)
get_strings_deform(bot::TensegrityRobot) = get_strings_deform(bot.tg)
get_strings_restlen(bot::TensegrityRobot) = get_strings_restlen(bot.tg)
get_strings_len_dot(bot::TensegrityRobot) = get_strings_len_dot(bot.tg)
get_strings_tension(bot::TensegrityRobot) = get_strings_tension(bot.tg)
get_strings_stiffness(bot::TensegrityRobot) = get_strings_stiffness(bot.tg)

function get_strings_len!(tg::AbstractTensegrity,q)
    distribute_q_to_rbs!(tg,q,zero(q))
    update_strings_apply_forces!(tg)
    get_strings_len(tg)
end

function get_strings_stiffness(tg::AbstractTensegrity)
    [s.k for s in tg.strings]
end

function get_strings_len(tg::AbstractTensegrity)
    [s.state.length for s in tg.strings]
end

function get_strings_len_dot(tg::AbstractTensegrity)
    [s.state.lengthdot for s in tg.strings]
end

function get_strings_deform(tg::AbstractTensegrity)
    [s.state.length - s.state.restlen for s in tg.strings]
end

function get_strings_restlen(tg::AbstractTensegrity)
    [s.state.restlen for s in tg.strings]
end

function get_strings_tension(tg::AbstractTensegrity)
    [s.state.tension for s in tg.strings]
end

function get_original_restlen(botinput::TensegrityRobot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_strings_restlen(bot.tg)
end

function force_densities_to_restlen(tg::AbstractTensegrity,γs)
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
    ConstrainedCoordinatesTrajectory([t0], [q0], [q̇0], [λ0])
end

function new_trajectory(tg::ClusterTensegrityStructure)
    t0 = zero(get_numbertype(tg))
    q0, q̇0  = get_initial(tg)
    λ0 = get_λ(tg)
    s̄0 = get_s̄(tg)
    SlidingConstrainedCoordinatesTrajectory([t0], [q0], [q̇0], [λ0], [s̄0])
end

function TensegrityRobot(tg,hub)
	reset_forces!(tg)
    # update_strings_apply_forces!(tg)
	# check_jacobian_singularity(tg)
	# check_stability(tg)
    TensegrityRobot(tg,hub,new_trajectory(tg))
end

function reset!(bot::TensegrityRobot)
    @unpack tg, traj = bot
    reset!(tg,traj)
    reset!(traj)
end

function reset!(tg::TensegrityStructure,traj)
    @unpack qs,q̇s = traj
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,qs[begin],q̇s[begin])
    update_strings!(tg)
end

function reset!(tg::ClusterTensegrityStructure,traj)
    @unpack qs,q̇s,s̄s = traj
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,qs[begin],q̇s[begin])
    distribute_s̄!(tg,s̄s[begin])
    update_strings!(tg)
end

function reset!(traj::ConstrainedCoordinatesTrajectory)
    @unpack ts, qs, q̇s, λs = traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
end

function reset!(traj::SlidingConstrainedCoordinatesTrajectory)
    @unpack ts,qs,q̇s,λs,s̄s= traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
    resize!(s̄s,1)
end

function set_new_initial!(bot::TensegrityRobot,q,q̇=zero(q))
    @unpack tg, traj = bot
    traj.qs[begin] .= q
    traj.q̇s[begin] .= q̇
    reset!(bot)
end
