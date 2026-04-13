
function activate_frictional_contacts!(structure,contact_env::ContactRigidBodies,q::AbstractVector;checkpersist=true,)
    (;  
        friction_coefficients,
        restitution_coefficients,
        persistent_bits,
    ) = structure.contacts_related
    activated_bits = copy(structure.contacts_related.activated_bits)

    (;bodies) = structure
    cnt = structure.connectivity
    (;bodyid2sys_locus_id) = cnt
    (;num_of_bodies) = cnt
    (;contact_bodies) = contact_env
    T = eltype(q)
    nq = length(q)
    na_ref = Ref(0)
    for contact in contact_bodies
        contact.pid = 0#contact.init_pid
        # INFO 3: the first locus stays fixed; start from the second one
    end
    update_bodies!(structure,q)

    # clear bits
    foreach(bodies) do body
        (;prop) = body
        bid = prop.id
        activated_bits[bodyid2sys_locus_id[bid]] .= false
        persistent_bits[bodyid2sys_locus_id[bid]] .= false
    end

    # activate!
    # loop all bodies (contactable) against all contact_bodies
    foreach(bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci_states) = state
        if body.prop.contactable
            # check each locus of the body
            for pid in eachindex(loci_states)
                locus_state = loci_states[pid]
                (;frame,contact_state) = locus_state
                gap, normal, cid, local_pos = contact_gap_and_normal(frame.position,contact_bodies,bodies)
                if !checkpersist
                    contact_state.active = false
                end
                activate!(contact_state,gap)
                if contact_state.active
                    # @show gap, frame.position
                    activated_bits[bodyid2sys_locus_id[bid][pid]] = true
                    # @show bid, pid
                    ibid = contact_bodies[cid].bid
                    ipid = contact_bodies[cid].pid
                    # @show ibid, ipid
                    contact_state.axes = Axes(normal)
                    if contact_state.persistent
                        persistent_bits[bodyid2sys_locus_id[bid][pid]] = true
                    end
                    foreach(bodies) do ibody
                        if ibody.prop.id == ibid
                            icontact_state = ibody.state.loci_states[ipid].contact_state
                            activate!(icontact_state, gap)
                            activated_bits[bodyid2sys_locus_id[ibid][ipid]] = true
                            icontact_state.axes = Axes(-normal)
                            (;axes, friction_coefficient, restitution_coefficient) = ibody.prop.loci[ipid]
                            ibody.prop.loci[ipid] = Locus(local_pos, axes, friction_coefficient, restitution_coefficient)
                            ibody.state.loci_states[ipid].frame.position = frame.position
                            coords = ibody.coords
                            ncoords = get_num_of_coords(coords)
                            iq = @MVector zeros(ncoords)
                            c(x) = to_local_coords(coords,x)
                            C(c) = to_position_jacobian(coords,iq,c)

                            ibody.cache.Cps[ipid] = C(c(local_pos))

                            if icontact_state.persistent
                                persistent_bits[bodyid2sys_locus_id[ibid][ipid]] = true
                            end
                        end
                    end
                    na_ref[] += 2
                    # INFO 6: one collision corresponds to two active contacts
                end
            end
        end
    end

    # else inactive
    foreach(bodies) do body
        (;prop) = body
        bid = prop.id
        for contact in contact_bodies
            if contact.bid == bid
                ipid = contact.pid
                for locus in body.state.loci_states[ipid+1:end]
                    locus.contact_state.active = false
                end
            end
        end
    end

    na = na_ref[]
    inv_friction_coefficients = ones(T,3na)
    for (i,μ) in enumerate(friction_coefficients[activated_bits])
        inv_friction_coefficients[3(i-1)+1] = 1/μ
    end
    H = Diagonal(inv_friction_coefficients)
    activated_restitution_coefficients = restitution_coefficients[activated_bits]
    # member's points idx to system's active points' idx
    bodyid2act_idx = deepcopy(bodyid2sys_locus_id)
    act_start = 0
    persistent_idx = Int[]
    for bid = 1:num_of_bodies
        bodyid2act_idx[bid] .= 0
        contacts_bits_body = findall(activated_bits[bodyid2sys_locus_id[bid]])
        nactive_body = length(contacts_bits_body)
        mem_idx = act_start+1:act_start+nactive_body
        bodyid2act_idx[bid][contacts_bits_body] .= mem_idx
        mem_per_idx = findall(persistent_bits[bodyid2sys_locus_id[bid]][contacts_bits_body])
        append!(persistent_idx,mem_idx[mem_per_idx])
        act_start += nactive_body
    end
    L = zeros(T,3na,3na)
    Lv = deepcopy(L)
    D = zeros(T,3na,nq)
    Dper = zero(D)
    Dimp = zero(D)
    ∂Dq̇∂q = zeros(T,3na,nq)
    ∂DᵀΛ∂q = zeros(T,nq,nq)
    ŕ = Vector{T}(undef,3na)
    cache = @eponymtuple(
        na, bodyid2act_idx, persistent_idx, activated_bits, H, activated_restitution_coefficients, D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, L, Lv
    )
end

function activate_frictional_contacts!(structure,contact_env::StaticEnvironment,q::AbstractVector;checkpersist=true,)
    (;  
        friction_coefficients,
        restitution_coefficients,
        persistent_bits,
    ) = structure.contacts_related
    activated_bits = copy(structure.contacts_related.activated_bits)
    cnt = structure.connectivity
    (;bodyid2sys_locus_id) = cnt
    (;num_of_bodies) = cnt
    surfaces = contact_env.geometry
    T = eltype(q)
    nq = length(q)
    na = 0
    update_bodies!(structure,CoordinatesState(nothing,q,zero(q),zero(q),zeros(T,nq),zeros(T,0),nothing,zeros(T,0),T[]))
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci_states) = state
        activated_bits[bodyid2sys_locus_id[bid]] .= false
        persistent_bits[bodyid2sys_locus_id[bid]] .= false
        if body.prop.contactable
            for pid in eachindex(loci_states)
                locus_state = loci_states[pid]
                (;frame,contact_state) = locus_state
                gap, normal = contact_gap_and_normal(frame.position,surfaces)
                if !checkpersist
                    contact_state.active = false
                end
                activate!(contact_state,gap)
                if contact_state.active
                    activated_bits[bodyid2sys_locus_id[bid][pid]] = true
                    contact_state.axes = Axes(normal)
                    if contact_state.persistent
                        persistent_bits[bodyid2sys_locus_id[bid][pid]] = true
                    end
                    na += 1
                end
            end
        end
    end
    # @show na, length(active_contacts)
    inv_friction_coefficients = ones(T,3na)
    for (i,μ) in enumerate(friction_coefficients[activated_bits])
        inv_friction_coefficients[3(i-1)+1] = 1/μ
    end
    H = Diagonal(inv_friction_coefficients)
    activated_restitution_coefficients = restitution_coefficients[activated_bits]
    # member's points idx to system's active points' idx
    bodyid2act_idx = deepcopy(bodyid2sys_locus_id)
    act_start = 0
    persistent_idx = Int[]
    for bid = 1:num_of_bodies
        bodyid2act_idx[bid] .= 0
        contacts_bits_body = findall(activated_bits[bodyid2sys_locus_id[bid]])
        nactive_body = length(contacts_bits_body)
        mem_idx = act_start+1:act_start+nactive_body
        bodyid2act_idx[bid][contacts_bits_body] .= mem_idx
        mem_per_idx = findall(persistent_bits[bodyid2sys_locus_id[bid]][contacts_bits_body])
        append!(persistent_idx,mem_idx[mem_per_idx])
        act_start += nactive_body
    end
    L = zeros(T,3na,3na)
    Lv = deepcopy(L)
    D = zeros(T,3na,nq)
    Dper = zero(D)
    Dimp = zero(D)
    ∂Dq̇∂q = zeros(T,3na,nq)
    ∂DᵀΛ∂q = zeros(T,nq,nq)
    ŕ = Vector{T}(undef,3na)
    Λ = zeros(T,3na)
    Γ  = zeros(T,3na)
    cache = @eponymtuple(
        na, 
        bodyid2act_idx, persistent_idx, activated_bits, 
        H, activated_restitution_coefficients, 
        D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, 
        ŕ, 
        L, Lv,
        Λ, Γ
    )
end

function activate_contacts!(structure,contact_env::ContactRigidBodies,q::AbstractVector;checkpersist=true,)
    (;  
        friction_coefficients,
        restitution_coefficients,
        persistent_bits,
    ) = structure.contacts_related
    activated_bits = copy(structure.contacts_related.activated_bits)
    (;bodies) = structure
    cnt = structure.connectivity
    (;bodyid2sys_locus_id) = cnt
    (;num_of_bodies) = cnt
    (;contact_bodies) = contact_env
    T = eltype(q)
    nq = length(q)
    na_ref = Ref(0)
    update_bodies!(structure,q)
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci_states) = state
        activated_bits[bodyid2sys_locus_id[bid]] .= false
        persistent_bits[bodyid2sys_locus_id[bid]] .= false
        if body.prop.contactable
            for pid in eachindex(loci_states)
                locus_state = loci_states[pid]
                (;frame,contact_state) = locus_state
                gap, normal = contact_gap_and_normal(frame.position,contact_bodies,bodies)
                if !checkpersist
                    contact_state.active = false
                end
                activate!(contact_state,gap)
                if contact_state.active
                    activated_bits[bodyid2sys_locus_id[bid][pid]] = true
                    contact_state.axes = Axes(normal)
                    if contact_state.persistent
                        persistent_bits[bodyid2sys_locus_id[bid][pid]] = true
                    end
                    na_ref[] += 1
                end
            end
        end
    end
    # @show na, length(active_contacts)
    na = na_ref[]
    inv_friction_coefficients = ones(T,na)
    for (i,μ) in enumerate(friction_coefficients[activated_bits])
        inv_friction_coefficients[(i-1)+1] = 1
    end
    H = Diagonal(inv_friction_coefficients)
    activated_restitution_coefficients = restitution_coefficients[activated_bits]
    # member's points idx to system's active points' idx
    bodyid2act_idx = deepcopy(bodyid2sys_locus_id)
    act_start = 0
    persistent_idx = Int[]
    for bid = 1:num_of_bodies
        bodyid2act_idx[bid] .= 0
        contacts_bits_body = findall(activated_bits[bodyid2sys_locus_id[bid]])
        nactive_body = length(contacts_bits_body)
        mem_idx = act_start+1:act_start+nactive_body
        bodyid2act_idx[bid][contacts_bits_body] .= mem_idx
        mem_per_idx = findall(persistent_bits[bodyid2sys_locus_id[bid]][contacts_bits_body])
        append!(persistent_idx,mem_idx[mem_per_idx])
        act_start += nactive_body
    end
    # @show bodyid2act_idx
    Ls = [
        begin 
            na_body = count(!iszero, mem)
            zeros(T,na_body,na_body)
        end
        for mem in bodyid2act_idx
    ]
    L = BlockDiagonal(Ls)
    Lv = deepcopy(L)
    D = zeros(T,na,nq)
    Dper = zero(D)
    Dimp = zero(D)
    ∂Dq̇∂q = zeros(T,na,nq)
    ∂DᵀΛ∂q = zeros(T,nq,nq)
    ŕ = Vector{T}(undef,na)
    cache = @eponymtuple(
        na, bodyid2act_idx, persistent_idx, activated_bits, H, activated_restitution_coefficients, D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, L, Lv
    )
end

function activate_contacts!(structure,contact_env::StaticEnvironment,q::AbstractVector;checkpersist=true,)
    (;  
        friction_coefficients,
        restitution_coefficients,
        persistent_bits,
    ) = structure.contacts_related
    activated_bits = copy(structure.contacts_related.activated_bits)
    cnt = structure.connectivity
    (;bodyid2sys_locus_id) = cnt
    (;num_of_bodies) = cnt
    surfaces = contact_env.geometry
    T = eltype(q)
    nq = length(q)
    na = 0
    update_bodies!(structure,q)
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci_states) = state
        activated_bits[bodyid2sys_locus_id[bid]] .= false
        persistent_bits[bodyid2sys_locus_id[bid]] .= false
        if body.prop.contactable
            for pid in eachindex(loci_states)
                locus_state = loci_states[pid]
                (;frame,contact_state) = locus_state
                gap, normal = contact_gap_and_normal(frame.position,surfaces)
                if !checkpersist
                    contact_state.active = false
                end
                activate!(contact_state,gap)
                if contact_state.active
                    activated_bits[bodyid2sys_locus_id[bid][pid]] = true
                    contact_state.axes = Axes(normal)
                    if contact_state.persistent
                        persistent_bits[bodyid2sys_locus_id[bid][pid]] = true
                    end
                    na += 1
                end
            end
        end
    end
    # @show na, length(active_contacts)
    inv_friction_coefficients = ones(T,na)
    for (i,μ) in enumerate(friction_coefficients[activated_bits])
        inv_friction_coefficients[(i-1)+1] = 1
    end
    H = Diagonal(inv_friction_coefficients)
    activated_restitution_coefficients = restitution_coefficients[activated_bits]
    # member's points idx to system's active points' idx
    bodyid2act_idx = deepcopy(bodyid2sys_locus_id)
    act_start = 0
    persistent_idx = Int[]
    for bid = 1:num_of_bodies
        bodyid2act_idx[bid] .= 0
        contacts_bits_body = findall(activated_bits[bodyid2sys_locus_id[bid]])
        nactive_body = length(contacts_bits_body)
        mem_idx = act_start+1:act_start+nactive_body
        bodyid2act_idx[bid][contacts_bits_body] .= mem_idx
        mem_per_idx = findall(persistent_bits[bodyid2sys_locus_id[bid]][contacts_bits_body])
        append!(persistent_idx,mem_idx[mem_per_idx])
        act_start += nactive_body
    end
    # @show bodyid2act_idx
    Ls = [
        begin 
            na_body = count(!iszero, mem)
            zeros(T,na_body,na_body)
        end
        for mem in bodyid2act_idx
    ]
    L = BlockDiagonal(Ls)
    Lv = deepcopy(L)
    D = zeros(T,na,nq)
    Dper = zero(D)
    Dimp = zero(D)
    ∂Dq̇∂q = zeros(T,na,nq)
    ∂DᵀΛ∂q = zeros(T,nq,nq)
    ŕ = Vector{T}(undef,na)
    cache = @eponymtuple(
        na, bodyid2act_idx, persistent_idx, activated_bits, H, activated_restitution_coefficients, D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, L, Lv
    )
end

function get_frictional_directions_and_positions!(structure,contact_cache, inst_state, Λ, )
    (;
        D, Dper,Dimp, 
        ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, 
        bodyid2act_idx,
    ) = contact_cache
    update_bodies!(structure,inst_state)
    ∂Dq̇∂q .= 0
    ∂DᵀΛ∂q .= 0
    foreach(structure.bodies) do body
        (;prop,state,coords,cache) = body
        if body isa AbstractRigidBody
            Cps = cache.Cps
        elseif body isa AbstractFlexibleBody
            Cps = cache.Sps
        end
        bid = prop.id
        (;loci_states) = state
        for (pid,locus_state) in enumerate(loci_states)
            (;frame,contact_state) = locus_state
            if contact_state.active
                (;normal,tangent,bitangent) = contact_state.axes
                C = Cps[pid]
                CT = C*build_T(structure,bid)
                dm = hcat(normal,tangent,bitangent) |> transpose
                ci = bodyid2act_idx[bid][pid]
                # @show bid,pid,ci
                # @show bodyid2act_idx[bid]
                epi = 3(ci-1)+1:3ci
                D[epi,:] = dm*CT
                ŕ[epi]   = dm*frame.position
                #TODO proper dispatch for different coordinates
                # if coords isa QCF.QC
                #     Tbody = build_T(structure,bid)
                #     locus = prop.loci[pid]
                #     ∂Cq̇∂q = QCF.to_velocity_jacobian(coords,Tbody*q,Tbody*q̇,locus.position)*Tbody
                #     ∂Dq̇∂q[epi,:] = dm*∂Cq̇∂q
                #     Λi = @view Λ[epi]
                #     fi = dm'*Λi
                #     ∂Cᵀfi∂q = QCF.to_force_jacobian(Tbody*q,fi,locus.position)
                #     ∂DᵀΛ∂q .+= transpose(Tbody)*∂Cᵀfi∂q*Tbody
                # end
                if contact_state.persistent
                    Dper[epi,:] .= D[epi,:]
                else
                    Dimp[epi,:] .= D[epi,:]
                end
            end
        end
    end
end

function get_directions_and_positions!(structure,contact_cache, q, q̇, Λ, )
    (;
        D, Dper,Dimp, 
        ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, 
        bodyid2act_idx,
    ) = contact_cache
    update_bodies!(structure,q)
    ∂Dq̇∂q .= 0
    ∂DᵀΛ∂q .= 0
    foreach(structure.bodies) do body
        (;prop,state,coords,cache) = body
        bid = prop.id
        (;loci_states) = state
        for (pid,locus_state) in enumerate(loci_states)
            (;contact_state) = locus_state
            if contact_state.active
                (;position) = loci_states[pid].frame
                (;normal,tangent,bitangent) = contact_state.axes
                C = cache.Cps[pid]
                CT = C*build_T(structure,bid)
                dm = normal |> transpose
                ci = bodyid2act_idx[bid][pid]
                epi = ci
                D[epi,:] = dm*CT
                ŕ[epi]   = dm*position
                if coords isa QCF.QC
                    Tbody = build_T(structure,bid)
                    locus = prop.loci[pid]
                    ∂Cq̇∂q = QCF.to_velocity_jacobian(coords,Tbody*q,Tbody*q̇,locus.position)*Tbody
                    ∂Dq̇∂q[epi,:] = dm*∂Cq̇∂q
                    fi = Λ[epi]*normal
                    ∂Cᵀfi∂q = QCF.to_force_jacobian(Tbody*q,fi,locus.position)
                    ∂DᵀΛ∂q .+= transpose(Tbody)*∂Cᵀfi∂q*Tbody
                end
                if contact_state.persistent
                    Dper[epi,:] .= D[epi,:]
                else
                    Dimp[epi,:] .= D[epi,:]
                end
            end
        end
    end
end

function get_distribution_law!(structure,contact_cache,inst_state)
    (;
        D,L,Lv,H,bodyid2act_idx
    ) = contact_cache
    N_in = intrinsic_nullspace(structure,inst_state)
    # N_ex = extrinsic_nullspace(structure,q)
    A = cstr_jacobian(structure,inst_state)
    N_ex = nullspace(A*N_in)
    N = N_in*N_ex
    #fixme intrinsic_nullspace not correct for shared coordinates or boundary conditions
    N = nullspace(A)
    # @show rank(N), (A*N |> norm)
    R = D*N
    L .= (I-pinv(R)'*R')*H
    Lv .= (I-R*R')
end
