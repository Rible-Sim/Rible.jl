material_properties = Table(
    [
        (name = "ABS",              density = 1045u"kg/m^3",     modulus_elas = 3.1u"GPa"),
        # (name = "Acetate", density = 2.0u"GPa"),
        (name = "Acrylic",          density = 1185u"kg/m^3",     modulus_elas = 3.5u"GPa"),
        # (name = "Cartilage", density = 2.0u"GPa"),
        (name = "Diamond",          density = 3300u"kg/m^3",     modulus_elas = 1200u"GPa"),
        (name = "Graphite fibre",   density = 2300u"kg/m^3",     modulus_elas = 700u"GPa"),
        (name = "Ceramics",         density = 3600u"kg/m^3",     modulus_elas = 276u"GPa"),
        # (name = "Carbides", density = 2.0u"GPa"),
        # (name = "CFR epoxy", density = 2.0u"GPa"),
        (name = "Concrete",         density = 2400u"kg/m^3",     modulus_elas = 14u"GPa"),
        # (name = "Douglas fir", density = 2.0u"GPa"),
        (name = "Epoxy",            density = 1150u"kg/m^3",     modulus_elas = 8u"GPa"),
        # (name = "GFR polyester", density = 2.0u"GPa"),
        (name = "Glasses",          density = 4100u"kg/m^3",     modulus_elas = 70u"GPa"),
        # (name = "HDPE", density = 2.0u"GPa"),
        # (name = "LDPE", density = 2.0u"GPa"),
        (name = "Kevlar",           density = 1435u"kg/m^3",     modulus_elas = 131u"GPa"),
        (name = "Nylon 66",         density = 1150u"kg/m^3",     modulus_elas = 3.3u"GPa"),
        (name = "Oak",              density =  650u"kg/m^3",     modulus_elas =  12u"GPa"),
        (name = "Teak",             density =  630u"kg/m^3",     modulus_elas =  12u"GPa"),
        (name = "Polyester",        density = 1350u"kg/m^3",     modulus_elas = 2.4u"GPa"),
        (name = "PC",               density = 1200u"kg/m^3",     modulus_elas = 2.4u"GPa"),
        (name = "PP",               density = 900u"kg/m^3",      modulus_elas = 1.6u"GPa"),
        (name = "PVC",              density = 1700u"kg/m^3",     modulus_elas = 4.1u"GPa"),
        (name = "Aluminium",        density = 2710u"kg/m^3",     modulus_elas = 71u"GPa"),
        (name = "Al-Cu alloy",      density = 2800u"kg/m^3",     modulus_elas = 75u"GPa"),
        (name = "Al-Mg alloy",      density = 2725u"kg/m^3",     modulus_elas = 71u"GPa"),
        (name = "Alloy steel",      density = 7900u"kg/m^3",     modulus_elas = 210u"GPa"),
        (name = "Brass",            density = 8500u"kg/m^3",     modulus_elas = 104u"GPa"),
        (name = "Bronze",           density = 8800u"kg/m^3",     modulus_elas = 117u"GPa"),
        (name = "Cast iron",        density = 7150u"kg/m^3",     modulus_elas = 97u"GPa"),
        (name = "Copper",           density = 8950u"kg/m^3",     modulus_elas = 117u"GPa"),
        (name = "Gold",             density = 19300u"kg/m^3",    modulus_elas = 79u"GPa"),
        (name = "Iron",             density = 7850u"kg/m^3",     modulus_elas = 206u"GPa"),
        (name = "Lead",             density = 11370u"kg/m^3",    modulus_elas = 18u"GPa"),
        (name = "Ni-Alloy",         density = 9000u"kg/m^3",     modulus_elas = 200u"GPa"),
        (name = "Platinum",         density = 21040u"kg/m^3",    modulus_elas = 164u"GPa"),
        (name = "Silver",           density = 10530u"kg/m^3",    modulus_elas = 78u"GPa"),
        (name = "Stainless steel",  density = 7930u"kg/m^3",     modulus_elas = 200u"GPa"),
        (name = "Titanium",         density = 4540u"kg/m^3",     modulus_elas = 118u"GPa"),
    ]
)

function generalize_force!(F,bot,q,q̇,t;actuate=false,gravity=true,(user_defined_force!)=(F,t)->nothing)
    if actuate
        actuate!(bot,[t])
    end
    (;structure) = bot
    clear_forces!(structure)
    lazy_update_bodies!(structure,q,q̇)
    update_tensiles!(structure)
    if gravity
        apply_gravity!(structure)
    end
    F .= assemble_force!(structure)
    user_defined_force!(F,t)
end

function generalize_force_jacobain!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    (;structure) = bot
    ∂F∂q̌ .= 0
    ∂F∂q̌̇ .= 0
    clear_forces!(structure)
    lazy_update_bodies!(structure,q,q̇)
    update_tensiles!(structure)
    build_∂Q̌∂q̌!(∂F∂q̌,structure)
    build_∂Q̌∂q̌̇!(∂F∂q̌̇,structure)
end

struct ContactCache{cacheType}
    cache::cacheType
end

function prepare_contacts(bot,env;
        checkpersist = true,
    )
    (;structure) = bot
    (;bodyid2sys_loci_idx) = structure.connectivity.numbered
    npoints = length.(bodyid2sys_loci_idx) |> sum
    contacts_bits = BitVector(undef,npoints)
    persistent_bits = BitVector(undef,npoints)
    T = get_numbertype(structure)
    μs_sys = ones(T,npoints)
    es_sys = zeros(T,npoints)
    gaps_sys = fill(typemax(T),npoints)

    # initilize
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci) = prop
        μs_sys[bodyid2sys_loci_idx[bid]] .= [locus.friction_coefficient for locus in loci]
        es_sys[bodyid2sys_loci_idx[bid]] .= [locus.restitution_coefficient for locus in loci]
    end

    prepared = @eponymtuple(
        contacts_bits,
        persistent_bits,
        μs_sys,
        es_sys,
        gaps_sys
    )
end

function activate_frictional_contacts!(structure,contact_env,solver_cache,q;checkpersist=true)
    (;  
        μs_sys,
        es_sys,
        contacts_bits,
        persistent_bits,
    ) = solver_cache.cache

    (;bodyid2sys_loci_idx) = structure.connectivity.numbered
    (;surfaces) = contact_env
    T = eltype(q)
    nq = length(q)
    na = 0
    update_bodies!(structure,q)
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci_states) = state
        contacts_bits[bodyid2sys_loci_idx[bid]] .= false
        persistent_bits[bodyid2sys_loci_idx[bid]] .= false
        if body isa AbstractRigidBody
            if body.prop.contactable
                for pid in eachindex(loci_states)
                    locus_state = loci_states[pid]
                    (;position,contact_state) = locus_state
                    gap, normal = contact_gap_and_normal(position,surfaces)
                    if !checkpersist
                        contact_state.active = false
                    end
                    activate!(contact_state,gap)
                    if contact_state.active
                        contacts_bits[bodyid2sys_loci_idx[bid][pid]] = true
                        contact_state.axes = spatial_axes(normal)
                        if contact_state.persistent
                            persistent_bits[bodyid2sys_loci_idx[bid][pid]] = true
                        end
                        na += 1
                    end
                end
            end
        end
    end
    # @show na, length(active_contacts)
    inv_friction_coefficients = ones(T,3na)
    for (i,μ) in enumerate(μs_sys[contacts_bits])
        inv_friction_coefficients[3(i-1)+1] = 1/μ
    end
    H = Diagonal(inv_friction_coefficients)
    restitution_coefficients = es_sys[contacts_bits]
    # member's points idx to system's active points' idx
    bodyid2act_idx = deepcopy(bodyid2sys_loci_idx)
    act_start = 0
    persistent_idx = Int[]
    for bid = 1:structure.num_of_bodies
        bodyid2act_idx[bid] .= 0
        contacts_bits_body = findall(contacts_bits[bodyid2sys_loci_idx[bid]])
        nactive_body = length(contacts_bits_body)
        mem_idx = act_start+1:act_start+nactive_body
        bodyid2act_idx[bid][contacts_bits_body] .= mem_idx
        mem_per_idx = findall(persistent_bits[bodyid2sys_loci_idx[bid]][contacts_bits_body])
        append!(persistent_idx,mem_idx[mem_per_idx])
        act_start += nactive_body
    end
    Ls = [
        begin 
            na_body = count(!iszero, mem)
            zeros(T,3na_body,3na_body)
        end
        for mem in bodyid2act_idx
    ]
    L = BlockDiagonal(Ls)
    D = zeros(T,3na,nq)
    Dper = zero(D)
    Dimp = zero(D)
    ∂Dq̇∂q = zeros(T,3na,nq)
    ∂DᵀΛ∂q = zeros(T,nq,nq)
    ŕ = Vector{T}(undef,3na)
    cache = @eponymtuple(
        na, bodyid2act_idx, persistent_idx, contacts_bits, H, restitution_coefficients, D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, L
    )
    ContactCache(cache)
end

function activate_contacts!(structure,contact_env,solver_cache,q;checkpersist=true)
    (;  
        μs_sys,
        es_sys,
        contacts_bits,
        persistent_bits,
    ) = solver_cache.cache

    (;bodyid2sys_loci_idx) = structure.connectivity.numbered
    (;surfaces) = contact_env
    T = eltype(q)
    nq = length(q)
    na = 0
    update_bodies!(structure,q)
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci_states) = state
        contacts_bits[bodyid2sys_loci_idx[bid]] .= false
        persistent_bits[bodyid2sys_loci_idx[bid]] .= false
        if body isa AbstractRigidBody
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
                        contacts_bits[bodyid2sys_loci_idx[bid][pid]] = true
                        contact_state.axes = spatial_axes(normal)
                        if contact_state.persistent
                            persistent_bits[bodyid2sys_loci_idx[bid][pid]] = true
                        end
                        na += 1
                    end
                end
            end
        end
    end
    # @show na, length(active_contacts)
    inv_friction_coefficients = ones(T,na)
    for (i,μ) in enumerate(μs_sys[contacts_bits])
        inv_friction_coefficients[(i-1)+1] = 1
    end
    H = Diagonal(inv_friction_coefficients)
    restitution_coefficients = es_sys[contacts_bits]
    # member's points idx to system's active points' idx
    bodyid2act_idx = deepcopy(bodyid2sys_loci_idx)
    act_start = 0
    persistent_idx = Int[]
    for bid = 1:structure.num_of_bodies
        bodyid2act_idx[bid] .= 0
        contacts_bits_body = findall(contacts_bits[bodyid2sys_loci_idx[bid]])
        nactive_body = length(contacts_bits_body)
        mem_idx = act_start+1:act_start+nactive_body
        bodyid2act_idx[bid][contacts_bits_body] .= mem_idx
        mem_per_idx = findall(persistent_bits[bodyid2sys_loci_idx[bid]][contacts_bits_body])
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
    D = zeros(T,na,nq)
    Dper = zero(D)
    Dimp = zero(D)
    ∂Dq̇∂q = zeros(T,na,nq)
    ∂DᵀΛ∂q = zeros(T,nq,nq)
    ŕ = Vector{T}(undef,na)
    cache = @eponymtuple(
        na, bodyid2act_idx, persistent_idx, contacts_bits, H, restitution_coefficients, D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, L
    )
    ContactCache(cache)
end

function get_frictional_directions_and_positions!(structure,cache, q, q̇, Λ, )
    (;
        D, Dper,Dimp, 
        ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, 
        bodyid2act_idx,
    ) = cache.cache
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
                (;position) = loci_states[pid]
                (;normal,tangent,bitangent) = contact_state.axes
                C = cache.Cps[pid]
                CT = C*build_T(structure,bid)
                dm = hcat(normal,tangent,bitangent) |> transpose
                ci = bodyid2act_idx[bid][pid]
                epi = 3(ci-1)+1:3ci
                D[epi,:] = dm*CT
                ŕ[epi]   = dm*position
                if coords.nmcs isa QCF.QC
                    Tbody = build_T(structure,bid)
                    locus = prop.loci[pid]
                    ∂Cq̇∂q = QCF.∂Cẋ∂x(Tbody*q,Tbody*q̇,locus.position)*Tbody
                    ∂Dq̇∂q[epi,:] = dm*∂Cq̇∂q
                    Λi = @view Λ[epi]
                    fi = dm'*Λi
                    ∂Cᵀfi∂q = QCF.∂Cᵀf∂x(Tbody*q,fi,locus.position)
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

function get_directions_and_positions!(structure,cache, q, q̇, Λ, )
    (;
        D, Dper,Dimp, 
        ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, 
        bodyid2act_idx,
    ) = cache.cache
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
                if coords.nmcs isa QCF.QC
                    Tbody = build_T(structure,bid)
                    locus = prop.loci[pid]
                    ∂Cq̇∂q = QCF.∂Cẋ∂x(Tbody*q,Tbody*q̇,locus.position)*Tbody
                    ∂Dq̇∂q[epi,:] = dm*∂Cq̇∂q
                    fi = Λ[epi]*normal
                    ∂Cᵀfi∂q = QCF.∂Cᵀf∂x(Tbody*q,fi,locus.position)
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

function get_frictional_distribution_law!(structure,cache,q)
    (;
        L,bodyid2act_idx
    ) = cache.cache
    T = eltype(q)
    update_bodies!(structure,q)
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci) = prop
        (;loci_states) = state
        active_idx = findall(!iszero,bodyid2act_idx[bid])
        na_body = length(active_idx)
        if na_body > 1
            R = zeros(T,3na_body,6)
            inv_μs_body = ones(T,3na_body)
            for (i,pid) in enumerate(active_idx)
                locus_state = loci_states[pid]
                locus = loci[pid]
                (;position,contact_state) = locus_state
                (;normal,tangent,bitangent) = contact_state.axes
                inv_μs_body[3(i-1)+1] = 1/locus.friction_coefficient
                dm = hcat(normal,tangent,bitangent) |> transpose
                R[3(i-1)+1:3(i-1)+3,1:3] = dm
                R[3(i-1)+1:3(i-1)+3,4:6] = dm*(-skew(position))
            end
            BlockDiagonals.blocks(L)[bid] .= (I-pinv(R)'*R')*Diagonal(inv_μs_body)
        end
    end
end


function get_distribution_law!(structure,cache,q)
    (;
        L,bodyid2act_idx
    ) = cache.cache
    T = eltype(q)
    update_bodies!(structure,q)
    foreach(structure.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci) = prop
        (;loci_states) = state
        active_idx = findall(!iszero,bodyid2act_idx[bid])
        na_body = length(active_idx)
        if na_body > 1
            R = zeros(T,na_body,6)
            inv_μs_body = ones(T,na_body)
            for (i,pid) in enumerate(active_idx)
                locus_state = loci_states[pid]
                locus = loci[pid]
                (;position,contact_state) = locus_state
                (;normal,tangent,bitangent) = contact_state.axes
                inv_μs_body[i] = 1/locus.friction_coefficient
                dm = hcat(normal,tangent,bitangent) |> transpose
                R[i,1:3] = dm
                R[i,4:6] = dm*(-skew(position))
            end
            BlockDiagonals.blocks(L)[bid] .= (I-pinv(R)'*R')*Diagonal(inv_μs_body)
        end
    end
end


function make_pres_actor(μ0,μ1,start,stop)
    nμ = length(μ0)

    function itp(t)
        scaled_itps = extrapolate(
            Interpolations.scale(
                interpolate(
                    hcat(μ0,μ1),
                    (NoInterp(),BSpline(Linear()))
                    # (NoInterp(),BSpline(Quadratic(Flat(OnGrid()))))
                ),
                1:nμ, start:stop-start:stop
            ),
            (Throw(),Flat())
        )
        [scaled_itps(j,t) for j in 1:nμ]
    end

    PrescribedActuator(
        1,
        ManualActuator(1,collect(1:nμ),zeros(nμ),Uncoupled()),
        itp
    )
end
