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

function dynfuncs(bot;actuate=false,gravity=false,(Fˣ!)=(F,t)->nothing)
    (;st) = bot
    function F!(F,q,q̇,t)
        if actuate
            RB.actuate!(bot,[t])
        end
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_tensiles!(st)
        if gravity
            RB.apply_gravity!(st)
        end
        F .= RB.generate_forces!(st)
        Fˣ!(F,t)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
        ∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_tensiles!(st)
        RB.build_∂Q̌∂q̌!(∂F∂q̌,st)
        RB.build_∂Q̌∂q̌̇!(∂F∂q̌̇,st)
    end
    @eponymtuple(F!,Jac_F!)
end

function contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0]),
        checkpersist = true,
    )
    (;st) = bot
    (;mem2num) = st.connectivity.numbered
    npoints = length.(mem2num) |> sum
    contacts_bits = BitVector(undef,npoints)
    persistent_bits = BitVector(undef,npoints)
    T = RB.get_numbertype(st)
    μs_sys = ones(T,npoints)
    es_sys = zeros(T,npoints)
    gaps_sys = fill(typemax(T),npoints)

    # initilize
    foreach(st.bodies) do body
        (;prop,state) = body
        bid = prop.id
        (;loci) = prop
        μs_sys[mem2num[bid]] .= [locus.friction_coefficient for locus in loci]
        es_sys[mem2num[bid]] .= [locus.restitution_coefficient for locus in loci]
    end

    function F!(F,q,q̇,t)
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_tensiles!(st)
        RB.apply_gravity!(st)
        F .= RB.generate_forces!(st)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
        ∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_tensiles!(st)
        RB.build_∂Q̌∂q̌!(∂F∂q̌,st)
        RB.build_∂Q̌∂q̌̇!(∂F∂q̌̇,st)
    end
    
    function prepare_contacts!(q)
        T = eltype(q)
        nq = length(q)
        na = 0
        RB.update_rigids!(st,q)
        foreach(st.bodies) do body
            (;prop,state) = body
            bid = prop.id
            (;contacts,loci_states) = state
            contacts_bits[mem2num[bid]] .= false
            persistent_bits[mem2num[bid]] .= false
            if body isa RB.AbstractRigidBody
                for pid in eachindex(loci_states)
                    (;position) = loci_states[pid]
                    contact = contacts[pid]
                    (;e,μ) = contact
                    gap = RB.signed_distance(position,flatplane)
                    if !checkpersist
                        contact.state.active = false
                    end
                    RB.activate!(contact,gap)
                    if contact.state.active
                        contacts_bits[mem2num[bid][pid]] = true
                        contact.state.frame = RB.spatial_frame(flatplane.n)
                        if contact.state.persistent
                            persistent_bits[mem2num[bid][pid]] = true
                        end
                        na += 1
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
        # member's points indices to system's active points' indices
        mem2act_idx = deepcopy(mem2num)
        act_start = 0
        persistent_indices = Int[]
        for bid = 1:st.nbodies
            mem2act_idx[bid] .= 0
            contacts_bits_body = findall(contacts_bits[mem2num[bid]])
            nactive_body = length(contacts_bits_body)
            mem_idx = act_start+1:act_start+nactive_body
            mem2act_idx[bid][contacts_bits_body] .= mem_idx
            mem_per_idx = findall(persistent_bits[mem2num[bid]][contacts_bits_body])
            append!(persistent_indices,mem_idx[mem_per_idx])
            act_start += nactive_body
        end
        Ls = [
            begin 
                na_body = count(!iszero, mem)
                zeros(T,3na_body,3na_body)
            end
            for mem in mem2act_idx
        ]
        L = BlockDiagonal(Ls)
        D = zeros(T,3na,nq)
        Dper = zero(D)
        Dimp = zero(D)
        ∂Dq̇∂q = zeros(T,3na,nq)
        ∂DᵀΛ∂q = zeros(T,nq,nq)
        ŕ = Vector{T}(undef,3na)
        na, mem2act_idx, persistent_indices, contacts_bits, H, restitution_coefficients, D, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, L
    end

    function get_directions_and_positions!(D, Dper,Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕ, q, q̇, Λ, mem2act_idx,)
        RB.update_rigids!(st,q)
        ∂Dq̇∂q .= 0
        ∂DᵀΛ∂q .= 0
        foreach(st.bodies) do body
            (;prop,state) = body
            bid = prop.id
            (;contacts,loci_states) = state
            for (pid,contact) in enumerate(contacts)
                if contact.state.active
                    (;position) = loci_states[pid]
                    (;normal,tangent,bitangent) = contact.state.frame
                    C = state.cache.Cps[pid]
                    CT = C*RB.build_T(st,bid)
                    dm = hcat(normal,tangent,bitangent) |> transpose
                    ci = mem2act_idx[bid][pid]
                    epi = 3(ci-1)+1:3ci
                    D[epi,:] = dm*CT
                    ŕ[epi]   = dm*position
                    if state.cache.funcs.nmcs isa RB.QBF.QC
                        Tbody = RB.build_T(st,bid)
                        locus = prop.loci[pid]
                        ∂Cẋ∂x = RB.QBF.make_∂Cẋ∂x(locus.position)
                        ∂Cq̇∂q = ∂Cẋ∂x(Tbody*q,Tbody*q̇)*Tbody
                        ∂Dq̇∂q[epi,:] = dm*∂Cq̇∂q
                        ∂Cᵀf∂x = RB.QBF.make_∂Cᵀf∂x(locus.position)
                        Λi = @view Λ[epi]
                        fi = dm'*Λi
                        ∂DᵀΛ∂q .+= transpose(Tbody)*∂Cᵀf∂x(Tbody*q,fi)*Tbody
                    end
                    if contact.state.persistent
                        Dper[epi,:] .= D[epi,:]
                    else
                        Dimp[epi,:] .= D[epi,:]
                    end
                end
            end
        end
    end

    function get_distribution_law!(L,mem2act_idx,q)
        T = eltype(q)
        RB.update_rigids!(st,q)
        foreach(st.bodies) do body
            (;prop,state) = body
            bid = prop.id
            (;contacts,loci_states) = state
            active_idx = findall(!iszero,mem2act_idx[bid])
            na_body = length(active_idx)
            if na_body > 1
                R = zeros(T,3na_body,6)
                inv_μs_body = ones(T,3na_body)
                for (i,pid) in enumerate(active_idx)
                    (;position) = loci_states[pid]
                    contact = contacts[pid]
                    (;e,μ) = contact
                    (;normal,tangent,bitangent) = contact.state.frame
                    inv_μs_body[3(i-1)+1] = 1/μ
                    dm = hcat(normal,tangent,bitangent) |> transpose
                    R[3(i-1)+1:3(i-1)+3,1:3] = dm
                    R[3(i-1)+1:3(i-1)+3,4:6] = dm*(-RB.skew(position))
                end
                blocks(L)[bid] .= (I-pinv(R)'*R')*Diagonal(inv_μs_body)
            end
        end
    end

    @eponymtuple(
        F!,Jac_F!,
        prepare_contacts!,
        get_directions_and_positions!,
        get_distribution_law!
    )
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

    RB.PrescribedActuator(
        1,
        RB.ManualActuator(1,collect(1:nμ),zeros(nμ),RB.Uncoupled()),
        itp
    )
end
