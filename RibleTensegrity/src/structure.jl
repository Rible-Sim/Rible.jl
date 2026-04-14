"""
TensegrityStructure Type.
$(TYPEDEF)
"""
struct TensegrityStructure{BodyType,TenType,CntType,StateType,CRType,CacheType} <: AbstractStructure{BodyType,TenType,CntType}
    num_of_dim::Int
    bodies::BodyType
    apparatuses::TenType
    connectivity::CntType
    state::StateType
    contacts_related::CRType
    cache::CacheType
end

"""
TensegrityStructure Constructor.
$(TYPEDSIGNATURES)
"""
function TensegrityStructure(bodies,apparatuses,cnt::AbstractConnectivity)
    num_of_dims = get_num_of_dims(bodies)
    state = StructureState(bodies,apparatuses,cnt)
    (;bodyid2sys_locus_id) = cnt
    num_of_sys_loci = length.(bodyid2sys_locus_id) |> sum
    activated_bits = BitVector(undef,num_of_sys_loci)
    persistent_bits = BitVector(undef,num_of_sys_loci)
    T = get_numbertype(bodies)
    friction_coefficients = ones(T,num_of_sys_loci)
    restitution_coefficients = zeros(T,num_of_sys_loci)
    gaps = fill(typemax(T),num_of_sys_loci)

    # initilize
    foreach(bodies) do body
        (;prop,) = body
        bid = prop.id
        (;loci) = prop
        friction_coefficients[bodyid2sys_locus_id[bid]] .= [locus.friction_coefficient for locus in loci]
        restitution_coefficients[bodyid2sys_locus_id[bid]] .= [locus.restitution_coefficient for locus in loci]
    end

    contacts_related = @eponymtuple(
        activated_bits,
        persistent_bits,
        friction_coefficients,
        restitution_coefficients,
        gaps
    )
    cache = StructureCache(bodies, cnt)
    structure = TensegrityStructure(
        num_of_dims,
        bodies,
        apparatuses,
        cnt,
        state,
        contacts_related,
        cache
    )
    check_jacobian_singularity(structure)
    check_constraints_consistency(structure)
    structure
end

function build_material_stiffness_matrix!(structure::AbstractStructure,num_of_full_coords,bodyid2sys_full_coords,q,k)
    (;num_of_dim) = structure
    structure.state.system.q .= q
    update!(structure)
    Jj = zeros(eltype(q),num_of_dim,num_of_full_coords)
    retǨm = zeros(eltype(q),num_of_full_coords,num_of_full_coords)
    foreach(structure.apparatuses) do appar
        if appar.joint isa CableJoint
            j = appar.id
            (;force,joint) = appar
            (;hen,egg) = joint.hen2egg
            rb1 = hen.body
            rb2 = egg.body
            ap1id = hen.pid
            ap2id = egg.pid
            C1 = rb1.cache.Cps[ap1id]
            C2 = rb2.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            (;state) = force
            (;length,) = state
            s = 1/length
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            Ūjq = Uj*q
            retǨm .+= k[j]*s^2*(Ūjq*transpose(Ūjq))
        end
    end
    retǨm
end

function build_material_stiffness_matrix!(structure::AbstractStructure,q,k)
    (;num_of_full_coords,bodyid2sys_full_coords) = structure.connectivity
    build_material_stiffness_matrix!(structure,num_of_full_coords,bodyid2sys_full_coords,q,k)
end

function build_material_stiffness_matrix!(structure::AbstractStructure{bType,aType,<:PresFreeConnectivity},q,k) where {bType,aType}
    (;num_of_full_coords,bodyid2sys_full_coords,sys_free_coords_idx) = structure.connectivity
    Km = build_material_stiffness_matrix!(structure,num_of_full_coords,bodyid2sys_full_coords,q,k)
    @view Km[sys_free_coords_idx,sys_free_coords_idx] 
end

function build_geometric_stiffness_matrix!(structure::AbstractStructure,num_of_full_coords,bodyid2sys_full_coords,q,f)
    (;num_of_dim) = structure
    structure.state.system.q .= q
    update!(structure)
    Jj = zeros(eltype(q),num_of_dim,num_of_full_coords)
    retǨg = zeros(eltype(q),num_of_full_coords,num_of_full_coords)
    foreach(structure.apparatuses) do appar
        if appar.joint isa CableJoint
            j = appar.id
            (;force,joint) = appar
            (;hen,egg) = joint.hen2egg
            rb1 = hen.body
            rb2 = egg.body
            ap1id = hen.pid
            ap2id = egg.pid
            C1 = rb1.cache.Cps[ap1id]
            C2 = rb2.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            (;state) = force
            (;length,) = state
            s = 1/length
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            Ǔj = Uj
            Ūjq = Uj*q
            retǨg .+= f[j]/length*(Ǔj-s^2*Ūjq*transpose(Ūjq))
        end
    end
    retǨg
end

function build_geometric_stiffness_matrix!(structure::AbstractStructure,q,f) 
    (;num_of_full_coords,bodyid2sys_full_coords) = structure.connectivity
    Kg = build_geometric_stiffness_matrix!(structure,num_of_full_coords,bodyid2sys_full_coords,q,f)
end

function build_geometric_stiffness_matrix!(structure::AbstractStructure{bType,aType,<:PresFreeConnectivity},q,f) where {bType,aType}
    (;num_of_full_coords,bodyid2sys_full_coords,sys_free_coords_idx) = structure.connectivity
    Kg = build_geometric_stiffness_matrix!(structure,num_of_full_coords,bodyid2sys_full_coords,q,f)
    @view Kg[sys_free_coords_idx,sys_free_coords_idx] 
end

"""
$(TYPEDSIGNATURES)
"""
function linearize(tginput,λ,u,q,q̇=zero(q))
    structure = deepcopy(tginput)
    set_restlen!(structure,u)
    reset_forces!(structure)
    distribute_q_to_rbs!(structure,q,q̇)
    update_cables_apply_forces!(structure)
    M = build_massmatrix(structure)
    A = build_A(structure)
    Q̃ = build_Q(structure)
    ∂L∂q,∂L∂q̇ = build_tangent(structure)
    (;ncoords,num_of_cstr) = structure
    nz = ncoords + num_of_cstr
    M̂ = zeros(eltype(q),nz,nz)
    Ĉ  = zeros(eltype(q),nz,nz)
    K̂ = zeros(eltype(q),nz,nz)
    M̂[1:ncoords,1:ncoords] .= M
    Ĉ[1:ncoords,1:ncoords] .= -Q̃*∂L∂q̇

    # fjac = test_fvector(structure,q)
    K̂[1:ncoords,1:ncoords] .= -Q̃*∂L∂q .+ cstr_forces_jacobian(structure,λ)
    Aq = A(q)
    c = maximum(abs.(K̂[1:ncoords,1:ncoords]))
    K̂[1:ncoords,ncoords+1:nz] .= c.*transpose(Aq)
    K̂[ncoords+1:nz,1:ncoords] .= c.*Aq
    M̂,Ĉ,K̂
end



function build_Y(bot)
    (;structure, hub) = bot
    (;actuators) = hub
    (;cables) = structure.apparatuses
    ncables = length(cables)
    nact = length(actuators)
    ret = spzeros(Int,ncables,nact)
    foreach(actuators) do actuator
        (;id,coupler,reg) = actuator
        if coupler isa Serial
            ret[actuator.reg.ids,id] .= 1
        elseif coupler isa Ganged
            is1,is2 = actuator.reg.ids
            ret[is1,id] =  1
            ret[is2,id] = -1
        else
            error("Unknown actuator type")
        end
    end
    ret
end