

"""
$(TYPEDEF)
"""
struct CableJoint{hen2eggType} <: AbstractJoint
    hen2egg::hen2eggType
    num_of_cstr::Int
end

function ApparatusCache(joint::CableJoint, force::DistanceSpringDamper, full_coords_idx, )
    T = get_numbertype(joint)
    nd = get_num_of_dims(joint)
    #note to be overwrite with prepare_cache!
    num_appar_full_idx = length(full_coords_idx)
    (; num_of_cstr) = joint

    φ = zeros(T, num_of_cstr)
    A = spzeros(T, num_of_cstr, num_appar_full_idx)
    ∂Aᵀλ∂q = zeros(T, num_appar_full_idx, num_appar_full_idx)
    ∂Aq̇∂q = zeros(T, num_of_cstr, num_appar_full_idx)
    Q = zeros(T, num_appar_full_idx)
    K = zeros(T, num_appar_full_idx, num_appar_full_idx)
    C = zeros(T, num_appar_full_idx, num_appar_full_idx)

    J = zeros(T, nd, num_appar_full_idx)
    DJ = similar(J)
    D = @MMatrix zeros(T, nd, nd)
    Im = SMatrix{nd,nd}(one(T)*I(nd))
    named = @eponymtuple(J, DJ, D, Im)

    ApparatusCache(
        φ,
        A,
        ∂Aᵀλ∂q,
        ∂Aq̇∂q,
        Q,
        K,
        C,
        named,
    )
end

function stretch!(appar::Apparatus{<:CableJoint,<:DistanceSpringDamper},c)
    (;joint,force) = appar
    κ, η, μ = c
    #todo κ stiffness, η damping
    force.state.restlen = μ
end

function get_numbertype(joint::CableJoint)
    get_numbertype(joint.hen2egg.hen.body)
end

function get_joint_idx(joint::CableJoint)
    (;hen,egg) = joint.hen2egg
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    ncoords_hen = get_num_of_coords(coords_hen)
    ncoords_egg = get_num_of_coords(coords_egg)
    full_idx = vcat(
        collect(1:ncoords_hen),
        collect(1:ncoords_egg) .+ ncoords_hen
    )
end

get_num_of_params(appar::Apparatus{<:CableJoint,<:DistanceSpringDamper},) = 3

function get_params!(params,appar::Apparatus{<:CableJoint,<:DistanceSpringDamper},)
    (;force) = appar
    params .= (k = force.k, c = force.c, restlen = force.state.restlen) |> collect
end

get_num_of_dims(joint::CableJoint) = get_num_of_dims(joint.hen2egg.hen.body)

function get_appar_idx(appar::Apparatus{<:CableJoint},bodyid2sys_full_coords)
    (;joint) = appar
    (;hen,egg) = joint.hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    full_hen = bodyid2sys_full_coords[id_hen]
    full_egg = bodyid2sys_full_coords[id_egg]
    appar_full_coords_idx = vcat(
        full_hen,
        full_egg
    )
    appar_full_coords_idx
end

function prepare_cache!(appar::Apparatus{<:CableJoint},cnt::AbstractConnectivity)
    T = get_numbertype(appar)
    nq = cnt.num_of_full_coords
    nd = get_num_of_dims(appar)
    (;cache, ) = appar

    if (size(cache.K, 1) != nq) || (size(cache.K, 2) != nq)
        cache.K = spzeros(T, nq, nq)
    else
        zero!(cache.K)
    end
    if (size(cache.C, 1) != nq) || (size(cache.C, 2) != nq)
        cache.C = spzeros(T, nq, nq)
    else
        zero!(cache.C)
    end

    named = cache.named
    J = named.J
    D = named.D
    if (size(J, 1) != nd) || (size(J, 2) != nq) || (size(D, 1) != nd) || (size(D, 2) != nd)
        cache.named = @eponymtuple(
            J = zeros(T, nd, nq),
            DJ = zeros(T, nd, nq),
            D = (@MMatrix zeros(T, nd, nd)),
            Im = SMatrix{nd,nd}(one(T)*I(nd))
        )
    end
end

@inline function fill_cable_J!(J, C_hen, C_egg, sys_full_idx_hen, sys_full_idx_egg, nd)
    fill!(J, zero(eltype(J)))
    @inbounds for (local_idx, global_idx) in pairs(sys_full_idx_hen)
        for row = 1:nd
            J[row, global_idx] -= C_hen[row, local_idx]
        end
    end
    @inbounds for (local_idx, global_idx) in pairs(sys_full_idx_egg)
        for row = 1:nd
            J[row, global_idx] += C_egg[row, local_idx]
        end
    end
end

function add_tangent_stiffness_matrix!(∂Q∂q,appar::Apparatus{<:CableJoint},st::AbstractStructure,q,s=nothing)
    
    cnt = st.connectivity
    (;
        bodyid2sys_full_coords,
    ) = cnt
    (;joint, force, cache) = appar
    (;J,DJ,D,Im) =  cache.named
    T = eltype(D)

    num_of_dim = get_num_of_dims(st)
    (;hen,egg) = joint.hen2egg
    body_hen = hen.body
    body_egg = egg.body
    C_hen = body_hen.cache.Cps[hen.pid]
    C_egg = body_egg.cache.Cps[egg.pid]

    sys_full_idx_hen = bodyid2sys_full_coords[hen.body.prop.id]
    sys_full_idx_egg = bodyid2sys_full_coords[egg.body.prop.id]

    (;k,c,state,slack) = force
    (;direction,tension) = state
    l = state.length
    l̇ = state.lengthdot
    if !(slack && (tension==0))
        density = tension/l
        β = c*l̇/l + density
        D .= (k - β) .* (direction * direction') .+ β .* Im
        fill_cable_J!(J, C_hen, C_egg, sys_full_idx_hen, sys_full_idx_egg, num_of_dim)
        mul!(DJ, D, J)
        mul!(∂Q∂q, transpose(J), DJ, -one(T), one(T))
    end
    nothing
end


function add_tangent_damping_matrix!(∂Q∂q̇,appar::Apparatus{<:CableJoint},st::AbstractStructure,q)
    
    cnt = st.connectivity
    (;
        bodyid2sys_full_coords,
    ) = cnt
    (;joint, force, cache) = appar
    (;J,DJ,D) = cache.named
    T = eltype(D)
    num_of_dim = get_num_of_dims(st)

    (;hen,egg) = joint.hen2egg
    body_hen = hen.body
    body_egg = egg.body
    C_hen = body_hen.cache.Cps[hen.pid]
    C_egg = body_egg.cache.Cps[egg.pid]
    sys_full_idx_hen = bodyid2sys_full_coords[body_hen.prop.id]
    sys_full_idx_egg = bodyid2sys_full_coords[body_egg.prop.id]
    (;c,state,slack) = force
    (;direction,tension) = state
    if !(slack && (tension == 0))
        D .= c .* (direction * direction')
        fill_cable_J!(J, C_hen, C_egg, sys_full_idx_hen, sys_full_idx_egg, num_of_dim)
        mul!(DJ, D, J)
        mul!(∂Q∂q̇, transpose(J), DJ, -one(T), one(T))
    end
    nothing
end

function update_apparatus!(st::AbstractStructure, appar::Apparatus{<:CableJoint}, s)
    (;id,joint,force) = appar
    (;hen,egg) = joint.hen2egg
    locus_state_hen = hen.body.state.loci_states[hen.pid]
    locus_state_egg = egg.body.state.loci_states[egg.pid]
    p_hen = locus_state_hen.frame.position
    ṗ_hen = locus_state_hen.frame.velocity
    p_egg = locus_state_egg.frame.position
    ṗ_egg = locus_state_egg.frame.velocity
    update!(force,p_hen,p_egg,ṗ_hen,ṗ_egg)
end

function apply_force_to_bodies!(structure::AbstractStructure, appar::Apparatus{<:CableJoint})
    (;joint,force) = appar
    (;hen,egg) = joint.hen2egg
    locus_state_hen = hen.body.state.loci_states[hen.pid]
    locus_state_egg = egg.body.state.loci_states[egg.pid]
    f_hen = locus_state_hen.force
    f_egg = locus_state_egg.force
    f_hen .+= force.state.force
    f_egg .-= force.state.force
end

function execute!(structure::AbstractStructure,actuator::RegisterActuator{<:Vector{<:NamedTuple{
        (:cable,:val_id),}},
        <:FunctionOperator}, u)

    (;id, signifier, operator,register) = actuator
    (;func!, func_vals) = operator
    (;values) = register
    func!(func_vals,values,u)
    foreach(signifier) do sig
        (;cable,val_id) = sig
        cable.force.state.restlen = func_vals[val_id]
    end
end

function gen_force_para_jacobian!(∂Q̌∂c,appar::Apparatus{<:CableJoint},st::AbstractStructure,q,q̇,t,s) 
    (;
        bodyid2sys_full_coords,
    ) = st.connectivity
    (;
        hen2egg,
    ) = appar.joint
    (;k,) = appar.force
    (;direction,length,restlen) =  appar.force.state
    (;hen,egg) = hen2egg
    C_hen = hen.body.cache.Cps[hen.pid]
    C_egg = egg.body.cache.Cps[egg.pid]
    sys_full_hen = bodyid2sys_full_coords[hen.body.prop.id]
    sys_full_egg = bodyid2sys_full_coords[egg.body.prop.id]
    # f = k*(l-μ)
    ∂f∂k = length-restlen
    ∂f∂η = 0.0
    ∂f∂μ = -k
    ∂f∂c = [∂f∂k ∂f∂η ∂f∂μ;]
    d∂f∂c = direction*∂f∂c
    ∂Q̌∂c[sys_full_egg,:] .-= transpose(C_egg)*d∂f∂c
    ∂Q̌∂c[sys_full_hen,:] .+= transpose(C_hen)*d∂f∂c
end

function gen_force_actu_jacobian!(∂F∂u, structure::AbstractStructure,actuator::RegisterActuator{<:Vector{<:NamedTuple{
        (:cable,:val_id),}},
        <:FunctionOperator},u)
    (;
        bodyid2sys_full_coords,
    ) = structure.connectivity
    (;id, signifier, operator,register) = actuator
    (;jac!, Jac_vals) = operator
    (;values) = register
    jac!(Jac_vals,values,u)
    ∂μ∂u = Jac_vals
    ndir = get_num_of_dims(structure)
    foreach(signifier) do sig
        (;cable,val_id) = sig
        (;
            hen2egg,
        ) = cable.joint
        (;k,) = cable.force
        (;direction,length,restlen) =  cable.force.state
        (;hen,egg) = hen2egg
        C_hen = hen.body.cache.Cps[hen.pid]
        C_egg = egg.body.cache.Cps[egg.pid]
        sys_full_hen = bodyid2sys_full_coords[hen.body.prop.id]
        sys_full_egg = bodyid2sys_full_coords[egg.body.prop.id]
        # f = k*(l-μ)
        ∂f∂μ = -k
        μrow = @view ∂μ∂u[val_id, :]
        @inbounds for j = axes(∂F∂u, 2)
            s = ∂f∂μ * μrow[j]
            if s == 0
                continue
            end
            for (local_idx, gidx) in enumerate(sys_full_egg)
                c = C_egg[1, local_idx] * direction[1]
                for d = 2:ndir
                    c += C_egg[d, local_idx] * direction[d]
                end
                ∂F∂u[gidx, j] -= c * s
            end
            for (local_idx, gidx) in enumerate(sys_full_hen)
                c = C_hen[1, local_idx] * direction[1]
                for d = 2:ndir
                    c += C_hen[d, local_idx] * direction[d]
                end
                ∂F∂u[gidx, j] += c * s
            end
        end
    end
end
