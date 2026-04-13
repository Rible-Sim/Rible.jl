
function cstr_forces_jacobian(st::AbstractStructure,inst_state::AbstractCoordinatesState)
    (;q,λ) = inst_state
    cstr_forces_jacobian(st,q,λ)
end

function cstr_forces_jacobian!(ret, st::AbstractStructure, inst_state::AbstractCoordinatesState)
    cstr_forces_jacobian!(ret, st, inst_state.q, inst_state.λ)
end

function cstr_forces_jacobian(st::AbstractStructure,q,λ)
    num_of_coords = get_num_of_coords(st)
    ret = zeros(eltype(λ),num_of_coords,num_of_coords)
    cstr_forces_jacobian!(ret, st,q,λ)
    ret
end

function cstr_forces_jacobian!(ret, st::AbstractStructure,q,λ)
    (;bodies,apparatuses) = st
    cnt = st.connectivity
    (;
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_full_coords,
        num_of_cstr,
        num_of_full_coords,
        num_of_extrinsic_cstr,
        apparid2sys_extrinsic_cstr_idx
    ) = cnt
    fill!(ret, zero(eltype(λ)))
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        add_cstr_forces_jacobian!((@view ret[memfull,memfull]),
            body.coords,
            (@view λ[memincst])
        )
    end
    #todo skip 2D for now
    if get_num_of_dims(st) == 3
        foreach(apparatuses) do appar
            cstr_idx = apparid2sys_extrinsic_cstr_idx[appar.id]
            appar_sys_full_idx = cnt.apparid2sys_full_coords_idx[appar.id]
            add_cstr_forces_jacobian!((@view ret[appar_sys_full_idx,appar_sys_full_idx]), 
                appar,
                cnt,q,
                (@view λ[cstr_idx])
            )
        end
    end
    ret
end

function cstr_velocity_jacobian!(ret, st::AbstractStructure, inst_state::AbstractCoordinatesState)
    (;bodies,apparatuses,connectivity) = st
    (;
        num_of_cstr,
        num_of_full_coords,
        bodyid2sys_full_coords,
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        apparid2sys_full_coords_idx,
        apparid2sys_extrinsic_cstr_idx
    ) = connectivity
    (;q,q̇) = inst_state
    @assert size(ret, 1) == num_of_cstr && size(ret, 2) == num_of_full_coords
    fill!(ret, zero(eltype(q̇)))
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        if !isempty(memincst)
            cstr_velocity_jacobian!(
                (@view ret[memincst,memfull]),
                body,
                @view q̇[memfull]
            )
        end
    end
    if get_num_of_dims(st) == 3
        foreach(apparatuses) do appar
            cstr_idx = apparid2sys_extrinsic_cstr_idx[appar.id]
            appar_sys_full_idx = apparid2sys_full_coords_idx[appar.id]
            cstr_velocity_jacobian!(
                (@view ret[cstr_idx,appar_sys_full_idx]),
                appar, connectivity, q, q̇
            )
        end
    end
    ret
end

function cstr_velocity_jacobian!(ret, st::AbstractStructure, q, q̇)
    inst_state = st.state.system
    inst_state.q .= q
    inst_state.q̇ .= q̇
    cstr_velocity_jacobian!(ret, st, inst_state)
end

function cstr_velocity_jacobian(st::AbstractStructure, inst_state::AbstractCoordinatesState)
    (;num_of_cstr,num_of_full_coords) = st.connectivity
    ret = zeros(eltype(inst_state.q̇), num_of_cstr, num_of_full_coords)
    cstr_velocity_jacobian!(ret, st, inst_state)
end

function cstr_velocity_jacobian(st::AbstractStructure, q, q̇)
    (;num_of_cstr,num_of_full_coords) = st.connectivity
    ret = zeros(eltype(q̇), num_of_cstr, num_of_full_coords)
    cstr_velocity_jacobian!(ret, st, q, q̇)
    ret
end

function intrinsic_nullspace(st::AbstractStructure,)
    intrinsic_nullspace(st,st.state.system)
end

function intrinsic_nullspace(st::AbstractStructure,inst_state)
    (;bodies,connectivity) = st
    (;q) = inst_state
    cnt = connectivity
    (;num_of_full_coords,bodyid2sys_full_coords,num_of_dof_unconstrained,bodyid2sys_dof_idx,) = cnt
    ret = zeros(eltype(q),num_of_full_coords,num_of_dof_unconstrained)
    foreach(bodies) do body
        bodyid = body.prop.id
        mem2full = bodyid2sys_full_coords[bodyid]
        ret[mem2full,bodyid2sys_dof_idx[bodyid]] = nullspace_mat(body.coords,q[mem2full])
    end
    ret
end

function extrinsic_nullspace(st,inst_state)
    (;bodies,connectivity,num_of_dof) = st
    cnt = connectivity
    (;q) = inst_state
    (;num_of_full_coords,bodyid2sys_full_coords,num_of_dof_unconstrained,) = cnt
    ret = zeros(eltype(q),num_of_dof_unconstrained,num_of_dof)
    bodyid2sys_dof_idx = deepcopy(cnt.bodyid2sys_dof_idx)
    nullspaces = Vector{Matrix{eltype(q)}}(undef,2)
    foreach(apparatuses) do appar
        apparid = appar.id
        if appar isa PrototypeJoint
            nullspaces[1] = zeros(eltype(q),num_of_dof_unconstrained,7)
            for i in bodyid2sys_dof_idx[1]
                nullspaces[1][i,i] = 1
            end
            (;hen,egg) = appar.hen2egg
            coords_hen = hen.body.coords
            coords_egg = egg.body.coords
            id_hen = hen.body.prop.id
            id_egg = egg.body.prop.id
            q_hen = q[bodyid2sys_full_coords[id_hen]]
            q_egg = q[bodyid2sys_full_coords[id_egg]]
            c_hen = to_local_coords(coords_hen,hen.position)
            c_egg = to_local_coords(coords_egg,egg.position)
            r_hen =  to_position_jacobian(coords_hen,q_hen,c_hen)*q_hen
            r_egg =  to_position_jacobian(coords_egg,q_egg,c_egg)*q_egg
            ro_hen = to_position_jacobian(coords_hen,q_hen,zero(c_hen))*q_hen
            ro_egg = to_position_jacobian(coords_egg,q_egg,zero(c_egg))*q_egg

            # hen.body.prop.loci[hen.trlid].axes
            # egg.body.prop.loci[egg.trlid].axes
            # hen.body.prop.loci[hen.rotid].axes
            invX̄_hen = coords_hen.data.invX̄
            invX̄_egg = coords_egg.data.invX̄
            axes_rot_egg = invX̄_egg*egg.rot_axes
            axes_idx = [
                (3,1), #normal * tangent
                (3,2), #normal * bitangent
                (1,2)  #tangent * bitangent
            ]
            id_axis_egg = 3
            X_hen = NCF.get_X(coords_hen,q_hen)
            X_egg = NCF.get_X(coords_egg,q_egg)
            axis_egg = X_egg*axes_rot_egg.X
            normal = axis_egg[:,id_axis_egg]
            I3 = I(3)
            O3 = zero(I3)
            o3 = O3[:,1]
            nullspaces[1][bodyid2sys_dof_idx[id_egg],vcat(bodyid2sys_dof_idx[id_hen],7)] = [
                I3 skew((r_egg-ro_egg)-(r_hen-ro_hen)) (r_egg-ro_egg)×normal;
                O3             I3                      normal
            ]
        else
            if apparid == 2
                nullspaces[2] = zeros(eltype(q),7,3)
                nullspaces[2][3,1] = 1
                nullspaces[2][4,2] = 1
                nullspaces[2][7,3] = 1
            end
        end
    end
    ret .= nullspaces[1]*nullspaces[2]
end


# In-place ∂Q∂q for cables and flexible bodies
function add_tangent_stiffness_matrix!(∂Q∂q, body::AbstractBody, st::AbstractStructure, q, s=nothing)
    add_tangent_stiffness_matrix!(∂Q∂q, bodyy, st, NoField(), q, s)
end

function add_tangent_stiffness_matrix!(∂Q∂q,body::AbstractBody, st::AbstractStructure, field::AbstractField, q,s=nothing)
    field_jacobian!(∂Q∂q, body, field, q)
end

function add_tangent_stiffness_matrix!(∂Q∂q,body::FlexibleBody,st::AbstractStructure,q,s=nothing)
    (;coords,cache) = body
    (;e,) = cache
    ancs = coords
    ∂Q∂e = ANCF.make_∂Q∂e(ancs)(e)
    mfull = bodyid2sys_full_coords[body.prop.id]
    (;full_idx) = body.coords
    ∂Q∂q[mfull,mfull] .-= ∂Q∂e[full_idx,full_idx]
end

function add_tangent_stiffness_matrix!(∂Q∂q,appar::Apparatus,st::AbstractStructure,q,s=nothing)
    #default do 
    nothing
end


function add_tangent_stiffness_matrix!(∂Q∂q,appar::Apparatus{JointType,<:TorsionalSpringDamper},st::AbstractStructure,q,s=nothing) where {JointType}
    (;
        bodyid2sys_full_coords,
    ) = st.connectivity
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;k,state) = appar.force
    (;
        relative_core
    ) = cache
    (;hen,egg) = hen2egg
    coords_egg = egg.body.coords
    id_egg = egg.body.prop.id
    egg_sys_full_idx = bodyid2sys_full_coords[id_egg]
    invX̄_egg = coords_egg.data.invX̄
    q_egg = @view q[egg_sys_full_idx]
    state.angle = s[1]
    deformation = state.angle - state.rest_angle
    state.torque = -k*deformation
    loci_axes_rot_egg = egg.rot_axes
    axes_rot_egg = invX̄_egg*loci_axes_rot_egg

    r̄_home = egg.position
    r̄_away = r̄_home + loci_axes_rot_egg.X[:,1]
    C_home = to_position_jacobian(coords_egg,q_egg,to_local_coords(coords_egg,r̄_home))
    C_away = to_position_jacobian(coords_egg,q_egg,to_local_coords(coords_egg,r̄_away))
    
    ∂F∂q = state.torque*(-transpose(C_home) + transpose(C_away))*kron(hcat(0,transpose(axes_rot_egg.X[:,2])),I(3))
    ∂Q∂q[egg_sys_full_idx,egg_sys_full_idx] .+= ∂F∂q
    nothing
end

function add_tangent_stiffness_matrix!(∂Q∂q,appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper},st::AbstractStructure,q,s=nothing)
    (;
        bodyid2sys_full_coords,
    ) = st.connectivity
    (;joint,force) = appar
    body = joint.body
    id_body = body.prop.id
    coords_body = body.coords
    body_sys_full_idx = bodyid2sys_full_coords[id_body]
    
    (;k,state) = force
    state.angle = s[1]
    deformation = state.angle - state.rest_angle
    state.torque = -k*deformation # torque on body
    
    u,v = NCF.get_uv(body.coords,q)

    # F = torque * (v[1] * ∂x∂q + v[2] * ∂y∂q)
    # ∂F/∂q = torque * (∂x∂q' ⊗ ∂v[1]/∂q + ∂y∂q' ⊗ ∂v[2]/∂q)
    # Note: ∂x∂q and ∂y∂q are constant vectors (rows of C)
    # For NC2D4C: v = [-uy, ux], so ∂v/∂q = [-∂y∂q; ∂x∂q]
    # For NC2D6C: v is explicit, so ∂v/∂q terms come from rows 5,6 of C

    num_of_dim = NCF.get_num_of_dims(body.coords)
    select_uvw = body.coords.conversion_to_X
    ∂x∂q = @view select_uvw[num_of_dim+1, :]
    ∂y∂q = @view select_uvw[num_of_dim+2, :]

    ∂vx∂q = -∂y∂q
    ∂vy∂q =  ∂x∂q

    # ∂Q∂q += state.torque * (∂x∂q * ∂vx∂q' + ∂y∂q * ∂vy∂q')
    # Using outer product addition
    @views ∂Q∂q[body_sys_full_idx,body_sys_full_idx] .+= state.torque .* (∂x∂q .* ∂vx∂q' .+ ∂y∂q .* ∂vy∂q')
end

function gen_force_auxi_jacobian!(∂Q∂s,appar::Apparatus,st::AbstractStructure,q,q̇,t, s=nothing)
    #default do nothing
end

function gen_force_auxi_jacobian!(∂Q∂s,appar::Apparatus{JointType,<:TorsionalSpringDamper},st::AbstractStructure,q,q̇,t,s) where {JointType}
    (;
        bodyid2sys_full_coords,
    ) = st.connectivity
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;k,state) = appar.force
    (;hen,egg) = hen2egg
    coords_egg = egg.body.coords
    id_egg = egg.body.prop.id
    egg_sys_full_idx = bodyid2sys_full_coords[id_egg]
    invX̄_egg = coords_egg.data.invX̄
    q_egg = @view q[egg_sys_full_idx]
    state.angle = s[1]
    deformation = state.angle - state.rest_angle
    state.torque = -k*deformation
    loci_axes_rot_egg = egg.rot_axes
    axes_rot_egg = invX̄_egg*loci_axes_rot_egg

    r̄_home = egg.position
    r̄_away = r̄_home + loci_axes_rot_egg.X[:,1]
    C_home = to_position_jacobian(coords_egg,q_egg,to_local_coords(coords_egg,r̄_home))
    C_away = to_position_jacobian(coords_egg,q_egg,to_local_coords(coords_egg,r̄_away))
    
    ∂F∂s = -k*(-transpose(C_home) + transpose(C_away))*axes_rot_egg.X[:,2]
    ## @show ∂F∂s
    ∂Q∂s[egg_sys_full_idx,:] .+= ∂F∂s
end

function gen_force_auxi_jacobian!(∂Q∂s,appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper},st::AbstractStructure,q,q̇,t,s)
    (;
        bodyid2sys_full_coords,
    ) = st.connectivity
    (;joint,force) = appar
    body = joint.body 
    (;k,state) = force
    deformation = state.angle - state.rest_angle
    
    num_of_dim = NCF.get_num_of_dims(body.coords)
    select_uvw = body.coords.conversion_to_X
    ∂x∂q = @view select_uvw[num_of_dim+1, :]
    ∂y∂q = @view select_uvw[num_of_dim+2, :]
    
    indices = bodyid2sys_full_coords[body.prop.id]
    q_body = @view q[indices]
    
    u,v = NCF.get_uv(body.coords, q_body)
    
    # F = torque * v (Geometric couple)
    # ∂F/∂s = ∂torque/∂s * v = -k * v (mapped to q)
    # This matches ProtoJoint implementation (see add_tangent_stiffness_matrix!)
    # No singularity at Ss=0.
    @views ∂Q∂s[indices, 1] .+= (-k) .* (v[1] .* ∂x∂q .+ v[2] .* ∂y∂q)
end

"""
    gen_force_para_jacobian!(∂Q∂c,appar::Apparatus,st,q,q̇,t, s=nothing)

Compute the Jacobian of the generalized force w.r.t. the structual parameters.
"""
function gen_force_para_jacobian!(∂Q∂c,appar::Apparatus,st,q,q̇,t, s=nothing)
    #default do nothing
end


# In-place ∂Q∂q for cables and flexible bodies
function add_tangent_stiffness_matrix!(∂Q∂q, st::AbstractStructure, field=NoField())
    (;bodies,apparatuses,connectivity) = st
    (;
        bodyid2sys_full_coords,
        apparid2sys_aux_var_idx
    ) = connectivity
    (;q,s) = st.state.system
    foreach(bodies) do body
        add_tangent_stiffness_matrix!(
            (@view ∂Q∂q[bodyid2sys_full_coords[body.prop.id], bodyid2sys_full_coords[body.prop.id]]), 
            body, st, field, 
            (@view q[bodyid2sys_full_coords[body.prop.id]])
        )
    end
    foreach(apparatuses) do appar
        add_tangent_stiffness_matrix!(
            ∂Q∂q,
            appar,
            st,
            q,
            (@view s[apparid2sys_aux_var_idx[appar.id]])
        )
    end
end

# In-place ∂Q∂q̇ for cables and flexible bodies
function add_tangent_damping_matrix!(∂Q∂q̇,body::AbstractBody,st::AbstractStructure,q)
    #default do nothing
end

function add_tangent_damping_matrix!(∂Q∂q̇,appar::Apparatus,st::AbstractStructure,q)
    #default do 
    nothing
end


# In-place ∂Q∂q̇ for cables
function add_tangent_damping_matrix!(∂Q∂q̇, st::AbstractStructure, field=NoField())
    (;apparatuses,) = st
    (;q,) = st.state.system
    foreach(apparatuses) do appar
        add_tangent_damping_matrix!(∂Q∂q̇, appar, st , q)
    end
    nothing
end


function build_Ǩ(st)
    _,λ = check_static_equilibrium_output_multipliers(st)
    q = get_coords(st)
    build_Ǩ(st,q,λ)
end

function build_Ǩ(st,q,λ)
    (;num_of_free_coords) = st.connectivity
    T = get_numbertype(st)
    ∂Q∂q = zeros(T,num_of_free_coords,num_of_free_coords)
    add_tangent_stiffness_matrix!(∂Q∂q,st)
    Ǩa = -cstr_forces_jacobian(st,q,λ)
    #note ∂Q∂q = -(Ǩm+Ǩg)
    Ǩ = - ∂Q∂q .+ Ǩa
end

"""
$(TYPEDSIGNATURES)
"""
function find_nullspace(c)
    Nc,Nu = size(c)
    P = VectorOfArray([ begin
                            p = zeros(eltype(c),Nu)
                            p[i] = 1
                            p
                        end for i = 1:Nu])
    for i = 1:Nc
        # @show i
        cPⁱ⁻¹ = c*P
        rowi = cPⁱ⁻¹[i,:]
        _ℓⁱ = findall((x)->!iszero(x),rowi)
        ℓⁱ = _ℓⁱ[sortperm(rowi[_ℓⁱ],order=Base.Order.Reverse)]
        if isempty(ℓⁱ)
            continue
        end
        # display(cPⁱ⁻¹)
        # @show ℓⁱ
        # for k = 1:(length(ℓ)-1)
        for k = 1:(length(ℓⁱ)-1)
            αⁱₖ = -cPⁱ⁻¹[i,ℓⁱ[k]]./cPⁱ⁻¹[i,ℓⁱ[k+1]]
            P[:,ℓⁱ[k]] .= P[:,ℓⁱ[k]] + αⁱₖ*P[:,ℓⁱ[k+1]]
        end
        deleteat!(P.u,ℓⁱ[end])
    end
    Array(P)
end

