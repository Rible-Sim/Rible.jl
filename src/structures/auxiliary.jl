
function prepare_cache!(appar::Apparatus,cnt::AbstractConnectivity)
end


"""
Return the auxiliary variables of the entire structure.
$(TYPEDSIGNATURES)
"""
function auxi_function(structure::AbstractStructure, inst_state::AbstractCoordinatesState)
    (; num_of_aux_var) = structure.connectivity
    ret = Vector{eltype(inst_state.q)}(undef,num_of_aux_var)
    auxi_function!(ret, structure, inst_state)
    ret
end

"""
Compute and store the auxiliary variables of the entire structure into `ret`.
$(TYPEDSIGNATURES)
"""
function auxi_function!(ret, structure::AbstractStructure, inst_state::AbstractCoordinatesState)
    (;q,s) = inst_state
    (;apparatuses,connectivity) = structure
    cnt = connectivity
    (;  
        apparid2sys_aux_var_idx,
    ) = cnt
    
    lazy_update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)

    foreach(apparatuses) do appar
        appar_aux_var_idx = apparid2sys_aux_var_idx[appar.id]
        auxi_function!(
            (@view ret[appar_aux_var_idx]),
            appar,cnt,
            q,
            (@view  s[appar_aux_var_idx]),
        )
    end
    nothing
end

"""
Compute and store the Jacobian matrices of the structure's auxiliary function.
$(TYPEDSIGNATURES)
"""
function auxi_jacobian!(∂S∂q,∂S∂s,structure::AbstractStructure,inst_state::AbstractCoordinatesState)
    (;q,s) = inst_state
    (;apparatuses,connectivity) = structure
    cnt = connectivity
    (; 
        apparid2sys_full_coords_idx,
        apparid2sys_aux_var_idx
    ) = cnt
    
    lazy_update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)
    
    foreach(apparatuses) do appar
        appar_aux_var_idx = apparid2sys_aux_var_idx[appar.id]
        appar_sys_full_coords_idx = apparid2sys_full_coords_idx[appar.id]
        auxi_jacobian!(
            (@view ∂S∂q[appar_aux_var_idx,:]),
            (@view ∂S∂s[appar_aux_var_idx,appar_aux_var_idx]),
            appar,cnt,
            q,
            (@view  s[appar_aux_var_idx]),
        )
    end
end

function auxi_function!(ret, appar::Apparatus,cnt::AbstractConnectivity,q, s)
    #default do nothing
end

function auxi_function!(ret, appar::Apparatus{<:PrototypeJoint,<:TorsionalSpringDamper},
        cnt::AbstractConnectivity, q, s
    )
    (;bodyid2sys_full_coords) = cnt
    (;joint,force) = appar
    (;
        hen2egg,
        violations,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = joint
    (;
        relative_core
    ) = cache
    ## (;

    ## ) = force
    num_of_aux_var = get_num_of_aux_var(force)
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    # cstr violations
    # @show joint.id
    
    X_hen = NCF.get_X(coords_hen,q_hen)
    X_egg = NCF.get_X(coords_egg,q_egg)
    invX̄_hen = coords_hen.data.invX̄
    invX̄_egg = coords_egg.data.invX̄
    axes_rot_hen = invX̄_hen*hen.rot_axes
    axes_rot_egg = invX̄_egg*egg.rot_axes
    R = (transpose((X_hen*axes_rot_hen).X)*X_egg*axes_rot_egg).X
    θ = s[1]
    
    auxi_functions = [
        (θ,R) -> sin(θ)-R[2,1],
        (θ,R) -> cos(θ)-R[1,1],
        (θ,R) -> sin(θ)-R[2,1],
        (θ,R) -> cos(θ)-R[1,1]
    ]

    ret[1] = auxi_functions[force.state.mode_id](θ,R)
    ## c1 = [1.0,0,0]×(R[:,1]/norm(R[:,1]))
    ## c2 = [0,1.0,0]×(R[:,2]/norm(R[:,2]))
    ## c3 = [0,0,1.0]×(R[:,3]/norm(R[:,3]))
    ## asin.( [0,0,c3[1]])
    ## @show R[2,2:3], rad2deg(θ)
end


"""
$(TYPEDSIGNATURES)

Compute and store the Jacobian matrices w.r.t. the free coordinates and auxiliary variables of the auxiliary function of an apparatus.

# Arguments
- `∂S∂q`: Output matrix for partial derivatives w.r.t. free coordinates
- `∂S∂s`: Output matrix for partial derivatives w.r.t. auxiliary variables
- `appar`: The apparatus 
- `structure`: The mechanical structure
- `q`: Current generalized coordinates
- `s`: Current auxiliary variables
"""
function auxi_jacobian!(∂S∂q,∂S∂s, appar::Apparatus, cnt::AbstractConnectivity, q, s)
    #default do nothing
end


function auxi_jacobian!(∂S∂q,∂S∂s, appar::Apparatus{<:PrototypeJoint,<:TorsionalSpringDamper},
        cnt::AbstractConnectivity, q, s
    )
    (;bodyid2sys_full_coords,apparid2sys_full_coords_idx) = cnt
    (;joint,force) = appar
    (;
        hen2egg,
        violations,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = joint
    (;
        relative_core
    ) = cache
    full_idx = appar.full_coords_idx
    num_of_aux_var = get_num_of_aux_var(force)
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q_appar = vcat(q_hen,q_egg)
    # cstr violations
    # @show joint.id

    X_hen = NCF.get_X(coords_hen,q_hen)
    X_egg = NCF.get_X(coords_egg,q_egg)
    invX̄_hen = coords_hen.data.invX̄
    invX̄_egg = coords_egg.data.invX̄
    axes_rot_hen = invX̄_hen*hen.rot_axes
    axes_rot_egg = invX̄_egg*egg.rot_axes
    select_uvw_hen = BlockDiagonal([coords_hen.conversion_to_X,zero(coords_egg.conversion_to_X)])
    select_uvw_egg = BlockDiagonal([zero(coords_hen.conversion_to_X),coords_egg.conversion_to_X])
    
    ## R = (transpose((X_hen*axes_rot_hen).X)*X_egg*axes_rot_egg).X

    axis_zero = zero(axes_rot_hen.X[:,1])
    
    half_11 =  select_uvw_hen'*
            kron(vcat(0,axes_rot_hen.X[:,1],0,axis_zero),I(3))*
            kron(vcat(0,axis_zero,0,axes_rot_egg.X[:,1]),I(3))'*
            select_uvw_egg
    half_21 = select_uvw_hen'*
            kron(vcat(0,axes_rot_hen.X[:,2],0,axis_zero),I(3))*
            kron(vcat(0,axis_zero,0,axes_rot_egg.X[:,1]),I(3))'*
            select_uvw_egg
    hess_11 = Symmetric(half_11 + half_11')
    hess_21 = Symmetric(half_21 + half_21')
    
    θ = s[1] 
    if force.state.mode_id == 1
        ## @show θ0, "sin"
        ∂S∂q[1,apparid2sys_full_coords_idx[appar.id]] .= -(q_appar'*hess_21)[full_idx] #R[2,1]
        ∂S∂s[1,1]  = cos(θ)
    elseif force.state.mode_id == 2
        ## @show θ0, "cos"
        ∂S∂q[1,apparid2sys_full_coords_idx[appar.id]] .= -(q_appar'*hess_11)[full_idx] #R[1,1]
        ∂S∂s[1,1]  = sin(θ)
    elseif force.state.mode_id == 3
        ## @show θ0, "-sin"
        ∂S∂q[1,apparid2sys_full_coords_idx[appar.id]] .= -(q_appar'*hess_21)[full_idx] #R[2,1]
        ∂S∂s[1,1]  = cos(θ)
    elseif force.state.mode_id == 4
        ## @show θ0, "-cos"
        ∂S∂q[1,apparid2sys_full_coords_idx[appar.id]] .= -(q_appar'*hess_11)[full_idx] #R[1,1]
        ∂S∂s[1,1]  = sin(θ)
    else
        error("Invalid case id: $(force.mode_id)")
    end

    ## c1 = [1.0,0,0]×(R[:,1]/norm(R[:,1]))
    ## c2 = [0,1.0,0]×(R[:,2]/norm(R[:,2]))
    ## c3 = [0,0,1.0]×(R[:,3]/norm(R[:,3]))
    ## asin.( [0,0,c3[1]])
end

"""
    switch!(st::AbstractStructure, inst_state::AbstractCoordinatesState)

Switch the representation of a structure's apparatuses based on the current configuration `q` and auxiliary variables `s`.

The representation include the auxilary functions and hyper parameters of the apparatuses, but does **not** include the state of the apparatuses.

This function is usually used to prepare the structure for the next time step.

# Arguments
- `st::AbstractStructure`: The structure to update
- `q`: Current configuration coordinates 
- `s`: Current auxiliary variables

"""
function switch!(st::AbstractStructure,inst_state::AbstractCoordinatesState)
    (;apparatuses,connectivity) = st
    cnt = connectivity
    (;  
        num_of_full_coords,
        num_of_aux_var,
        apparid2sys_aux_var_idx
    ) = cnt
    
    clear_forces!(st)
    lazy_update_bodies!(st, inst_state)
    update_apparatuses!(st, inst_state)
    (;q, s) = inst_state

    foreach(apparatuses) do appar
        appar_aux_var_idx = apparid2sys_aux_var_idx[appar.id]
        switch!(
            appar,st,
            q, (@view  s[appar_aux_var_idx]),
        )
    end
end

"""
    switch!(appar::Apparatus, structure::AbstractStructure, q, s)

Switch the representation of an apparatus based on the current configuration `q` and auxiliary variables `s`.

The representation include the auxilary functions and hyper parameters of the apparatuses, but does **not** include the state of the apparatuses.

This function is usually used to prepare the structure for the next time step.

# Arguments
- `appar::Apparatus`: The apparatus to update
- `structure::AbstractStructure`: The parent structure
- `q`: Current configuration coordinates
- `s`: Current auxiliary variables for this apparatus
"""
function switch!(appar::Apparatus,structure::AbstractStructure, q, s)
    #default do 
    nothing
end


function switch!(appar::Apparatus{<:PrototypeJoint,<:TorsionalSpringDamper},
        structure::AbstractStructure,
        q, s
    )
    (;joint, force) = appar
    
    ## @show R
    θ = s[1]
    force.state.mode_id = angle2mode(θ)
end

function get_auxilary(appar::Apparatus) 
    T = get_numbertype(appar)
    zeros(T,appar.num_of_aux_var)
end

function get_auxilary(appar::Apparatus{<:AbstractJoint,<:TorsionalSpringDamper})
    (;force) = appar
    get_auxilary(force)
end

function auxi_function!(ret, appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper},
        cnt::AbstractConnectivity, q, s
    )
    (;bodyid2sys_full_coords) = cnt
    (;joint,force) = appar
    body = joint.body
    q_body = @view q[bodyid2sys_full_coords[body.prop.id]]
    u, v = NCF.get_uv(body.coords,q_body)
    x, y = u[1:2]
    θ = s[1]
    
    mode = force.state.mode_id

    l = hypot(x, y)
    if mode == 1
        ret[1] = sin(θ) - y/l
    elseif mode == 2
        ret[1] = cos(θ) - x/l
    elseif mode == 3
        ret[1] = sin(θ) - y/l
    elseif mode == 4
        ret[1] = cos(θ) - x/l
    end
end


function auxi_jacobian!(∂S∂q,∂S∂s, appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper},
        cnt::AbstractConnectivity, q, s
    )
    (;bodyid2sys_full_coords, apparid2sys_full_coords_idx) = cnt
    (;joint,force) = appar
    body = joint.body
    q_body = @view q[bodyid2sys_full_coords[body.prop.id]]
    u, v = NCF.get_uv(body.coords,q_body)
    
    num_of_dim = NCF.get_num_of_dims(body.coords)
    select_uvw = body.coords.conversion_to_X
    ∂x∂q = @view select_uvw[num_of_dim+1, :]
    ∂y∂q = @view select_uvw[num_of_dim+2, :]
    
    θ = s[1]
    mode = force.state.mode_id
    
    if mode == 1
        # res = sin(θ) - y/l
        ∂S∂s[1,1] = cos(θ)
        ∂S∂x, ∂S∂y = ForwardDiff.gradient((u) -> begin l = hypot(u[1],u[2]); -u[2]/l end ,u)
        ∂S∂q_body = ∂S∂x * ∂x∂q + ∂S∂y * ∂y∂q
    elseif mode == 2
        # res = cos(θ) - x/l
        ∂S∂s[1,1] = -sin(θ)
        ∂S∂x, ∂S∂y = ForwardDiff.gradient((u) -> begin l = hypot(u[1],u[2]); -u[1]/l end ,u)
        ∂S∂q_body = ∂S∂x * ∂x∂q + ∂S∂y * ∂y∂q
    elseif mode == 3
        # res = sin(θ) - y/l
        ∂S∂s[1,1] = cos(θ)
        ∂S∂x, ∂S∂y = ForwardDiff.gradient((u) -> begin l = hypot(u[1],u[2]); -u[2]/l end ,u)
        ∂S∂q_body = ∂S∂x * ∂x∂q + ∂S∂y * ∂y∂q
    elseif mode == 4
        # res = cos(θ) - x/l
        ∂S∂s[1,1] = -sin(θ)
        ∂S∂x, ∂S∂y = ForwardDiff.gradient((u) -> begin l = hypot(u[1],u[2]); -u[1]/l end ,u)
        ∂S∂q_body = ∂S∂x * ∂x∂q + ∂S∂y * ∂y∂q
    else
         error("Invalid mode_id: $mode")
    end
    ∂S∂q[1, bodyid2sys_full_coords[body.prop.id]] .= ∂S∂q_body
end

function switch!(appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper},
        structure::AbstractStructure,
        q, s
    )
    (;force) = appar
    θ = s[1]
    force.state.mode_id = angle2mode(θ)
end


function get_auxilary(appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper})
    (;joint,force) = appar
    get_auxilary(force)
end
