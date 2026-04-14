
"""
Prescribed-Free Rigid Body Connectivity Type.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct PresFreeConnectivity{idxType,id2idxType} <: AbstractConnectivity
    num_of_bodies::Int
    num_of_apparatuses::Int
    num_of_joint_apparatuses::Int
    num_of_force_apparatuses::Int
    num_of_full_coords::Int
    num_of_free_coords::Int
    num_of_pres_coords::Int
    num_of_intrinsic_cstr::Int
    num_of_extrinsic_cstr::Int
    num_of_cstr::Int
    num_of_dof_unconstrained::Int
    num_of_dof::Int
    num_of_aux_var::Int
    sys_free_coords_idx::idxType
    sys_pres_coords_idx::idxType
    bodyid2sys_full_coords::id2idxType
    bodyid2sys_free_coords::id2idxType
    bodyid2sys_pres_coords::id2idxType
    bodyid2sys_intrinsic_cstr_idx::id2idxType
    bodyid2sys_dof_idx::id2idxType
    apparid2sys_extrinsic_cstr_idx::id2idxType
    apparid2sys_full_coords_idx::id2idxType
    apparid2sys_aux_var_idx::id2idxType
    "body's loci' idx to System's loci' idx"
    bodyid2sys_locus_id::Vector{Vector{Int}}
    "System's loci' idx to Signifier"
    sys_locus_id2sig::Vector{Signifier{Int64}}
    "System's loci' idx to System's loci' coords' idx"
    sys_locus_id2coords_idx::Vector{Vector{Int}}
    "body's loci' idx to System's loci' coords' idx"
    bodyid2sys_loci_coords_idx::Vector{Vector{Int}}
    "Number of the System's loci"
    num_of_sys_loci::Int
    "Number of the System's loci' coords"
    num_of_sys_loci_coords::Int
    num_of_body_params::Int
    num_of_appar_params::Int
    num_of_params::Int
    apparid2params_idx::Vector{Vector{Int}}
end


function presfree_process_coords_idx(bodies; sharing_matrix=Int[;;])
    bodies_ids,num_of_bodies = check_id_sanity(bodies)
    if size(sharing_matrix,2) > num_of_bodies
        @warn "Cropping the sharing matrix."
        sharing = sharing_matrix[:,1:num_of_bodies]
    else
        sharing = sharing_matrix[:,:]
    end
    sys_full_coords_idx = Int[]
    sys_pres_coords_idx = Int[]
    sys_free_coords_idx = Int[]
    bodyid2sys_full_coords = Vector{Int}[]
    bodyid2sys_pres_coords = Vector{Int}[]
    bodyid2sys_free_coords = Vector{Int}[]
    ntotal_by_body = zeros(Int,num_of_bodies)
    pres_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    free_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    foreach(bodies) do body
        bodyid = body.prop.id
        ntotal_by_body[bodyid] = get_num_of_coords(body)
        pres_idx_by_body[bodyid] = body.coords.pres_idx
        free_idx_by_body[bodyid] = body.coords.free_idx
    end
    for bodyid = 1:num_of_bodies
        ntotal = ntotal_by_body[bodyid]
        pres = pres_idx_by_body[bodyid]
        free = free_idx_by_body[bodyid]
        push!(bodyid2sys_full_coords,fill(-1,ntotal))
        push!(bodyid2sys_pres_coords,Int[])
        push!(bodyid2sys_free_coords,Int[])
        unshareds = collect(1:ntotal)
        shared_idx = Int[]
        for row in eachrow(sharing)
            rbids = findall(!iszero,row)
            if bodyid in rbids[begin+1:end]
                myindex = row[bodyid]
                formerid = first(rbids)
                formerindex = row[formerid]
                bodyid2sys_full_coords[bodyid][myindex] = bodyid2sys_full_coords[formerid][formerindex]
                push!(shared_idx,myindex)
            end
        end
        deleteat!(unshareds,sort(shared_idx))
        nusi = length(unshareds)
        bodyid2sys_full_coords[bodyid][unshareds] = collect(length(sys_full_coords_idx)+1:length(sys_full_coords_idx)+nusi)
        append!(sys_full_coords_idx,bodyid2sys_full_coords[bodyid][unshareds])
        for i in unshareds
            if i in pres
                # pres
                push!(sys_pres_coords_idx,bodyid2sys_full_coords[bodyid][i])
            else
                # free
                push!(sys_free_coords_idx,bodyid2sys_full_coords[bodyid][i])
            end
        end
        for i in free
            free_idx = findfirst((x)->x==bodyid2sys_full_coords[bodyid][i],sys_free_coords_idx)
            push!(bodyid2sys_free_coords[bodyid],free_idx)
        end
        for i in pres
            pres_idx = findfirst((x)->x==bodyid2sys_full_coords[bodyid][i],sys_pres_coords_idx)
            push!(bodyid2sys_pres_coords[bodyid],pres_idx)
        end
    end
    num_of_free_coords = length(sys_free_coords_idx)
    num_of_pres_coords = length(sys_pres_coords_idx)
    num_of_full_coords = length(sys_full_coords_idx)
    @eponymtuple(
        num_of_bodies, 
        sys_free_coords_idx, sys_pres_coords_idx, sys_full_coords_idx, 
        num_of_free_coords,num_of_pres_coords,num_of_full_coords,
        bodyid2sys_free_coords, bodyid2sys_pres_coords, bodyid2sys_full_coords
    )
end


"""
Constructor for the `PresFreeConnectivity` type.
$(TYPEDSIGNATURES)
"""
function PresFreeConnectivity(bodies,apparatuses=Int[];sharing_matrix::AbstractMatrix=Int[;;])
    apparatuses_ids,num_of_apparatuses = check_id_sanity(apparatuses)
    num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,num_of_dof_unconstrained,bodyid2sys_dof_idx = index_incstr(bodies)
    (;
        num_of_bodies, 
        sys_free_coords_idx, sys_pres_coords_idx, sys_full_coords_idx, 
        num_of_free_coords,num_of_pres_coords,num_of_full_coords,
        bodyid2sys_free_coords, bodyid2sys_pres_coords, bodyid2sys_full_coords 
    ) = presfree_process_coords_idx(bodies;sharing_matrix)

    (;
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        num_of_extrinsic_cstr,
        num_of_aux_var,
        apparid2sys_full_coords_idx,
        apparid2sys_aux_var_idx,
        apparid2sys_extrinsic_cstr_idx
    ) = process_apparatus_indices(apparatuses, bodyid2sys_full_coords, num_of_intrinsic_cstr)

    num_of_cstr = num_of_intrinsic_cstr + num_of_extrinsic_cstr
    num_of_dof = num_of_free_coords - num_of_cstr
    if num_of_dof <= 0
        @warn "Non positive degree of freedom: $num_of_dof."
    end

    (;
        bodyid2sys_locus_id,
        sys_locus_id2sig,
        sys_locus_id2coords_idx,
        bodyid2sys_loci_coords_idx,
        num_of_sys_loci,
        num_of_sys_loci_coords
    ) = process_loci_indices(bodies)

    num_of_body_params = num_of_sys_loci_coords

    (;
        num_of_appar_params,
        num_of_params,
        apparid2params_idx
    ) = process_param_indices(apparatuses, num_of_body_params)

    PresFreeConnectivity(
        num_of_bodies,
        num_of_apparatuses,
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        num_of_full_coords,
        num_of_free_coords,
        num_of_pres_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        num_of_dof_unconstrained,
        num_of_dof,
        num_of_aux_var,
        sys_free_coords_idx,
        sys_pres_coords_idx,
        bodyid2sys_full_coords,
        bodyid2sys_free_coords,
        bodyid2sys_pres_coords,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_dof_idx,
        apparid2sys_extrinsic_cstr_idx,
        apparid2sys_full_coords_idx,
        apparid2sys_aux_var_idx,
        bodyid2sys_locus_id,
        sys_locus_id2sig,
        sys_locus_id2coords_idx,
        bodyid2sys_loci_coords_idx,
        num_of_sys_loci,
        num_of_sys_loci_coords,
        num_of_body_params,
        num_of_appar_params,
        num_of_params,
        apparid2params_idx,
    )
end



"""
Build a `StructureState` from bodies, apparatuses, and a `PresFreeConnectivity`.
$(TYPEDSIGNATURES)
"""
function StructureState(bodies,apparatuses,cnt::PresFreeConnectivity)
    (;
        num_of_bodies,
        num_of_full_coords,
        num_of_free_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        num_of_aux_var,
        sys_free_coords_idx,
        sys_pres_coords_idx,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_full_coords,
        bodyid2sys_free_coords,
        apparid2sys_aux_var_idx,
        bodyid2sys_loci_coords_idx,
    ) = cnt
    
    pres_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    free_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    foreach(bodies) do body
        bodyid = body.prop.id
        pres_idx_by_body[bodyid] = body.coords.pres_idx
        free_idx_by_body[bodyid] = body.coords.free_idx
    end
    T = get_numbertype(bodies)
    t = zero(T)
    q = zeros(T,num_of_full_coords)
    q̇ = zero(q)
    q̈ = zero(q)
    F = zero(q)
    λ = zeros(T,num_of_cstr)
    s = zeros(T,num_of_aux_var)
    p = zero(q)
    p̌ = zeros(T,num_of_free_coords)
    c = get_params(bodies,apparatuses,cnt)
    system = PresFreeCoordinatesState(t,q,q̇,q̈,F,p,p̌,λ,s,sys_free_coords_idx,sys_pres_coords_idx,c)
    members = [
        begin
            qmem = @view q[bodyid2sys_full_coords[bodyid]]
            q̇mem = @view q̇[bodyid2sys_full_coords[bodyid]]
            q̈mem = @view q̈[bodyid2sys_full_coords[bodyid]]
            Fmem = @view F[bodyid2sys_full_coords[bodyid]]
            λmem = @view λ[bodyid2sys_intrinsic_cstr_idx[bodyid]]
            smem = @view s[Int[]]
            cmem = @view c[bodyid2sys_loci_coords_idx[bodyid]]
            pmem = zero(p[bodyid2sys_full_coords[bodyid]])
            p̌mem = zero(p[bodyid2sys_free_coords[bodyid]])
            PresFreeCoordinatesState(
                t,
                qmem,q̇mem,q̈mem,Fmem,
                pmem,p̌mem,λmem,
                smem,
                free_idx_by_body[bodyid],
                pres_idx_by_body[bodyid],
                cmem
            )
        end
        for bodyid = 1:num_of_bodies
    ]
    M = spzeros(T,num_of_full_coords,num_of_full_coords)
    foreach(bodies) do body
        q, q̇ = body_state2coords_state(body)
        members[body.prop.id].q .= q
        members[body.prop.id].q̇ .= q̇
        memfull = bodyid2sys_full_coords[body.prop.id]
        M[memfull,memfull] .+= body.cache.inertia_cache.M
    end
    system.p = M*system.q̇
    system.p̌ .= system.p[sys_free_coords_idx]
    
    foreach(apparatuses) do appar
        prepare_cache!(appar,cnt)
        s[apparid2sys_aux_var_idx[appar.id]] = get_auxilary(appar)
    end
    StructureState(system,members)
end

# out of space for convinience only
function cstr_jacobian(structure::AbstractStructure, inst_state::PresFreeCoordinatesState)
    (;num_of_cstr,num_of_free_coords) = structure.connectivity
    ret = zeros(eltype(inst_state.q),num_of_cstr,num_of_free_coords)
    cstr_jacobian!(ret,structure,inst_state)
    ret
end

function cstr_jacobian!(ret, structure::AbstractStructure,inst_state::PresFreeCoordinatesState)
    (;
        sys_free_coords_idx,num_of_cstr,num_of_full_coords
    ) = structure.connectivity
    ret_full = zeros(eltype(inst_state.q),num_of_cstr,num_of_full_coords)
    cstr_jacobian!(ret_full, structure,CoordinatesState(inst_state))
    ret .= ret_full[:,sys_free_coords_idx] 
end


get_num_of_coords(cnt::PresFreeConnectivity) = cnt.num_of_full_coords
get_num_of_free_coords(cnt::PresFreeConnectivity) = cnt.num_of_free_coords
get_num_of_pres_coords(cnt::PresFreeConnectivity) = cnt.num_of_pres_coords
get_free_coords_idx(cnt::PresFreeConnectivity) = cnt.sys_free_coords_idx
get_pres_coords_idx(cnt::PresFreeConnectivity) = cnt.sys_pres_coords_idx


function make_nullspace(st::AbstractStructure{bType,sType,<:PresFreeConnectivity},q0::AbstractVector) where {bType,sType}
	(;bodies,connectivity) = st
    (;num_of_free_coords,num_of_full_coords,sys_pres_coords_idx,sys_free_coords_idx,bodyid2sys_free_coords,bodyid2sys_intrinsic_cstr_idx,num_of_intrinsic_cstr) = connectivity
    function inner_nullspace(q̌)
        T = eltype(q̌)
		q = Vector{T}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        ret = zeros(T,num_of_free_coords,num_of_free_coords-num_of_intrinsic_cstr)
        foreach(bodies) do body
            bodyid = body.prop.id
			memfree = bodyid2sys_free_coords[bodyid]
            if !isempty(bodyid2sys_intrinsic_cstr_idx[bodyid])
                if body.coords isa NCF.NC3D12C
                        u,v,w = NCF.get_uvw(body.coords,q̌[memfree])
                        N = @view ret[bodyid2sys_free_coords[bodyid],bodyid2sys_intrinsic_cstr_idx[bodyid]]
                        N[1:3,1:3]   .= Matrix(1I,3,3)
                        N[4:6,4:6]   .= -skew(u)
                        N[7:9,4:6]   .= -skew(v)
                        N[10:12,4:6] .= -skew(w)
                elseif body.coords isa NCF.NC2D6C                    
                        u,v = NCF.get_uv(body.coords,q̌[memfree])
                        N = @view ret[bodyid2sys_free_coords[bodyid],bodyid2sys_intrinsic_cstr_idx[bodyid]]
                        N[1:2,1:2] .= Matrix(1I,2,2)
                        N[3:4,3] .= -skew(u)
                        N[5:6,3] .= -skew(v)
                end
            end
        end
        ret
    end
end

function gen_force_state_jacobian!(∂F∂q̌, ∂F∂q̌̇, ∂F∂u, bot::Robot, field::AbstractField, policy::AbstractPolicy, inst_state::PresFreeCoordinatesState, ∂F∂s=nothing)
    (;structure) = bot
    nq = get_num_of_full_coords(structure)
    nu = get_num_of_actions(bot)
    ns = get_num_of_aux_var(structure)

    ∂F∂q_full = zeros(eltype(∂F∂q̌), nq, nq)
    ∂F∂q̇_full = zeros(eltype(∂F∂q̌̇), nq, nq)
    ∂F∂u_full = zeros(eltype(∂F∂u), nq, nu)
    ∂F∂s_full = zeros(eltype(∂F∂s), nq, ns)
    
    (;q,q̇,t,s) = inst_state
    (;structure,hub) = bot
    ∂F∂q̌ .= 0
    ∂F∂q̌̇ .= 0
    ∂F∂u .= 0
    clear_forces!(structure)
    actuate!(bot,policy,inst_state)
    lazy_update_bodies!(structure,inst_state)
    update_apparatuses!(structure, s)
    if s !== nothing
        ∂F∂s_full .= 0
        gen_force_auxi_jacobian!(∂F∂s_full,structure,inst_state,)
    end
    assemble_forces!(inst_state, structure)
    gen_force_state_jacobian!(∂F∂q_full, ∂F∂q̇_full, ∂F∂s_full, ∂F∂u_full, policy, bot, inst_state)
    add_tangent_stiffness_matrix!(∂F∂q_full, structure, field)
    add_tangent_damping_matrix!(∂F∂q̇_full, structure, field)

    free_idx = structure.connectivity.sys_free_coords_idx
    ∂F∂q̌ .= @view ∂F∂q_full[free_idx, free_idx]
    ∂F∂q̌̇ .= @view ∂F∂q̇_full[free_idx, free_idx]
    ∂F∂u .= @view ∂F∂u_full[free_idx, :]
    ∂F∂s .= @view ∂F∂s_full[free_idx, :]
end

function cstr_forces_jacobian(st::AbstractStructure{bType,sType,<:PresFreeConnectivity},q,λ) where {bType,sType}
    (;num_of_free_coords) = st.connectivity
    ret = zeros(eltype(λ),num_of_free_coords,num_of_free_coords)
    cstr_forces_jacobian!(ret,st,q,λ)
    ret
end

#todo avoid duplicated code
function cstr_forces_jacobian!(ret, st::AbstractStructure{bType,sType,<:PresFreeConnectivity},q,λ) where {bType,sType}
    (;bodies,apparatuses) = st
    cnt = st.connectivity
    (;
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_full_coords,
        num_of_full_coords,
        apparid2sys_extrinsic_cstr_idx,
        sys_free_coords_idx
    ) = cnt
    ret_full = zeros(eltype(λ),num_of_full_coords,num_of_full_coords)
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        add_cstr_forces_jacobian!((@view ret_full[memfull,memfull]),
            body.coords,
            λ[memincst]
        )
    end
    #todo skip 2D for now
    if get_num_of_dims(st) == 3
        foreach(apparatuses) do appar
            cstr_idx = apparid2sys_extrinsic_cstr_idx[appar.id]
            appar_sys_full_idx = cnt.apparid2sys_full_coords_idx[appar.id]
            add_cstr_forces_jacobian!((@view ret_full[appar_sys_full_idx, appar_sys_full_idx]),
                appar,
                cnt,q,λ[cstr_idx]
            )
        end
    end
    ret .= ret_full[sys_free_coords_idx,sys_free_coords_idx]
end

function auxi_jacobian!(∂S∂q̌,∂S∂s,structure::AbstractStructure,inst_state::PresFreeCoordinatesState)
    (;q,s) = inst_state
    nq = structure.connectivity.num_of_full_coords
    ns = structure.connectivity.num_of_aux_var
    ∂S∂q_full = zeros(eltype(∂S∂q̌), ns, nq)
    
    auxi_jacobian!(∂S∂q_full,∂S∂s,structure,CoordinatesState(inst_state))
    
    free_idx = structure.connectivity.sys_free_coords_idx
    ∂S∂q̌ .= ∂S∂q_full[:, free_idx]
end


function build_Ǩ(st::AbstractStructure{bType,aType,<:PresFreeConnectivity},q,λ) where {bType, aType}
    (;num_of_full_coords, sys_free_coords_idx) = st.connectivity
    T = get_numbertype(st)
    ∂Q∂q = zeros(T,num_of_full_coords,num_of_full_coords)
    add_tangent_stiffness_matrix!(∂Q∂q,st)
    Ǩa = -cstr_forces_jacobian(st,q,λ)
    #note ∂Q∂q = -(Ǩm+Ǩg)
    Ǩ = - ∂Q∂q[sys_free_coords_idx,sys_free_coords_idx] .+ Ǩa
end