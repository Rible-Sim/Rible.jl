
function cstr_forces_jacobian(st::AbstractStructure,q,λ)
    (;bodies,apparatuses) = st
    (;numbered,indexed,) = st.connectivity
    (;
        num_of_free_coords,
        num_of_intrinsic_cstr,
        bodyid2sys_free_coords,
        bodyid2sys_intrinsic_cstr_idx,
        num_of_cstr,
        num_of_extrinsic_cstr,
        apparid2sys_extrinsic_cstr_idx
    ) = indexed
    ret = zeros(eltype(λ),num_of_free_coords,num_of_free_coords)
    foreach(bodies) do body
        bodyid = body.prop.id
        memfree = bodyid2sys_free_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        free_idx = body.coords.free_idx
        if !isempty(memincst)
            ret[memfree,memfree] .+= cstr_forces_jacobian(
                body.coords,
                λ[memincst]
            )
        end
    end
    #todo skip 2D for now
    if get_num_of_dims(st) == 3
        foreach(apparatuses) do appar
            joint_cstr_idx = num_of_intrinsic_cstr .+ apparid2sys_extrinsic_cstr_idx[appar.id]
            jointed_sys_free_idx = indexed.apparid2sys_free_coords_idx[appar.id]
            ret[jointed_sys_free_idx,jointed_sys_free_idx] .+= cstr_forces_jacobian(appar,st,q,λ[joint_cstr_idx])
        end
    end
    ret
end

function cstr_velocity_jacobian(st::AbstractStructure,q,q̇)
    (;num_of_free_coords) = st.connectivity.indexed
    (;bodies,apparatuses) = st
    (;indexed,) = st.connectivity
    (;num_of_cstr,num_of_extrinsic_cstr,apparid2sys_extrinsic_cstr_idx,
    num_of_intrinsic_cstr,bodyid2sys_free_coords,bodyid2sys_intrinsic_cstr_idx) = indexed
    ret = zeros(eltype(q̇),num_of_cstr,num_of_free_coords)
    foreach(bodies) do body
        bodyid = body.prop.id
        memfree = bodyid2sys_free_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        if !isempty(memincst)
            ret[memincst,memfree] .= cstr_velocity_jacobian(
                body.coords,
                q̇[memfree]
            )
        end
    end
    #todo skip 2D for now
    if get_num_of_dims(st) == 3
        foreach(apparatuses) do appar
            joint_cstr_idx = num_of_intrinsic_cstr .+ apparid2sys_extrinsic_cstr_idx[appar.id]
            jointed_sys_free_idx = indexed.apparid2sys_free_coords_idx[appar.id]
            ret[joint_cstr_idx,jointed_sys_free_idx] .= cstr_velocity_jacobian(appar,st,q,q̇)
        end
    end
    ret
end

function test_fvector(st,q0)
    function L(q)
        reset_forces!(st)
        distribute_q_to_rbs!(st,q,zero(q))
        update_cables_apply_forces!(st)
        fvector(st)
        [st.cables[i].state.length for i = 1:2]
    end
    FiniteDiff.finite_difference_jacobian(L,q0)
end

"""
$(TYPEDSIGNATURES)
"""
function linearize(tginput,λ,u,q,q̇=zero(q))
    st = deepcopy(tginput)
    set_restlen!(st,u)
    reset_forces!(st)
    distribute_q_to_rbs!(st,q,q̇)
    update_cables_apply_forces!(st)
    M = build_massmatrix(st)
    A = build_A(st)
    Q̃ = build_Q̃(st)
    ∂L∂q,∂L∂q̇ = build_tangent(st)
    (;ncoords,num_of_cstr) = st
    nz = ncoords + num_of_cstr
    M̂ = zeros(eltype(q),nz,nz)
    Ĉ  = zeros(eltype(q),nz,nz)
    K̂ = zeros(eltype(q),nz,nz)
    M̂[1:ncoords,1:ncoords] .= M
    Ĉ[1:ncoords,1:ncoords] .= -Q̃*∂L∂q̇

    # fjac = test_fvector(st,q)
    K̂[1:ncoords,1:ncoords] .= -Q̃*∂L∂q .+ cstr_forces_jacobian(st,λ)
    Aq = A(q)
    c = maximum(abs.(K̂[1:ncoords,1:ncoords]))
    K̂[1:ncoords,ncoords+1:nz] .= c.*transpose(Aq)
    K̂[ncoords+1:nz,1:ncoords] .= c.*Aq
    M̂,Ĉ,K̂
end

function intrinsic_nullspace(st,q)
    (;bodies,connectivity) = st
    (;indexed,) = connectivity
    (;num_of_full_coords,bodyid2sys_full_coords,num_of_dof_unconstrained,bodyid2sys_dof_idx,) = indexed
    ret = zeros(eltype(q),num_of_full_coords,num_of_dof_unconstrained)
    foreach(bodies) do body
        bodyid = body.prop.id
        (;nmcs) = body.coords
        mem2full = bodyid2sys_full_coords[bodyid]
        ret[mem2full,bodyid2sys_dof_idx[bodyid]] = nullspace_mat(nmcs,q[mem2full])
    end
    ret
end

function extrinsic_nullspace(st,q)
    (;bodies,connectivity,num_of_dof) = st
    (;indexed,jointed) = connectivity
    (;num_of_full_coords,bodyid2sys_full_coords,num_of_dof_unconstrained,) = indexed
    ret = zeros(eltype(q),num_of_dof_unconstrained,num_of_dof)
    bodyid2sys_dof_idx = deepcopy(indexed.bodyid2sys_dof_idx)
    nullspaces = Vector{Matrix{eltype(q)}}(undef,2)
    foreach(jointed.joints) do joint
        apparid = joint.id
        if joint isa PrototypeJoint
            nullspaces[1] = zeros(eltype(q),num_of_dof_unconstrained,7)
            for i in bodyid2sys_dof_idx[1]
                nullspaces[1][i,i] = 1
            end
            (;hen,egg) = joint.hen2egg
            nmcs_hen = hen.body.coords.nmcs
            nmcs_egg = egg.body.coords.nmcs
            id_hen = hen.body.prop.id
            id_egg = egg.body.prop.id
            q_hen = q[bodyid2sys_full_coords[id_hen]]
            q_egg = q[bodyid2sys_full_coords[id_egg]]
            c_hen = to_local_coords(nmcs_hen,hen.body.prop.loci[hen.pid].position)
            c_egg = to_local_coords(nmcs_egg,egg.body.prop.loci[egg.pid].position)
            r_hen =  to_position_jacobian(nmcs_hen,q_hen,c_hen)*q_hen
            r_egg =  to_position_jacobian(nmcs_egg,q_egg,c_egg)*q_egg
            ro_hen = to_position_jacobian(nmcs_hen,q_hen,zero(c_hen))*q_hen
            ro_egg = to_position_jacobian(nmcs_egg,q_egg,zero(c_egg))*q_egg

            # hen.body.prop.loci[hen.trlid].axes
            # egg.body.prop.loci[egg.trlid].axes
            # hen.body.prop.loci[hen.rotid].axes
            invX̄_hen = nmcs_hen.data.invX̄
            invX̄_egg = nmcs_egg.data.invX̄
            axes_rot_egg = invX̄_egg*egg.body.prop.loci[egg.rotid].axes
            axes_idx = [
                (3,1), #normal * tangent
                (3,2), #normal * bitangent
                (1,2)  #tangent * bitangent
            ]
            id_axis_egg = 3
            X_hen = NCF.get_X(nmcs_hen,q_hen)
            X_egg = NCF.get_X(nmcs_egg,q_egg)
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

function frequencyshift(M̂,Ĉ,K̂,α::Real)
    M̄ = M̂
    C̄ = 2α*M̂ + Ĉ
    K̄ = α^2*M̂ + α*Ĉ + K̂
    M̄,C̄,K̄
end

function frequencyshift(M̂,K̂,α::Real)
    M̄ = M̂
    K̄ = α*M̂ + K̂
    M̄,K̄
end

function enlarge(M̄,C̄,K̄)
    T = eltype(M̄)
    nz = size(M̄)[1]
    M̃ = zeros(T,2nz,2nz)
    M̃[1:nz,1:nz] .= -C̄
    M̃[1:nz,nz+1:2nz] .= -M̄
    M̃[nz+1:2nz,1:nz] .= Matrix(one(T)*I,nz,nz)
    K̃ = zeros(T,2nz,2nz)
    K̃[1:nz,1:nz] .= K̄
    K̃[nz+1:2nz,nz+1:2nz] .= Matrix(one(T)*I,nz,nz)
    M̃,K̃
end

function find_finite(ω2,Z,num_of_dof)
    first_frequency_index = findfirst((x)->x>0,ω2)
    finite_ω2 = ω2[first_frequency_index:first_frequency_index+num_of_dof-1]
    finite_Z = Z[:,first_frequency_index:first_frequency_index+num_of_dof-1]
    finite_ω2,finite_Z
end

function build_material_stiffness_matrix!(st::Structure,q,k)
    (;num_of_dim) = st
    (;indexed,) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_coords_idx,bodyid2sys_full_coords) = indexed
    update!(st,q)
    Jj = zeros(eltype(q),num_of_dim,num_of_full_coords)
    retǨm = zeros(eltype(q),num_of_free_coords,num_of_free_coords)
    foreach(st.apparatuses) do appar
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
            Ūjq = Uj[sys_free_coords_idx,:]*q
            retǨm .+= k[j]*s^2*(Ūjq*transpose(Ūjq))
        end
    end
    retǨm
end

function build_geometric_stiffness_matrix!(st::Structure,q,f)
    (;num_of_dim) = st
    (;indexed,) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_coords_idx,bodyid2sys_full_coords) = indexed
    update!(st,q)
    Jj = zeros(eltype(q),num_of_dim,num_of_full_coords)
    retǨg = zeros(eltype(q),num_of_free_coords,num_of_free_coords)
    foreach(st.apparatuses) do appar
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
            Ǔj = @view Uj[sys_free_coords_idx,sys_free_coords_idx]
            Ūjq = Uj[sys_free_coords_idx,:]*q
            retǨg .+= f[j]/length*(Ǔj-s^2*Ūjq*transpose(Ūjq))
        end
    end
    retǨg
end

function make_Ǩm_Ǩg(st,q0)
    (;num_of_dim) = st
    (;numbered,indexed,) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_pres_coords_idx,sys_free_coords_idx,bodyid2sys_full_coords) = indexed
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = numbered
    function inner_Ǩm_Ǩg(q̌,s,μ,k,c)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        retǨm = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        retǨg = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        foreach(st.apparatuses) do appar
            j = appar.id
            (;hen,egg) = appar.joint.hen2egg
            rb1 = hen.body
            rb2 = egg.body
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = hen.pid
            ap2id = egg.pid
            c1 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb1id][ap1id]]]
            c2 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb2id][ap2id]]]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            Ǔj = @view Uj[sys_free_coords_idx,sys_free_coords_idx]
            Ūjq = Uj[sys_free_coords_idx,:]*q
            retǨm .+= k[j]*s[j]^2*(Ūjq*transpose(Ūjq))
            retǨg .+= k[j]*(1-μ[j]*s[j])*(Ǔj-s[j]^2*Ūjq*transpose(Ūjq))
        end
        retǨm,retǨg
    end
    function inner_Ǩm_Ǩg(q̌)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        retǨm = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        retǨg = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        foreach(st.apparatuses) do appar
            j = appar.id
            (;force,joint) = appar
            (;hen,egg) = joint.hen2egg
            rb1 = hen.body
            rb2 = egg.body
            ap1id = hen.pid
            ap2id = egg.pid
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            (;k,c,state,slack) = force
            (;direction,tension,length,lengthdot) = state
            s = 1/length
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            Ǔj = @view Uj[sys_free_coords_idx,sys_free_coords_idx]
            Ūjq = Uj[sys_free_coords_idx,:]*q
            retǨm .+= k*s^2*(Ūjq*transpose(Ūjq))
            retǨg .+= tension/length*(Ǔj-s^2*Ūjq*transpose(Ūjq))
        end
        retǨm,retǨg
    end
end

function make_S(st,q0)
    (;num_of_dim) = st
    (;numbered,indexed,) = st.connectivity
    (;sys_pres_coords_idx,sys_free_coords_idx,num_of_full_coords,bodyid2sys_full_coords) = indexed
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = numbered
    (;cables) = st.apparatuses
    ncables = length(cables)
    function inner_S(q̌,s)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        ret = zeros(eltype(q̌),ncables)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(st.apparatuses) do appar
            j = appar.id
            (;force,joint) = appar
            (;hen,egg) = joint.hen2egg
            rb1 = hen.body
            rb2 = egg.body
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = hen.pid
            ap2id = egg.pid
            # c1 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb1id][ap1id]]]
            # c2 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb2id][ap2id]]]
            # C1 = rb1.state.cache.funcs.C(c1)
            # C2 = rb2.state.cache.funcs.C(c2)
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            ret[j] = transpose(q)*Uj*q*s[j]^2 - 1
        end
        ret
    end
    function inner_S(q̌,s,c)
        q = Vector{eltype(q̌)}(undef,num_of_full_coords)
        q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
        q[sys_free_coords_idx] .= q̌
        ret = zeros(eltype(q̌),ncables)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(st.apparatuses) do appar
            j = appar.id
            (;hen,egg) = appar.joint.hen2egg
            rb1 = hen.body
            rb2 = egg.body
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = hen.pid
            ap2id = egg.pid
            c1 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb1id][ap1id]]]
            c2 = c[sys_loci2coords_idx[bodyid2sys_loci_idx[rb2id][ap2id]]]
            C1 = rb1.state.cache.funcs.C(c1)
            C2 = rb2.state.cache.funcs.C(c2)
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            ret[j] = transpose(q)*Uj*q*s[j]^2 - 1
        end
        ret
    end
    inner_S
end

# In-place ∂Q̌∂q̌ for cables and flexible bodies
function build_tangent_stiffness_matrix!(∂Q̌∂q̌,st)
    (;bodies,apparatuses,connectivity) = st
    (;indexed,numbered) = connectivity
    (;
        num_of_free_coords,
        bodyid2sys_free_coords,
        bodyid2sys_full_coords,
        apparid2sys_free_coords_idx
    ) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    # ∂Q̌∂q̌ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    q = get_coords(st)

    foreach(bodies) do body
        if body isa FlexibleBody
            (;coords,cache) = body
            (;e,) = cache
            ancs = coords.nmcs
            ∂Q∂e = ANCF.make_∂Q∂e(ancs)(e)
            mfree = bodyid2sys_free_coords[body.prop.id]
            (;free_idx) = body.coords
            ∂Q̌∂q̌[mfree,mfree] .-= ∂Q∂e[free_idx,free_idx]
        end
    end

    foreach(apparatuses) do appar
        if appar.joint isa CableJoint
            (;hen,egg) = appar.joint.hen2egg
            body_hen = hen.body
            body_egg = egg.body
            C_hen = body_hen.cache.Cps[hen.pid]
            C_egg = body_egg.cache.Cps[egg.pid]
            free_idx_hen = body_hen.coords.free_idx
            free_idx_egg = body_egg.coords.free_idx
            mfree_hen = bodyid2sys_free_coords[body_hen.prop.id]
            mfree_egg = bodyid2sys_free_coords[body_egg.prop.id]
            (;k,c,state,slack) = appar.force
            (;direction,tension) = state
            l = state.length
            l̇ = state.lengthdot
            if slack && (tension==0)
                ∂Q̌∂q̌ .-= 0
            else
                D .= direction*transpose(direction)
                density = tension/l
                β = c*l̇/l + density
                D .*= k-β
                D .+= β.*Im
                J̌ .= 0
                J̌[:,mfree_egg] .+= C_egg[:,free_idx_egg]
                J̌[:,mfree_hen] .-= C_hen[:,free_idx_hen]
                ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
            end
        elseif appar.joint isa ClusterJoint
            foreach(appar.force) do seg
                (;hen, egg) = seg.joint.hen2egg
                body_hen = hen.body
                body_egg = egg.body
                C_hen = body_hen.cache.Cps[hen.pid]
                C_egg = body_egg.cache.Cps[egg.pid]
                free_idx_hen = body_hen.coords.free_idx
                free_idx_egg = body_egg.coords.free_idx
                mfree_hen = bodyid2sys_free_coords[body_hen.prop.id]
                mfree_egg = bodyid2sys_free_coords[body_egg.prop.id]
                (;k, c, state) = seg.force
                (;direction, tension) = state
                l = state.length
                l̇ = state.lengthdot
                D .= direction*transpose(direction)
                density = tension / l
                β = c*l̇/l + density
                D .*= k - β
                D .+= β .* Im
                J̌ .= 0
                J̌[:,mfree_egg] .+= C_egg[:,free_idx_egg]
                J̌[:,mfree_hen] .-= C_hen[:,free_idx_hen]
                ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
            end
        elseif appar.force isa RotationalSpringDamper
            (;
                num_of_cstr,
                hen2egg,
                cache,
                mask_1st,mask_2nd,mask_3rd,mask_4th
            ) = appar.joint
            (;
                relative_core
            ) = cache
            full_idx = appar.full_coords_idx
            free_idx = appar.free_coords_idx
            sys_free_coords_idx = apparid2sys_free_coords_idx[appar.id]
            (;mask,k) = appar.force
            (;hen,egg) = hen2egg
            nmcs_hen = hen.body.coords.nmcs
            nmcs_egg = egg.body.coords.nmcs
            num_of_coords_hen = get_num_of_coords(nmcs_hen)
            num_of_coords_egg = get_num_of_coords(nmcs_egg)
            id_hen = hen.body.prop.id
            id_egg = egg.body.prop.id
            q_hen = @view q[bodyid2sys_full_coords[id_hen]]
            q_egg = @view q[bodyid2sys_full_coords[id_egg]]
            q_jointed = vcat(
                q_hen,
                q_egg
            )
            jointed2angles = make_jointed2angles(hen2egg,relative_core)
            nq = length(q_jointed)
            ## angles = jointed2angles(q_jointed)
	        ## torques = k.*angles
            ## angles_jacobian = ForwardDiff.jacobian(jointed2angles,q_jointed)
            ## angles_hessians = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(jointed2angles, x), q_jointed)
            # angles_hessians = FiniteDiff.finite_difference_jacobian(x -> ForwardDiff.jacobian(jointed2angles, x), q_jointed)
            ## reshaped_angles_hessians = reshape(angles_hessians,3,nq,nq)
            # @show sys_free_coords_idx, free_idx
            for i in mask
                ## angle = angles[i]
                ## torque = torques[i]
                ## generalized_force_jacobian = torque.*reshaped_angles_hessians[i,:,:] .+ k.*angles_jacobian[i,:]*angles_jacobian[[i],:]
                # @show generalized_force_jacobian
                ## ∂Q̌∂q̌[sys_free_coords_idx,sys_free_coords_idx] .-= generalized_force_jacobian[free_idx,free_idx]
            end
        end
    end

    return ∂Q̌∂q̌
end

# In-place ∂Q̌∂q̌̇ for cables
function build_tangent_damping_matrix!(∂Q̌∂q̌̇,st)
    (;apparatuses,connectivity) = st
    (;numbered,indexed) = connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_coords_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    # ∂Q̌∂q̌̇ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    foreach(apparatuses) do appar
        if appar.joint isa CableJoint
            spring_damper = appar.force
            (;hen,egg) = appar.joint.hen2egg
            body_hen = hen.body
            body_egg = egg.body
            C_hen = body_hen.cache.Cps[hen.pid]
            C_egg = body_egg.cache.Cps[egg.pid]
            free_idx_hen = body_hen.coords.free_idx
            free_idx_egg = body_egg.coords.free_idx
            mfree_hen = bodyid2sys_free_coords[body_hen.prop.id]
            mfree_egg = bodyid2sys_free_coords[body_egg.prop.id]
            (;k,c,state,slack) = spring_damper
            (;direction,tension) = state
            if slack && (tension == 0)
                ∂Q̌∂q̌̇ .-= 0
            else
                D .= direction*transpose(direction)
                D .*= c
                J̌ .= 0
                J̌[:,mfree_egg] .+= C_egg[:,free_idx_egg]
                J̌[:,mfree_hen] .-= C_hen[:,free_idx_hen]

                ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
            end
            # ∂Q̌∂q̌_full[mfree_egg,mfree_egg] .+= transpose(C_egg)*D*C_egg
            # ∂Q̌∂q̌_full[mfree_hen,mfree_egg] .-= transpose(C_hen)*D*C_egg
            # ∂Q̌∂q̌_full[mfree_egg,mfree_hen] .-= transpose(C_egg)*D*C_hen
            # ∂Q̌∂q̌_full[mfree_hen,mfree_hen] .+= transpose(C_hen)*D*C_hen
        end
    end
end

function build_∂Q̌∂s̄(st)
    (; bodies, apparatuses, connectivity) = st
    (; indexed, numbered) = connectivity
    (;
        num_of_free_coords,
        bodyid2sys_free_coords,
        bodyid2sys_full_coords,
        apparid2sys_free_coords_idx
    ) = indexed
    ns = 0
    nc = 0
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            ns += length(appar.joint.sps)
            nc += 1
        end
    end
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ∂Q̌∂s̄ = zeros(T, 2ns, num_of_free_coords)
    D = zeros(T, num_of_dim)
    lkn = zeros(T, 2ns, num_of_dim)
    J̌ = zeros(T, num_of_dim, num_of_free_coords)

    N_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    kc = zeros(T, ns + nc)
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force) = appar
            nsegs = length(force)
            for seg in force
                (; id, num_of_add_var) = seg
                kc[id+num_of_add_var] = seg.force.k
            end
            N = zeros(Float64, nsegs, 2nsegs - 2)
            N[1, 1:2] = [1 -1]
            N[end, end-1:end] = [-1 1]
            for i in 2:nsegs-1
                N[i, 2i-3:2i] = [-1 1 1 -1]
            end
            push!(N_list, N)
        end
    end
    N = reduce(blockdiag, N_list)
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force) = appar
            foreach(force) do seg
                idx = seg.id + seg.num_of_add_var
                (; hen, egg) = seg.joint.hen2egg
                body_hen = hen.body
                body_egg = egg.body
                C_hen = body_hen.cache.Cps[hen.pid]
                C_egg = body_egg.cache.Cps[egg.pid]
                free_idx_hen = body_hen.coords.free_idx
                free_idx_egg = body_egg.coords.free_idx
                mfree_hen = bodyid2sys_free_coords[body_hen.prop.id]
                mfree_egg = bodyid2sys_free_coords[body_egg.prop.id]
                (; k, c, state) = seg.force
                (; direction, tension) = state
                if tension == 0
                    ∂Q̌∂s̄ .-= 0
                else
                    D .= direction
                    J̌ .= 0
                    J̌[:, mfree_egg] .+= C_egg[:, free_idx_egg]
                    J̌[:, mfree_hen] .-= C_hen[:, free_idx_hen]
                    kN = kc[idx] .* N[idx, :]
                    @tullio lkn[k, l] = D[l] * kN[k]
                    ∂Q̌∂s̄ .-= lkn * J̌
                end
            end
        end
    end
    return ∂Q̌∂s̄'
end

function build_ζ(st)
    (; apparatuses, connectivity) = st
    (; apparid2sys_add_var_idx) = connectivity.indexed
    ns = 0
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            ns += length(appar.joint.sps)
        end
    end
    ζ = Vector{Float64}(undef, 2ns)
    ζ⁺ = @view ζ[begin:2:end]
    ζ⁻ = @view ζ[begin+1:2:end]
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force, joint, id) = appar
            (; sps) = joint
            nsp = length(sps)
            idx_begin = Int((apparid2sys_add_var_idx[id][1]-1) / 2)
            for i in 1:nsp
                ζ⁺[idx_begin+i] = force[i+1].force.state.tension / sps[i].α - force[i].force.state.tension
                ζ⁻[idx_begin+i] = force[i].force.state.tension - sps[i].α * force[i+1].force.state.tension
            end
        end
    end
    return ζ
end

function get_clusterA(st)
    (; apparatuses) = st
    nc = 0
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            nc += 1
        end
    end
    A_list = [Matrix{Float64}(undef, 1, 1) for _ in 1:nc]
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force, joint, id) = appar
            (; sps) = joint
            n = length(sps)
            ap = [-force[i+1].force.k for i in 1:n-1]
            ap0 = [force[i].force.k + force[i+1].force.k/sps[i].α for i in 1:n]
            ap1 = [-force[i+1].force.k/sps[i].α for i in 1:n-1]
            an0 = [force[i].force.k + force[i+1].force.k*sps[i].α for i in 1:n]
            an1 = [-force[i+1].force.k*sps[i].α for i in 1:n-1]
            A⁺ = diagm(-1=>ap, 0=>ap0, 1=>ap1)
            A⁻ = diagm(-1=>ap, 0=>an0, 1=>an1)
            A = [A⁺ -A⁺; -A⁻ A⁻]
            A_list[id] = A
        end
        return A_list
    end
    return A_list
end

function get_TransMatrix(n)
    T = zeros(Float64, 2n, 2n)
    for (j, i) in enumerate(1:2:2n)
        T[i, j] = 1
    end
    for (j, i) in enumerate(2:2:2n)
        T[i, j+n] = 1
    end
    return T
end

function build_∂ζ∂s̄(st)
    (; apparatuses) = st
    nc = 0
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            nc += 1
        end
    end
    A_list = [SparseMatrixCSC{Float64,Int64}(undef, 1, 1) for _ in 1:nc]
    clusterA = get_clusterA(st)
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force, id) = appar
            nsegs = length(force)
            T = get_TransMatrix(nsegs - 1)
            A = T*clusterA[id]*T'
            A_list[id] = A
        end
    end
    return reduce(blockdiag, A_list)
end

function build_∂ζ∂q(st)
    (; apparatuses, connectivity) = st
    (; indexed) = connectivity
    (; num_of_free_coords, bodyid2sys_free_coords, apparid2sys_add_var_idx) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ns = 0
    nc = 0
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            ns += length(appar.joint.sps)
            nc += 1
        end
    end
    ∂l∂q = zeros(T, ns+nc, num_of_free_coords)
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force) = appar
            foreach(force) do seg
                idx = seg.id + seg.num_of_add_var
                J̌ = zeros(T, num_of_dim, num_of_free_coords)
                (; hen, egg) = seg.joint.hen2egg
                body_hen = hen.body
                body_egg = egg.body
                C_hen = body_hen.cache.Cps[hen.pid]
                C_egg = body_egg.cache.Cps[egg.pid]
                free_idx_hen = body_hen.coords.free_idx
                free_idx_egg = body_egg.coords.free_idx
                mfree_hen = bodyid2sys_free_coords[body_hen.prop.id]
                mfree_egg = bodyid2sys_free_coords[body_egg.prop.id]
                (; direction) = seg.force.state
                J̌[:, mfree_egg] .+= C_egg[:, free_idx_egg]
                J̌[:, mfree_hen] .-= C_hen[:, free_idx_hen]
                ∂l∂q[idx, :] = direction'*J̌
            end
        end
    end
    kc = zeros(T, ns + nc)
    αc = zeros(T, ns)
    b_list = [SparseMatrixCSC{Float64,Int64}(undef, 1, 1) for _ in 1:nc]
    T_list = [SparseMatrixCSC{Float64,Int64}(undef, 1, 1) for _ in 1:nc]
    foreach(apparatuses) do appar
        if isa(appar.joint, ClusterJoint)
            (; force, joint, id) = appar
            (; sps) = joint
            nsegs = length(force)
            nsp = length(sps)
            sp_idx = apparid2sys_add_var_idx[id]
            sp_idx_begin = Int((sp_idx[1]-1) / 2)
            foreach(force) do seg
                kc_idx = seg.id + seg.num_of_add_var
                kc[kc_idx] = seg.force.k
            end
            for i in 1:nsp
                αc[sp_idx_begin+i] = sps[i].α
            end
            b⁺ = zeros(Float64, nsegs-1, nsegs)
            b⁻ = zeros(Float64, nsegs-1, nsegs)
            for i in 1:nsegs-1
                b⁺[i, i:i+1] = [-kc[i] kc[i+1]/αc[i]]
                b⁻[i, i:i+1] = [kc[i] -kc[i+1]*αc[i]]
            end
            b_list[id] = vcat(b⁺, b⁻)
            T_list[id] = get_TransMatrix(nsegs - 1)
        end
    end
    b = reduce(blockdiag, b_list)
    T = reduce(blockdiag, T_list)
    return T*b*∂l∂q
end

function build_Ǩ(st)
    _,λ = check_static_equilibrium_output_multipliers(st)
    q = get_coords(st)
    build_Ǩ(st,q,λ)
end

function build_Ǩ(st,q,λ)
    (;num_of_free_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    ǨmǨg = zeros(T,num_of_free_coords,num_of_free_coords)
    build_tangent_stiffness_matrix!(ǨmǨg,st)
    Ǩa = cstr_forces_jacobian(st,q,λ)
    Ǩ = - ǨmǨg .- Ǩa
    # Ǩ .= Ǩ
    Ǩ
end

function norm_wrt!(Z,M)
    n = size(Z,2)
    for i = 1:n
        z = @view Z[:,i]
        zmz = transpose(z)*M*z
        z ./= sqrt(zmz)
    end
    Z
end

function undamped_eigen(st;gravity=false)
    _,λ = check_static_equilibrium_output_multipliers(st;gravity)
    q = get_coords(st)
    q̌ = get_free_coords(st)
    M̌ = assemble_M̌(st)
    Ǩ = build_Ǩ(st,q,λ)
    Ǎ = make_cstr_jacobian(st)(q)
    Ň = nullspace(Ǎ)
    ℳ = transpose(Ň)*M̌*Ň
    𝒦 = transpose(Ň)*Ǩ*Ň
    # @show ℳ, 𝒦
    ω²,ξ = eigen(Symmetric(𝒦),Symmetric(ℳ))
    # @show transpose(ξ)*ℳ*ξ
    Ňξ = Ň*ξ
    # @show transpose(Ňξ)*M̌*Ňξ
    norm_wrt!(Ňξ,M̌)
    δq̌ = [v for v in eachcol(Ňξ)]
    ω²,δq̌
    # nq = length(q̌)
    # nλ = length(λ)
    # nx = nq + nλ
    # M̂ = zeros(eltype(q),nx,nx)
    # K̂ = zeros(eltype(q),nx,nx)
    # M̂[1:nq,1:nq] .= M̌
    # K̂[1:nq,1:nq] .= Ǩ
    # c = maximum(abs.(K̂[1:nq,1:nq]))
    # K̂[1:nq,nq+1:nx] .= c.*transpose(Ǎ)
    # K̂[nq+1:nx,1:nq] .= c.*Ǎ
    #
    # eigen(K̂,M̂)
end

function old_undamped_eigen(st)
    λ0 = check_static_equilibrium_output_multipliers(st)
    M̂,Ĉ,K̂ = linearize(st,q0,λ0)
    α = 10
    M̄,K̄ = frequencyshift(M̂,K̂,α)
    # @show size(K̄),rank(K̄),cond(K̄),rank(M̄)
    d,aug_Z = eigen(K̄,M̄)
    aug_ω2 = d .- α
    (;ncoords, num_of_dof) = st
    # @show aug_ω2
    ω2,Z = find_finite(aug_ω2,aug_Z,num_of_dof)
    ω = sqrt.(ω2)
    Zq = Z[1:ncoords,:]
    M = build_massmatrix(st)
    normalize_wrt_mass!(Zq,M)
    ω, Zq#, Z
end

function undamped_modal_solve!(st,q0,q̇0,λ0,tf,dt)
    M̂,Ĉ,K̂ = linearize(st,q0,λ0)
    # show(stdout,"text/plain",K̂)
    # showtable(K̂)
    # M̄,C̄,K̄ = RB.frequencyshift(M̂,Ĉ,K̂,0.1)
    # M̃,K̃ = RB.enlarge(M̄,C̄,K̄)
    aug_ω2,aug_Z = eigen(K̂,M̂)
    ω2,Z = find_finite(aug_ω2,aug_Z)
    # @show aug_ω2,ω2
    normalize_wrt_mass!(Z,M̂)
    # @show transpose(Z)*M̂*Z
    # @show transpose(Z)*K̂*Z
    ω = sqrt.(ω2)
    z0 = vcat(zero(q0),λ0)
    ż0 = vcat(q̇0,zero(λ0))
    ζ0 = transpose(Z)*M̂*z0
    ζd0 = transpose(Z)*M̂*ż0

    d = length(ζ0)
    step = Integer(tf/dt)
    ζ = Matrix{eltype(q0)}(undef,d,step+1)
    for it in 0:step
        t = dt*it
        ζ[:,it+1] .= ζ0.*cos.(ω.*t) .+ ζd0./ω.*sin.(ω.*t)
    end
    z = Z*ζ
    q = z[1:length(q0),:]
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

function make_nullspace(st::Structure,q0::AbstractVector)
	(;bodies,connectivity) = st
    (;num_of_free_coords,num_of_full_coords,sys_pres_coords_idx,sys_free_coords_idx,bodyid2sys_free_coords,bodyid2sys_intrinsic_cstr_idx,num_of_intrinsic_cstr) = connectivity.indexed
    function inner_nullspace(q̌)
        T = eltype(q̌)
		q = Vector{T}(undef,num_of_full_coords)
		q[sys_pres_coords_idx] .= q0[sys_pres_coords_idx]
		q[sys_free_coords_idx] .= q̌
        ret = zeros(T,num_of_free_coords,num_of_free_coords-num_of_intrinsic_cstr)
        foreach(bodies) do body
            bodyid = body.prop.id
            (;nmcs) = body.state.cache.funcs
			memfree = bodyid2sys_free_coords[bodyid]
            if !isempty(bodyid2sys_intrinsic_cstr_idx[bodyid])
                if nmcs isa NCF.NC3D12C
                        u,v,w = NCF.get_uvw(nmcs,q̌[memfree])
                        N = @view ret[bodyid2sys_free_coords[bodyid],bodyid2sys_intrinsic_cstr_idx[bodyid]]
                        N[1:3,1:3]   .= Matrix(1I,3,3)
                        N[4:6,4:6]   .= -skew(u)
                        N[7:9,4:6]   .= -skew(v)
                        N[10:12,4:6] .= -skew(w)
                elseif nmcs isa NCF.NC2D6C                    
                        u,v = NCF.get_uv(nmcs,q̌[memfree])
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

"""
$(TYPEDSIGNATURES)
"""
function check_stability(st::Structure;F̌=nothing,verbose=false)
    static_equilibrium,λ = check_static_equilibrium_output_multipliers(st;F=F̌)
    @assert static_equilibrium
    check_stability(st,λ;verbose)
end

function check_stability(st::Structure,λ;verbose=false)
    q = get_coords(st)
    c = get_local_coords(st)
    A = make_cstr_jacobian(st,q)
    Ň(q̌,c) = nullspace(A(q̌))
    check_stability(st,λ,Ň;verbose)
end

function check_stability(st::Structure,λ,Ň;verbose=false)
    q̌ = get_free_coords(st)
    c = get_local_coords(st)
    Ǩ0 = build_Ǩ(st,λ)
    Ň0 = Ň(q̌,c)
    𝒦0 = transpose(Ň0)*Ǩ0*Ň0
    eigen_result = eigen(𝒦0)
    nn = count(x -> x < 0, eigen_result.values)
    if nn > 1
        @warn "Instability detected! Number of negative eigenvalues: $nn"
        isstable = false
    else
        isstable = true
    end
    isstable, Ň0, eigen_result
end

