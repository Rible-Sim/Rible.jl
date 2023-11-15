
function make_cstr_function(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;bodyid2sys_full_coords) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,
        mask_1st,
        violations,
        halves,
        transformations,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    function _inner_cstr_function(q,c)
        T = eltype(q)
        ret = zeros(T,num_of_cstr)
        q_hen = @view q[bodyid2sys_full_coords[id_hen]]
        q_egg = @view q[bodyid2sys_full_coords[id_egg]]
        q = vcat(q_hen,q_egg)
        # cstr violations
        if (nmcs_hen isa NCF.NC) && (nmcs_egg isa NCF.NC)
            for icstr = 1:num_of_cstr
                ret[icstr] = q'*halves[icstr]*q
            end
            ret[mask_1st] .+= transformations*q
            ret .-= violations
        elseif (nmcs_hen isa QCF.QC) && (nmcs_egg isa QCF.QC)
            quat_hen = Quaternion(q_hen[4:7]...)
            quat_egg = Quaternion(q_egg[4:7]...)
            quat_trl_rel_hen = QuatRotation(hen.rbsig.prop.loci[hen.rotid].axes.X[:,[2,3,1]]).q
            quat_trl_rel_egg = QuatRotation(egg.rbsig.prop.loci[egg.trlid].axes.X[:,[2,3,1]]).q
            quat_rot_rel_hen = QuatRotation(hen.rbsig.prop.loci[hen.rotid].axes.X[:,[2,3,1]]).q
            quat_rot_rel_egg = QuatRotation(egg.rbsig.prop.loci[egg.rotid].axes.X[:,[2,3,1]]).q
            quat_trl_hen = quat_hen*quat_trl_rel_hen
            quat_trl_egg = quat_egg*quat_trl_rel_egg
            quat_rot_egg = quat_egg*quat_rot_rel_egg
            c_hen = to_local_coords(nmcs_hen,hen.rbsig.prop.loci[hen.pid].position)
            c_egg = to_local_coords(nmcs_egg,egg.rbsig.prop.loci[egg.pid].position)
            r_hen = QCF.to_position(nmcs_hen,q_hen,c_hen)
            r_egg = QCF.to_position(nmcs_egg,q_egg,c_egg)
            d = r_egg - r_hen
            # cstr 1st
            vio_1st = d
            vio_4th = [d'*d]
            vio_3rd = vcat(
                QCF.Rmat(quat_trl_hen)'*d,
                QCF.Rmat(quat_trl_egg)'*d
            )
            vio_2nd = QCF.vec(conj(quat_rot_egg)*quat_hen)[2:4]
            ret .= vcat(
                vio_1st[mask_1st],
                vio_4th[mask_4th],
                vio_3rd[mask_3rd],
                vio_2nd[mask_2nd];
            ) .- violations

            # @show ret
            ret
        end
    end
    function inner_cstr_function(q)
        c = get_local_coords(st)
        _inner_cstr_function(q,c)
    end
    inner_cstr_function
end

function make_cstr_jacobian(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
    ) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,
        mask_1st,
        hessians,
        transformations,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    free_idx_hen =  hen.rbsig.coords.free_idx
    free_idx_egg =  egg.rbsig.coords.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    function _inner_cstr_jacobian(q,c)
        T = eltype(q)
        ret = zeros(T,num_of_cstr,num_of_free_coords)
        q_hen = @view q[bodyid2sys_full_coords[id_hen]]
        q_egg = @view q[bodyid2sys_full_coords[id_egg]]
        q = vcat(q_hen,q_egg)
        # translate
        if (nmcs_hen isa NCF.NC) && (nmcs_egg isa NCF.NC)
            for icstr in 1:num_of_cstr
                A = q'*hessians[icstr]
                ret[icstr,free_hen] .= A[1,free_idx_hen]
                ret[icstr,free_egg] .= A[1,(free_idx_egg.+num_of_coords_hen)]
            end
            ret[mask_1st,free_hen] .+= transformations[mask_1st,free_idx_hen]
            ret[mask_1st,free_egg] .+= transformations[mask_1st,(free_idx_egg.+num_of_coords_hen)]
            ret
        elseif (nmcs_hen isa QCF.QC) && (nmcs_egg isa QCF.QC)
            quat_hen = Quaternion(q_hen[4:7]...)
            quat_egg = Quaternion(q_egg[4:7]...)
            quat_trl_rel_hen = QuatRotation(hen.rbsig.prop.loci[hen.rotid].axes.X[:,[2,3,1]]).q
            quat_trl_rel_egg = QuatRotation(egg.rbsig.prop.loci[egg.trlid].axes.X[:,[2,3,1]]).q
            quat_rot_rel_hen = QuatRotation(hen.rbsig.prop.loci[hen.rotid].axes.X[:,[2,3,1]]).q
            quat_rot_rel_egg = QuatRotation(egg.rbsig.prop.loci[egg.rotid].axes.X[:,[2,3,1]]).q
            quat_trl_hen = quat_hen*quat_trl_rel_hen
            quat_trl_egg = quat_egg*quat_trl_rel_egg
            quat_rot_egg = quat_egg*quat_rot_rel_egg
            # @show quat_trl_hen
            # @show quat_trl_egg
            # @show quat_rot_egg
            c_hen = to_local_coords(nmcs_hen,hen.rbsig.prop.loci[hen.pid].position)
            c_egg = to_local_coords(nmcs_egg,egg.rbsig.prop.loci[egg.pid].position)
            C_hen = to_transformation(nmcs_hen,q_hen,c_hen)
            C_egg = to_transformation(nmcs_egg,q_egg,c_egg)
            J = [-C_hen C_egg;] 
            r_hen = QCF.to_position(nmcs_hen,q_hen,c_hen)
            r_egg = QCF.to_position(nmcs_egg,q_egg,c_egg)
            d = r_egg - r_hen
            # hes 1st
            jac_1st = zeros(T,3,num_of_jointed_coords)
            # hes 4th
            jac_4th = zeros(T,1,num_of_jointed_coords)
            # hes 3rd
            jac_3rd = zeros(T,6,num_of_jointed_coords)
            # rotate
            jac_2nd = zeros(T,3,num_of_jointed_coords)
            # jac 1st 
            jac_1st[1:3,:] .= J
            # jac 4th 
            jac_4th[:,:]   .= 2d'*J
            # jac 3rd on hen
            # translate on hen
            jac_3rd[1:3,:]                         .= QCF.Rmat(quat_trl_hen)'*J
            jac_3rd[1:3,4:7]                      .+= QCF.∂Rᵀη∂q(quat_trl_hen,d)*QCF.Mmat(quat_trl_rel_hen)
            # jac 3rd on egg
            # translate on egg
            jac_3rd[4:6,:]                         .= QCF.Rmat(quat_trl_egg)'*J
            jac_3rd[4:6,num_of_coords_hen.+(4:7)] .+= QCF.∂Rᵀη∂q(quat_trl_egg,d)*QCF.Mmat(quat_trl_rel_egg)
            # jac 2nd
            # rotate of egg
            O43 = @SMatrix zeros(T,4,3)
            jac_2nd[1:3,:] = (
                [
                    O43 QCF.Pmat(conj(quat_rot_egg))*QCF.Pmat(conj(quat_egg)) O43 QCF.Pmat(conj(quat_rot_egg))*QCF.Mmat(quat_hen)*QCF.Inv_mat
                ]
            )[2:4,:]
            jac = vcat(
                jac_1st[mask_1st,:],
                jac_4th[mask_4th,:],
                jac_3rd[mask_3rd,:],
                jac_2nd[mask_2nd,:];
            )
            # @show jac
            ret[:,free_hen] .= jac[:,free_idx_hen]
            ret[:,free_egg] .= jac[:,(free_idx_egg.+num_of_coords_hen)]
            ret
        end
    end
    function inner_cstr_jacobian(q,c)
        _inner_cstr_jacobian(q,c)
    end
    function inner_cstr_jacobian(q)
        c = get_local_coords(st)
        _inner_cstr_jacobian(q,c)
    end
    inner_cstr_jacobian
end

function get_jointed_free_idx(cst)
    (;num_of_cstr,hen2egg,) = cst
    (;hen,egg) = hen2egg
    free_idx_hen = hen.rbsig.coords.free_idx
    free_idx_egg = egg.rbsig.coords.free_idx
    ncoords_hen = NCF.get_num_of_coords(hen.rbsig.coords.nmcs)
    # ncoords_egg = NCF.get_num_of_coords(egg.rbsig.coords.nmcs)
    free_idx = vcat(
        free_idx_hen,
        free_idx_egg .+ ncoords_hen
    )
end

function get_jointed_free(cst,indexed)
    (;num_of_cstr,hen2egg,) = cst
    (;bodyid2sys_free_coords,bodyid2sys_full_coords,num_of_free_coords) = indexed
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    mem2sysfree1 = bodyid2sys_free_coords[id_hen]
    mem2sysfree2 = bodyid2sys_free_coords[id_egg]
    cst2sysfree = vcat(
        mem2sysfree1,
        mem2sysfree2
    )
end

function make_cstr_forces_jacobian(cst::PrototypeJoint,st)
    (;
        num_of_cstr,
        hessians
    ) = cst
    free_idx = get_jointed_free_idx(cst)
    function cstr_forces_jacobian(λ)
        ret = [
            begin
                a = -λ[i] .* hessians[i][free_idx,free_idx]
                # display(a)
                a 
            end
            for i = 1:num_of_cstr
        ]
        sum(ret)
    end
end

function get_jointed_free_idx(cst::LinearJoint)
    free_idx = collect(1:size(cst.A,2))
end

function get_jointed_free(cst::LinearJoint,indexed)
    cst2sysfree = collect(1:indexed.num_of_free_coords)
end

function make_cstr_forces_jacobian(cst::LinearJoint,st)
    free_idx = get_jointed_free_idx(cst)
    function cstr_forces_jacobian(λ)
        zeros(eltype(λ),length(free_idx),length(free_idx))
    end
end

function get_jointed_free(cst::FixedIndicesConstraint,indexed)
    cst2sysfree = collect(1:indexed.num_of_free_coords)
end

function make_cstr_forces_jacobian(cst::FixedIndicesConstraint,st)
    (;num_of_free_coords) = st.connectivity.indexed
    function cstr_forces_jacobian(λ)
        zeros(eltype(λ),num_of_free_coords,num_of_free_coords)
    end
end