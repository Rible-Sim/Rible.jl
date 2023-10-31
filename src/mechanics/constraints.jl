
function make_cstr_function(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;bodyid2sys_full_coords) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,values,
        axes_trl_hen,axes_trl_egg,
        axes_rot_hen,axes_rot_egg,
        mask_1st, mask_2nd, 
        mask_3rd, mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    function _inner_cstr_function(q,d,c)
        T = eltype(q)
        ret = zeros(T,num_of_cstr)
        q_hen = @view q[bodyid2sys_full_coords[id_hen]]
        q_egg = @view q[bodyid2sys_full_coords[id_egg]]
        X_hen = NCF.get_X(nmcs_hen,q_hen)
        X_egg = NCF.get_X(nmcs_egg,q_egg)
        # translate
        c_hen = c[sys_loci2coords_idx[bodyid2sys_loci_idx[id_hen]][hen.pid]]
        c_egg = c[sys_loci2coords_idx[bodyid2sys_loci_idx[id_egg]][egg.pid]]
        C_hen = to_transformation(hen.rbsig.state.cache.funcs.nmcs,c_hen)
        C_egg = to_transformation(egg.rbsig.state.cache.funcs.nmcs,c_egg)
        # @show C_egg,q_egg
        J = [-C_hen C_egg;]
        q = vcat(q_hen,q_egg)
        r_hen2egg = J*q
        # first    
        Φ_1st = r_hen2egg
        # fourth    
        Φ_4th = [r_hen2egg'*r_hen2egg]
        # third
        Φ_3rd = zeros(T,6)
        # second
        Φ_2nd = zeros(T,3)
        if (nmcs_hen isa NCF.NC3D12C) && (nmcs_egg isa NCF.NC3D12C)
            # translate on hen
            trl_hen = (X_hen*axes_trl_hen).X
            # translate on egg
            trl_egg = (X_egg*axes_trl_egg).X
            # rotate of egg
            rot_hen_t = X_hen*axes_rot_hen.tangent
            rot_hen_b = X_hen*axes_rot_hen.bitangent
            rot_egg_n = X_egg*axes_rot_egg.normal
            rot_egg_b = X_egg*axes_rot_egg.bitangent
            # third
            Φ_3rd .= [
                trl_hen'*r_hen2egg;
                trl_egg'*r_hen2egg;
            ]
            # second
            Φ_2nd .= [
                rot_hen_t'*rot_egg_b,
                rot_hen_t'*rot_egg_n,
                rot_hen_b'*rot_egg_n
            ]
        end
        # cstr values
        ret = vcat(
            Φ_1st[mask_1st],
            Φ_4th[mask_4th],
            Φ_3rd[mask_3rd], 
            Φ_2nd[mask_2nd], 
        ) .- values
        ret
    end
    function inner_cstr_function(q,d,c)
        _inner_cstr_function(q,d,c)
    end
    function inner_cstr_function(q)
        d = get_d(st)
        c = get_local_coords(st)
        _inner_cstr_function(q,d,c)
    end
    inner_cstr_function
end

function make_cstr_jacobian(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;bodyid2sys_free_coords,bodyid2sys_full_coords,num_of_free_coords) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,
        axes_trl_hen,axes_trl_egg,
        axes_rot_hen,axes_rot_egg, 
        mask_1st, mask_2nd, 
        mask_3rd, mask_3rd, mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    Y_hen = nmcs_hen.conversion_to_std
    Y_egg = nmcs_egg.conversion_to_std
    function _inner_cstr_jacobian(q,c)
        T = eltype(q)
        ret = zeros(T,num_of_cstr,num_of_free_coords)
        q_hen = @view q[bodyid2sys_full_coords[id_hen]]
        q_egg = @view q[bodyid2sys_full_coords[id_egg]]
        X_hen = NCF.get_X(nmcs_hen,q_hen)
        X_egg = NCF.get_X(nmcs_egg,q_egg)
        # translate 
        c_hen = c[sys_loci2coords_idx[bodyid2sys_loci_idx[id_hen][hen.pid]]]
        c_egg = c[sys_loci2coords_idx[bodyid2sys_loci_idx[id_egg][egg.pid]]]
        C_hen = to_transformation(hen.rbsig.state.cache.funcs.nmcs,c_hen)
        C_egg = to_transformation(egg.rbsig.state.cache.funcs.nmcs,c_egg)
        J = [-C_hen C_egg;]
        q = vcat(q_hen,q_egg)
        r_hen2egg = J*q
        o3 = zeros(T,3)
        # jac 
        A_1st = zeros(T,3,num_of_free_coords)
        A_4th = zeros(T,1,num_of_free_coords)
        A_3rd = zeros(T,6,num_of_free_coords)
        A_2nd = zeros(T,3,num_of_free_coords)
        # first
        A_1st[:,free_hen] .= -C_hen[:,uci_hen]
        A_1st[:,free_egg] .=  C_egg[:,uci_egg]
        # jac 4th
        A_4th[1,free_hen] =  2r_hen2egg'*(-C_hen[:,uci_hen])
        A_4th[1,free_egg] =  2r_hen2egg'*( C_egg[:,uci_egg])

        if (nmcs_hen isa NCF.NC3D12C) && (nmcs_egg isa NCF.NC3D12C) 
            # translate on hen
            trl_hen = (X_hen*axes_trl_hen).X
            # translate on egg
            trl_egg = (X_egg*axes_trl_egg).X
            # rotate of egg
            rot_hen_t = X_hen*axes_rot_hen.tangent
            rot_hen_b = X_hen*axes_rot_hen.bitangent
            rot_egg_n = X_egg*axes_rot_egg.normal
            rot_egg_b = X_egg*axes_rot_egg.bitangent
            # jac third on hen
            A_3rd[1:3,free_hen] .= kron(
                    hcat(
                        o3, transpose(trl_hen)
                    ),
                    transpose(r_hen2egg)
                )[:,uci_hen]
            A_3rd[1:3,free_hen] .-= transpose(trl_hen)*C_hen[:,uci_hen]
            A_3rd[1:3,free_egg] .+= transpose(trl_hen)*C_egg[:,uci_egg]
            # jac third on egg
            A_3rd[4:6,free_egg] .= kron(
                    hcat(
                        o3, transpose(trl_egg)
                    ),
                    transpose(r_hen2egg)
                )[:,uci_egg]
            A_3rd[4:6,free_hen] .-= transpose(trl_egg)*C_hen[:,uci_hen]
            A_3rd[4:6,free_egg] .+= transpose(trl_egg)*C_egg[:,uci_egg]
            # jac 2nd
            A_2nd[1,free_hen] = (
                    transpose(
                        kron(vcat(0,axes_rot_hen.tangent),rot_egg_b)
                    )*Y_hen
                )[:,uci_hen]
            A_2nd[1,free_egg] = (
                    transpose(
                        kron(vcat(0,axes_rot_egg.bitangent),rot_hen_t)
                    )*Y_egg
                )[:,uci_egg]
            A_2nd[2,free_hen] = (
                    transpose(
                        kron(vcat(0,axes_rot_hen.tangent),rot_egg_n)
                    )*Y_hen
                )[:,uci_hen]
            A_2nd[2,free_egg] = (
                    transpose(
                        kron(vcat(0,axes_rot_egg.normal),rot_hen_t)
                    )*Y_egg
                )[:,uci_egg]
            A_2nd[3,free_hen] = (
                    transpose(
                        kron(vcat(0,axes_rot_hen.bitangent),rot_egg_n)
                    )*Y_hen
                )[:,uci_hen]
            A_2nd[3,free_egg] = (
                    transpose(
                        kron(vcat(0,axes_rot_egg.normal),rot_hen_b)
                    )*Y_egg
                )[:,uci_egg]
        end
        # cstr values
        ret = vcat(
            A_1st[mask_1st,:],
            A_4th[mask_4th,:],
            A_3rd[mask_3rd,:], 
            A_2nd[mask_2nd,:],
        )
        ret
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

function make_cstr_hessians(cst::PrototypeJoint,st::Structure)
    (;numbered) = st.connectivity
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = numbered
    (;
        hen2egg,
        axes_trl_hen,axes_trl_egg,
        axes_rot_hen,axes_rot_egg,
        mask_1st,mask_2nd,
        mask_3rd,mask_4th,
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    nld_hen = get_num_of_local_dims(nmcs_hen)
    nld_egg = get_num_of_local_dims(nmcs_egg)
    Y_hen = nmcs_hen.conversion_to_std
    Y_egg = nmcs_egg.conversion_to_std
    cv = BlockDiagonal([Y_hen,Y_egg])
    ndim = get_num_of_dims(nmcs_hen)
    T = get_numbertype(cst)
    I_Int = I(ndim)
    # translate
    c = get_local_coords(st)
    c_hen = c[sys_loci2coords_idx[bodyid2sys_loci_idx[id_hen][hen.pid]]]
    c_egg = c[sys_loci2coords_idx[bodyid2sys_loci_idx[id_egg][egg.pid]]]
    C_hen = to_transformation(hen.rbsig.state.cache.funcs.nmcs,c_hen)
    C_egg = to_transformation(egg.rbsig.state.cache.funcs.nmcs,c_egg)
    # first
    cstr_hessians_1st = [
        begin
            ret = kron(
                zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                ),
                I_Int
            )
        end
        for i = 1:3
    ]
    cstr_hessians_4th = [
        begin
            ret = kron(
                zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                ),
                I_Int
            )
            ret[1:(1+nld_hen)*ndim,1:(1+nld_hen)*ndim] = 2C_hen'C_hen
            ret[(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim,1:(1+nld_hen)*ndim] = -2C_egg'C_hen
            ret[1:(1+nld_hen)*ndim,(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim] = -2C_hen'C_egg
            ret[(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim,(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim] = 2C_egg'C_egg
            ret
        end
    ]
    cstr_hessians_3rd = fill(zero(cstr_hessians_4th[1]),6)
    cstr_hessians_2nd = fill(zero(cstr_hessians_4th[1]),3)
    if (nmcs_hen isa NCF.NC3D12C) && (nmcs_egg isa NCF.NC3D12C)
        cstr_hessians_3rd[1:3] = [
            begin
                d̄ = axes_trl_hen.X[:,i]
                ret_raw = zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                )
                ret_raw[1,2:1+nld_hen] = -d̄
                ret_raw[2:1+nld_hen,1] = -d̄
                ret_raw[2:1+nld_hen,2:1+nld_hen] = -c_hen*d̄'-d̄*c_hen'
                ret_raw[(1+nld_hen)+1,2:1+nld_egg] = d̄
                ret_raw[2:1+nld_egg,(1+nld_hen)+1] = d̄
                ret_raw[(1+nld_hen)+2:end,2:1+nld_egg] = c_egg*d̄'
                ret_raw[2:1+nld_egg,(1+nld_hen)+2:end] = d̄*c_egg'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
                ret                
            end
            for i = 1:3
        ]
        cstr_hessians_3rd[4:6] = [
            begin
                ret_raw = zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                )
                d̄ = axes_trl_egg.X[:,i]
                ret_raw[1,2:1+nld_hen] = -d̄
                ret_raw[2:1+nld_hen,1] = -d̄
                ret_raw[(1+nld_hen)+1,2:1+nld_egg] = d̄
                ret_raw[2:1+nld_egg,(1+nld_hen)+1] = d̄
                ret_raw[(1+nld_hen)+2:end,2:1+nld_egg] = -c_hen*d̄'
                ret_raw[2:1+nld_egg,(1+nld_hen)+2:end] = -d̄*c_hen'
                ret_raw[(1+nld_hen)+2:end,(1+nld_hen)+2:end] = c_egg*d̄'+d̄*c_hen'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
                ret
            end
            for i = 1:3
        ]
        cstr_hessians_2nd .= [
            begin
                ret_raw = zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                )
                ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.tangent*axes_rot_egg.bitangent'
                ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.bitangent*axes_rot_hen.tangent'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
                ret
            end,
            begin
                ret_raw = zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                )
                ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.tangent*axes_rot_egg.normal'
                ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.normal*axes_rot_hen.tangent'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
                ret
            end,
            begin
                ret_raw = zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                )
                ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.bitangent*axes_rot_egg.normal'
                ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.normal*axes_rot_hen.bitangent'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
                ret
            end,
        ]
    end
    cstr_hessians = vcat(
        cstr_hessians_1st[mask_1st],
        cstr_hessians_4th[mask_4th],
        cstr_hessians_3rd[mask_3rd],
        cstr_hessians_2nd[mask_2nd],
    )
    cstr_hessians
end

function get_jointed_free_idx(cst)
    (;num_of_cstr,hen2egg,) = cst
    (;hen,egg) = hen2egg
    uci_hen = hen.rbsig.state.cache.free_idx
    uci_egg = egg.rbsig.state.cache.free_idx
    ncoords1 = NCF.get_num_of_coords(hen.rbsig.state.cache.funcs.nmcs)
    # ncoords2 = NCF.get_num_of_coords(egg.rbsig.state.cache.nmcs)
    free_coords_idx = vcat(
        uci_hen,
        uci_egg .+ ncoords1
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

function make_cstr_forces_jacobian(cst,st)
    (;num_of_cstr) = cst
    free_coords_idx = get_jointed_free_idx(cst)
    cstr_hessians = make_cstr_hessians(cst,st)
    function cstr_forces_jacobian(λ)
        ret = [
            begin
                a = cstr_hessians[i][free_coords_idx,free_coords_idx] .* λ[i]
                # display(a)
                a 
            end
            for i = 1:num_of_cstr
        ]
        sum(ret)
    end
end

function get_jointed_free_idx(cst::LinearJoint)
    free_coords_idx = collect(1:size(cst.A,2))
end

function get_jointed_free(cst::LinearJoint,indexed)
    cst2sysfree = collect(1:indexed.num_of_free_coords)
end

function make_cstr_forces_jacobian(cst::LinearJoint,st)
    free_coords_idx = get_jointed_free_idx(cst)
    function cstr_forces_jacobian(λ)
        zeros(eltype(λ),length(free_coords_idx),length(free_coords_idx))
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