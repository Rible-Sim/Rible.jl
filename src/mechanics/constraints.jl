
function make_constraints_function(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;mem2sysfull) = indexed
    (;num2sys,mem2num) = numbered
    (;
        nconstraints,e2e,values,
        axes_trl_hen,axes_trl_egg,axes_rot_hen,axes_rot_egg,
        mask_1st, mask_2nd, mask_3rd_hen, mask_3rd_egg, mask_4th
    ) = cst
    (;hen,egg) = e2e
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    function _inner_constraints_function(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        R_hen = NCF.find_rotation(nmcs_hen,q_hen)
        R_egg = NCF.find_rotation(nmcs_egg,q_egg)
        # translate
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        # @show C_egg,q_egg
        r_hen2egg = C_egg*q_egg .- C_hen*q_hen
        # translate on hen
        trl_hen_n = R_hen*axes_trl_hen.n
        trl_hen_t1 = R_hen*axes_trl_hen.t1
        trl_hen_t2 = R_hen*axes_trl_hen.t2
        # translate on egg
        trl_egg_n = R_egg*axes_trl_egg.n
        trl_egg_t1 = R_egg*axes_trl_egg.t1
        trl_egg_t2 = R_egg*axes_trl_egg.t2
        # rotate of egg
        rot_hen_t1 = R_hen*axes_rot_hen.t1
        rot_hen_t2 = R_hen*axes_rot_hen.t2
        rot_egg_n = R_egg*axes_rot_egg.n
        rot_egg_t2 = R_egg*axes_rot_egg.t2
        # first    
        Φ_1st = r_hen2egg
        # second
        Φ_2nd = [
            rot_hen_t1'*rot_egg_t2,
            rot_hen_t1'*rot_egg_n,
            rot_hen_t2'*rot_egg_n
        ]
        # third
        Φ_3rd_hen = [
            trl_hen_n'*r_hen2egg,
            trl_hen_t1'*r_hen2egg,
            trl_hen_t2'*r_hen2egg,
        ]
        Φ_3rd_egg = [
            trl_egg_n'*r_hen2egg,
            trl_egg_t1'*r_hen2egg,
            trl_egg_t2'*r_hen2egg,
        ]
        # fourth    
        Φ_4th = [r_hen2egg'*r_hen2egg]
        # constraint values
        ret = vcat(
            Φ_4th[mask_4th],
            Φ_3rd_hen[mask_3rd_hen], # translate prior to 
            Φ_3rd_egg[mask_3rd_egg], 
            Φ_1st[mask_1st],
            Φ_2nd[mask_2nd], # rotate
        ) .- values
        ret
    end
    function inner_constraints_function(q,d,c)
        _inner_constraints_function(q,d,c)
    end
    function inner_constraints_function(q)
        d = get_d(st)
        c = get_c(st)
        _inner_constraints_function(q,d,c)
    end
    inner_constraints_function
end

function make_constraints_jacobian(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    (;
        nconstraints,e2e,
        axes_trl_hen,axes_trl_egg,axes_rot_hen,axes_rot_egg, 
        mask_1st, mask_2nd, mask_3rd_hen, mask_3rd_egg, mask_4th
    ) = cst
    (;hen,egg) = e2e
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    Y_hen = NCF.get_conversion(nmcs_hen)
    Y_egg = NCF.get_conversion(nmcs_egg)
    function _inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        R_hen = NCF.find_rotation(nmcs_hen,q_hen)
        R_egg = NCF.find_rotation(nmcs_egg,q_egg)
        # translate 
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        r_hen2egg = C_egg*q_egg .- C_hen*q_hen
        # translate on hen
        trl_hen_n = R_hen*axes_trl_hen.n
        trl_hen_t1 = R_hen*axes_trl_hen.t1
        trl_hen_t2 = R_hen*axes_trl_hen.t2
        trl_hen = [trl_hen_n trl_hen_t1 trl_hen_t2]
        # translate on egg
        trl_egg_n = R_egg*axes_trl_egg.n
        trl_egg_t1 = R_egg*axes_trl_egg.t1
        trl_egg_t2 = R_egg*axes_trl_egg.t2
        trl_egg = [trl_egg_n trl_egg_t1 trl_egg_t2]
        # rotate of egg
        rot_hen_t1 = R_hen*axes_rot_hen.t1
        rot_hen_t2 = R_hen*axes_rot_hen.t2
        rot_egg_n = R_egg*axes_rot_egg.n
        rot_egg_t2 = R_egg*axes_rot_egg.t2
        # jac
        o3 = zeros(eltype(q),3)
        # jac first
        A_1st = zeros(eltype(q),3,nfree)
        A_1st[:,free_hen] .= -C_hen[:,uci_hen]
        A_1st[:,free_egg] .=  C_egg[:,uci_egg]
        # jac third on hen
        A_3rd_hen = zeros(eltype(q),3,nfree)
        A_3rd_hen[:,free_hen] .= kron(
                hcat(
                    o3, transpose(trl_hen)
                ),
                transpose(r_hen2egg)
            )[:,uci_hen]
        A_3rd_hen[:,free_hen] .-= transpose(trl_hen)*C_hen[:,uci_hen]
        A_3rd_hen[:,free_egg] .+= transpose(trl_hen)*C_egg[:,uci_egg]
        # jac third on egg
        A_3rd_egg = zeros(eltype(q),3,nfree)
        A_3rd_egg[:,free_egg] .= kron(
                hcat(
                    o3, transpose(trl_egg)
                ),
                transpose(r_hen2egg)
            )[:,uci_egg]
        A_3rd_egg[:,free_hen] .-= transpose(trl_egg)*C_hen[:,uci_hen]
        A_3rd_egg[:,free_egg] .+= transpose(trl_egg)*C_egg[:,uci_egg]
        # jac 2nd
        A_2nd = zeros(eltype(q),3,nfree)
        if (nmcs_hen isa NCF.LNC3D12C) && (nmcs_egg isa NCF.LNC3D12C) 
            A_2nd[1,free_hen] = (
                    transpose(
                        kron(vcat(0,axes_rot_hen.t1),rot_egg_t2)
                    )*Y_hen
                )[:,uci_hen]
            A_2nd[1,free_egg] = (
                    transpose(
                        kron(vcat(0,axes_rot_egg.t2),rot_hen_t1)
                    )*Y_egg
                )[:,uci_egg]
            A_2nd[2,free_hen] = (
                    transpose(
                        kron(vcat(0,axes_rot_hen.t1),rot_egg_n)
                    )*Y_hen
                )[:,uci_hen]
            A_2nd[2,free_egg] = (
                    transpose(
                        kron(vcat(0,axes_rot_egg.n),rot_hen_t1)
                    )*Y_egg
                )[:,uci_egg]
            A_2nd[3,free_hen] = (
                    transpose(
                        kron(vcat(0,axes_rot_hen.t2),rot_egg_n)
                    )*Y_hen
                )[:,uci_hen]
            A_2nd[3,free_egg] = (
                    transpose(
                        kron(vcat(0,axes_rot_egg.n),rot_hen_t2)
                    )*Y_egg
                )[:,uci_egg]
        end
        # jac 4th
        A_4th = zeros(eltype(q),1,nfree)
        A_4th[1,free_hen] = -2r_hen2egg'*C_hen[:,uci_hen]
        A_4th[1,free_egg] =  2r_hen2egg'*C_egg[:,uci_egg]
        # constraint values
        ret = vcat(
            A_4th[mask_4th,:],
            A_3rd_hen[mask_3rd_hen,:], # translate prior to 
            A_3rd_egg[mask_3rd_egg,:], 
            A_1st[mask_1st,:],
            A_2nd[mask_2nd,:], # rotate
        )
        ret
    end
    function inner_A(q,c)
        _inner_A(q,c)
    end
    function inner_A(q)
        c = get_c(st)
        _inner_A(q,c)
    end
    inner_A
end

function make_constraints_functionqᵀq(cst::PrototypeJoint,st::Structure)
    (;numbered) = st.connectivity
    (;mem2num,num2sys) = numbered
    (;
        e2e,
        axes_trl_hen,axes_trl_egg,axes_rot_hen,axes_rot_egg,
        mask_1st,mask_2nd,mask_3rd_hen,mask_3rd_egg,mask_4th,
    ) = cst
    (;hen,egg) = e2e
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    nld_hen = NCF.get_num_of_local_dims(nmcs_hen)
    nld_egg = NCF.get_num_of_local_dims(nmcs_egg)
    Y_hen = NCF.get_conversion(nmcs_hen)
    Y_egg = NCF.get_conversion(nmcs_egg)
    cv = BlockDiagonal([Y_hen,Y_egg])
    ndim = NCF.get_num_of_dims(nmcs_hen)
    T = get_numbertype(cst)
    I_Int = NCF.make_I(Int,ndim)
    # translate
    c = get_c(st)
    c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
    c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
    C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
    C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
    # first
    constraints_hessian_1st =[
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
    constraints_hessian_3rd_hen = [
        begin
            d̄ = axes_trl_hen.X[:,i]
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            if (nmcs_hen isa NCF.LNC3D12C) && (nmcs_egg isa NCF.LNC3D12C)
                ret_raw[1,2:1+nld_hen] = -d̄
                ret_raw[2:1+nld_hen,1] = -d̄
                ret_raw[2:1+nld_hen,2:1+nld_hen] = -c_hen*d̄'-d̄*c_hen'
                ret_raw[(1+nld_hen)+1,2:1+nld_egg] = d̄
                ret_raw[2:1+nld_egg,(1+nld_hen)+1] = d̄
                ret_raw[(1+nld_hen)+2:end,2:1+nld_egg] = c_egg*d̄'
                ret_raw[2:1+nld_egg,(1+nld_hen)+2:end] = d̄*c_egg'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
            end
            ret                
        end
        for i = 1:3
    ]
    constraints_hessian_3rd_egg = [
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            if (nmcs_hen isa NCF.LNC3D12C) && (nmcs_egg isa NCF.LNC3D12C)
                d̄ = axes_trl_egg.X[:,i]
                ret_raw[1,2:1+nld_hen] = -d̄
                ret_raw[2:1+nld_hen,1] = -d̄
                ret_raw[(1+nld_hen)+1,2:1+nld_egg] = d̄
                ret_raw[2:1+nld_egg,(1+nld_hen)+1] = d̄
                ret_raw[(1+nld_hen)+2:end,2:1+nld_egg] = -c_hen*d̄'
                ret_raw[2:1+nld_egg,(1+nld_hen)+2:end] = -d̄*c_hen'
                ret_raw[(1+nld_hen)+2:end,(1+nld_hen)+2:end] = c_egg*d̄'+d̄*c_hen'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
            end
            ret
        end
        for i = 1:3
    ]
    constraints_hessian_2nd = [
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            if (nmcs_hen isa NCF.LNC3D12C) && (nmcs_egg isa NCF.LNC3D12C)
                ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.t1*axes_rot_egg.t2'
                ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.t2*axes_rot_hen.t1'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
            end
            ret
        end,
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            if (nmcs_hen isa NCF.LNC3D12C) && (nmcs_egg isa NCF.LNC3D12C)
                ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.t1*axes_rot_egg.n'
                ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.n*axes_rot_hen.t1'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
            end
            ret
        end,
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            if (nmcs_hen isa NCF.LNC3D12C) && (nmcs_egg isa NCF.LNC3D12C)
                ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.t2*axes_rot_egg.n'
                ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.n*axes_rot_hen.t2'
                ret = transpose(cv)*kron(ret_raw,I_Int)*cv
            end
            ret
        end,
    ]    
    constraints_hessian_4th = [
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
    constraints_hessian = vcat(
        constraints_hessian_4th[mask_4th],
        constraints_hessian_3rd_hen[mask_3rd_hen],
        constraints_hessian_3rd_egg[mask_3rd_egg],
        constraints_hessian_1st[mask_1st],
        constraints_hessian_2nd[mask_2nd],
    )
    constraints_hessian
end

function get_jointed_free_idx(cst)
    (;nconstraints,e2e,) = cst
    (;hen,egg) = e2e
    uci_hen = hen.rbsig.state.cache.free_idx
    uci_egg = egg.rbsig.state.cache.free_idx
    ncoords1 = NCF.get_num_of_coordinates(hen.rbsig.state.cache.funcs.nmcs)
    # ncoords2 = NCF.get_num_of_coordinates(egg.rbsig.state.cache.nmcs)
    unconstrained_indices = vcat(
        uci_hen,
        uci_egg .+ ncoords1
    )
end

function get_jointed_free(cst,indexed)
    (;nconstraints,e2e,) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;hen,egg) = e2e
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    mem2sysfree1 = mem2sysfree[id_hen]
    mem2sysfree2 = mem2sysfree[id_egg]
    cst2sysfree = vcat(
        mem2sysfree1,
        mem2sysfree2
    )
end

function make_∂Aᵀλ∂q(cst,st)
    (;nconstraints) = cst
    unconstrained_indices = get_jointed_free_idx(cst)
    constraints_hessian = make_constraints_functionqᵀq(cst,st)
    function ∂Aᵀλ∂q(λ)
        ret = [
            begin
                a = constraints_hessian[i][unconstrained_indices,unconstrained_indices] .* λ[i]
                # display(a)
                a 
            end
            for i = 1:nconstraints
        ]
        sum(ret)
    end
end

function get_jointed_free_idx(cst::LinearJoint)
    unconstrained_indices = collect(1:size(cst.A,2))
end

function get_jointed_free(cst::LinearJoint,indexed)
    cst2sysfree = collect(1:indexed.nfree)
end

function make_∂Aᵀλ∂q(cst::LinearJoint,st)
    unconstrained_indices = get_jointed_free_idx(cst)
    function ∂Aᵀλ∂q(λ)
        zeros(eltype(λ),length(unconstrained_indices),length(unconstrained_indices))
    end
end


function get_jointed_free(cst::FixedIndicesConstraint,indexed)
    cst2sysfree = collect(1:indexed.nfree)
end

function make_∂Aᵀλ∂q(cst::FixedIndicesConstraint,st)
    (;nfree) = st.connectivity.indexed
    function ∂Aᵀλ∂q(λ)
        zeros(eltype(λ),nfree,nfree)
    end
end