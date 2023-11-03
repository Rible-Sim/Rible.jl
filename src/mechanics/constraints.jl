
function make_cstr_function(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;bodyid2sys_full_coords) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,values,
        half_1st, half_2nd, 
        half_3rd, half_4th,
        transformations,
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    function _inner_cstr_function(q,c)
        T = eltype(q)
        ret = zeros(T,num_of_cstr)
        q_hen = @view q[bodyid2sys_full_coords[id_hen]]
        q_egg = @view q[bodyid2sys_full_coords[id_egg]]
        q = vcat(q_hen,q_egg)
        # cstr values
        Refq = Ref(q)
        RefqT = Ref(q')
        Φ_1st = transformations*q
        Φ_4th = RefqT.*half_4th.*Refq
        Φ_3rd = RefqT.*half_3rd.*Refq
        Φ_2nd = RefqT.*half_2nd.*Refq
        ret = vcat(
            Φ_1st,
            Φ_4th,
            Φ_3rd, 
            Φ_2nd, 
        ) .- values
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
        hess_1st, hess_2nd, 
        hess_3rd, hess_4th,
        transformations,
    ) = cst
    (;hen,egg) = hen2egg
    free_idx_hen =  hen.rbsig.state.cache.free_idx
    free_idx_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
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
        J = transformations
        RefqT = Ref(q')
        A_zeros = spzeros(T,0,num_of_jointed_coords)
        A_1st = J
        A_4th = reduce(vcat,RefqT.*hess_4th;init = A_zeros)
        A_3rd = reduce(vcat,RefqT.*hess_3rd;init = A_zeros)
        A_2nd = reduce(vcat,RefqT.*hess_2nd;init = A_zeros)
        A = vcat(
            A_1st,
            A_4th,
            A_3rd, 
            A_2nd,
        )
        ret[:,free_hen] .= A[:, free_idx_hen]
        ret[:,free_egg] .= A[:,(free_idx_egg.+num_of_coords_hen)]
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

function get_jointed_free_idx(cst)
    (;num_of_cstr,hen2egg,) = cst
    (;hen,egg) = hen2egg
    uci_hen = hen.rbsig.state.cache.free_idx
    uci_egg = egg.rbsig.state.cache.free_idx
    ncoords1 = NCF.get_num_of_coords(hen.rbsig.state.cache.funcs.nmcs)
    # ncoords2 = NCF.get_num_of_coords(egg.rbsig.state.cache.nmcs)
    free_idx = vcat(
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
    (;
        num_of_cstr,
        hess_1st, hess_2nd, 
        hess_3rd, hess_4th,
    ) = cst
    free_idx = get_jointed_free_idx(cst)
    hessians = vcat(
        hess_1st,
        hess_4th,
        hess_3rd,
        hess_2nd
    )
    function cstr_forces_jacobian(λ)
        ret = [
            begin
                a = hessians[i][free_idx,free_idx] .* λ[i]
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