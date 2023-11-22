
function cstr_forces_jacobian(st::AbstractStructure,q,λ)
    (;numbered,indexed,jointed) = st.connectivity
    (;num_of_free_coords,
    num_of_intrinsic_cstr,
    bodyid2sys_free_coords,
    bodyid2sys_intrinsic_cstr_idx) = indexed
    (;njoints,num_of_extrinsic_cstr,joints,jointid2sys_extrinsic_cstr_idx) = jointed
    ret = zeros(eltype(λ),num_of_free_coords,num_of_free_coords)
    (;bodies,num_of_cstr) = st
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
        foreach(joints) do joint
            joint_cstr_idx = num_of_intrinsic_cstr .+ jointid2sys_extrinsic_cstr_idx[joint.id]
            jointed_sys_free_idx = joint.sys_free_idx
            ret[jointed_sys_free_idx,jointed_sys_free_idx] .+= cstr_forces_jacobian(joint,st,q,λ[joint_cstr_idx])
        end
    end
    ret
end

function ∂Aq̇∂q(st,q̇)
    (;num_of_free_coords) = st.connectivity.indexed
    (;bodies,num_of_cstr) = st
    (;indexed,jointed) = st.connectivity
    (;num_of_intrinsic_cstr,bodyid2sys_free_coords,bodyid2sys_intrinsic_cstr_idx) = indexed
    ret = zeros(eltype(q̇),num_of_cstr,num_of_free_coords)
    foreach(bodies) do body
        bodyid = body.prop.id
        memfree = bodyid2sys_free_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        free_idx = body.state.cache.free_idx
        cstr_idx = body.state.cache.cstr_idx
        if !isempty(memincst)
            ret[memincst,memfree] .+= body.state.cache.funcs.∂Aq̇∂q(q̇[memfree])[cstr_idx,free_idx]
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
线性化。
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
    (;num_of_full_coords,bodyid2sys_full_coords,sys_num_of_dof,bodyid2sys_dof_idx,) = indexed
    ret = zeros(eltype(q),num_of_full_coords,sys_num_of_dof)
    foreach(bodies) do body
        bodyid = body.prop.id
        (;nmcs) = body.coords
        mem2full = bodyid2sys_full_coords[bodyid]
        ret[mem2full,bodyid2sys_dof_idx[bodyid]] = nullspace_mat(nmcs,q[mem2full])
    end
    ret
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

function build_Ǩ(st)
    _,λ = check_static_equilibrium_output_multipliers(st)
    build_Ǩ(st,λ)
end

function build_material_stiffness_matrix!(st::Structure,q,k)
    (;num_of_dim) = st
    (;indexed,tensioned) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_full_coords) = indexed
    (;connected) = tensioned
    (;cables) = st.tensiles
    update!(st,q)
    Jj = zeros(eltype(q),num_of_dim,num_of_full_coords)
    retǨm = zeros(eltype(q),num_of_free_coords,num_of_free_coords)
    foreach(connected) do scnt
        j = scnt.id
        rb1 = scnt.hen.rbsig
        rb2 = scnt.egg.rbsig
        ap1id = scnt.hen.pid
        ap2id = scnt.egg.pid
        C1 = rb1.cache.Cps[ap1id]
        C2 = rb2.cache.Cps[ap2id]
        mfull1 = bodyid2sys_full_coords[rb1.prop.id]
        mfull2 = bodyid2sys_full_coords[rb2.prop.id]
        cable = cables[j]
        (;state) = cable
        (;length,) = state
        s = 1/length
        Jj .= 0
        Jj[:,mfull2] .+= C2
        Jj[:,mfull1] .-= C1
        Uj = transpose(Jj)*Jj
        Ūjq = Uj[sys_free_idx,:]*q
        retǨm .+= k[j]*s^2*(Ūjq*transpose(Ūjq))
    end
    retǨm
end

function build_geometric_stiffness_matrix!(st::Structure,q,f)
    (;num_of_dim) = st
    (;indexed,tensioned) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_full_coords) = indexed
    (;connected) = tensioned
    (;cables) = st.tensiles
    update!(st,q)
    Jj = zeros(eltype(q),num_of_dim,num_of_full_coords)
    retǨg = zeros(eltype(q),num_of_free_coords,num_of_free_coords)
    foreach(connected) do scnt
        j = scnt.id
        rb1 = scnt.hen.rbsig
        rb2 = scnt.egg.rbsig
        ap1id = scnt.hen.pid
        ap2id = scnt.egg.pid
        C1 = rb1.cache.Cps[ap1id]
        C2 = rb2.cache.Cps[ap2id]
        mfull1 = bodyid2sys_full_coords[rb1.prop.id]
        mfull2 = bodyid2sys_full_coords[rb2.prop.id]
        cable = cables[j]
        (;state) = cable
        (;length,) = state
        s = 1/length
        Jj .= 0
        Jj[:,mfull2] .+= C2
        Jj[:,mfull1] .-= C1
        Uj = transpose(Jj)*Jj
        Ǔj = @view Uj[sys_free_idx,sys_free_idx]
        Ūjq = Uj[sys_free_idx,:]*q
        retǨg .+= f[j]/length*(Ǔj-s^2*Ūjq*transpose(Ūjq))
    end
    retǨg
end

function make_Ǩm_Ǩg(st,q0)
    (;num_of_dim) = st
    (;numbered,indexed,tensioned) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_pres_idx,sys_free_idx,bodyid2sys_full_coords) = indexed
    (;connected) = tensioned
    (;cables) = st.tensiles
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = numbered
    function inner_Ǩm_Ǩg(q̌,s,μ,k,c)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_idx] .= q0[sys_pres_idx]
		q[sys_free_idx] .= q̌
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        retǨm = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        retǨg = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.hen.rbsig
            rb2 = scnt.egg.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
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
            Ǔj = @view Uj[sys_free_idx,sys_free_idx]
            Ūjq = Uj[sys_free_idx,:]*q
            retǨm .+= k[j]*s[j]^2*(Ūjq*transpose(Ūjq))
            retǨg .+= k[j]*(1-μ[j]*s[j])*(Ǔj-s[j]^2*Ūjq*transpose(Ūjq))
        end
        retǨm,retǨg
    end
    function inner_Ǩm_Ǩg(q̌)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_idx] .= q0[sys_pres_idx]
		q[sys_free_idx] .= q̌
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        retǨm = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        retǨg = zeros(eltype(q̌),num_of_free_coords,num_of_free_coords)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.hen.rbsig
            rb2 = scnt.egg.rbsig
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
            C1 = rb1.state.cache.Cps[ap1id]
            C2 = rb2.state.cache.Cps[ap2id]
            mfull1 = bodyid2sys_full_coords[rb1.prop.id]
            mfull2 = bodyid2sys_full_coords[rb2.prop.id]
            cable = cables[j]
            (;k,c,state,slack) = cable
            (;direction,tension,length,lengthdot) = state
            s = 1/length
            Jj .= 0
            Jj[:,mfull2] .+= C2
            Jj[:,mfull1] .-= C1
            Uj = transpose(Jj)*Jj
            Ǔj = @view Uj[sys_free_idx,sys_free_idx]
            Ūjq = Uj[sys_free_idx,:]*q
            retǨm .+= k*s^2*(Ūjq*transpose(Ūjq))
            retǨg .+= tension/length*(Ǔj-s^2*Ūjq*transpose(Ūjq))
        end
        retǨm,retǨg
    end
end

function make_S(st,q0)
    (;num_of_dim) = st
    (;numbered,indexed,tensioned) = st.connectivity
    (;sys_pres_idx,sys_free_idx,num_of_full_coords,bodyid2sys_full_coords) = indexed
    (;bodyid2sys_loci_idx,sys_loci2coords_idx) = numbered
    (;connected) = tensioned
    (;cables) = st.tensiles
    ncables = length(cables)
    function inner_S(q̌,s)
		q = Vector{eltype(q̌)}(undef,num_of_full_coords)
		q[sys_pres_idx] .= q0[sys_pres_idx]
		q[sys_free_idx] .= q̌
        ret = zeros(eltype(q̌),ncables)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.hen.rbsig
            rb2 = scnt.egg.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
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
        q[sys_pres_idx] .= q0[sys_pres_idx]
        q[sys_free_idx] .= q̌
        ret = zeros(eltype(q̌),ncables)
        Jj = zeros(eltype(q̌),num_of_dim,num_of_full_coords)
        foreach(connected) do scnt
            j = scnt.id
            rb1 = scnt.hen.rbsig
            rb2 = scnt.egg.rbsig
            rb1id = rb1.prop.id
            rb2id = rb2.prop.id
            ap1id = scnt.hen.pid
            ap2id = scnt.egg.pid
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

# Out-of-place ∂Q̌∂q̌ (dispatch)
function build_∂Q̌∂q̌(st)
    build_∂Q̌∂q̌(st, st.connectivity.tensioned)
end

# Out-of-place ∂Q̌∂q̌ for cables and clustered cables
function build_∂Q̌∂q̌(st, @eponymargs(connected, clustered))
    ∂Q̌∂q̌1 = build_∂Q̌∂q̌(st, @eponymtuple(connected))
    ∂Q̌∂q̌2 = build_∂Q̌∂q̌(st, @eponymtuple(clustered))
    return ∂Q̌∂q̌1 + ∂Q̌∂q̌2
end

# Out-of-place ∂Q̌∂q̌ for cables
function build_∂Q̌∂q̌(st,@eponymargs(connected,))
    (;cables) = st.tensiles
    (;indexed) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ∂Q̌∂q̌ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;hen,egg) = cc
        rb1 = hen.rbsig
        rb2 = egg.rbsig
        C1 = rb1.state.cache.Cps[hen.pid]
        C2 = rb2.state.cache.Cps[egg.pid]
        free_idx1 = rb1.state.cache.free_idx
        free_idx2 = rb2.state.cache.free_idx
        mfree1 = bodyid2sys_free_coords[rb1.prop.id]
        mfree2 = bodyid2sys_free_coords[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension,length,lengthdot) = state
        if slack && (tension==0)
            ∂Q̌∂q̌ .= 0
        else
            D .= direction*transpose(direction)
            density = tension/length
            β = c*lengthdot/length + density
            D .*= k-β
            D .+= β.*Im
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]
            ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌
end

# Out-of-place ∂Q̌∂q̌ for cluster cables
function build_∂Q̌∂q̌(st,@eponymargs(clustered))
    (;clustercables) = st.tensiles
    (;indexed) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ∂Q̌∂q̌ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    i = 0
    foreach(clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            cable = clustercables[i].segs[cc.id]
            (;hen,egg) = cc
            rb1 = hen.rbsig
            rb2 = egg.rbsig
            C1 = rb1.state.cache.Cps[hen.pid]
            C2 = rb2.state.cache.Cps[egg.pid]
            free_idx1 = rb1.state.cache.free_idx
            free_idx2 = rb2.state.cache.free_idx
            mfree1 = bodyid2sys_free_coords[rb1.prop.id]
            mfree2 = bodyid2sys_free_coords[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension,length,lengthdot) = state
            if tension==0
                ∂Q̌∂q̌ .-= 0
            else
                D .= direction*transpose(direction)
                density = tension/length
                β = c*lengthdot/length + density
                D .*= k-β
                D .+= β.*Im
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,free_idx2]
                J̌[:,mfree1] .-= C1[:,free_idx1]
                ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
            end
        end
    end
    return ∂Q̌∂q̌
end

# In-place ∂Q̌∂q̌ for cables and flexible bodies
function build_∂Q̌∂q̌!(∂Q̌∂q̌,st)
    (;bodies,connectivity) = st
    (;tensioned,indexed) = connectivity
    (;cables) = st.tensiles
    (;connected) = tensioned
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    # ∂Q̌∂q̌ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;hen,egg) = cc
        rb1 = hen.rbsig
        rb2 = egg.rbsig
        C1 = rb1.cache.Cps[hen.pid]
        C2 = rb2.cache.Cps[egg.pid]
        free_idx1 = rb1.coords.free_idx
        free_idx2 = rb2.coords.free_idx
        mfree1 = bodyid2sys_free_coords[rb1.prop.id]
        mfree2 = bodyid2sys_free_coords[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension,length,lengthdot) = state
        if slack && (tension==0)
            ∂Q̌∂q̌ .-= 0
        else
            D .= direction*transpose(direction)
            density = tension/length
            β = c*lengthdot/length + density
            D .*= k-β
            D .+= β.*Im
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]
            ∂Q̌∂q̌ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end

    foreach(bodies) do body
        if body isa FlexibleBody
            (;cache) = body.state
            (;e,funcs) = cache
            (;ancs) = funcs
            ∂Q∂e = ANCF.make_∂Q∂e(ancs)(e)
            mfree = bodyid2sys_free_coords[body.prop.id]
            free_idx = body.state.cache.free_idx
            ∂Q̌∂q̌[mfree,mfree] .-= ∂Q∂e[free_idx,free_idx]
        end
    end

    return ∂Q̌∂q̌
end

# Out-of-place ∂Q̌∂q̌̇ (dispatch)
function build_∂Q̌∂q̌̇(st)
    build_∂Q̌∂q̌̇(st, st.connectivity.tensioned)
end

# Out-of-place ∂Q̌∂q̌̇ for cables and clustered cables
function build_∂Q̌∂q̌̇(st, @eponymargs(connected, clustered))
    ∂Q̌∂q̌̇1 = build_∂Q̌∂q̌̇(st, @eponymtuple(connected))
    ∂Q̌∂q̌̇2 = build_∂Q̌∂q̌̇(st, @eponymtuple(clustered))
    return ∂Q̌∂q̌̇1 + ∂Q̌∂q̌̇2
end

# Out-of-place ∂Q̌∂q̌̇ for cables
function build_∂Q̌∂q̌̇(st, @eponymargs(connected, ))
    (;cables) = st.tensiles
    (;indexed) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ∂Q̌∂q̌̇ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;hen,egg) = cc
        rb1 = hen.rbsig
        rb2 = egg.rbsig
        C1 = rb1.state.cache.Cps[hen.pid]
        C2 = rb2.state.cache.Cps[egg.pid]
        free_idx1 = rb1.state.cache.free_idx
        free_idx2 = rb2.state.cache.free_idx
        mfree1 = bodyid2sys_free_coords[rb1.prop.id]
        mfree2 = bodyid2sys_free_coords[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension) = state
        if slack && (tension == 0)
            ∂Q̌∂q̌̇ .-= 0
        else
            D .= direction*transpose(direction)
            D .*= c
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]

            ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfrjiexee2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
    return ∂Q̌∂q̌̇
end

# Out-of-place ∂Q̌∂q̌̇ for clustered cables
function build_∂Q̌∂q̌̇(st, @eponymargs(clustered, ))
    (;clustercables) = st.tensiles
    (;indexed) = st.connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ∂Q̌∂q̌̇ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    i = 0
    foreach(clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            cable = clustercables[i].segs[cc.id]
            (;hen,egg) = cc
            rb1 = hen.rbsig
            rb2 = egg.rbsig
            C1 = rb1.state.cache.Cps[hen.pid]
            C2 = rb2.state.cache.Cps[egg.pid]
            free_idx1 = rb1.state.cache.free_idx
            free_idx2 = rb2.state.cache.free_idx
            mfree1 = bodyid2sys_free_coords[rb1.prop.id]
            mfree2 = bodyid2sys_free_coords[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension) = state
            if tension == 0
                ∂Q̌∂q̌̇ .-= 0
            else
                D .= direction*transpose(direction)
                D .*= c
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,free_idx2]
                J̌[:,mfree1] .-= C1[:,free_idx1]

                ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
            end
        end
    end
    return ∂Q̌∂q̌̇
end

# In-place ∂Q̌∂q̌̇ for cables
function build_∂Q̌∂q̌̇!(∂Q̌∂q̌̇,st)
    (;tensioned,indexed) = st.connectivity
    (;connected) = tensioned
    (;cables) = st.tensiles
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    # ∂Q̌∂q̌̇ = zeros(T,num_of_free_coords,num_of_free_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)
    foreach(connected) do cc
        cable = cables[cc.id]
        (;hen,egg) = cc
        rb1 = hen.rbsig
        rb2 = egg.rbsig
        C1 = rb1.cache.Cps[hen.pid]
        C2 = rb2.cache.Cps[egg.pid]
        free_idx1 = rb1.coords.free_idx
        free_idx2 = rb2.coords.free_idx
        mfree1 = bodyid2sys_free_coords[rb1.prop.id]
        mfree2 = bodyid2sys_free_coords[rb2.prop.id]
        (;k,c,state,slack) = cable
        (;direction,tension) = state
        if slack && (tension == 0)
            ∂Q̌∂q̌̇ .-= 0
        else
            D .= direction*transpose(direction)
            D .*= c
            J̌ .= 0
            J̌[:,mfree2] .+= C2[:,free_idx2]
            J̌[:,mfree1] .-= C1[:,free_idx1]

            ∂Q̌∂q̌̇ .-= transpose(J̌)*D*J̌
        end
        # ∂Q̌∂q̌_full[mfree2,mfree2] .+= transpose(C2)*D*C2
        # ∂Q̌∂q̌_full[mfree1,mfree2] .-= transpose(C1)*D*C2
        # ∂Q̌∂q̌_full[mfree2,mfree1] .-= transpose(C2)*D*C1
        # ∂Q̌∂q̌_full[mfree1,mfree1] .+= transpose(C1)*D*C1
    end
end

function build_∂Q̌∂s̄(st)
    (;connectivity) = st
    (;cables,clustercables) = st.tensiles
    nclustercables = length(clustercables)
    (;tensioned,indexed) = connectivity
    (;num_of_full_coords,num_of_free_coords,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_full_coords) = indexed
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    ∂Q̌∂s̄ = zeros(T,2ns,num_of_free_coords)
    D = zeros(T, num_of_dim)
    lkn = zeros(T, 2ns, num_of_dim)
    # Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J̌ = zeros(T,num_of_dim,num_of_free_coords)

    N_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    kc = Vector{Float64}()
    for (cid,clustercable) in enumerate(clustercables)
        (;segs) = clustercable
        nsegs = length(segs)
        for (sid, seg) in enumerate(segs)
            push!(kc, seg.k)
        end
        N = sparse(zeros(Float64, nsegs, 2nsegs-2))
        N[1,1:2]=[1 -1]; N[end,end-1:end]=[-1 1]
        for i in 2:nsegs-1
            N[i, 2i-3:2i] = [-1 1 1 -1]
        end
        push!(N_list, N)
    end
    N = reduce(blockdiag,N_list)
    i = 0; j = 0
    foreach(tensioned.clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            j += 1
            cable = clustercables[i].segs[cc.id]
            (;hen,egg) = cc
            rb1 = hen.rbsig
            rb2 = egg.rbsig
            C1 = rb1.state.cache.Cps[hen.pid]
            C2 = rb2.state.cache.Cps[egg.pid]
            free_idx1 = rb1.state.cache.free_idx
            free_idx2 = rb2.state.cache.free_idx
            mfree1 = bodyid2sys_free_coords[rb1.prop.id]
            mfree2 = bodyid2sys_free_coords[rb2.prop.id]
            (;k,c,state) = cable
            (;direction,tension) = state
            if tension == 0
                ∂Q̌∂s̄ .-= 0
            else
                D .= direction
                J̌ .= 0
                J̌[:,mfree2] .+= C2[:,free_idx2]
                J̌[:,mfree1] .-= C1[:,free_idx1]
                kN = kc[j] .* N[j,:]
                @tullio lkn[k, l] = D[l] * kN[k]
                ∂Q̌∂s̄ .-= lkn * J̌
            end
        end
    end
    return ∂Q̌∂s̄'
end

function build_Ǩ(st,λ)
    (;num_of_free_coords) = st.connectivity.indexed
    T = get_numbertype(st)
    # Ǩ = zeros(T,num_of_free_coords,num_of_free_coords)
    Ǩ = -build_∂Q̌∂q̌(st) .- cstr_forces_jacobian(st,λ)
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
    Ǩ = build_Ǩ(st,λ)
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

function undamped_eigen!(bot::Robot;gravity=false,scaling=0.01)
    (;st,traj) = bot
    q̌ = get_free_coords(st)
    ω²,δq̌ = undamped_eigen(st;gravity)
    neg_idx = findall(ω².<=0)
    if !isempty(neg_idx)
        @warn "Negative ω² occurs, idx $neg_idx, zeroing."
        ω²[neg_idx] .= 0
    end
    ω = sqrt.(ω²)
    resize!(traj,1)
    nω = length(ω)
    for i = 1:nω
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = ω[i]
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    bot
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
返回零空间。
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

"""
校核稳定性。
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

function check_stability!(bot::Robot,Ň;
        gravity=false,
        scaling=0.01,
        scalings=nothing
    )
    (;st,traj) = bot
    static_equilibrium,λ = check_static_equilibrium_output_multipliers(st)
    @assert static_equilibrium
    q̌ = get_free_coords(st)
    _, Ň0, er = check_stability(bot.structure,λ,Ň;verbose=true)
    resize!(traj,1)
    for i in 1:length(er.values)
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = er.values[i]
        δq̌i = Ň0*er.vectors[:,i]
        # @show δq̌i, er.vectors[:,i]
        if scalings isa Nothing
            si = scaling
        else
            si = scalings[i]
        end
        ratio = norm(δq̌i) / norm(q̌) 
        traj.q̌[end] .= q̌ .+ si.*δq̌i/ratio
    end
    bot
end

function make_nullspace(st::Structure,q0::AbstractVector)
	(;bodies,connectivity) = st
    (;num_of_free_coords,num_of_full_coords,sys_pres_idx,sys_free_idx,bodyid2sys_free_coords,bodyid2sys_intrinsic_cstr_idx,num_of_intrinsic_cstr) = connectivity.indexed
    function inner_nullspace(q̌)
        T = eltype(q̌)
		q = Vector{T}(undef,num_of_full_coords)
		q[sys_pres_idx] .= q0[sys_pres_idx]
		q[sys_free_idx] .= q̌
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

function get_poly(bot_input;
        Ň
    )
    bot = deepcopy(bot_input)
    (;st) = bot
    # (;num_of_dof,num_of_cstr,connectivity) = bot.structure
    # (;cables) = st.tensiles
    # (;num_of_full_coords,num_of_free_coords) = connectivity.indexed
    # ncables = length(cables)
    # nλ = num_of_cstr
    gue = get_initial(st)
    Φ = make_cstr_function(st,gue.q)
    A = make_cstr_jacobian(st,gue.q)
    Q̌ = make_Q̌(st,gue.q)
    S = make_S(st,gue.q)
    Ǩm_Ǩg = make_Ǩm_Ǩg(st,gue.q)

    pv = get_polyvar(st)

    pnq̌ = 1.0pv.q̌ .+ 0.0
    pns = 1.0pv.s .+ 0.0
    pnλ = 1.0pv.λ .+ 0.0
    pnd = 1.0pv.d .+ 0.0
    pnc = 1.0pv.c .+ 0.0
    pnk = 1.0pv.k .+ 0.0
    pnμ = 1.0pv.μ .+ 0.0
    polyΦ = Φ(pnq̌,pnd,pnc)
    polyA = A(pnq̌,pnc)
    polyQ̌ = Q̌(pnq̌,pns,pnμ,pnk,pnc)
    polyS = S(pnq̌,pns,pnc)
    polyQ̌a = transpose(polyA)*pnλ
    polyǨa = reduce(hcat,differentiate.(-polyQ̌a,Ref(pv.q̌))) |> transpose
    polyǨm, polyǨg = Ǩm_Ǩg(pnq̌,pns,pnμ,pnk,pnc)
    polyǨ = polyǨm .+ polyǨg .+ polyǨa
    polyŇ = Ň(pnq̌,pnc)
    poly𝒦 = transpose(polyŇ)*polyǨ*polyŇ

    polyP = [
        - polyQ̌ .- transpose(polyA)*pnλ ;
        polyS;
        polyΦ;
        # poly𝒦*pnξ.-pnζ.*pnξ;
        # transpose(pnξ)*pnξ-1;
    ]

    # Ǩ0 = RB.build_Ǩ(bot.structure,gue.λ)
    # Ǩx = map(polyǨ) do z
    # 		z(
    # 			pv.q̌=>gue.q̌,
    # 			pv.s=>gue.s,
    # 			pv.λ=>gue.λ,
    # 			pv.μ=>gue.μ,
    # 			pv.k=>gue.k,
    # 			pv.d=>gue.d,
    # 			pv.c=>gue.c
    # 		)
    # 	end
    # # @show Ǩ0
    # @show Ǩ0.- Ǩx |> norm

    # P0 = map(polyP) do z
    # 	z(
    # 		pvq̌=>q̌0,
    # 		pvs=>s0,
    # 		pvλ=>λ0,
    # 		# pvξ=>ξ0,
    # 		pvμ=>μ0,
    # 		pvk=>k0,
    # 	    pvd=>d0,
    # 		pvc=>c0,
    # 		# pv.ζ=>ζ0
    # 	)
    # end
    # @show P0[                 1:num_of_free_coords] |> norm
    # @show P0[           num_of_free_coords+1:num_of_free_coords+ncables] |> norm
    # @show P0[   num_of_free_coords+ncables+1:num_of_free_coords+ncables+nλ] |> norm
    # @show P0[num_of_free_coords+ncables+nλ+1:num_of_free_coords+ncables+nλ+num_of_dof]
    # @show P0[end]
    polyP,poly𝒦,gue,pv
end

function pinpoint(bot_input;
        Ň
    )
	polyP, poly𝒦, gue, pv = get_poly(bot_input;Ň)
	ň = length(pv.q̌)
	ns = length(pv.s)
	nλ = length(pv.λ)
    function make_bf()
        function inner_pp!(f,x)
            q̌x = @view x[        1:ň]
            sx = @view x[      ň+1:ň+ns]
            λx = @view x[   ň+ns+1:ň+ns+nλ]
			Px = map(polyP) do z
                z(
                    pv.q̌=>q̌x,
                    pv.s=>sx,
                    pv.λ=>λx,
                    pv.μ=>gue.μ,
                    pv.k=>gue.k,
                    pv.d=>gue.d,
                    pv.c=>gue.c,
                )
			end

			f .= Px
        end
    end
    f_holder = zeros(ň+ns+nλ)
    x_initial = vcat(gue.q̌,gue.s,gue.λ)
    pp! = make_bf()

    pp = nlsolve(pp!,x_initial,ftol=1e-10,iterations=100,method=:newton)
    # @show
    pp!(f_holder,pp.zero)
	# @show f_holder |> norm
    # @show f_holder[                 1:ň+ns+nλ] |> norm
    # @show f_holder[ň+ns+nλ+1:ň+ns+nλ+num_of_dof] |> norm
    # @show f_holder[end]
    # @show  pp.zero[ň+ns+nλ+1:ň+ns+nλ+num_of_dof]
    # @show  pp.zero[end]
    q̌ = pp.zero[        1:ň]
    s = pp.zero[      ň+1:ň+ns]
    λ = pp.zero[   ň+ns+1:ň+ns+nλ]
    ini = @eponymtuple(
			q̌,s,λ,
			isconverged=converged(pp),
			d=gue.d, c=gue.c, μ=gue.μ, k=gue.k
	)
	polyP, poly𝒦, ini, pv
end

function path_follow(bot_input;Ň)
	polyP, poly𝒦, ini, pv = pinpoint(bot_input;Ň)
	variable_groups = [pv.q̌,pv.s,pv.λ]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ]]
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
	Psys = System(polyP;parameters)
	result = HomotopyContinuation.solve(
			Psys,
			startsols;
			start_parameters,
			target_parameters,
			threading = false
	)
	path_results = results(result)
	if length(path_results) != 1
		@show failed(result)
		error("Tracking failed.")
	end
	path_result1 = path_results[1]
	sol = real(solution(path_result1))
	q̌,s,λ = split_by_lengths(sol,length.(variable_groups))
	@eponymtuple(q̌,s,λ)
end

function path_follow_critical(bot_input)
	polyP, ini, pv = pinpoint_critical(bot_input)
	variable_groups = [pv.q̌,pv.s,pv.λ,pv.ξ,[pv.ζ]]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ;ini.ξ;ini.ζ]]
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
	Psys = System(polyP;parameters)
	result = HomotopyContinuation.solve(
			Psys,
			startsols;
			start_parameters,
			target_parameters,
			threading = false
	)
	path_results = results(result)
	if length(path_results) != 1
		@show failed(result)
		error("Tracking failed.")
	end
	path_result1 = path_results[1]
	sol = real(solution(path_result1))
	q̌,s,λ,ξ,ζ = RB.split_by_lengths(sol,length.(variable_groups))
	@eponymtuple(q̌,s,λ,ξ,ζ)
end


