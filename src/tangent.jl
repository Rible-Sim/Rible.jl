struct TensionTangent{T}
    J::Vector{Array{T,2}}
    U::Vector{Array{T,2}}
    l::Vector{T}
    l̇::Vector{T}
    ∂l∂q::Vector{Array{T,2}}
    ∂l̇∂q::Vector{Array{T,2}}
    ∂l̇∂q̇::Vector{Array{T,2}}
    ∂f∂q::Vector{Array{T,2}}
    ∂f∂q̇::Vector{Array{T,2}}
    ∂l̂∂q::Vector{Array{T,2}}
end

function build_tangent(tg)
    q, q̇ = get_q(tg)
    J = [build_Ji(tg,i) for i = 1:tg.nstrings]
    U = [transpose(Ji)*Ji for Ji in J]
    l = [sqrt(transpose(q)*Ui*q) for Ui in U]
    # display(l)
    # display([s.state.length for s in tg.strings])
    # @show l.-[s.state.length for s in tg.strings]
    l̇ = [(transpose(q)*Ui*q̇)/li for (Ui,li) in zip(U,l)]
    ∂l∂q = [(transpose(q)*Ui)./li for (Ui,li) in zip(U,l)]
    ∂l̇∂q̇ = ∂l∂q
    ∂l̇∂q = [(li*transpose(q̇)-l̇i*transpose(q))/li^2*Ui for (Ui,li,l̇i) in zip(U,l,l̇)]
    k = [s.k for s in tg.strings]
    c = [s.c for s in tg.strings]
    ∂f∂q = [ki*∂li∂q+ci*∂l̇i∂q for (ki,ci,∂li∂q,∂l̇i∂q) in zip(k,c,∂l∂q,∂l̇∂q)]
    ∂f∂q̇ = [ci*∂l̇i∂q̇ for (ci,∂l̇i∂q̇) in zip(c,∂l̇∂q̇)]
    l̂ = [Ji*q/li for (Ji,li) in zip(J,l)]
    ∂l̂∂q = [(Ji-l̂i*∂li∂q)/li for (li,l̂i,Ji,∂li∂q) in zip(l,l̂,J,∂l∂q)]
    # @show l̂.-[s.state.direction for s in tg.strings]
    f = [s.state.tension for s in tg.strings]
    ∂𝐟∂q = [l̂i*∂fi∂q+fi*∂l̂i∂q for (l̂i,fi,∂fi∂q,∂l̂i∂q) in zip(l̂,f,∂f∂q,∂l̂∂q)]
    ∂𝐟∂q̇ = [l̂i*∂fi∂q̇ for (l̂i,∂fi∂q̇) in zip(l̂,∂f∂q̇)]
    @unpack ndim,ncoords,nstrings = tg
    ∂Γ∂q = zeros(eltype(q),ndim*nstrings,ncoords)
    ∂Γ∂q̇ = zeros(eltype(q),ndim*nstrings,ncoords)
    for i in 1:nstrings
        ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] = ∂𝐟∂q[i]
        ∂Γ∂q̇[(i-1)*ndim+1:i*ndim,:] = ∂𝐟∂q̇[i]
    end
    ∂Γ∂q,∂Γ∂q̇
end

function build_Jac_Γ(tg::TensegrityStructure)
    ns = tg.nstrings
    @unpack ncoords,ndim = tg
    J = [build_Ji(tg,i) for i = 1:ns]
    U = [transpose(Ji)*Ji for Ji in J]
    k = [s.k for s in tg.strings]
    c = [s.c for s in tg.strings]
    function inner_Jac_Γ(q,q̇)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg,q,q̇)
        update_strings!(tg)
        f = [s.state.tension for s in tg.strings]
        l = [s.state.length for s in tg.strings]
        u = [s.state.restlen for s in tg.strings]
        qᵀU = [transpose(q)*U[i] for i = 1:ns]
        # l = [sqrt(qᵀU[i]*q) for i = 1:ns]
        l̇ = [(qᵀU[i]*q̇)/l[i] for i = 1:ns]
        l̂ = [J[i]*q/l[i] for i = 1:ns]
        ∂Γ∂q = zeros(eltype(q),ndim*ns,ncoords)
        ∂Γ∂q̇ = zeros(eltype(q),ndim*ns,ncoords)
        for i in 1:ns
            l̂qᵀUi = l̂[i]*qᵀU[i]
            # ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] .= l̂[i]*(k[i]./l[i].*qᵀU[i] .+ c[i]./l[i].*transpose(q̇).-c[i]*l̇[i]./l[i]^2 .*qᵀU[i])
            # ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] .+= f[i].*(J[i]./l[i].-l̂[i]*qᵀU[i]./l[i]^2)
            ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] .= (k[i]*u[i]-c[i]*l̇[i])/l[i]^2 .*l̂qᵀUi .+ c[i]/l[i].*(l̂[i]*transpose(q̇)*U[i]) .+ f[i]/l[i] .*J[i]
            ∂Γ∂q̇[(i-1)*ndim+1:i*ndim,:] .= c[i]/l[i].*l̂qᵀUi
        end
        ∂Γ∂q,∂Γ∂q̇
    end
end

function build_Jac_Γ(tg::ClusterTensegrityStructure)
    @unpack nstrings, nclusterstrings, nslidings = tg
    @unpack ncoords,ndim,clusterstrings = tg
    nclustersegs = nclusterstrings + nslidings
    J = [build_Ji(tg,i) for i = 1:nstrings]
    k = [s.k for s in tg.strings]
    c = [s.c for s in tg.strings]
    Jc = Vector{SparseMatrixCSC{Float64,Int64}}()
    kc = Vector{Float64}()
    cc = Vector{Float64}()
    for (cid,clusterstring) in enumerate(clusterstrings)
        @unpack segs = clusterstring
        for (sid, seg) in enumerate(segs)
            push!(Jc, build_Ji(tg, cid, sid))
            push!(kc, seg.k)
            push!(cc, seg.c)
        end
    end
    U = [transpose(Ji)*Ji for Ji in J]
    Uc = [transpose(Jci)*Jci for Jci in Jc]

    function _inner_Jac_Γ(q,q̇,s̄)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg,q,q̇)
        distribute_s̄!(tg,s̄)
        update_strings!(tg)
        #@unpack clusterstrings = tg
        f = [s.state.tension for s in tg.strings]
        l = [s.state.length for s in tg.strings]
        u = [s.state.restlen for s in tg.strings]

        uc = Vector{Float64}()
        fc = Vector{Float64}()
        lc = Vector{Float64}()
        for clusterstring in clusterstrings
            @unpack s = clusterstring.sps
            @unpack state = clusterstring.segs
            n = length(clusterstring.segs)
            res = [state[i].restlen for i in 1:n]
            res[1:n-1] += s
            res[2:n] -= s
            append!(uc,res)
            @unpack segs = clusterstring
            for seg in segs
                push!(fc, seg.state.tension)
                push!(lc, seg.state.length)
            end
        end
        qᵀU = [transpose(q)*U[i] for i = 1:nstrings]
        l̇ = [(qᵀU[i]*q̇)/l[i] for i = 1:nstrings]
        l̂ = [J[i]*q/l[i] for i = 1:nstrings]

        qᵀUc = [transpose(q)*Uc[i] for i = 1:nclustersegs]
        l̇c = [(qᵀUc[i]*q̇)/lc[i] for i = 1:nclustersegs]
        l̂c = [Jc[i]*q/lc[i] for i = 1:nclustersegs]

        ∂Γ∂qs = zeros(eltype(q),ndim*nstrings,ncoords)
        ∂Γ∂q̇s = zeros(eltype(q),ndim*nstrings,ncoords)
        ∂Γ∂qc = zeros(eltype(q),ndim*nclustersegs,ncoords)
        ∂Γ∂q̇c = zeros(eltype(q),ndim*nclustersegs,ncoords)
        for i in 1:nstrings
            l̂qᵀUi = l̂[i]*qᵀU[i]
            ∂Γ∂qs[(i-1)*ndim+1:i*ndim,:] .= (k[i]*u[i]-c[i]*l̇[i])/l[i]^2 .*l̂qᵀUi .+ c[i]/l[i].*(l̂[i]*transpose(q̇)*U[i]) .+ f[i]/l[i] .*J[i]
            ∂Γ∂q̇s[(i-1)*ndim+1:i*ndim,:] .= c[i]/l[i].*l̂qᵀUi
        end
        for i in 1:nclustersegs
            l̂qᵀUic = l̂c[i]*qᵀUc[i]
            ∂Γ∂qc[(i-1)*ndim+1:i*ndim,:] .+= (kc[i]*uc[i]-cc[i]*l̇c[i])/lc[i]^2 .*l̂qᵀUic .+ cc[i]/lc[i].*(l̂c[i]*transpose(q̇)*Uc[i]) .+ fc[i]/lc[i] .*Jc[i]
            ∂Γ∂q̇c[(i-1)*ndim+1:i*ndim,:] .+= cc[i]/lc[i].*l̂qᵀUic
        end
        vcat(∂Γ∂qs,∂Γ∂qc), vcat(∂Γ∂q̇s,∂Γ∂q̇c)
    end
    inner_Jac_Γ(q,q̇,s̄) = _inner_Jac_Γ(q,q̇,s̄)
    function inner_Jac_Γ(q,q̇)
        s̄ = get_s̄(tg)
        _inner_Jac_Γ(q,q̇,s̄)
    end
    inner_Jac_Γ
end

function build_∂Γ∂s̄(tg::ClusterTensegrityStructure)
    @unpack nstrings, nslidings, nclusterstrings = tg
    @unpack clusterstrings,ncoords,ndim = tg
    nclustersegs = nclusterstrings + nslidings
    Jc = Vector{SparseMatrixCSC{Float64,Int64}}()
    kc = Vector{Float64}()
    N_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    for (cid,clusterstring) in enumerate(clusterstrings)
        @unpack segs = clusterstring
        nsegs = length(segs)
        for (sid, seg) in enumerate(segs)
            push!(Jc, build_Ji(tg, cid, sid))
            push!(kc, seg.k)
        end
        N = sparse(zeros(Float64, nsegs, 2nsegs-2))
        N[1,1:2]=[1 -1]; N[end,end-1:end]=[-1 1]
        for i in 2:nsegs-1
            N[i, 2i-3:2i] = [-1 1 1 -1]
        end
        push!(N_list, N)
    end
    #N = repeat(reduce(blockdiag,N_list),inner=(ndim,1))
    N = reduce(blockdiag,N_list)
    function _inner_∂Γ∂s̄_(q)
        lc = Vector{Float64}()
        for clusterstring in clusterstrings
            @unpack segs = clusterstring
            for seg in segs
                push!(lc, seg.state.length)
            end
        end

        l̂c = [Jc[i]*q/lc[i] for i = 1:nslidings+nclusterstrings]
        ∂Γ∂s̄s = zeros(eltype(q),ndim*nstrings,2nslidings)
        ∂Γ∂s̄c = zeros(eltype(q),ndim*(nclustersegs),2nslidings)
        for i in 1:nclustersegs
            #∂Γ∂s̄c[(i-1)*ndim+1,:] .= kc[i].*l̂c[i][1].*N[ndim*(i-1)+1,:]
            #∂Γ∂s̄c[i*ndim,:] .= kc[i].*l̂c[i][2].*N[ndim*i,:]
            #@show kc[i],l̂c[i],N[(i-1)*ndim+1:i*ndim,:]
            for j in 1:ndim
                ∂Γ∂s̄c[(i-1)*ndim+j,:] .= -kc[i]*l̂c[i][j]*N[i,:]
            end
            #@show size(∂Γ∂s̄c[(i-1)*ndim+1:i*ndim,:])
        end
        return vcat(∂Γ∂s̄s, ∂Γ∂s̄c)
    end

    function inner_∂Γ∂s̄(q,q̇,s̄)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg,q,q̇)
        distribute_s̄!(tg, s̄)
        update_strings!(tg)
        _inner_∂Γ∂s̄_(q)
    end

    function inner_∂Γ∂s̄(q)
        _inner_∂Γ∂s̄_(q)
    end
    return inner_∂Γ∂s̄
end

function build_∂ζ∂q(tg)
    @unpack nslidings, nclusterstrings = tg
    @unpack clusterstrings,ncoords,ndim = tg
    nclustersegs = nslidings+nclusterstrings

    Jc = Vector{SparseMatrixCSC{Float64,Int64}}()
    kc = Vector{Float64}()
    αc = Vector{Float64}()
    b_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    T_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    for (cid,clusterstring) in enumerate(clusterstrings)
        @unpack segs,sps = clusterstring
        nsegs = length(segs)
        for (sid, seg) in enumerate(segs)
            push!(Jc, build_Ji(tg, cid, sid))
            push!(kc, seg.k)
        end
        for sp in sps
            append!(αc,sp.α)
        end
        b⁺ = sparse(zeros(Float64, nsegs-1, nsegs))
        b⁻ = sparse(zeros(Float64, nsegs-1, nsegs))
        for i in 1:nsegs-1
            b⁺[i,i:i+1] = [-kc[i] kc[i+1]/αc[i]]
            b⁻[i,i:i+1] = [kc[i] -αc[i]*kc[i+1]]
        end
        push!(b_list, vcat(b⁺,b⁻))
        push!(T_list, get_TransMatrix(nsegs-1))
    end
    Uc = [transpose(Jci)*Jci for Jci in Jc]
    b = reduce(blockdiag, b_list)
    T = reduce(blockdiag, T_list)

    function inner_∂ζ∂q(q)
        lc = Vector{Float64}()
        for clusterstring in clusterstrings
            for seg in clusterstring.segs
                push!(lc, seg.state.length)
            end
        end
        qᵀUc = [transpose(q)*Uc[i] for i = 1:nclustersegs]
        dldq = [(qᵀUc[i])/lc[i] for i = 1:nclustersegs]
        #∂ζ∂q = 1/2*T*b*reduce(vcat,dldq)
        ∂ζ∂q = T*b*reduce(vcat,dldq)
        return ∂ζ∂q
    end
end

function get_TransMatrix(n)
    T = zeros(Float64, 2n, 2n)
    for (j,i) in enumerate(1:2:2n)
        T[i,j] = 1
    end
    for (j,i) in enumerate(2:2:2n)
        T[i,j+n] = 1
    end
    return sparse(T)
end

function build_∂ζ∂s̄(tg)
    @unpack clusterstrings = tg
    A_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    clusterA = get_clusterA(tg)
    for (cid,clusterstring) in enumerate(clusterstrings)
        nseg = length(clusterstring.segs)
        T = get_TransMatrix(nseg-1)
        A = sparse(T*clusterA[cid]*T')
        push!(A_list,A)
    end
    return reduce(blockdiag,A_list)
end

function make_testtangent(tgstruct)
    @unpack ncoords,nstrings = tgstruct
    function build_𝐟(x)
        q = x[1:ncoords]
        q̇ = x[ncoords+1:2ncoords]
        reset_forces!(tgstruct)
        distribute_q_to_rbs!(tgstruct,q,q̇)
        update_strings_apply_forces!(tgstruct)
        vcat([tgstruct.strings[i].state.direction*tgstruct.strings[i].state.tension
            for i = 1:nstrings]...)
    end
end
