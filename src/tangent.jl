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

function build_Jac_Γ(tg::TensegrityRobots.TensegrityStructure)
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

function build_Jac_Γ(tg::TensegrityRobots.ClusterTensegrityStructure,strings_index::Vector{Int64},cluster_index::Vector{Int64})
    ns = tg.nstrings
    ncs = tg.nclusterstrings
    nsss = Int64((tg.nclusterstrings + tg.nslidings)/ncs)
    nsnc  = maximum(vcat(strings_index,cluster_index))
    @unpack ncoords,ndim = tg
    J = [build_Ji(tg,i) for i = 1:ns]
    k = [s.k for s in tg.strings]
    c = [s.c for s in tg.strings]

    Jc = [build_Ji(tg,i,j) for i in 1:ncs for j in 1:nsss]
    kc = [tg.clusterstrings[i].segs[j].k for i in 1:ncs for j in 1:nsss]
    cc = [tg.clusterstrings[i].segs[j].c for i in 1:ncs for j in 1:nsss]

    U = [transpose(Ji)*Ji for Ji in J]
    Uc = [transpose(Jci)*Jci for Jci in Jc]
    
    function inner_Jac_Γ(q,q̇,s̄)
        #reset_forces!(tg)
        #distribute_q_to_rbs!(tg,q,q̇)
        #distribute_s̄!(tg,s̄)
        #update_strings!(tg)
        @unpack clusterstrings = tg
        f = [s.state.tension for s in tg.strings]
        l = [s.state.length for s in tg.strings]
        u = [s.state.restlen for s in tg.strings]

        fc = [tg.clusterstrings[i].segs[j].state.tension for i in 1:ncs for j in 1:nsss]
        lc = [tg.clusterstrings[i].segs[j].state.length for i in 1:ncs for j in 1:nsss]
        uc = Vector{Float64}()
        for clusterstring in clusterstrings
            @unpack s = clusterstring.sps
            @unpack state = clusterstring.segs
            n = length(clusterstring.segs)
            res = [state[i].restlen for i in 1:n]
            res[1:n-1] += s
            res[2:n] -= s
            append!(uc,res)
        end
        qᵀU = [transpose(q)*U[i] for i = 1:ns]
        l̇ = [(qᵀU[i]*q̇)/l[i] for i = 1:ns]
        l̂ = [J[i]*q/l[i] for i = 1:ns]

        qᵀUc = [transpose(q)*Uc[i] for i = 1:nsss*ncs]
        l̇c = [(qᵀUc[i]*q̇)/lc[i] for i = 1:nsss*ncs]
        l̂c = [Jc[i]*q/lc[i] for i = 1:nsss*ncs]

        ∂Γ∂q = zeros(eltype(q),ndim*nsnc,ncoords)
        ∂Γ∂q̇ = zeros(eltype(q),ndim*nsnc,ncoords)
        for (i,j) in enumerate(strings_index)
            l̂qᵀUi = l̂[i]*qᵀU[i]
            ∂Γ∂q[(j-1)*ndim+1:j*ndim,:] .= (k[i]*u[i]-c[i]*l̇[i])/l[i]^2 .*l̂qᵀUi .+ c[i]/l[i].*(l̂[i]*transpose(q̇)*U[i]) .+ f[i]/l[i] .*J[i]
            ∂Γ∂q̇[(j-1)*ndim+1:j*ndim,:] .= c[i]/l[i].*l̂qᵀUi
        end
        for (i,j) in enumerate(cluster_index)
            l̂qᵀUic = l̂c[i]*qᵀUc[i]
            ∂Γ∂q[(j-1)*ndim+1:j*ndim,:] .+= (kc[i]*uc[i]-cc[i]*l̇c[i])/lc[i]^2 .*l̂qᵀUic .+ cc[i]/lc[i].*(l̂c[i]*transpose(q̇)*Uc[i]) .+ fc[i]/lc[i] .*Jc[i]
            ∂Γ∂q̇[(j-1)*ndim+1:j*ndim,:] .+= cc[i]/lc[i].*l̂qᵀUic
        end
        ∂Γ∂q,∂Γ∂q̇
    end
end

function build_∂Γ∂s̄(tg::ClusterTensegrityStructure,strings_index::Vector{Int64},cluster_index::Vector{Int64})
    nsnc = maximum(vcat(strings_index, cluster_index))
    nsnc = tg.nstrings
    nss = tg.nslidings
    ncs = tg.nclusterstrings
    nsss = Int64((tg.nclusterstrings + tg.nslidings)/ncs)
    @unpack ncoords,ndim = tg
    Jc = [build_Ji(tg,i,j) for i in 1:ncs for j in 1:nsss]
    kc = [1/tg.clusterstrings[i].segs[j].k for i in 1:ncs for j in 1:nsss]
    nk = Int64(length(kc)/2)
    N = zeros(Float64, nk, nss)
    N0 = zeros(Float64, nk, nss)
    N[1,1] = 1; N[1,2] = -1; N[nk,nss] = 1; N[nk,nss-1] = -1
    for i in 2:nk-1
        N[i,2i-3:2i] = [-1 1 1 -1]
    end
    N = repeat([N N0;N0 N],inner=(ndim,1))

    function inner_∂Γ∂s̄(q, q̇, s̄)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg,q,q̇)
        distribute_s̄!(tg, s̄)
        update_strings!(tg)
        lc = [tg.clusterstrings[i].segs[j].state.length for i in 1:ncs for j in 1:nsss]
        
        l̂c = [Jc[i]*q/lc[i] for i = 1:nsss*ncs]

        ∂Γ∂s̄ = zeros(eltype(q),ndim*nsnc,nss*2)
        
        for (i,j) in enumerate(cluster_index)
            #∂Γ∂s̄[(j-1)*ndim+1:j*ndim,:] .= -kc[i].*l̂c[i]'*N[ndim*(i-1)+1:ndim*i,:]
            ∂Γ∂s̄[(j-1)*ndim+1,:] .= 1/2*kc[i].*l̂c[i][1].*N[ndim*(i-1)+1,:]
            ∂Γ∂s̄[j*ndim,:] .= 1/2*kc[i].*l̂c[i][2].*N[ndim*i,:]
            #∂Γ∂s̄[(j-1)*ndim+1,:] .= KN[i,:].*l̂c[i][1]
            #∂Γ∂s̄[j*ndim,:] .= KN[i,:].*l̂c[i][2]
        end
        return ∂Γ∂s̄
    end
    
end

function build_∂ζ∂q(tg)    
    ncs = tg.nclusterstrings
    nsss = Int64((tg.nclusterstrings + tg.nslidings)/ncs)
    @unpack ncoords,ndim = tg

    Jc = [build_Ji(tg,i,j) for i in 1:ncs for j in 1:nsss]
    kc = [tg.clusterstrings[i].segs[j].k for i in 1:ncs for j in 1:nsss]
    αc = [tg.clusterstrings[i].sps[j].α for i in 1:ncs for j in 1:nsss-1]
    Uc = [transpose(Jci)*Jci for Jci in Jc]
    nk = Int64(length(kc)/2)

    b0 = zeros(Float64, nk-1,nk)
    b1⁺ = zeros(Float64, nk-1, nk)
    b1⁻ = deepcopy(b1⁺)
    b2⁺ = deepcopy(b1⁺)
    b2⁻ = deepcopy(b1⁺)

    for i in 1:nk-1
        b1⁺[i,i:i+1] = [-kc[i] kc[i+1]/αc[i]]
        b1⁻[i,i:i+1] = [kc[i] -αc[i]*kc[i+1]]
        b2⁺[i,i:i+1] = [-kc[i+nk] kc[i+1+nk]/αc[i+nk-1]]
        b2⁻[i,i:i+1] = [kc[i+nk] -αc[i+nk-1]*kc[i+1+nk]]
    end
    b = [b1⁺ b0;b1⁻ b0;b0 b2⁺;b0 b2⁻]

    T0 = zeros(Float64, 2*(nsss-1), 2*(nsss-1))
    T = get_TransMatrix(2*(nsss-1))    
    T = [T T0;T0 T]

    function inner_∂W∂q(q)
        lc = [tg.clusterstrings[i].segs[j].state.length for i in 1:ncs for j in 1:nsss]
        qᵀUc = [transpose(q)*Uc[i] for i = 1:nsss*ncs]
        dldq = [(qᵀUc[i])/lc[i] for i = 1:nsss*ncs] 
        a = zeros(Float64,length(dldq),length(dldq[1]))
        for i in 1:length(dldq)
            a[i,:] = dldq[i]
        end
        ∂W∂q = 1/2*T*b*a
        return ∂W∂q
    end
end

function get_TransMatrix(n)
    T = zeros(Float64, n, n)
    for (j,i) in enumerate(1:2:n)
        T[i,j] = 1
    end
    for (j,i) in enumerate(2:2:n)
        T[i,j+Int64(n/2)] = 1
    end
    return T
end

function build_∂ζ∂s̄(tg)
    ns = tg.nslidings
    A = get_clusterA(tg)
    A0 = zeros(Float64, ns, ns)
    T = get_TransMatrix(ns)
    AA = [T*A[1]*T' A0; A0 T*A[2]*T']
    return AA
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
