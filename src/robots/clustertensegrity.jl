function distribute_s̄!(tg::AbstractTensegrityStructure, s̄)
    s⁺ = @view s̄[begin:2:end]
    s⁻ = @view s̄[begin+1:2:end]
    is = 0
    for cs in tg.tensiles.clustercables
        nsi = length(cs.sps)
        cs.sps.s⁺ .= s⁺[is+1:is+nsi]
        cs.sps.s⁻ .= s⁻[is+1:is+nsi]
        cs.sps.s  .= s⁺[is+1:is+nsi] - s⁻[is+1:is+nsi]
        is += nsi
    end
end

function get_s̄(tg::AbstractTensegrityStructure)
    ns = sum([length(tg.tensiles.clustercables[i].sps) for i in 1:tg.nclustercables])
    s̄ = zeros(get_numbertype(tg),2ns)
    s⁺ = @view s̄[begin:2:end]
    s⁻ = @view s̄[begin+1:2:end]
    is = 0
    for cs in tg.clusterstrings
        nsi = length(cs.sps)
        s⁺[is+1:is+nsi] .= cs.sps.s⁺
        s⁻[is+1:is+nsi] .= cs.sps.s⁻
        is += nsi
    end
    s̄
end

struct FischerBurmeister{T}
    ϵ::T
    X₁::T
    X₂::T
end
FischerBurmeister() = FischerBurmeister(1e-14,1.0,1.0)
function (f::FischerBurmeister)(x,y)
    √((x/f.X₁)^2+(y*f.X₂)^2+2f.ϵ^2) - (x/f.X₁+y*f.X₂)
end

function make_Ψ(tg::AbstractTensegrityStructure)
    (;clustercables) = tg.tensiles
    nclustercables = length(clustercables)
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    FB = FischerBurmeister(1e-14,10.,10.)
    function _inner_Ψ(s̄)
        ψ = Vector{eltype(s̄)}(undef,2ns)
        ψ⁺ = @view ψ[begin:2:end]
        ψ⁻ = @view ψ[begin+1:2:end]
        s⁺ = @view s̄[begin:2:end]
        s⁻ = @view s̄[begin+1:2:end]
        is = 0
        for (csid,cs) in enumerate(clustercables)
            @unpack sps,segs = cs
            nsi = length(cs.sps)
            for (i,sp) in enumerate(sps)
                ζ⁺ = segs[i+1].state.tension/sp.α - segs[i].state.tension
                ζ⁻ = segs[i].state.tension - sp.α*segs[i+1].state.tension
                #@show ζ⁺,s⁺[is+i]
                ψ⁺[is+i] = FB(ζ⁺,s⁺[is+i])
                ψ⁻[is+i] = FB(ζ⁻,s⁻[is+i])
            end
            is += nsi
        end
        ψ
    end

    function inner_Ψ(s̄)
        _inner_Ψ(s̄)
    end

    function inner_Ψ(q,s̄)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg, q)
        distribute_s̄!(tg,s̄)
        update_tensiles!(tg)
        _inner_Ψ(s̄)
    end
    inner_Ψ
end

function build_ζ(tg::AbstractTensegrityStructure)
    (;clustercables) = tg.tensiles
    nclustercables = length(clustercables)
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    ζ = Vector{Float64}(undef,2ns)
    ζ⁺ = @view ζ[begin:2:end]
    ζ⁻ = @view ζ[begin+1:2:end]
    is = 0
    for cs in clustercables
        @unpack sps,segs = cs
        nsi = length(cs.sps)
        for (i,sp) in enumerate(sps)
            ζ⁺[is+i] = segs[i+1].state.tension/sp.α - segs[i].state.tension
            ζ⁻[is+i] = segs[i].state.tension - sp.α*segs[i+1].state.tension
        end
        is += nsi
    end
    ζ
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

function build_∂ζ∂q(tg::AbstractTensegrityStructure,q̌)
    (;ndim, connectivity) = tg
    (;clustercables) = tg.tensiles
    nclustercables = length(clustercables)
    (;tensioned, indexed) = connectivity
    # (;q̌) = tg.state.system
    (;nfull, nfree, sysfree, mem2sysfree, mem2sysfull) = indexed
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    nclustersegs = ns + nclustercables
    Type = get_numbertype(tg)
    ∂l∂q = zeros(Type, nclustersegs, nfree)
    i = 0; j = 0
    foreach(tensioned.clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            J̌ = zeros(Type,ndim,nfree)
            j += 1
            cable = clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            uci1 = rb1.state.cache.free_idx
            uci2 = rb2.state.cache.free_idx
            mfree1 = mem2sysfree[rb1.prop.id]
            mfree2 = mem2sysfree[rb2.prop.id]
            (;length, direction) = cable.state
            J̌[:,mfree2] .+= C2[:,uci2]
            J̌[:,mfree1] .-= C1[:,uci1]
            # ∂l∂q[j,:] = q̌'*transpose(J̌)*J̌/length
            ∂l∂q[j,:] = direction' * J̌
            # @show J̌*q̌/length, direction
        end
    end
    kc = Vector{Float64}()
    αc = Vector{Float64}()
    b_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    T_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    for (cid,clustercable) in enumerate(clustercables)
        @unpack segs,sps = clustercable
        nsegs = length(segs)
        for (sid, seg) in enumerate(segs)
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
    b = reduce(blockdiag, b_list)
    T = reduce(blockdiag, T_list)
    return T*b*∂l∂q
end

function Record_build_∂ζ∂q(tg::AbstractTensegrityStructure,q̌, xlsxname, sheetname)
    (;nclustercables, clustercables, ndim, connectivity) = tg
    (;tensioned, indexed) = connectivity
    # (;q̌) = tg.state.system
    (;nfull, nfree, sysfree, mem2sysfree, mem2sysfull) = indexed
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    nclustersegs = ns + nclustercables
    Type = get_numbertype(tg)
    ∂l∂q = zeros(Type, nclustersegs, nfree)
    i = 0; j = 0
    foreach(tensioned.clustercables) do clustercable
        i += 1
        foreach(clustercable) do cc
            J̌ = zeros(Type,ndim,nfree)
            j += 1
            cable = tg.clustercables[i].segs[cc.id]
            (;end1,end2) = cc
            rb1 = end1.rbsig
            rb2 = end2.rbsig
            C1 = rb1.state.cache.Cps[end1.pid]
            C2 = rb2.state.cache.Cps[end2.pid]
            uci1 = rb1.state.cache.free_idx
            uci2 = rb2.state.cache.free_idx
            mfree1 = mem2sysfree[rb1.prop.id]
            mfree2 = mem2sysfree[rb2.prop.id]
            (;length, direction) = cable.state
            J̌[:,mfree2] .+= C2[:,uci2]
            J̌[:,mfree1] .-= C1[:,uci1]
            # ∂l∂q[j,:] = q̌'*transpose(J̌)*J̌/length
            if (i==1) & (j==2)
                display(C1)
                display(C2)
                display(J̌)
                display(C1[:, uci1])
                @show uci1, uci2, mfree1, mfree2
            end
            ∂l∂q[j,:] = direction' * J̌
        end
    end
    kc = Vector{Float64}()
    αc = Vector{Float64}()
    b_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    T_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    for (cid,clustercable) in enumerate(clustercables)
        @unpack segs,sps = clustercable
        nsegs = length(segs)
        for (sid, seg) in enumerate(segs)
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
    b = reduce(blockdiag, b_list)
    T = reduce(blockdiag, T_list)
    h1,w1 = size(b); h2, w2 = size(T); h3, w3 = size(∂l∂q)
    XLSX.openxlsx(xlsxname, mode="rw") do xf
        sh = xf[sheetname]
        for i in 1:h1
            sh[i, 1:w1] = b[i,:]
        end
        for i in h1+1:h1+h2
            sh[i, 1:w2] = T[i-h1, :]
        end
        for i in h1+h2+1:h1+h2+h3
            sh[i, 1:w3] = ∂l∂q[i-h1-h2,:]
        end
    end
    return T*b*∂l∂q
end

function get_clusterA(tg)
    (;clustercables) = tg.tensiles
    A_list = [Matrix{Float64}(undef,1,1) for i in 1:length(clustercables)]
    for (csid, cs) in enumerate(clustercables)
        (;k) = cs.segs
        (;α) = cs.sps
        n = length(α)
        ap = [-k[i+1] for i in 1:n-1]
        ap0 = [k[i]+k[i+1]/α[i] for i in 1:n]
        ap1 = [-k[i+1]/α[i] for i in 1:n-1]
        an0 = [k[i]+α[i]*k[i+1] for i in 1:n]
        an1 = [-α[i]*k[i+1] for i in 1:n-1]
        A⁺ = diagm(-1=>ap, 0=>ap0, 1=>ap1)
        A⁻ = diagm(-1=>ap, 0=>an0, 1=>an1)
        A = [A⁺ -A⁺;-A⁻ A⁻]
        A_list[csid] = A
    end
    return A_list
end

function build_∂ζ∂s̄(tg)
    (;clustercables) = tg.tensiles
    A_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    clusterA = get_clusterA(tg)
    for (cid,clustercable) in enumerate(clustercables)
        nseg = length(clustercable.segs)
        T = get_TransMatrix(nseg-1)
        A = sparse(T*clusterA[cid]*T')
        push!(A_list,A)
    end
    return reduce(blockdiag,A_list)
end
