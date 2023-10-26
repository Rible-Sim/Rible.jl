
mutable struct ClusterNonminimalCoordinatesState{T,qT,qviewT}
    t::T
    q::qT
    q̇::qT
    q̈::qT
    F::qT
    λ::qT
    s::qT
    q̌::qviewT
    q̌̇::qviewT
    q̌̈::qviewT
    q̃::qviewT
    q̃̇::qviewT
    q̃̈::qviewT
    F̌::qviewT
end

function ClusterNonminimalCoordinatesState(t,q,q̇,q̈,F,λ,freei,presi,s)
    q̌ = @view q[freei]
    q̌̇ = @view q̇[freei]
    q̌̈ = @view q̈[freei]
    q̃ = @view q[presi]
    q̃̇ = @view q̇[presi]
    q̃̈ = @view q̈[presi]
    F̌ = @view F[freei]
    ClusterNonminimalCoordinatesState(t,q,q̇,q̈,F,λ,s,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌)
end



function StructureState(bodies,tensiles,cnt::Connectivity{<:Any,<:Any,<:NamedTuple{(:connected, :clustered)},<:Any})
    (;clustercables) = tensiles
    (;indexed,jointed) = cnt
    (;nfull,ninconstraints,sysfree,syspres) = indexed
    (;mem2sysincst,mem2sysfull,mem2sysfree,mem2syspres) = indexed
    (;nexconstraints) = jointed
    nclustercables = length(clustercables)
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    nconstraints = ninconstraints + nexconstraints
    nb = length(bodies)
    pres_idx_by_mem = Vector{Vector{Int}}(undef,nb)
    free_idx_by_mem = Vector{Vector{Int}}(undef,nb)
    foreach(bodies) do rb
        rbid = rb.prop.id
        pres_idx_by_mem[rbid] = rb.state.cache.pres_idx
        free_idx_by_mem[rbid] = rb.state.cache.free_idx
    end
    T = get_numbertype(bodies)
    t = zero(T)
    q = Vector{T}(undef,nfull)
    q̇ = zero(q)
    q̈ = zero(q)
    F = zero(q)
    s = zeros(T, 2ns)
    λ = Vector{T}(undef,nconstraints)
    system = ClusterNonminimalCoordinatesState(t,q,q̇,q̈,F,λ,sysfree,syspres,s)
    parts = [
        begin
            qmem = @view q[mem2sysfull[rbid]]
            q̇mem = @view q̇[mem2sysfull[rbid]]
            q̈mem = @view q̈[mem2sysfull[rbid]]
            Fmem = @view F[mem2sysfull[rbid]]
            λmem = @view λ[mem2sysincst[rbid]]
            NonminimalCoordinatesState(t,qmem,q̇mem,q̈mem,Fmem,λmem,
                                    free_idx_by_mem[rbid],pres_idx_by_mem[rbid])
        end
        for rbid = 1:nb
    ]
    foreach(bodies) do rb
        (;ro,R,ṙo,ω,cache) = rb.state
        q,q̇ = NCF.rigidstate2naturalcoords(cache.funcs.lncs,ro,R,ṙo,ω)
        parts[rb.prop.id].q .= q
        parts[rb.prop.id].q̇ .= q̇
    end
    StructureState(system,parts)
end

function update_tensiles!(st, @eponymargs(clustered,))
    (;clustercables) = st.tensiles
    id = 0
    foreach(clustered) do scnt
        id += 1
        clustercable = clustercables[id]
        (;s) = clustercable.sps
        for (segid, seg) in enumerate(clustercable.segs)
            (;state, k, c, original_restlen, prestress) = seg
            (;restlen) = state
            state1 = scnt[segid].hen.rbsig.state
            state2 = scnt[segid].egg.rbsig.state
            pid1 = scnt[segid].hen.pid
            pid2 = scnt[segid].egg.pid
            p1 = state1.rps[pid1]
            ṗ1 = state1.ṙps[pid1]
            f1 = state1.fps[pid1]
            p2 = state2.rps[pid2]
            ṗ2 = state2.ṙps[pid2]
            f2 = state2.fps[pid2]
            Δr = p2 - p1
            Δṙ = ṗ2 - ṗ1
            state.length,state.direction = lengthdir(p2-p1)
            l = state.length
            τ = state.direction
            state.lengthdot = (transpose(Δr)*Δṙ)/l
            if segid == 1
                u = restlen + s[segid]
            elseif segid == length(clustercable.segs)
                u = restlen - s[segid-1]
            else
                u = restlen + s[segid] - s[segid-1]
            end
            state.tension = k*(l-u) + prestress
            state.tension = max(state.tension, 0)
            f1 .+=  τ*state.tension
            f2 .-=  τ*state.tension
            # state.force = τ*state.tension
        end
    end
end

function update_tensiles!(st, @eponymargs(connected, clustered))
    update_tensiles!(st, @eponymtuple(connected))
    update_tensiles!(st, @eponymtuple(clustered))
end


function distribute_s̄!(st::AbstractStructure, s̄)
    s⁺ = @view s̄[begin:2:end]
    s⁻ = @view s̄[begin+1:2:end]
    is = 0
    for cs in st.tensiles.clustercables
        nsi = length(cs.sps)
        cs.sps.s⁺ .= s⁺[is+1:is+nsi]
        cs.sps.s⁻ .= s⁻[is+1:is+nsi]
        cs.sps.s  .= s⁺[is+1:is+nsi] - s⁻[is+1:is+nsi]
        is += nsi
    end
end

function get_s̄(st::AbstractStructure)
    ns = sum([length(st.tensiles.clustercables[i].sps) for i in 1:st.nclustercables])
    s̄ = zeros(get_numbertype(st),2ns)
    s⁺ = @view s̄[begin:2:end]
    s⁻ = @view s̄[begin+1:2:end]
    is = 0
    for cs in st.clusterstrings
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

function make_Ψ(st::AbstractStructure)
    (;clustercables) = st.tensiles
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
        reset_forces!(st)
        distribute_q_to_rbs!(st, q)
        distribute_s̄!(st,s̄)
        update_tensiles!(st)
        _inner_Ψ(s̄)
    end
    inner_Ψ
end

function build_ζ(st::AbstractStructure)
    (;clustercables) = st.tensiles
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

function build_∂ζ∂q(st::AbstractStructure,q̌)
    (;ndim, connectivity) = st
    (;clustercables) = st.tensiles
    nclustercables = length(clustercables)
    (;tensioned, indexed) = connectivity
    # (;q̌) = st.state.system
    (;nfull, nfree, sysfree, mem2sysfree, mem2sysfull) = indexed
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    nclustersegs = ns + nclustercables
    Type = get_numbertype(st)
    ∂l∂q = zeros(Type, nclustersegs, nfree)
    i = 0; j = 0
    foreach(tensioned.clustered) do clustercable
        i += 1
        foreach(clustercable) do cc
            J̌ = zeros(Type,ndim,nfree)
            j += 1
            cable = clustercables[i].segs[cc.id]
            (;hen,egg) = cc
            rb1 = hen.rbsig
            rb2 = egg.rbsig
            C1 = rb1.state.cache.Cps[hen.pid]
            C2 = rb2.state.cache.Cps[egg.pid]
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

function Record_build_∂ζ∂q(st::AbstractStructure,q̌, xlsxname, sheetname)
    (;nclustercables, clustercables, ndim, connectivity) = st
    (;tensioned, indexed) = connectivity
    # (;q̌) = st.state.system
    (;nfull, nfree, sysfree, mem2sysfree, mem2sysfull) = indexed
    ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
    nclustersegs = ns + nclustercables
    Type = get_numbertype(st)
    ∂l∂q = zeros(Type, nclustersegs, nfree)
    i = 0; j = 0
    foreach(tensioned.clustercables) do clustercable
        i += 1
        foreach(clustercable) do cc
            J̌ = zeros(Type,ndim,nfree)
            j += 1
            cable = st.clustercables[i].segs[cc.id]
            (;hen,egg) = cc
            rb1 = hen.rbsig
            rb2 = egg.rbsig
            C1 = rb1.state.cache.Cps[hen.pid]
            C2 = rb2.state.cache.Cps[egg.pid]
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

function get_clusterA(st)
    (;clustercables) = st.tensiles
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

function build_∂ζ∂s̄(st)
    (;clustercables) = st.tensiles
    A_list = Vector{SparseMatrixCSC{Float64,Int64}}()
    clusterA = get_clusterA(st)
    for (cid,clustercable) in enumerate(clustercables)
        nseg = length(clustercable.segs)
        T = get_TransMatrix(nseg-1)
        A = sparse(T*clusterA[cid]*T')
        push!(A_list,A)
    end
    return reduce(blockdiag,A_list)
end
