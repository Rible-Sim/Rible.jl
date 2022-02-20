function distribute_s̄!(tg::ClusterTensegrityStructure,s̄)
    css = tg.clusterstrings
    s⁺ = @view s̄[begin:2:end]
    s⁻ = @view s̄[begin+1:2:end]
    is = 0
    for cs in tg.clusterstrings
        nsi = length(cs.sps)
        cs.sps.s⁺ .= s⁺[is+1:is+nsi]
        cs.sps.s⁻ .= s⁻[is+1:is+nsi]
        cs.sps.s  .= s⁺[is+1:is+nsi] - s⁻[is+1:is+nsi]
        is += nsi
    end
end

function get_s̄(tg::ClusterTensegrityStructure)
    ns = tg.nslidings
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

function build_Ψ(tg::ClusterTensegrityStructure)
    @unpack clusterstrings = tg
    FB = FischerBurmeister(1e-14,1000.,1000.)
    function _inner_Ψ(s̄)
        ψ = Vector{eltype(s̄)}(undef,2tg.nslidings)
        ψ⁺ = @view ψ[begin:2:end]
        ψ⁻ = @view ψ[begin+1:2:end]
        s⁺ = @view s̄[begin:2:end]
        s⁻ = @view s̄[begin+1:2:end]
        is = 0
        for (csid,cs) in enumerate(clusterstrings)
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
        update_strings!(tg)
        _inner_Ψ(s̄)
    end
    inner_Ψ
end

function build_ζ(tg::ClusterTensegrityStructure)
    @unpack clusterstrings = tg
    ζ = Vector{Float64}(undef,2tg.nslidings)
    ζ⁺ = @view ζ[begin:2:end]
    ζ⁻ = @view ζ[begin+1:2:end]
    is = 0
    for cs in clusterstrings
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
function get_clusterA(tg)
    @unpack clusterstrings  = tg
    A_list = [Matrix{Float64}(undef,1,1) for i in 1:length(clusterstrings)]
    for (csid, cs) in enumerate(clusterstrings)
        @unpack k = cs.segs
        @unpack α = cs.sps
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

function get_clusterB(tg)
    @unpack clusterstrings = tg
    B_list = [Vector{Float64}() for i in 1:length(clusterstrings)]
    for (csid, cs) in enumerate(clusterstrings)
        @unpack k,state = cs.segs
        @unpack α,s = cs.sps
        n = length(α)
        restlen = [state[i].restlen for i in 1:n+1]
        for rid in 1:length(restlen)
            if rid == 1
                restlen[rid] += s[rid]
            elseif rid == length(restlen)
                restlen[rid] += -s[rid-1]
            else
                restlen[rid] += s[rid] - s[rid-1]
            end
        end
        len = [state[i].length for i in 1:n+1]
        B⁺ = [-k[i]*(len[i]-restlen[i])+k[i+1]*(len[i+1]-restlen[i+1])/α[i] for i in 1:n]
        B⁻ = [k[i]*(len[i]-restlen[i])-α[i]*k[i+1]*(len[i+1]-restlen[i+1]) for i in 1:n]
        B = [B⁺;B⁻]
        B_list[csid] = B
    end
    return B_list
end

function get_clusterZ(tg)
    @unpack clusterstrings = tg
    Z_list = [Vector{Float64}() for i in 1:length(clusterstrings)]
    for (csid, cs) in enumerate(clusterstrings)
        @unpack sps = cs
        Z_list[csid] = vcat(sps.s⁺,sps.s⁻)
    end
    return Z_list
end

function get_clusterζ(tg::ClusterTensegrityStructure)
    @unpack clusterstrings = tg
    ζ_list = [Vector{Float64}() for i in 1:length(clusterstrings)]
    for (csid,cs) in enumerate(clusterstrings)
        @unpack sps,segs = cs
        n = length(sps)
        ζ⁺ = [segs[i+1].state.tension/sps[i].α - segs[i].state.tension for i in 1:n]
        ζ⁻ = [segs[i].state.tension - sps[i].α*segs[i+1].state.tension for i in 1:n]
        ζ = vcat(ζ⁺,ζ⁻)
        ζ_list[csid] = ζ
    end
    return ζ_list
end
