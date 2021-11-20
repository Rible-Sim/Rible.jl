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
    FB = FischerBurmeister()
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
