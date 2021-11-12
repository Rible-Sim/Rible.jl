function build_Ψ(tg)
    @unpack clusterstrings = tg.tensiles
    function inner_Ψ(s̄)
        s⁺ = @view s̄[begin:2:end]
        s⁻ = @view s̄[begin+1:2:end]
        
    end
    function inner_Ψ(q,s̄)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg, q)
        update_strings_forces!(tg)
        update_clusterstrings_apply_forces!(tg, s)
        inner_Ψ(s̄)
    end
end
