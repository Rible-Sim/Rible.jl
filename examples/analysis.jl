function analyse_energy(tgstruct,sol;gravity=false,elasticity=false,factor=1)
    kes = [TR.kinetic_energy_coords(tgstruct,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
    epes = factor.*[TR.elastic_potential_energy(tgstruct,q) for q in sol.qs]
    gpes = factor.*[TR.gravity_potential_energy(tgstruct,q) for q in sol.qs]

    es = kes
    if elasticity
        es += epes
    end
    if gravity
        es += gpes
    end
    es1 = es[1]
    es_err = abs.((es.-es1)./es1)

    kes,epes,gpes,es,es_err
end

function string_potential(tgstruct,sol)
    [
    [begin
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q)
        TR.update_strings_apply_forces!(tgstruct)
        TR.potential_energy(s)
    end
    for q in sol.qs] for s in tgstruct.strings]
end

function analyse_slackness(tgstruct,q)
    TR.distribute_q_to_rbs!(tgstruct,q)
    for (i,s) in enumerate(tgstruct.strings)
        len = s.state.length
        restlen = s.state.restlen
        if len < restlen
            @info "String $i is slack."
        end
    end
end
