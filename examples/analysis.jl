function analyse_energy(tgstruct,sol;gravity=false,elasticity=false)
    kes = [TR.kinetic_energy_coords(tgstruct,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]

    epes = [TR.elastic_potential_energy(tgstruct,q) for q in sol.qs]
    gpes = [TR.gravity_potential_energy(tgstruct,q) for q in sol.qs]


    es = kes
    if elasticity
        es += epes
    end
    if gravity
        es += gpes
    end
    es_err = (abs.(es.-es[1]))./es

    kes,epes,gpes,es,es_err
end
