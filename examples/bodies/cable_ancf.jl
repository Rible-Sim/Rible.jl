function cable_ancf(pres_idx, 𝐞, L = 1.0) 
    radius = 2.0e-3
    # ancs = RB.ANCF.ANC3DRURU(8.96e3;E=110e6,L,radius)
    # mat_cable = filter(
    #     row->row.name == "Nylon 66",
    #     material_properties
    # )[1]
    # ρ = ustrip(Unitful.kg/Unitful.m^3,mat_cable.density)
    # E = ustrip(Unitful.Pa,mat_cable.modulus_elas)
    ρ = 1.03e3
    E = 0.2e9
    @myshow ρ, E
    ancs = RB.ANCF.ANC3DRURU(ρ;E,L,radius)
    mass = RB.ANCF.build_mass(ancs)
    @show mass
    T = typeof(L)
    mass_locus = SVector(L/2,T(0),T(0))
    r̄p1 = SVector(T(0),T(0),T(0))
    # r̄p2 = SVector(L,T(0),T(0))
    loci = [
        r̄p1,
    ]
    prop = RB.FlexibleBodyProperty(
        1,
        :cable,
        mass,
        mass_locus,
        # length(loci),
        loci
    )
    # cache = RB.BodyCache(prop,ancs,𝐞)
    state, coords, cache = RB.FlexibleBodyState(prop,ancs,𝐞;pres_idx)
    fb = RB.FlexibleBody(prop,state, coords, cache)
end