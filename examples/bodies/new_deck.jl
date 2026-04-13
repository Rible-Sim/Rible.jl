function new_deck(id,loci,ro,R,ri,box;
        contactable = true,
        visible = true,
        pres_idx = Int[],
        cstr_idx = collect(1:6),
    )
    mat = filter(
        row -> row.name == "Teak", 
        RB.material_properties
    )[1]
    density = mat.density |> ustrip
    box_volume = Meshes.volume(box) |> ustrip
    mass = density*box_volume
    # bar_inertia_exp = :(1/12*$mass*$bar_length^2)
    wid,len,hei = Meshes.sides(box) .|> ustrip
    # m = 0.2999233976
    Īg = SMatrix{3,3}(
        Matrix(
            Diagonal(
                [
                    mass*(hei^2+len^2)/12, #Ix
                    mass*(hei^2+wid^2)/12, #Iy
                    mass*(len^2+wid^2)/12, #Iz
                ]
            )
        )
    )
    mass_locus = SVector(0,0,0.0)
    @debug(
        SortedDict(
            [
                ("id", id),
                ("mass", mass),
                ("density", density),
                ("inertia", Īg),
                ("wid,len,hei", (wid,len,hei)),
                ("box_volume", box_volume)
            ]
        )
    )
    axes_normals = [
        SVector(1.0,0,0)
        for _ in loci
    ]
    prop = RB.RigidBodyProperty(id, contactable, mass, Īg, RB.Locus(mass_locus), [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(loci, axes_normals, zeros(length(loci)), zeros(length(loci)))]; visible=visible)
    ṙo = zero(ro)
    ω = zero(ro)
    nmcs = RB.NCF.NC1P3V(ri,ro,R)
    # cf = RB.NCF.CoordinateFunctions(nmcs,q0,pres_idx,free_idx)
    # @show typeof(nmcs)
    boxmesh = Meshes.boundary(box) |> RB.simple2mesh
    state = RB.RigidBodyState(prop, ro, R, ṙo, ω)
    coords = RB.PresFreeCoordinates(
        RB.NCF.NC(nmcs; cstr_idx=cstr_idx),
        pres_idx
    )
    RB.RigidBody(prop,state,coords,boxmesh)
end
