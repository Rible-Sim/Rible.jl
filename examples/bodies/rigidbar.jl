function rigidbar(i,
        ri,rj;
        ṙi,ṙj,
        m = 0.05,
        μ = 0.9,
        e = 0.0,
        contactable = true,
        visible = true,
        ci = Int[],
        cstr_idx = [1],
        isbody = false,
        loadmesh = true,
    )
    ro = (ri + rj)/2
    ṙo = (ṙi + ṙj)/2
    u = rj - ri
    u̇ = ṙj - ṙi
    u /= norm(u)
    v,w = RB.HouseholderOrthogonalization(u)
    R = SMatrix{3,3}(hcat(u,v,w))
    # if isodd(pos)
    b = norm(rj-ri)
    mass_locus  = SVector{3}(   0,0.0,0.0)
    r̄p1 = SVector{3}(-b/2,0.0,0.0)
    r̄p2 = SVector{3}( b/2,0.0,0.0)
    r̄p3 = SVector{3}(   0,0.0,0.0)
    inertia = m*b^2/12
    Ī = SMatrix{3,3}([
        inertia 0          0;
        0       inertia    0;
        0       0     inertia
    ])
    loci = [r̄p1,r̄p2]
    axes = [SVector(1.0,0,0),SVector(1.0,0,0),SVector(1.0,0,0)]
    friction_coefficients = fill(μ,length(loci))
    restitution_coefficients = fill(e,length(loci))
    prop = RB.RigidBodyProperty(
        i,
        contactable,
        m,
        Ī,
        mass_locus,
        loci,
        axes,
        friction_coefficients,
        restitution_coefficients;
        visible = visible,
    )
    pretty_table(
        SortedDict(
            [
                ("id", i),
                ("bar length", b),
            ]
        )
    )
    if isbody
        @show isbody
        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        ω = zero(ṙo)
    else
        nmcs = RB.NCF.NC3D1P1V(ri, u, ro, R)
        ω = RB.NCF.find_angular_velocity(nmcs,vcat(ri,u),vcat(ṙi,u̇))
        # @show ω, ṙi, u̇
    end
    state = RB.RigidBodyState(prop, ro, R, ṙo, ω)
    coords = RB.NonminimalCoordinates(
        nmcs, ci, cstr_idx
    )
    # leg_mesh = load("400杆.STL")
    if loadmesh
        barmesh = load(joinpath(assetdir,"BZ.STL")) |> make_patch(;
            scale=1/1000,
            rot = RotY(π/2),
        )
    else
        barmesh = RB.endpoints2mesh(r̄p1,r̄p2;radius = norm(r̄p2-r̄p1)/40)
    end
    body = RB.RigidBody(prop, state, coords, barmesh)
end