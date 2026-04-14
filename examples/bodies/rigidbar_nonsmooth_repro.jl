function rigidbar(i,
        ri,rj;
        ṙi,ṙj,
        m = 0.05,
        μ = 0.9,
        e = 0.0,
        contactable = true,
        visible = true,
        loadmesh = true,
        isbody = false, #note not use
    )
    ro = (ri + rj)/2
    ṙo = (ṙi + ṙj)/2
    u = rj - ri
    u̇ = ṙj - ṙi
    dist = norm(u)
    u /= dist
    #note u̇ /= dist
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
    loci = [r̄p1,r̄p2,r̄p3]
    axes = [SVector(1.0,0,0),SVector(1.0,0,0),SVector(1.0,0,0)]
    friction_coefficients = fill(μ,length(loci))
    restitution_coefficients = fill(e,length(loci))
    prop = RB.RigidBodyProperty(i, contactable, m, Ī, RB.Locus(mass_locus), [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(loci, axes, friction_coefficients, restitution_coefficients)]; visible = visible)

    nmcs = RB.NCF.NC3D1P1V(ri, u, ro, R)
    #note ω = RB.NCF.find_angular_velocity(nmcs,vcat(ri,u),vcat(ṙi,u̇))
    q = vcat(ri,u)
    q̇ = vcat(ṙi,u̇)

    Ẋ =  RB.NCF.get_X(nmcs,q̇)
    X =  RB.NCF.get_X(nmcs,q)
    Ω = Ẋ*pinv(X)
    ω = SVector{3}(Ω[3,2],Ω[1,3],Ω[2,1])

    # @show i u̇ ω×u
    @debug "rigidbar" id=i bar_length=b isbody ri ṙi ω u̇ ω×u
    state = RB.RigidBodyState(prop, ro, R, ṙo, ω)
    coords = nmcs
    # leg_mesh = load("400杆.STL")
    if loadmesh
        barmesh = load(RB.assetpath("BZ.STL")) |> RB.make_patch(;
            scale=1/1000,
            rot = RotY(π/2),
        )
    else
        barmesh = RB.endpoints2mesh(r̄p1,r̄p2;radius = norm(r̄p2-r̄p1)/40)
    end
    body = RB.RigidBody(prop, state, coords, barmesh)

end