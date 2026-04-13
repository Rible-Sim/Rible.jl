function default_R(ri,rj,)
    u = rj - ri
    bar_length = norm(u)	
    û = u./bar_length
    v̂,ŵ = RB.NCF.HouseholderOrthogonalization(û)
    R = SMatrix{3,3}(hcat(û,v̂,ŵ))
end
function make_3d_bar(
        id,
        ri,rj;
        pres_idx = Int[],
        m = 0.080,
        e = 0.0,
        μ = 0.05,
        radius_ratio = 1/30,
        R = default_R(ri,rj),
        mat_name = nothing, #"Teak",
        loadmesh = false,
        loadmesh2 = false,
        coordsType = RB.NCF.NC3D2P,
        barcolor = :slategrey,
        loci_mode = :ends,
        Īg_override = nothing,
    )
    # @show id,ri,rj
    contactable = true
    visible = true
    u = rj - ri
    bar_length = norm(u)	
    û = u./bar_length
    radius = bar_length*radius_ratio

    mass_locus = SVector{3}([ bar_length/2,0,0])
    r̄p1  = SVector{3}([        0.0,  0,  0])
    r̄p2  = SVector{3}([ bar_length,  0,  0])
    if loci_mode == :ends
        loci = [r̄p1,r̄p2]
    else
        loci = [mass_locus]
    end
    density = NaN
    if mat_name isa Nothing
        mass = m
    else
        mat = filter(
            row -> row.name == mat_name, 
            RB.material_properties
        )[1]
        density = mat.density |> ustrip
        sec_area = π*radius^2
        density_per_unit_len = density*sec_area
        mass = density_per_unit_len*bar_length
    end
    bar_inertia_exp = :(1/12*$mass*$bar_length^2)
    bar_inertia = eval(bar_inertia_exp)
    if Īg_override isa Nothing
        Īg = SMatrix{3,3}(
            [	
                bar_inertia 0.0 0.0;
                0.0 0.0 0.0;
                0.0 0.0 0.0

            ]
        )
    else
        Īg = Īg_override
    end
    @debug "make_3d_bar" id radius density bar_length mass Īg bar_inertia_exp

    axes = [SVector(1.0,0.0,0.0) for _ in eachindex(loci)]
    prop = RB.RigidBodyProperty(id, contactable, mass, Īg, RB.Locus(mass_locus), [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(loci, axes, fill(μ, length(loci)), fill(e, length(loci)))]; visible=visible)
    # @show prop.inertia
    ro = ri
    ṙo = zero(ri)
    ω = zero(ri)
    # @show ri,rj,q0
    state = RB.RigidBodyState(prop,ri,R,ṙo,ω)
    if coordsType isa typeof(RB.NCF.NC3D2P)
        nmcs = RB.NCF.NC3D2P(ri,rj,ri,R)
        coords = RB.PresFreeCoordinates(nmcs,pres_idx)
    elseif coordsType isa typeof(RB.NCF.NC3D1P1V)
        nmcs = RB.NCF.NC3D1P1V(ri,û,ri,R)
        coords = RB.PresFreeCoordinates(nmcs,pres_idx)
    elseif coordsType isa typeof(RB.NCF.NC2P2V)
        v = R*[0.0,1.0,0.0]
        w = R*[0.0,0.0,1.0]
        nmcs = RB.NCF.NC2P2V(ri,rj,v,w,ro,R)
        coords = RB.PresFreeCoordinates(nmcs,pres_idx)
    else
        qcs = QCF.QC(mass,Īg[1]*I(3))
        coords = RB.PresFreeCoordinates(qcs,pres_idx)
    end
    if loadmesh
        barmesh = load(RB.assetpath("assembly3.STL")) |> RB.make_patch(;
            # trans=[0,0,0.025],
            scale=1/500,
            color=:palegreen3,
        )
    elseif loadmesh2
        barmesh = load(RB.assetpath("ball_joint/part2.STL")) |> RB.make_patch(; 
            scale=1/200,
            color=:palegreen3,
        )
    else
        barmesh = RB.endpoints2mesh(r̄p1,r̄p2;radius,color=barcolor)
    end
    body = RB.RigidBody(prop,state,coords,barmesh)
end
