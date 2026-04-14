function make_3d_tri(
        id,
        loci_positions,
        ro,
        R,
        ri,
        rj = nothing,
        rk = nothing,
        rl = nothing;
        contactable = true,
        visible = true,
        pres_idx = Int[],
        cstr_idx = collect(1:6),
        radius = 0.005,
        height = 1e-2,
        color = :darkorchid4,
        loadmesh = false,
        has_plate = true,
        m = 3,
        mass_center = [0,0,0.1562442983-0.1*√2]
    )
    # free_idx = collect(1:6)
    mass = 0.2999233976
    Īg = SMatrix{3,3}(
        Matrix(Diagonal([
            7.6639282053E-04,
            7.6638139752E-04,
            1.2464720496E-03
        ])
        )
    )

    axes = [
        SVector(1.0,0,0)
        for i = eachindex(loci_positions)
    ]
    
    prop = RB.RigidBodyProperty(id, contactable, mass, Īg, 
            RB.Locus(SVector{3}(mass_center)), 
            [RB.Locus(SVector{3}(pos), ax, mu, e) for (pos,ax,mu,e) in zip(
                loci_positions, 
                axes, 
                zeros(length(loci_positions)), 
                zeros(length(loci_positions))
            )]; 
        visible=visible
    )
    ṙo = zero(ro)
    ω = zero(ro)
    # u = R*(loci_positions[2] - loci_positions[1])
    # v = R*(loci_positions[3] - loci_positions[1])
    # w = R*(loci_positions[4] - loci_positions[1])
    if rj isa Nothing
        nmcs = RB.NCF.NC1P3V(ri,ro,R)
    elseif rk isa Nothing
        nmcs = RB.NCF.NC2P2V(ri,rj,ro,R)
    elseif rl isa Nothing
        nmcs = RB.NCF.NC3P1V(ri,rj,rk,ro,R)
    else
        nmcs = RB.NCF.NC4P(ri,rj,rk,rl,ro,R)
        @show Ref(R).*(loci_positions .+ Ref(ro))
        @show ri,rj,rk,rl
    end
    @debug "make_3d_tri" id mass Īg mass_center loci_positions
   
    # cf = RB.NCF.CoordinateFunctions(nmcs,q0,pres_idx,free_idx)
    # @show typeof(nmcs)
    # radius = norm(loci_positions[2]-loci_positions[1])/32
    # @show radius

    if loadmesh
        trimesh = load(
            RB.assetpath("ball_joint/part1.STL")
        ) |> RB.make_patch(;
            scale = 1/200, 
            color = :tomato1,
            trans = [0,0,-0.39]
        )
    else
        trimesh1 = GB.merge(
            [
                RB.endpoints2mesh(
                    loci_positions[i],loci_positions[j];
                    radius,color
                )
                for (i,j) in [
                    [1,2],[1,3],[1,4],
                    [2,3],[3,4],[4,2]
                ]
            ]
        ) |> RB.make_patch(;color = :darkslategrey)
        platemesh = RB.endpoints2mesh(
            SVector(0.0,0.0,-height/2),
            SVector(0.0,0.0, height/2);
            radius=norm(loci_positions[2]),
            n1=m,
            n2=2
        ) |> RB.make_patch(;color = :darkslategrey)
        if has_plate
            trimesh = GB.merge(
                [
                    trimesh1,
                    platemesh
                ]
            )
        else
            trimesh = trimesh1
        end
    end

    state = RB.RigidBodyState(prop,ro,R,ṙo,ω)
    coords = RB.PresFreeCoordinates(RB.NCF.NC(nmcs; cstr_idx=cstr_idx), pres_idx)
    body = RB.RigidBody(prop,state,coords,trimesh)
end
