
function make_3d_plate(
        id,
        loci_positions,
        ro,
        R,
        ri,
        rj=nothing,
        rk=nothing,
        rl=nothing;
        contactable = true,
        m = 3,
        height=1e-2,
        radius=1e-3,
        pres_idx = Int[],
        cstr_idx = collect(1:6),
        loadmesh = true,
        meshvisible = true,
    )
    # free_idx = collect(1:6)	
    visible = pres_idx != Int[]
    mat = filter(
        row -> row.name == "Teak", 
        RB.material_properties
    )[1]
    density = mat.density |> ustrip
    if meshvisible
        if loadmesh
            platemesh = load(RB.assetpath("wood_plate3.STL")) |> RB.make_patch(;rot=RotZ(π/4))
        else
            platemesh = RB.endpoints2mesh(
                SVector(0.0,0.0,-height/2),
                SVector(0.0,0.0,height/2);
                radius,n1=m,n2=2)
        end
    else
        platemesh = nothing
    end
    # mass = density*GB.volume(platemesh)
    if m == 4
        mass =  0.11008275 #kg
        Īg = SMatrix{3,3}(
            Matrix(Diagonal([
                0.00025403,
                0.00025403,
                0.00050623
            ])
            )
        )
    else
        mass = 0.23569851 #kg
        Īg = SMatrix{3,3}(
            Matrix(Diagonal([
                0.00085048,
                0.00085048,
                0.00169704				
            ])
            )
        )
    end
    mass_center = SVector(0.0,0.0,0.0)
    @debug "make_3d_plate" id mass Īg m radius height density 

    prop = RB.RigidBodyProperty(id, contactable, mass, Īg, 
        RB.Locus(mass_center), 
        [RB.Locus(SVector{3}(pos), ax, mu, e) for (pos,ax,mu,e) in zip(
            loci_positions, 
            [SVector(1.0,0.0,0.0) for _ in eachindex(loci_positions)], 
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
    nmcs = RB.NCF.NC1P3V(ri,ro,R)
    state = RB.RigidBodyState(prop,ro,R,ṙo,ω)
    coords = RB.PresFreeCoordinates(RB.NCF.NC(nmcs; cstr_idx=cstr_idx), pres_idx)
    body = RB.RigidBody(prop,state,coords,platemesh)
end
