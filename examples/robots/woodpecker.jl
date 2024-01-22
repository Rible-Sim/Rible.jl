function woodpecker(;
        coordsType = RB.NCF.NC,
        case=2, 
        frictional_coefﬁcient = 0.3
    )
    SVo3 = SVector{3}([0.0,0.0,0.0])
    pole_radius = 0.0025#m
    sleeve_radius = 0.0031 #m
    half_sleeve_height = 0.0058 #m
    distance_spring_COM_sleeve = 0.0100 #m
    distance_spring_COM_woodpecker = 0.0150 #m
    distance_beak_COM_woodpecker_z = 0.0200 #m
    distance_beak_COM_woodpecker_y = 0.0201 #m
    # inertial_properties
    sleeve_mass = 0.0003*(1e3) #g 
    woodpecker_mass = 0.0045*(1e3) #g 
    sleeve_moment_of_inertia = 5.0e-9*(1e3) #g*m2 
    woodpecker_moment_of_inertia = 7.0e-7*(1e3) #g*m2 
    force_elements_angular_stiffness = 0.0056*(1e3) #(gm/s^2)m/rad 
    gravity = 9.81 #m/s2
    # contact_parameters
    COR_in_normal_direction = 0.5
    COR_in_tangential_direction = 0
    # initial conditions 
    sleeve_angular_position = -0.1036 #rad 
    woodpecker_angular_position = -0.2788 #rad 
    sleeve_vertical_velocity = -0.3411 #m/s
    sleeve_angular_velocity = 0.0 #rad/s 
    woodpecker_angular_velocity = -7.4583 #rad/s 
    sleeve_position_X = 0 #m

    function make_sleeve(i;
            ro=SVo3,
            R=RotX(sleeve_angular_position),
            ω = SVo3,
            ṙo=SVector(0,0,sleeve_vertical_velocity)
        )
        contactable = true
        visible = true
        r̄g = SVo3
        r̄p1 = SVector(0.0,-sleeve_radius,-half_sleeve_height)
        r̄p2 = SVector(0.0, sleeve_radius,-half_sleeve_height)
        r̄p3 = SVector(0.0, sleeve_radius, half_sleeve_height)
        r̄p4 = SVector(0.0,-sleeve_radius, half_sleeve_height)
        r̄p5 = SVector(0.0, distance_spring_COM_sleeve, 0)
        loci_positions = ifelse(
            case==2,
            [r̄p2,r̄p3,r̄p5],
            [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]
        )
        num_of_loci = length(loci_positions)
        axes_normals = fill(
            SVector{3}([1.0,0.0,0.0]),
            num_of_loci
        )
        m = sleeve_mass
        Ī = SMatrix{3,3}(
            [
                sleeve_moment_of_inertia 0 0;
                0 sleeve_moment_of_inertia 0;
                0 0 sleeve_moment_of_inertia;
            ]
        )
        ri = ro

        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals,
            fill(frictional_coefﬁcient,num_of_loci),
            fill(0.0,num_of_loci);
            visible,
        )
        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            coords = RB.NonminimalCoordinates(nmcs,)
        else
            qcs = RB.QCF.QC(m,Ī)
            coords = RB.NonminimalCoordinates(qcs,)
        end
        # sleevemesh = load(RB.assetpath("crank_sleeve/sleeve.STL")) |> RB.make_patch(;
        #     trans=[0.0,0,0],
        #     scale=1/1000,
        #     color=:mediumpurple4
        # )
        sleevemesh = nothing
        RB.RigidBody(prop,state,coords,sleevemesh)
    end
    function make_link(i;
            ω = SVector{3}(woodpecker_angular_velocity,0,0),
            R = RotX(woodpecker_angular_position),
            spring_position = SVo3,
            spring_velocity = SVo3,
        )
        contactable = true
        visible = true
        r̄g  = SVo3
        r̄p1 = SVector{3}([0.0, -distance_spring_COM_woodpecker, 0])
        r̄p2 = SVector{3}([0.0, -distance_beak_COM_woodpecker_y, distance_beak_COM_woodpecker_z])
        loci_positions = [r̄p1,r̄p2,]
        num_of_loci = length(loci_positions)
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = woodpecker_mass
        Ī = SMatrix{3,3}(
            [
                woodpecker_moment_of_inertia 0 0;
                0 woodpecker_moment_of_inertia 0;
                0 0 woodpecker_moment_of_inertia;
            ]
        )
        ro = spring_position+R*(-r̄p1)
        ri = ro
        ṙo = spring_velocity+ω×(R*(-r̄p1))
        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals,
            fill(frictional_coefﬁcient,num_of_loci),
            fill(COR_in_normal_direction,num_of_loci);
            visible
        )

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            pres_idx = Int[] 
            cstr_idx = collect(1:6)
            coords = RB.NonminimalCoordinates(nmcs, pres_idx, cstr_idx)
        else
            qcs = RB.QCF.QC(m,Ī)
            pres_idx = Int[]
            cstr_idx = [1]
            coords = RB.NonminimalCoordinates(qcs, pres_idx, cstr_idx)
        end

        # linkmesh = load(RB.assetpath("crank_sleeve/crank$(i-1).STL")) |> RB.make_patch(;
        #     trans=[
        #         ifelse(
        #             i == 2,
        #             0.01*2,
        #             0.01
        #         ),
        #         -l[i-1]/2,
        #         0
        #     ],
        #     # rot = RotZ(π),
        #     scale=1/1000,
        #     color = :slategrey,
        # )
        linkmesh = nothing
        RB.RigidBody(prop,state,coords,linkmesh)
    end
    # pole = make_pole(1)
    sleeve1 = make_sleeve(1;)
    spring_loci_id = ifelse(case==2,3,5)
    spring_position = sleeve1.state.loci_states[spring_loci_id].frame.position
    spring_velocity = sleeve1.state.loci_states[spring_loci_id].frame.velocity
    link1 = make_link(2;spring_position,spring_velocity)
    @show spring_velocity
    rbs = [sleeve1,link1]
    rigdibodies = TypeSortedCollection(rbs)
    q_sleeve1,_ = RB.body_state2coords_state(sleeve1)
    if case == 2
        cstr_idx = [1,2,5,6]
    else
        cstr_idx = [1,5,6]
    end
    j1 = RB.RevoluteJoint(
        1,
        RB.Hen2Egg(RB.Signifier(sleeve1,spring_loci_id,1),RB.Signifier(link1,1,1)),
        RB.RotationalSpringDamper3D(
            [0.0,0,0],
            [3],
            force_elements_angular_stiffness,
        )
    )
    j2 = RB.FixedIndicesConstraint(2,sleeve1,cstr_idx,q_sleeve1[cstr_idx])
    apparatuses = TypeSortedCollection([j1,j2])
    indexed = RB.index(rigdibodies,apparatuses)
    numbered = RB.number(rigdibodies,apparatuses)
    cnt = RB.Connectivity(indexed,numbered)
    st = RB.Structure(rigdibodies,apparatuses,cnt)
    RB.Robot(st,)
end