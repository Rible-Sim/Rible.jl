function woodpecker(;coordsType = RB.NCF.NC)
    SVo3 = SVector{3}([0.0,0.0,0.0])
    pole_radius = 0.0025
    pole_height = 0.04
    sleeve_radius = 0.0031
    half_sleeve_height = 0.0058
    distance_spring_COM_sleeve = 0.0100
    distance_spring_COM_woodpecker = 0.0150
    distance_beak_COM_woodpecker_z = 0.0200
    distance_beak_COM_woodpecker_y = 0.0201
    # inertial_properties
    sleeve_mass = 0.0003 #kg 
    woodpecker_mass = 0.0045 #kg 
    sleeve_moment_of_inertia = 5.0e-9 #kg/m2 
    woodpecker_moment_of_inertia = 7.0e-7 #kg/m2 
    force_elements_angular_stiffness = 0.0056 #Nm/rad 
    gravity = 9.81 #m/s2 
    # contact_parameters
    COR_in_normal_direction = 0.5 
    COR_in_tangential_direction = 0
    frictional_coefﬁcient = 0.3
    # initial conditions 
    sleeve_angular_position = -0.1036 #rad 
    woodpecker_angular_position = -0.2788 #rad 
    woodpecker_vertical_velocity = -0.3411 #m/s
    sleeve_angular_velocity = 0.0 #rad/s 
    woodpecker_angular_velocity = -7.4583 #rad/s 
    sleeve_position_X = 0 #m
    
    function make_pole(i)
        contactable = false
        visible = true
        r̄g = SVo3
        r̄p1 = SVector{3}([0, -pole_radius, pole_height])
        r̄p2 = SVector{3}([0, -pole_radius, pole_height])
        r̄p3 = SVector{3}([0,  pole_radius, pole_height])
        r̄p4 = SVector{3}([0,  pole_radius, pole_height])
        loci_positions = [r̄p1,r̄p2,r̄p3,r̄p4,]
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = sleeve_mass
        Ī = SMatrix{3,3}(sleeve_moment_of_inertia*I(3))
        R = RotX(0.0)
        ω = SVo3
        ri = SVo3
        ro = ri
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals;
            visible
        )

        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
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

        # polemesh = load(RB.assetpath("crank_sleeve/pole.STL")) |> RB.make_patch(;
        #     # trans=[-1.0,0,0],
        #     # rot = RotZ(π),
        #     scale=1/1000,
        #     color = :silver,
        # )
        polemesh = nothing
        RB.RigidBody(prop,state,coords,polemesh)
    end
    function make_sleeve(i;
            ro=SVo3,
            R=RotX(sleeve_angular_position),
            ω = SVo3,
            spring_velocity = SVo3,
        )
        contactable = true
        visible = true
        r̄g = SVo3
        r̄p1 = SVector(0.0,-sleeve_radius,-half_sleeve_height)
        r̄p2 = SVector(0.0, sleeve_radius,-half_sleeve_height)
        r̄p3 = SVector(0.0, sleeve_radius, half_sleeve_height)
        r̄p4 = SVector(0.0,-sleeve_radius, half_sleeve_height)
        r̄p5 = SVector(0.0, distance_spring_COM_sleeve, 0)
        loci_positions = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]
        num_of_loci = length(loci_positions)
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = sleeve_mass
        Ī = SMatrix{3,3}(sleeve_moment_of_inertia*I(3))
        ri = ro
        
        ṙo=spring_velocity

        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals,
            fill(frictional_coefﬁcient,num_of_loci),
            fill(COR_in_normal_direction,num_of_loci);
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
            ṙo = SVector(0,0,woodpecker_vertical_velocity),
        )
        contactable = true
        visible = true
        r̄g  = SVector{3}([0.0, 0,     0])
        r̄p1 = SVector{3}([0.0, -distance_spring_COM_woodpecker, 0])
        r̄p2 = SVector{3}([0.0, -distance_beak_COM_woodpecker_y, distance_beak_COM_woodpecker_z])
        loci_positions = [r̄p1,r̄p2,]
        num_of_loci = length(loci_positions)
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = woodpecker_mass
        Ī = SMatrix{3,3}(woodpecker_moment_of_inertia*I(3))
        ro = spring_position+R*(-r̄p1)
        ri = ro

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
    pole = make_pole(1)
    spring_position = make_sleeve(2;).state.loci_states[5].frame.position
    link1 = make_link(3;spring_position)
    spring_velocity = link1.state.loci_states[1].frame.velocity
    @show spring_velocity
    sleeve1 = make_sleeve(2;spring_velocity)
    rbs = [pole,sleeve1,link1]
    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
    indexed = RB.index(rigdibodies,)

    ss = [
        RB.RotationalSpringDamper3D(
            1,
            [0.0,0,0],
            Int[],
            0.0,
        ),
        RB.RotationalSpringDamper3D(
            2,
            [0.0,0,0],
            [3],
            force_elements_angular_stiffness,
        )
    ]
    tensiles = (spring_dampers = ss, cables = Int[])
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)

    j1 = RB.FixedBodyConstraint(1,indexed,pole)
    j2 = RB.RevoluteJoint(2,indexed,RB.Hen2Egg(2,RB.ID(sleeve1 ,5,1),RB.ID(link1,1,1)))
 
    js = [j1,j2]

    jointed = RB.join(js,indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,)
end