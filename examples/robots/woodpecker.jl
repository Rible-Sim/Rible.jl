function woodpecker(;coordsType = RB.NCF.NC)
    SVo3 = SVector{3}([0.0,0.0,0.0])
    pole_radius = 0.0025e3#mm
    sleeve_radius = 0.0031e3 #mm
    half_sleeve_height = 0.0058e3 #mm
    distance_spring_COM_sleeve = 0.0100e3 #mm
    distance_spring_COM_woodpecker = 0.0150e3 #mm
    distance_beak_COM_woodpecker_z = 0.0200e3 #mm
    distance_beak_COM_woodpecker_y = 0.0201e3 #mm
    # inertial_properties
    sleeve_mass = 0.0003e3 #g 
    woodpecker_mass = 0.0045e3 #g 
    sleeve_moment_of_inertia = 5.0e-9*(1e9) #g*mm2 
    woodpecker_moment_of_inertia = 7.0e-7*(1e9) #g*mm2 
    force_elements_angular_stiffness = 0.0056e9 #(gmm/s^2)mm/rad 
    gravity = 9.81e3 #mm/s2 
    # contact_parameters
    COR_in_normal_direction = 0.5 
    COR_in_tangential_direction = 0
    frictional_coefﬁcient = 0.3
    # initial conditions 
    sleeve_angular_position = -0.1036 #rad 
    woodpecker_angular_position = -0.2788 #rad 
    sleeve_vertical_velocity = -0.3411e3 #mm/s
    sleeve_angular_velocity = 0.0 #rad/s 
    woodpecker_angular_velocity = -7.4583 #rad/s 
    sleeve_position_X = 0 #mm

    function make_sleeve(i;
            ro=SVo3,
            R=RotX(sleeve_angular_position),
            ω = SVo3,
            ṙo=SVector(0,0,sleeve_vertical_velocity)
        )
        contactable = true
        visible = true
        r̄g = SVo3
        # r̄p1 = SVector(0.0,-sleeve_radius,-half_sleeve_height)
        r̄p2 = SVector(0.0, sleeve_radius,-half_sleeve_height)
        r̄p3 = SVector(0.0, sleeve_radius, half_sleeve_height)
        # r̄p4 = SVector(0.0,-sleeve_radius, half_sleeve_height)
        r̄p5 = SVector(0.0, distance_spring_COM_sleeve, 0)
        loci_positions = [r̄p2,r̄p3,r̄p5]
        num_of_loci = length(loci_positions)
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            # SVector{3}([1.0,0.0,0.0]),
            # SVector{3}([1.0,0.0,0.0]),
        ]
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
    spring_position = sleeve1.state.loci_states[3].frame.position
    spring_velocity = sleeve1.state.loci_states[3].frame.velocity
    link1 = make_link(2;spring_position,spring_velocity)
    @show spring_velocity
    rbs = [sleeve1,link1]
    rigdibodies = TypeSortedCollection(rbs)
    q_sleeve1,_ = RB.body_state2coords_state(sleeve1)
    cstr_idx = [1,2,5,6]
    A_linearjoint = zeros(length(cstr_idx),12)
    for (i,j) in enumerate(cstr_idx)
        A_linearjoint[i,j] = 1
    end
    js = [
        RB.RevoluteJoint(1,RB.Hen2Egg(1,RB.ID(sleeve1,3,1),RB.ID(link1,1,1))),
        RB.LinearJoint(2,sleeve1,A_linearjoint,q_sleeve1[cstr_idx]),
        # RB.FixedBodyConstraint(2,indexed,sleeve1)
    ]
    jointed = RB.join(js,)
    ss = [
        RB.RotationalSpringDamper3D(
            1,
            [0.0,0,0],
            [3],
            force_elements_angular_stiffness,
        ),
        RB.RotationalSpringDamper3D(
            2,
            [0.0,0,0],
            Int[],
            force_elements_angular_stiffness,
        )
    ]
    force_elements = (spring_dampers = ss, cables = Int[])
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)
    indexed = RB.index(rigdibodies,jointed)
    numbered = RB.number(rigdibodies)
    cnt = RB.Connectivity(jointed,tensioned,indexed,numbered)
    st = RB.Structure(rigdibodies,force_elements,cnt)
    RB.Robot(st,)
end