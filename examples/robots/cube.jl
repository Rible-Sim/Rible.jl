function make_cube(origin_position = [0.0,0.0,0.5],
        R = one(RotMatrix{3}),
        origin_velocity = [2.0,0.0,0.0],
        ω = [0.0,0.0,5.0];
        μ = 0.9,
        e = 0.0
    )
    m = 1.0
    mass_locus = zeros(3)
    inertia = SMatrix{3,3}(Matrix(1.0I,3,3))
    h = 5.0
    loci = [h.*[x,y,z] for z = [-1,1] for y = [-1,1] for x = [-1,1]]
    axes = [SVector(1.0,0,0) for _ = 1:length(loci)]
    friction_coefficients = [μ for _ = 1:length(loci)]
    restitution_coefficients = [e for _ = 1:length(loci)]
    pts = Point3.(loci)
    fcs = GB.QuadFace.([
        [1,3,4,2],
        [3,7,8,4],
        [7,5,6,8],
        [5,1,2,6],
        [2,4,8,6],
        [1,5,7,3]
    ])
    nls = GB.normals(pts,fcs)
    cube_mesh = GB.Mesh(GB.meta(pts,normals=nls),fcs) |> RB.make_patch(;color=:snow)
    prop = RB.RigidBodyProperty(
        1,true,m,inertia,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients,
        )
    ri = copy(origin_position)
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)

    nmcs = RB.NCF.NC1P3V(ri,origin_position,R)
    coords = RB.NonminimalCoordinates(nmcs,)
    rb1 = RB.RigidBody(prop,state,coords,cube_mesh)
    rbs = TypeSortedCollection((rb1,))

    apparatuses = Int[]
    sharing_matrix = zeros(Int,0,0)
    indexed = RB.index(rbs,apparatuses;sharing_matrix)
    numbered = RB.number(rbs,apparatuses)

    cnt = RB.Connectivity(numbered,indexed,)

    st = RB.Structure(rbs,apparatuses,cnt,)
    
    gauges = Int[]
    actuators = Int[]
    ## actuators = [
    ##     RB.RegisterActuator(
    ##         1,
    ##         collect(1:ncables_prism),
    ##         zeros(ncables_prism)
    ##     ),
    ## ]
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )
    bot = RB.Robot(st,hub)
end