function RigidData(r̄g,m,i,r̄l,mesh=nothing)
    @eponymtuple(r̄g,m,i,r̄l,mesh)
end

# Define the second-version modular arm with 3 axial springs and 3 radial springs
function build_1_man(n, ks1, restlenc, delta, bodycolor=:lightcyan3)# delta is the initial overlap distance between two units: 69.613
    restlens1 = 25.0
    cs1 = 2.0 #1195
    ks2 = 0.037
    restlens2 = 55.0
    cs2 = 2.0
    kc = 10.0
    c_c = 0.0
    pre = 3.0
    μ = 0.0001

    num_Y_modul = 3
    row_outer_spring = 3
    row_inter_spring = 3
    nb = (num_Y_modul + 1) * n
    # Cable-driven module: CoM-to-bottom distance; CoM-to-top distance
    len1_b = -30.54
    len1_t = -49.96
    # r1 outer spring, r2 cable hole, r4 inner spring
    r1 = 37.3
    r2 = 36.5
    r4 = 31
    # One Y-shaped part: CoM-to-top distance; CoM-to-bottom distance; bottom spring hole; r42 inner spring
    len2_t = -15.57
    len2_b = -60.64
    r3 = 5.0
    r42 = 30

    # CoM coordinates
    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]
    # r̄g3 = [0.0, 0.0, 0.0]
    # r̄g4 = [0.0, 0.0, 0.0]
    # Mass properties
    m1 = 0.12
    m2 = 0.027
    m3 = 0.1 #m4 = 0.15
    # Inertia tensor
    i1 =  [4575.03 0 0; 0 5114.84 0; 0 0 4575.03]## base first section
    i2 = [19.333 0 0; 0 19.343 0; 0 0 15.492]# .* 0.01## single module

    # Key connection points
    ## Cable-driven module
    r̄l1 = [ 
        # Top outer spring holes
        [r1 * cos(pi / 3),  r1 * sin(pi / 3), len1_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len1_t],
        [-r1, 0.0, len1_t],
        # Top cable holes
        [ r2 * sqrt(2) / 2,  r2 * sqrt(2) / 2, len1_t],
        [ r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2,  r2 * sqrt(2) / 2, len1_t],
        # Top inner spring holes
        [ r4 * cos(pi / 3),  r4 * sin(pi / 3), len1_t],
        [ r4 * cos(pi / 3), -r4 * sin(pi / 3), len1_t],
        [-r4, 0.0, len1_t],
    ]
    # Y-shaped part
    r̄l2 = [ # Top outer springs
        [ r1 * cos(pi / 3),  r1 * sin(pi / 3), len2_t],
        [ r1 * cos(pi / 3), -r1 * sin(pi / 3), len2_t],
        [-r1, 0.0, len2_t],
        # Top cable holes
        [ r2 * sqrt(2) / 2,  r2 * sqrt(2) / 2, len2_t],
        [ r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2,  r2 * sqrt(2) / 2, len2_t],
        # Top inner spring holes
        [r42 * cos(pi / 3),  r42 * sin(pi / 3), len2_t],
        [r42 * cos(pi / 3), -r42 * sin(pi / 3), len2_t],
        [-r42, 0.0, len2_t],
        # Bottom inner springs
        [ r3 * cos(pi / 3),  r3 * sin(pi / 3), -len2_b],
        [ r3 * cos(pi / 3), -r3 * sin(pi / 3), -len2_b],
        [-r3, 0.0, -len2_b],
        [0.0,0.0,len2_t]
    ]

    np_rl1 = length(r̄l1)
    np_rl2 = length(r̄l2)

    mesh1 = load(RB.assetpath("octopus/actuator-v2.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=:lightcyan3)#
    mesh2 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor)#

    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)
    rd3 = RigidData(r̄g2, m1, i1, r̄l2, mesh2)
    
    ## Define the full arm; current issue: it is not defined by segments and has no rotational DOF
    
    if n == 1
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul-1], [rd2]])
    else
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul], [rd1], [rd2 for _ in 1:num_Y_modul2] for _ in 1:n-1])
    end
    # @show rds

    function get_r̄s(rds)# Define each module position in the initial configuration
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][3] = len1_b
        for i in 1:nb-1
            if length(rds[i].r̄l) == np_rl2 # 13 is the number of key node coordinates in the Y-shaped module
                r̄s[i+1][3] += r̄s[i][3]  + len2_t+len2_b - delta
            else
                length(rds[i].r̄l) == np_rl1 # 10 is the number of key node coordinates in the actuator module
                r̄s[i+1][3] += r̄s[i][3] + len2_b + len1_t  - delta
            end
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        movable = true
        ci = Vector{Int}() # What is this? pres_idx?
        Φi = collect(1:6) # What is this?

        ri = r̄si # Initial global CoM coordinates of each module
        # @show ri
        aps = rdsi.r̄l # Elastic cable connection positions

        theta = deg2rad(90) # Initial deflection of each section; set to zero here if needed
        R = [
             cos(theta) sin(theta) 0;
            -sin(theta) cos(theta) 0;
            0 0 1
        ]

        r̄g = rdsi.r̄g # CoM coordinates
        m = rdsi.m # Mass
        I = rdsi.i # Inertia tensor
        nrp = length(aps) # Number of nodes
        ṙo = zeros(3)
        ω = zeros(3)
        prop = RB.RigidBodyProperty(i, movable, m, SMatrix{3,3}(I), 
            RB.Locus(SVector{3}(r̄g)), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
                [SVector{3}(aps[i]) for i in 1:nrp], 
                [SVector{3}([0.0, 0.0, 1.0]) for i in 1:nrp], #normals, 
                zeros(nrp), 
                zeros(nrp)
            )]
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω, )
        nmcs = RB.NCF.NC1P3V(ri, ro,)# Initial natural coordinates?
        # coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        coords = nmcs
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)
    ## Shared coordinates?
    # matrix_sharing_raw = Vector{Matrix{Int}}()
    # for i = 1:n-1
    #     s = zeros(6,nb)
    #         s[1:6,6i] = 1:6
    #         s[1:6,6i+1] = 1:6
    #     push!(matrix_sharing_raw, s)
    # end
    # sharing_matrix = reduce(vcat,matrix_sharing_raw)

    # Springs##########################################################
    ncables = (nb - 1) * (row_outer_spring + row_inter_spring) # This refers to inner and outer springs
    vks1 = zeros((nb - 1) * row_inter_spring,1)
    vks1 .= ks1
    for i in 1:3:(nb - 1) * row_inter_spring
        vks1[i] = ks1     # vks1[i] = ks1
    end

    restlens1 = restlens1 # inner springs
    ks2 = ks2
    restlens2 = restlens2 # outer springs
    spring_dampers = [(i <= (nb - 1) * row_inter_spring) ? RB.DistanceSpringDamper3D(restlens1, vks1[i], cs1) : RB.DistanceSpringDamper3D(restlens2, ks2, cs2) for i in 1:ncables] # First half are inner springs, second half are outer springs

    matrix_cnt_raw = Vector{Matrix{Int}}()
    # Define inner/outer spring connections
    for j in 1:num_Y_modul # connect inner springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
        # end
    end
    # for j in num_Y_modul+3:nb-1 # connect inner springs
    #     s = zeros(3, nb)
    #     s[1, j] = 8
    #     s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
    #     s[2, j] = 9
    #     s[2, j+1] = -12
    #     s[3, j] = 10
    #     s[3, j+1] = -13
    #     push!(matrix_cnt_raw, s)
    # end
    for j in 1:num_Y_modul # connect outer springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        # s[4, j] = 4
        # s[4, j+1] = -4
        push!(matrix_cnt_raw, s)
        # end
    end
    # for j in num_Y_modul+3:nb-1
    #     s = zeros(4, nb)
    #     s[1, j] = 1
    #     s[1, j+1] = -1
    #     s[2, j] = 2
    #     s[2, j+1] = -2
    #     s[3, j] = 3
    #     s[3, j+1] = -3
    #     s[4, j] = 4
    #     s[4, j+1] = -4
    #     push!(matrix_cnt_raw, s)
    # end
    
    connecting_matrix = RT.change_connecting_format(rigidbodies, reduce(vcat, matrix_cnt_raw))
    
    # Sliding cable - begin #############################################
    # restlencs = zeros(nb-1) # Each segment has five sliding cable sections, but the first cable section is actually absent
    restlencs = zeros(nb - n)
    # restlenc1 = 63.77 # The middle rotating section differs
    for i in 1:nb-n #nb-1
        # if (mod(i,num_Y_modul+2) == 0) # The fifth section is the rotating-section length; it does not exist in the actual model
        #     restlencs[i] = restlenc1 
        # if (mod(i, num_Y_modul + 1) == 0) && (i != nb - n)
            # restlencs[i] = restlenc2
        # else
            restlencs[i] = restlenc
        # end
    end

    nclusters_per_seg = 4*n
    
    cluster_segs = [
        [RT.DistanceSpringDamperSegment(restlencs[i], kc,c=c_c, prestress=pre) for i in 1:nb-n]
        for _ = 1:nclusters_per_seg
    ]

    cluster_sps = [
        [RT.SlidingPoint(μ) for _ in 1:num_Y_modul-1]
        for _ = 1:nclusters_per_seg
    ]
    
  
    cluster_matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:nclusters_per_seg]
    @show nb
    s = zeros(Int, nb-1, 4, nclusters_per_seg) # segments; 4[bpbp]; cable count
    for k in 1:nclusters_per_seg
        for l = 1:nb-1
            s[l, :, k] = [l    , 3 + k, l + 1,  3 + k]
        end
    end
    for k in 1:nclusters_per_seg
        push!(cluster_matrix_cnt_raw[k], s[:, :, k])
    end
    # Sliding cable - end ###############################################
    display(cluster_matrix_cnt_raw[1])
    connecting_cluster_matrix = [reduce(vcat, cluster_matrix_cnt_raw[k]) for k in 1:nclusters_per_seg]
    cables, clusters = RT.connect_spring_and_clusters(
        rigidbodies, spring_dampers, 
        cluster_sps, cluster_segs, 
        connecting_matrix, connecting_cluster_matrix, 
        istart= nclusters_per_seg 
    )
    
    # cst2 = RB.FixedBodyApparatus(77, rbs[14])
    # bjs = [
    # 	RB.PinJoint(RB.End2End(i,RB.ID(rbs[i],2),RB.ID(rbs[i+1],1)))
    # 	for i = 1:nbodies-1
    # ]

    # Fix_joint = RB.FixedJoint(RB.End2End(1,RB.ID(rbs[6],13),RB.ID(rbs[7],13)))
    # joint = RB.join( Fix_joint, indexedcoords)

    # Rev_joint = RB.RevoluteJoint(RB.End2End(1,RB.ID(rbs[6],13,1),RB.ID(rbs[7],13,1)))
    cst1 = RB.FixedBodyApparatus(nclusters_per_seg + ncables + 1, rbs[1])# front
    apparatuses = TypeSortedCollection(vcat([cst1], cables, clusters))
    cnt = RB.Connectivity(rigidbodies, apparatuses; ) #？
    
    # indexed = RB.index(rigidbodies, apparatuses;) #？
    # numbered = RB.number(rigidbodies, apparatuses)

    # cnt = RB.Connectivity(numbered, indexed)
    structure = RB.Structure(rigidbodies, apparatuses, cnt)

    capta_gauges = Int[]
    error_gauges = Int[] # Gauges/sensors used to measure errors for cost functions, optimization, or feedback control; currently unused.
    actuators = reduce(vcat,[ # Actuators attached to one or more bodies or apparatuses; their types dispatch to execute! for different driving forces
        # Try whether ExternalForceActuator is sufficient; if not, create a new actuator
        [RB.ExternalForceActuator(
            i,
            cluster,
            # Operator, 1 control input
            RB.NaiveOperator(1),
            # External force, may not be needed
            0.0,
            # Control input values
            [0.0],
        )
        for (i, cluster) in enumerate(clusters)]
        ,
        RB.ExternalForceActuator(
            nclusters_per_seg +1,
            RB.Signifier(rbs[end], np_rl2),
            RB.NaiveOperator(1),
            [0;-1.0;0;;],
            [0.0],
            )]
    )
    hub = RB.ControlHub(
        structure,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(structure, gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)

end

function build_1_man_flCTC(n, ks1, restlenc, delta, bodycolor=:lightcyan3)# delta is the initial overlap distance between two units: 69.613
    restlens1 = 25.0
    cs1 = 2.0 #1195
    ks2 = 0.037
    restlens2 = 55.0
    cs2 = 2.0
    kc = 10.0
    @show c_c = 0.0
    pre = 3.0
    @show μ = 0.0001

    num_Y_modul = 6
    row_outer_spring = 3
    row_inter_spring = 3
    nb = (num_Y_modul + 1) * n
    # Cable-driven module: CoM-to-bottom distance; CoM-to-top distance
    len1_b = -30.54
    len1_t = -49.96
    # r1 outer spring, r2 cable hole, r4 inner spring
    r1 = 37.3
    r2 = 36.5
    r4 = 31
    # One Y-shaped part: CoM-to-top distance; CoM-to-bottom distance; bottom spring hole; r42 inner spring
    len2_t = -15.57
    len2_b = -60.64
    r3 = 5.0
    r42 = 30

    # CoM coordinates
    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]
    # r̄g3 = [0.0, 0.0, 0.0]
    # r̄g4 = [0.0, 0.0, 0.0]
    # Mass properties
    m1 = 0.12
    @show m2 = 0.027
    m3 = 0.1 #m4 = 0.15
    # Inertia tensor
    i1 = [4575.03 0 0; 0 5114.84 0; 0 0 4575.03]## base first section
    i2 = [19.333 0 0; 0 19.343 0; 0 0 15.492]# .* 0.01## single module

    # Key connection points
    ## Cable-driven module
    r̄l1 = [
        # Top outer spring holes
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len1_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len1_t],
        [-r1, 0.0, len1_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        # Top inner spring holes
        [r4 * cos(pi / 3), r4 * sin(pi / 3), len1_t],
        [r4 * cos(pi / 3), -r4 * sin(pi / 3), len1_t],
        [-r4, 0.0, len1_t],
    ]
    # Y-shaped part
    r̄l2 = [ # Top outer springs
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len2_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len2_t],
        [-r1, 0.0, len2_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        # Top inner spring holes
        [r42 * cos(pi / 3), r42 * sin(pi / 3), len2_t],
        [r42 * cos(pi / 3), -r42 * sin(pi / 3), len2_t],
        [-r42, 0.0, len2_t],
        # Bottom inner springs
        [r3 * cos(pi / 3), r3 * sin(pi / 3), -len2_b],
        [r3 * cos(pi / 3), -r3 * sin(pi / 3), -len2_b],
        [-r3, 0.0, -len2_b],
        [0.0, 0.0, len2_t]
    ]

    np_rl1 = length(r̄l1)
    np_rl2 = length(r̄l2)

    mesh1 = load(RB.assetpath("octopus/actuator-v2.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=:lightcyan3)#
    mesh2 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor)#

    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)
    rd3 = RigidData(r̄g2, m1, i1, r̄l2, mesh2)

    ## Define the full arm; current issue: it is not defined by segments and has no rotational DOF

    if n == 1
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul-1], [rd2]])
    else
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul], [rd1], [rd2 for _ in 1:num_Y_modul2] for _ in 1:n-1])
    end
    # @show rds

    function get_r̄s(rds)# Define each module position in the initial configuration
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][3] = len1_b
        for i in 1:nb-1
            if length(rds[i].r̄l) == np_rl2 # 13 is the number of key node coordinates in the Y-shaped module
                r̄s[i+1][3] += r̄s[i][3] + len2_t + len2_b - delta
            else
                length(rds[i].r̄l) == np_rl1 # 10 is the number of key node coordinates in the actuator module
                r̄s[i+1][3] += r̄s[i][3] + len2_b + len1_t - delta
            end
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        movable = true
        ci = Vector{Int}() # What is this? pres_idx?
        Φi = collect(1:6) # What is this?

        ri = r̄si # Initial global CoM coordinates of each module
        # @show ri
        aps = rdsi.r̄l # Elastic cable connection positions

        theta = deg2rad(90) # Initial deflection of each section; set to zero here if needed
        R = [
            cos(theta) sin(theta) 0;
            -sin(theta) cos(theta) 0;
            0 0 1
        ]

        r̄g = rdsi.r̄g # CoM coordinates
        m = rdsi.m # Mass
        I = rdsi.i # Inertia tensor
        nrp = length(aps) # Number of nodes
        ṙo = zeros(3)
        ω = zeros(3)
        prop = RB.RigidBodyProperty(i, movable, m, SMatrix{3,3}(I), 
            RB.Locus(SVector{3}(r̄g)), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
                [SVector{3}(aps[i]) for i in 1:nrp], 
                [Axes(SVector{3}(0.0,0.0,1.0)) for _ in 1:nrp], 
                zeros(nrp), 
                zeros(nrp)
            )]
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω,)
        nmcs = RB.NCF.NC1P3V(ri, ro,)# Initial natural coordinates?
        # coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        coords = nmcs
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)
    ## Shared coordinates?
    # matrix_sharing_raw = Vector{Matrix{Int}}()
    # for i = 1:n-1
    #     s = zeros(6,nb)
    #         s[1:6,6i] = 1:6
    #         s[1:6,6i+1] = 1:6
    #     push!(matrix_sharing_raw, s)
    # end
    # sharing_matrix = reduce(vcat,matrix_sharing_raw)

    # Springs##########################################################
    ncables = (nb - 1) * (row_outer_spring + row_inter_spring) # This refers to inner and outer springs
    vks1 = zeros((nb - 1) * row_inter_spring, 1)
    vks1 .= ks1
    for i in 1:3:(nb-1)*row_inter_spring
        vks1[i] = ks1     # vks1[i] = ks1
    end

    restlens1 = restlens1 # inner springs
    ks2 = ks2
    restlens2 = restlens2 # outer springs
    spring_dampers = [(i <= (nb - 1) * row_inter_spring) ? RB.DistanceSpringDamper3D(restlens1, vks1[i], cs1) : RB.DistanceSpringDamper3D(restlens2, ks2, cs2) for i in 1:ncables] # First half are inner springs, second half are outer springs

    matrix_cnt_raw = Vector{Matrix{Int}}()
    # Define inner/outer spring connections
    for j in 1:num_Y_modul # connect inner springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
        # end
    end
    # for j in num_Y_modul+3:nb-1 # connect inner springs
    #     s = zeros(3, nb)
    #     s[1, j] = 8
    #     s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
    #     s[2, j] = 9
    #     s[2, j+1] = -12
    #     s[3, j] = 10
    #     s[3, j+1] = -13
    #     push!(matrix_cnt_raw, s)
    # end
    for j in 1:num_Y_modul # connect outer springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        # s[4, j] = 4
        # s[4, j+1] = -4
        push!(matrix_cnt_raw, s)
        # end
    end
    # for j in num_Y_modul+3:nb-1
    #     s = zeros(4, nb)
    #     s[1, j] = 1
    #     s[1, j+1] = -1
    #     s[2, j] = 2
    #     s[2, j+1] = -2
    #     s[3, j] = 3
    #     s[3, j+1] = -3
    #     s[4, j] = 4
    #     s[4, j+1] = -4
    #     push!(matrix_cnt_raw, s)
    # end

    connecting_matrix = RB.change_connecting_format(rigidbodies, reduce(vcat, matrix_cnt_raw))

    # Sliding cable - begin #############################################
    # restlencs = zeros(nb-1) # Each segment has five sliding cable sections, but the first cable section is actually absent
    restlencs = zeros(nb - n)
    # restlenc1 = 63.77 # The middle rotating section differs
    for i in 1:nb-n #nb-1
        # if (mod(i,num_Y_modul+2) == 0) # The fifth section is the rotating-section length; it does not exist in the actual model
        #     restlencs[i] = restlenc1 
        # if (mod(i, num_Y_modul + 1) == 0) && (i != nb - n)
        # restlencs[i] = restlenc2
        # else
        restlencs[i] = restlenc
        # end
    end

    nclusters_per_seg = 4 * n

    cluster_segs = [
        [RB.DistanceSpringDamperSegment(restlencs[i], kc, c=c_c, prestress=pre) for i in 1:nb-n]
        for _ = 1:nclusters_per_seg
    ]

    cluster_sps = [
        [RB.SlidingPoint(μ) for _ in 1:num_Y_modul-1]
        for _ = 1:nclusters_per_seg
    ]


    cluster_matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:nclusters_per_seg]
    
    s = zeros(Int, nb - 1, 4, nclusters_per_seg) # segments; 4[bpbp]; cable count
    for k in 1:nclusters_per_seg
        for l = 1:nb-1
            s[l, :, k] = [l, 3 + k, l + 1, 3 + k]
        end
    end
    for k in 1:nclusters_per_seg
        push!(cluster_matrix_cnt_raw[k], s[:, :, k])
    end
    # Sliding cable - end ###############################################
    
    connecting_cluster_matrix = [reduce(vcat, cluster_matrix_cnt_raw[k]) for k in 1:nclusters_per_seg]
    cables, clusters = RB.connect_spring_and_frictionlessclusters(
        rigidbodies, spring_dampers,
        cluster_sps, cluster_segs,
        connecting_matrix, connecting_cluster_matrix,
        istart=nclusters_per_seg
    )

    # cst2 = RB.FixedBodyApparatus(77, rbs[14])
    # bjs = [
    # 	RB.PinJoint(RB.End2End(i,RB.ID(rbs[i],2),RB.ID(rbs[i+1],1)))
    # 	for i = 1:nbodies-1
    # ]

    # Fix_joint = RB.FixedJoint(RB.End2End(1,RB.ID(rbs[6],13),RB.ID(rbs[7],13)))
    # joint = RB.join( Fix_joint, indexedcoords)

    # Rev_joint = RB.RevoluteJoint(RB.End2End(1,RB.ID(rbs[6],13,1),RB.ID(rbs[7],13,1)))
    cst1 = RB.FixedBodyApparatus(nclusters_per_seg + ncables + 1, rbs[1])# front
    apparatuses = TypeSortedCollection(vcat([cst1], cables, clusters))
    cnt = RB.Connectivity(rigidbodies, apparatuses;) #？

    # indexed = RB.index(rigidbodies, apparatuses;) #？
    # numbered = RB.number(rigidbodies, apparatuses)

    # cnt = RB.Connectivity(numbered, indexed)
    structure = RB.Structure(rigidbodies, apparatuses, cnt)
    gauges = Int[] # Gauges/sensors used to measure errors for cost functions, optimization, or feedback control; currently unused.
    actuators = reduce(vcat, [ # Actuators attached to one or more bodies or apparatuses; their types dispatch to execute! for different driving forces
        # Try whether ExternalForceActuator is sufficient; if not, create a new actuator
        [RB.ExternalForceActuator(
            i,
            cluster,
            # Operator, 1 control input
            RB.NaiveOperator(1),
            # External force, may not be needed
            0.0,
            # Control input values
            [0.0],
        )
         for (i, cluster) in enumerate(clusters)],
        RB.ExternalForceActuator(
            nclusters_per_seg + 1,
            RB.Signifier(rbs[end], np_rl2),
            # Operator, 1 control input
            RB.NaiveOperator(1),
            # External force, may not be needed
            SVector(0.0, -1.0, 0.0),
            # Control input values
            [0.0],
        )]
    )
    hub = RB.ControlHub(
        structure,
        gauges,
        actuators,
        RB.Coalition(structure, gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)

end

function build_2_seg(n, restlenc, rot_ang)# delta is the initial overlap distance between two units: 69.613
    ks1 = 1.195; restlens1 = 25.0; cs1 = 1.0
    ks2 = 0.037; restlens2 = 55.0; cs2 = 1.0
    @show kc = 10.0
    @show c_c = 0.0
    pre = 3.0
    μ = 0.00001
    deltas = restlenc.-76.0

    num_Y_modul = 3
    num_Y_modul2 = 3
    row_outer_spring = 3
    row_inter_spring = 3
    nb = num_Y_modul  + n+ num_Y_modul2
    # Cable-driven module: CoM-to-bottom distance; CoM-to-top distance
    len1_b = -30.54;  len1_t = -49.96
    # r1 outer spring, r2 cable hole, r4 inner spring
    r1 = 37.3;  r2 = 36.5;  r4 = 31
    # One Y-shaped part: CoM-to-top distance; CoM-to-bottom distance; bottom spring hole; r42 inner spring
    len2_t = -15.57; len2_b = -60.64;  r3 = 5.0;  r42 = 30

    # CoM coordinates
    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]

    # Mass properties
    m1 = 0.12
    m2 = 0.027

    # Inertia tensor
    i1 = [123.78 0 0; 0 113.959 0; 0 0 82.522]## base first section
    i2 = [19.333 0 0; 0 19.343 0; 0 0 15.492] #.* 0.1## single module
    # Key connection points
    ## Cable-driven module
    r̄l1 = [
        # Top outer spring holes
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len1_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len1_t],
        [-r1, 0.0, len1_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        # Top inner spring holes
        [r4 * cos(pi / 3), r4 * sin(pi / 3), len1_t],
        [r4 * cos(pi / 3), -r4 * sin(pi / 3), len1_t],
        [-r4, 0.0, len1_t],
        # Rotation-axis points
        [0.0, 0.0, -len1_b],
        [0.0, 0.0, len1_t-deltas[2]]
    ]
    # Y-shaped part
    r̄l2 = [ # Top outer springs
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len2_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len2_t],
        [-r1, 0.0, len2_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        # Top inner spring holes
        [r42 * cos(pi / 3), r42 * sin(pi / 3), len2_t],
        [r42 * cos(pi / 3), -r42 * sin(pi / 3), len2_t],
        [-r42, 0.0, len2_t],
        # Bottom inner springs
        [r3 * cos(pi / 3), r3 * sin(pi / 3), -len2_b],
        [r3 * cos(pi / 3), -r3 * sin(pi / 3), -len2_b],
        [-r3, 0.0, -len2_b],
        # Rotation-axis points
        [0.0,0.0,len2_t],
        [0.0, 0.0, -len2_b]
    ]

    mesh1 = load(RB.assetpath("octopus/actuator-v2.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi))
    mesh2 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi))

    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)

    ## Define the full arm; current issue: it is not defined by segments and has no rotational DOF

    if n == 1
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul]])
    else
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul], [rd1], [rd2 for _ in 1:num_Y_modul2]])
    end

    function get_r̄s(rds)# Define each module position in the initial configuration
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][3] = len1_b
        for i in 1:nb-1
            i <= num_Y_modul+1 ? delta=deltas[1] : delta=deltas[2]
            if length(rds[i].r̄l) == 15 && i == num_Y_modul+1 # If the previous section is Y-shaped and is the end of a segment, then...
                r̄s[i+1][3] += r̄s[i][3] + len1_b + len2_t
            elseif length(rds[i].r̄l) == 15 && i != num_Y_modul + 1# 13 is the number of key node coordinates in the Y-shaped module
                r̄s[i+1][3] += r̄s[i][3] + len2_t + len2_b - delta
            else
                length(rds[i].r̄l) == 12 # 11 is the number of key node coordinates in the actuator module
                r̄s[i+1][3] += r̄s[i][3] + len2_b + len1_t - delta
                
            end
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        movable = true
        ci = Vector{Int}()
        Φi = collect(1:6)

        ri = r̄si # Initial CoM x-axis coordinates of each module
        # @show ri
        aps = rdsi.r̄l # Connection point positions

        if i >= num_Y_modul+2
            theta = deg2rad(rot_ang[2]) # Initial deflection of each section; set to zero here if needed
        else
            theta = deg2rad(rot_ang[1]) # Initial deflection of each section; set to zero here if needed
        end
        R = [
            cos(theta) sin(theta) 0 ;
            -sin(theta) cos(theta) 0 ;
            0 0 1
        ]

        r̄g = rdsi.r̄g # CoM coordinates
        m = rdsi.m # Mass
        I = rdsi.i # Inertia tensor
        nrp = length(aps) # Number of nodes
        ṙo = zeros(3)
        ω = zeros(3)
        prop = RB.RigidBodyProperty(i, movable, m, SMatrix{3,3}(I), 
            RB.Locus(SVector{3}(r̄g)), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
                [SVector{3}(aps[i]) for i in 1:nrp], 
                [SVector{3}([0.0, 0.0, 1.0]) for i in 1:nrp], #normals, 
                zeros(nrp), 
                zeros(nrp)
            )]
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω,)

        nmcs = RB.NCF.NC1P3V(ri, ro,)# Initial natural coordinates?
        # coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        coords = nmcs
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)

    #---Springs----##########################################################
    ncables = (num_Y_modul + num_Y_modul2) * (row_outer_spring + row_inter_spring) # This refers to inner and outer springs
    ks1 = ks1
    restlens1 = restlens1 # inner springs
    ks2 = ks2
    restlens2 = restlens2 # outer springs
    spring_dampers = [(i <= (num_Y_modul + num_Y_modul2) * row_inter_spring) ? RB.DistanceSpringDamper3D(restlens1, ks1, cs1) : RB.DistanceSpringDamper3D(restlens2, ks2, cs2) for i in 1:ncables] # First half are inner springs, second half are outer springs

    matrix_cnt_raw = Vector{Matrix{Int}}()
    # Define inner/outer spring connections
    for j in 1:num_Y_modul # connect inner springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
        # end
    end
    for j in num_Y_modul+2:nb-1 # connect inner springs
        s = zeros(3, nb)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
    end
    for j in 1:num_Y_modul # connect outer springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        push!(matrix_cnt_raw, s)
        # end
    end
    for j in num_Y_modul+2:nb-1
        s = zeros(3, nb)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        push!(matrix_cnt_raw, s)
    end
    connecting_matrix = RT.change_connecting_format(rigidbodies, reduce(vcat, matrix_cnt_raw))
    # Sliding cable - begin #############################################

    restlencs = zeros(nb - n)
    for i in 1:num_Y_modul #nb-1
        restlencs[i] = restlenc[1]
    end
    for i in num_Y_modul+1:nb-n
        restlencs[i] = restlenc[2]
    end
    # @show restlencs
    nclusters_per_seg = 4
    cluster_segs = vcat(
        [[RT.DistanceSpringDamperSegment(restlencs[i], kc, c=c_c, prestress=pre) for i in 1:num_Y_modul]
        for _ = 1:nclusters_per_seg],
        [[RT.DistanceSpringDamperSegment(restlencs[num_Y_modul+i], kc, c=c_c, prestress=pre) for i in 1:num_Y_modul2]
        for _ = 1:nclusters_per_seg]
    )

    cluster_sps = vcat(
        [[RT.SlidingPoint(μ) for _ in 1:num_Y_modul-1]
        for _ = 1:nclusters_per_seg],
        [[RT.SlidingPoint(μ) for _ in 1:num_Y_modul2-1]
        for _ = 1:nclusters_per_seg]
    )
    
    cluster_matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:nclusters_per_seg*n]
    s1 = zeros(Int, num_Y_modul, 4, nclusters_per_seg)
    for k in 1:nclusters_per_seg
        for l = 1:num_Y_modul
            s1[l, :, k] = [l, 3 + k, l + 1, 3 + k]
        end
    end
    for k in 1:nclusters_per_seg
        push!(cluster_matrix_cnt_raw[k], s1[:, :, k])
    end
    
    s2 = zeros(Int, num_Y_modul2, 4, nclusters_per_seg)
    for k in 1:nclusters_per_seg
        for l = 1:num_Y_modul2
            s2[l, :, k] = [l + num_Y_modul + 1, 3 + k, l + num_Y_modul + 2, 3 + k]
        end
    end
    for k in nclusters_per_seg+1:nclusters_per_seg*2
        push!(cluster_matrix_cnt_raw[k], s2[:, :, k-nclusters_per_seg])
    end

    # Sliding cable - end ###############################################
    
    connecting_cluster_matrix = [reduce(vcat, cluster_matrix_cnt_raw[k]) for k in 1:nclusters_per_seg*n]
    cables, clusters = RT.connect_spring_and_clusters(
        rigidbodies, spring_dampers,
        cluster_sps, cluster_segs,
        connecting_matrix, connecting_cluster_matrix,
        istart=nclusters_per_seg*2
    )
    cst1 = RB.FixedBodyApparatus(nclusters_per_seg*2 + ncables + 1, rbs[1])

    rj = RB.proto_joint_apparatus(
        nclusters_per_seg*2 + ncables + 2, 
        RB.Hen2Egg(
            RB.Anchor(rbs[num_Y_modul+1], 14, 1), 
            RB.Anchor(rbs[num_Y_modul+2], 11, 1)
        ),
        RB.NoForce(),
        :Revolute;
    )
    rjc = RB.proto_joint_apparatus(
        nclusters_per_seg*2 + ncables + 3,
        RB.Hen2Egg(
            RB.Anchor(rbs[num_Y_modul+1], 14, 1), 
            RB.Anchor(rbs[num_Y_modul+2], 11, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :FloatingUniversal;
    )
    # rj = RB.FixedJoint(nclusters_per_seg * 2 + ncables + 2, RB.Hen2Egg(RB.Anchor(rbs[num_Y_modul+1], 14, 1), RB.Anchor(rbs[num_Y_modul+2], 11, 1)))
    # rj = RB.FixedJoint(nclusters_per_seg * 2 + ncables + 2, RB.Hen2Egg(RB.Anchor(rbs[num_Y_modul+1], 14, 1), RB.Anchor(rbs[num_Y_modul+2], 11, 1)))
    # rj1 = RB.SphericalJoint(nclusters_per_seg * 2 + ncables + 3, RB.Hen2Egg(RB.Anchor(rbs[num_Y_modul+2], 12, 1), RB.Anchor(rbs[num_Y_modul+3], 15, 1)))
    apparatuses = TypeSortedCollection(vcat([cst1,rj,rjc], cables, clusters)) # middle revolute joint
    # apparatuses = TypeSortedCollection(vcat([cst1, rj], cables, clusters)) # use this when rigidly connecting the middle joint

    cnt = RB.Connectivity(rigidbodies, apparatuses;)
    structure = RB.Structure(rigidbodies, apparatuses, cnt)
    capta_gauges = Int[]
    error_gauges = Int[] # Gauges/sensors used to measure errors for cost functions, optimization, or feedback control; currently unused.
    actuators = reduce(vcat,[ # Actuators attached to one or more bodies or apparatuses; their types dispatch to execute! for different driving forces
        # Try whether ExternalForceActuator is sufficient; if not, create a new actuator
        
        [
            RB.ExternalForceActuator(
                i,
                cluster,
                # Operator, 1 control input
                RB.NaiveOperator(1),
                # External force, may not be needed
                0.0,
                # Control input values
                [0.0],
            )
            for (i, cluster) in enumerate(clusters)
        ],
        RB.RegisterActuator(
            n*nclusters_per_seg + 1,
            rjc, # Create apparatus collection
            # Operator, 1 control input
            RB.FunctionOperator(
                (f, values, u) -> begin
                    @inbounds @. f = values + u
                end,
                (jac, values, u) -> nothing,
                copy([0.0]), # func_vals
                Matrix{eltype([0.0])}(I(length([0.0]))),  # Jac_vals
            ), # operator
            (values = [0.0],)
        ) # Use this to drive the revolute joint
        # [RB.ExternalForceActuator(
        #     n*nclusters_per_seg+1,
        #     RB.Anchor(rbs[num_Y_modul+2], 11),
        #     # Operator, 1 control input
        #     RB.NaiveOperator(1),
        #     # External force, may not be needed
        #     [1.0,1.0,0.0],
        #     # Control input values
        #     [0.0],
        # )]
    ])

    hub = RB.ControlHub(
        structure,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(structure, capta_gauges, error_gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)

end

function build_cone_man(n, ks1, restlenc, delta, bodycolor=:lightcyan3)# delta is the initial overlap distance between two units: 69.613

    scale = [1.0,0.9,0.8,0.7,0.7,0.6,0.6,0.5,0.5]
    restlens1 = 16.0
    cs1 = 2.0 #1195
    ks2 = 0.037
    restlens2 = 40.0
    cs2 = 2.0
    kc = 1.0
    @show c_c = 0.0
    pre = 3.0
    @show μ = 0.00001

    num_Y_modul = 8
    row_outer_spring = 3
    row_inter_spring = 3
    nb = (num_Y_modul + 1) * n
    # Cable-driven module: CoM-to-bottom distance; CoM-to-top distance
    len1_b = -30.54
    len1_t = -49.96
    # r1 outer spring, r2 cable hole, r4 inner spring
    r1 = 37.3
    r2 = 36.5
    r4 = 31
    # One Y-shaped part: CoM-to-top distance; CoM-to-bottom distance; bottom spring hole; r42 inner spring
    len2_t = -15.57
    len2_b = -60.64
    r3 = 5.0
    r42 = 30

    # CoM coordinates
    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]
    # r̄g3 = [0.0, 0.0, 0.0]
    # r̄g4 = [0.0, 0.0, 0.0]
    # Mass properties
    m1 = 0.12
    @show m2 = 0.027
    m3 = 0.1 #m4 = 0.15
    # Inertia tensor
    i1 = [4575.03 0 0; 0 5114.84 0; 0 0 4575.03]## base first section
    i2 = [19.333 0 0; 0 19.343 0; 0 0 15.492]# .* 0.01## single module

    # Key connection points
    ## Cable-driven module
    r̄l1 = [
        # Top outer spring holes
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len1_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len1_t],
        [-r1, 0.0, len1_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        # Top inner spring holes
        [r4 * cos(pi / 3), r4 * sin(pi / 3), len1_t],
        [r4 * cos(pi / 3), -r4 * sin(pi / 3), len1_t],
        [-r4, 0.0, len1_t],
    ]
    # Y-shaped part
    r̄l2 = [ # Top outer springs
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len2_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len2_t],
        [-r1, 0.0, len2_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        # Top inner spring holes
        [r42 * cos(pi / 3), r42 * sin(pi / 3), len2_t],
        [r42 * cos(pi / 3), -r42 * sin(pi / 3), len2_t],
        [-r42, 0.0, len2_t],
        # Bottom inner springs
        [r3 * cos(pi / 3), r3 * sin(pi / 3), -len2_b],
        [r3 * cos(pi / 3), -r3 * sin(pi / 3), -len2_b],
        [-r3, 0.0, -len2_b],
    ]

    mesh1 = load(RB.assetpath("octopus/actuator-v2.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=:lightcyan3)#
    mesh2 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor)#
    mesh21 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor,scale=scale[2])#
    mesh22 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor, scale=scale[3])#
    mesh23 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor, scale=scale[4])#
    mesh24 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor, scale=scale[5])#
    mesh25 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor, scale=scale[6])#

    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)
    rd3 = RigidData(r̄g2, m2, i2, r̄l2 .* scale[2], mesh21)
    rd4 = RigidData(r̄g2, m2, i2, r̄l2 .* scale[3], mesh22)
    rd5 = RigidData(r̄g2, m2, i2, r̄l2 .* scale[4], mesh23)
    rd6 = RigidData(r̄g2, m2, i2, r̄l2 .* scale[5], mesh24)
    rd7 = RigidData(r̄g2, m2, i2, r̄l2 .* scale[6], mesh25)

    ## Define the full arm; current issue: it is not defined by segments and has no rotational DOF

    if n == 1
        # rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul]])
        rds = reduce(vcat, [[rd1], [rd2], [rd3], [rd4], [rd5], [rd6,rd6], [rd7,rd7]])
    else
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul], [rd1], [rd2 for _ in 1:num_Y_modul2] for _ in 1:n-1])
    end

    # function get_r̄s(rds)# Define each module position in the initial configuration
    #     r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
    #     r̄s[1][3] = len1_b
    #     for i in 1:nb-1
    #         if length(rds[i].r̄l) == 13 # 13 is the number of key node coordinates in the Y-shaped module
    #             r̄s[i+1][3] += r̄s[i][3] + len2_t + len2_b - delta
    #         else
    #             length(rds[i].r̄l) == 10 # 10 is the number of key node coordinates in the actuator module
    #             r̄s[i+1][3] += r̄s[i][3] + len2_b + len1_t - delta
    #         end
    #     end
    #     return r̄s
    # end
    function get_r̄s(rds)# Define each module position in the initial configuration
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][3] = len1_b
        r̄s[2][3] += r̄s[1][3] + len2_b * scale[2] + len1_t #- delta
        for i in 1:nb-2
            r̄s[i+2][3] += r̄s[i+1][3] + len2_t * scale[i+1] + len2_b * scale[i+2] #- delta * scale[i+2]
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        movable = true
        ci = Vector{Int}()
        Φi = collect(1:6)

        ri = r̄si # Initial CoM x-axis coordinates of each module
        # @show ri
        aps = rdsi.r̄l # Elastic cable connection positions

        theta = deg2rad(90) # Initial deflection of each section; set to zero here if needed
        R = [
            cos(theta) sin(theta) 0;
            -sin(theta) cos(theta) 0;
            0 0 1
        ]

        r̄g = rdsi.r̄g # CoM coordinates
        m = rdsi.m # Mass
        I = rdsi.i # Inertia tensor
        nrp = length(aps) # Number of nodes
        ṙo = zeros(3)
        ω = zeros(3)
        prop = RB.RigidBodyProperty(i, movable, m, SMatrix{3,3}(I), 
            RB.Locus(SVector{3}(r̄g)), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
                [SVector{3}(aps[i]) for i in 1:nrp], 
                [Axes(SVector{3}(0.0,0.0,1.0)) for _ in 1:nrp], 
                zeros(nrp), 
                zeros(nrp)
            )]
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω,)

        nmcs = RB.NCF.NC1P3V(ri, ro,)# Initial natural coordinates?
        # coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        coords = nmcs
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)
    ## Shared coordinates?
    # matrix_sharing_raw = Vector{Matrix{Int}}()
    # for i = 1:n-1
    #     s = zeros(6,nb)
    #         s[1:6,6i] = 1:6
    #         s[1:6,6i+1] = 1:6
    #     push!(matrix_sharing_raw, s)
    # end
    # sharing_matrix = reduce(vcat,matrix_sharing_raw)

    # Springs##########################################################
    ncables = (nb - 1) * (row_outer_spring + row_inter_spring) # This refers to inner and outer springs
    # vks1 = zeros((nb - 1) * row_inter_spring, 1)
    # vks1 .= ks1
    # for i in 1:3:(nb-1)*row_inter_spring
    #     vks1[i] = ks1 / 2
    # end
    # ks1 = ks1
    # restlens1 = restlens1 # inner springs
    # ks2 = ks2
    # restlens2 = restlens2 # outer springs
    vks1 = zeros((nb - 1) * row_inter_spring, 1)
    vrestlens1 = zeros((nb - 1) * row_inter_spring, 1)
    vks2 = zeros((nb - 1) * row_outer_spring, 1)
    vrestlens2 = zeros((nb - 1) * row_outer_spring, 1)

    for i in 1:nb-1
        # for j in 1:3
            vks1[3i-2:3i] .= ks1
            vrestlens1[3i-2:3i] .= restlens1*scale[i]
            vks2[3i-2:3i] .= ks2 
            vrestlens2[3i-2:3i] .= restlens2 * scale[i]
        # end
    end

    spring_dampers = [(i <= (nb - 1) * row_inter_spring) ? RB.DistanceSpringDamper3D(restlens1, ks1, cs1) : RB.DistanceSpringDamper3D(restlens2, ks2, cs2) for i in 1:ncables] # First half are inner springs, second half are outer springs

    matrix_cnt_raw = Vector{Matrix{Int}}()
    # Define inner/outer spring connections
    for j in 1:num_Y_modul # connect inner springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
        # end
    end
    # for j in num_Y_modul+3:nb-1 # connect inner springs
    #     s = zeros(3, nb)
    #     s[1, j] = 8
    #     s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
    #     s[2, j] = 9
    #     s[2, j+1] = -12
    #     s[3, j] = 10
    #     s[3, j+1] = -13
    #     push!(matrix_cnt_raw, s)
    # end
    for j in 1:num_Y_modul # connect outer springs
        # for j in 4i-3:4i
        s = zeros(3, nb)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        # s[4, j] = 4
        # s[4, j+1] = -4
        push!(matrix_cnt_raw, s)
        # end
    end
    # for j in num_Y_modul+3:nb-1
    #     s = zeros(4, nb)
    #     s[1, j] = 1
    #     s[1, j+1] = -1
    #     s[2, j] = 2
    #     s[2, j+1] = -2
    #     s[3, j] = 3
    #     s[3, j+1] = -3
    #     s[4, j] = 4
    #     s[4, j+1] = -4
    #     push!(matrix_cnt_raw, s)
    # end
    connecting_matrix = RB.change_connecting_format(rigidbodies, reduce(vcat, matrix_cnt_raw))
    # Sliding cable - begin #############################################
    # restlencs = zeros(nb-1) # Each segment has five sliding cable sections, but the first cable section is actually absent
    restlencs = zeros(nb - n)
    # restlenc1 = 63.77 # The middle rotating section differs
    for i in 1:nb-n #nb-1
    # if (mod(i,num_Y_modul+2) == 0) # The fifth section is the rotating-section length; it does not exist in the actual model
        #     restlencs[i] = restlenc1 
        # if (mod(i, num_Y_modul + 1) == 0) && (i != nb - n)
        # restlencs[i] = restlenc2
        # else
        restlencs[i] = restlenc*scale[i]
        # end
    end

    nclusters_per_seg = 4 * n

    cluster_segs = [
        [RB.DistanceSpringDamperSegment(restlencs[i], kc, c=c_c, prestress=pre) for i in 1:nb-n]
        for _ = 1:nclusters_per_seg
    ]

    cluster_sps = [
        [RB.SlidingPoint(μ) for _ in 1:num_Y_modul-1]
        for _ = 1:nclusters_per_seg
    ]


    cluster_matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:nclusters_per_seg]
    @show nb
    s = zeros(Int, nb - 1, 4, nclusters_per_seg)
    for k in 1:nclusters_per_seg
        for l = 1:nb-1
            s[l, :, k] = [l, 3 + k, l + 1, 3 + k]
        end
    end
    for k in 1:nclusters_per_seg
        push!(cluster_matrix_cnt_raw[k], s[:, :, k])
    end
    # Sliding cable - end ###############################################
    display(cluster_matrix_cnt_raw[1])
    connecting_cluster_matrix = [reduce(vcat, cluster_matrix_cnt_raw[k]) for k in 1:nclusters_per_seg]
    cables, clusters = RB.connect_spring_and_clusters(
        rigidbodies, spring_dampers,
        cluster_sps, cluster_segs,
        connecting_matrix, connecting_cluster_matrix,
        istart=nclusters_per_seg
    )

    # cst2 = RB.FixedBodyApparatus(77, rbs[14])
    # bjs = [
    # 	RB.PinJoint(RB.End2End(i,RB.ID(rbs[i],2),RB.ID(rbs[i+1],1)))
    # 	for i = 1:nbodies-1
    # ]

    # Fix_joint = RB.FixedJoint(RB.End2End(1,RB.ID(rbs[6],13),RB.ID(rbs[7],13)))
    # joint = RB.join( Fix_joint, indexedcoords)

    # Rev_joint = RB.RevoluteJoint(RB.End2End(1,RB.ID(rbs[6],13,1),RB.ID(rbs[7],13,1)))
    cst1 = RB.FixedBodyApparatus(nclusters_per_seg + ncables + 1, rbs[1])# front
    apparatuses = TypeSortedCollection(vcat([cst1], cables, clusters))
    indexed = RB.index(rigidbodies, apparatuses;) #？
    numbered = RB.number(rigidbodies, apparatuses)
    cnt = RB.Connectivity(numbered, indexed)
    structure = RB.Structure(rigidbodies, apparatuses, cnt)
    gauges = Int[] # Gauges/sensors used to measure errors for cost functions, optimization, or feedback control; currently unused.
    actuators = [ # Actuators attached to one or more bodies or apparatuses; their types dispatch to execute! for different driving forces
        # Try whether ExternalForceActuator is sufficient; if not, create a new actuator
        RB.ExternalForceActuator(
            i,
            cluster,
            # Operator, 1 control input
            RB.NaiveOperator(1),
            # External force, may not be needed
            0.0,
            # Control input values
            [0.0],
        )
        for (i, cluster) in enumerate(clusters)
    ]
    hub = RB.ControlHub(
        structure,
        gauges,
        actuators,
        RB.Coalition(structure, gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)

end

function build_1_man_Prismatic(n, ks1, restlenc, delta, bodycolor=:lightcyan3)# delta is the initial overlap distance between two units: 69.613
    restlens1 = 25.0
    cs1 = 2.0 #1195
    ks2 = 0.037
    restlens2 = 55.0
    cs2 = 2.0
    kc = 10.0
    c_c = 0.0
    pre = 3.0
    μ = 0.0001
    prelenth = 60.0

    num_Y_modul = 6
    row_outer_spring = 3
    row_inter_spring = 3
    nb = (num_Y_modul + 1) * n +2
    # Cable-driven module: CoM-to-bottom distance; CoM-to-top distance
    len1_b = -30.54
    len1_t = -49.96
    # r1 outer spring, r2 cable hole, r4 inner spring
    r1 = 37.3
    r2 = 36.5
    r4 = 31
    # One Y-shaped part: CoM-to-top distance; CoM-to-bottom distance; bottom spring hole; r42 inner spring
    len2_t = -15.57
    len2_b = -60.64
    r3 = 5.0
    r42 = 30

    # CoM coordinates
    r̄g0 = [0.0, 0.0, 0.0]
    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]
    # r̄g3 = [0.0, 0.0, 0.0]
    # r̄g4 = [0.0, 0.0, 0.0]
    # Mass properties
    m0 = 0.0001
    m1 = 0.12
    m2 = 0.027
    m3 = 0.1 #m4 = 0.15
    # Inertia tensor
    i0 = [1.68e6 0 0; 0 1.68e6 0; 0 0 1.68e6]#
    i1 = [4575.03 0 0; 0 5114.84 0; 0 0 4575.03]## base first section
    i2 = [19.333 0 0; 0 19.343 0; 0 0 15.492]# .* 0.01## single module

    # Key connection points
    r̄l0 = [
        [0.0, 0., 0.],
        [0., r2 * sqrt(2) / 2,  0.],
        [0.,-r2 * sqrt(2) / 2,  0.]
        ]
    ## Cable-driven module
    r̄l1 = [
        # Top outer spring holes
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len1_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len1_t],
        [-r1, 0.0, len1_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        # Top inner spring holes
        [r4 * cos(pi / 3), r4 * sin(pi / 3), len1_t],
        [r4 * cos(pi / 3), -r4 * sin(pi / 3), len1_t],
        [-r4, 0.0, len1_t],
        [r2 * sqrt(2) / 2, 0., -len1_b + prelenth],
        [-r2 * sqrt(2) / 2, 0., -len1_b + prelenth],
    ]
    # Y-shaped part
    r̄l2 = [ # Top outer springs
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len2_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len2_t],
        [-r1, 0.0, len2_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        # Top inner spring holes
        [r42 * cos(pi / 3), r42 * sin(pi / 3), len2_t],
        [r42 * cos(pi / 3), -r42 * sin(pi / 3), len2_t],
        [-r42, 0.0, len2_t],
        # Bottom inner springs
        [r3 * cos(pi / 3), r3 * sin(pi / 3), -len2_b],
        [r3 * cos(pi / 3), -r3 * sin(pi / 3), -len2_b],
        [-r3, 0.0, -len2_b],
        [0.0, 0.0, len2_t]
    ]
    np_rl0 = length(r̄l0)
    np_rl1 = length(r̄l1)
    np_rl2 = length(r̄l2)

    mesh0 = load(RB.assetpath("octopus/CTC_endpoint.STL"))#
    mesh1 = load(RB.assetpath("octopus/actuator-v2.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=:lightcyan3)#
    mesh2 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi), color=bodycolor)#

    rd0 = RigidData(r̄g0, m0, i0, r̄l0, mesh0)
    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)
    rd3 = RigidData(r̄g2, m1, i1, r̄l2, mesh2)

    ## Define the full arm; current issue: it is not defined by segments and has no rotational DOF

    if n == 1
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul-1], [rd2],[rd0,rd0]])
    else
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul], [rd1], [rd2 for _ in 1:num_Y_modul2] for _ in 1:n-1])
    end
    # @show rds

    function get_r̄s(rds)# Define each module position in the initial configuration
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][3] = len1_b
        r̄s[end-1][2] = -r2 * sqrt(2) / 2
        r̄s[end-1][3] = prelenth
        r̄s[end][2] = r2 * sqrt(2) / 2
        r̄s[end][3] = prelenth
        
        for i in 1:nb-3 
            if length(rds[i].r̄l) == np_rl2 # 13 is the number of key node coordinates in the Y-shaped module
                r̄s[i+1][3] += r̄s[i][3] + len2_t + len2_b - delta
            else
                length(rds[i].r̄l) == np_rl1 # 10 is the number of key node coordinates in the actuator module
                r̄s[i+1][3] += r̄s[i][3] + len2_b + len1_t - delta
            end
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        movable = true
        ci = Vector{Int}() # What is this? pres_idx?
        Φi = collect(1:6) # What is this?

        ri = r̄si # Initial global CoM coordinates of each module
        # @show ri
        aps = rdsi.r̄l # Elastic cable connection positions

        theta = deg2rad(90) # Initial deflection of each section; set to zero here if needed
        R = [
            cos(theta) sin(theta) 0;
            -sin(theta) cos(theta) 0;
            0 0 1
        ]

        r̄g = rdsi.r̄g # CoM coordinates
        m = rdsi.m # Mass
        I = rdsi.i # Inertia tensor
        nrp = length(aps) # Number of nodes
        ṙo = zeros(3)
        ω = zeros(3)
        prop = RB.RigidBodyProperty(i, movable, m, SMatrix{3,3}(I), 
            RB.Locus(SVector{3}(r̄g)), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
                [SVector{3}(aps[i]) for i in 1:nrp], 
                [SVector{3}([0.0, 0.0, 1.0]) for i in 1:nrp], #normals, 
                zeros(nrp), 
                zeros(nrp)
            )]
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω,)
        nmcs = RB.NCF.NC1P3V(ri, ro,)# Initial natural coordinates?
        coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)

    # Springs##########################################################
    ncables = (nb - 3) * (row_outer_spring + row_inter_spring) # This refers to inner and outer springs
    vks1 = zeros((nb - 3) * row_inter_spring, 1)
    vks1 .= ks1
    for i in 1:3:(nb-3)*row_inter_spring
        vks1[i] = ks1     # vks1[i] = ks1
    end

    restlens1 = restlens1 # inner springs
    ks2 = ks2
    restlens2 = restlens2 # outer springs
    spring_dampers = [(i <= (nb - 3) * row_inter_spring) ? RB.DistanceSpringDamper3D(restlens1, vks1[i], cs1) : RB.DistanceSpringDamper3D(restlens2, ks2, cs2) for i in 1:ncables] # First half are inner springs, second half are outer springs

    matrix_cnt_raw = Vector{Matrix{Int}}()
    # Define inner/outer spring connections
    for j in 1:num_Y_modul # connect inner springs
        # for j in 4i-3:4i
        s = zeros(3, nb-2)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
        # end
    end
    
    for j in 1:num_Y_modul # connect outer springs
        # for j in 4i-3:4i
        s = zeros(3, nb-2)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        # s[4, j] = 4
        # s[4, j+1] = -4
        push!(matrix_cnt_raw, s)
        # end
    end

    connecting_matrix = RB.change_connecting_format(rigidbodies, reduce(vcat, matrix_cnt_raw))

    # Sliding cable - begin #############################################
    # restlencs = zeros(nb-1) # Each segment has five sliding cable sections, but the first cable section is actually absent
    restlencs = zeros(nb - n - 1)
    # restlenc1 = 63.77 # The middle rotating section differs
    for i in 1:nb-n-1 #nb-1
        if i ==1
            restlencs[i] = prelenth - len1_b - len1_t - 0.5
        else
            restlencs[i] = restlenc
        end
    end

    nclusters_per_seg = 4 * n

    cluster_segs = [
        [RB.DistanceSpringDamperSegment(restlencs[i], kc, c=c_c, prestress=pre) for i in 1:nb-n-1]
        for _ = 1:nclusters_per_seg
    ]

    cluster_sps = [
        [RB.SlidingPoint(μ) for _ in 1:num_Y_modul]
        for _ = 1:nclusters_per_seg
    ]


    cluster_matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:nclusters_per_seg]
    # @show nb
    s = zeros(Int, nb - 2, 4, nclusters_per_seg) # segments; 4[bpbp]; cable count
    for k in 1:nclusters_per_seg # k-th CTC
        if k == 1
            s[1, :, k] = [nb-1, 2, 1, 3 + k] # first CTC segment
        elseif k==2
            s[1, :, k] = [nb-1, 3, 1, 3 + k] # first CTC segment
        elseif k==3
            s[1, :, k] = [nb, 3, 1, 3 + k] # first CTC segment
        elseif k==4
            s[1, :, k] = [nb, 2, 1, 3 + k] # first CTC segment
        end
        for l = 2:nb-2 # l-th CTC segment
            s[l, :, k] = [l-1, 3 + k, l, 3 + k]
        end
    end
    for k in 1:nclusters_per_seg
        push!(cluster_matrix_cnt_raw[k], s[:, :, k])
    end
    # Sliding cable - end ###############################################
    
    connecting_cluster_matrix = [reduce(vcat, cluster_matrix_cnt_raw[k]) for k in 1:nclusters_per_seg]
    cables, clusters = RB.connect_spring_and_frictionlessclusters(#
        rigidbodies, spring_dampers,
        cluster_sps, cluster_segs,
        connecting_matrix, connecting_cluster_matrix,
        istart=nclusters_per_seg
    )

    cst1 = RB.FixedBodyApparatus(nclusters_per_seg + ncables + 1, rbs[1])# front
    # rjctc1 = RB.FixedBodyApparatus(nclusters_per_seg + ncables + 2, rbs[end-1])
    # rjctc2 = RB.FixedBodyApparatus(nclusters_per_seg + ncables + 3, rbs[end])
    PrJoint1 = RB.proto_joint_apparatus(
        nclusters_per_seg + ncables + 2,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb], 1, 1),
            RB.Anchor(rbs[1], 12, 1)
        ),
        RB.NoForce(),
        :Prismatic;
    )
    PrJoint2 = RB.proto_joint_apparatus(
        nclusters_per_seg + ncables + 3,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-1], 1, 1),
            RB.Anchor(rbs[1], 11, 1)
        ),
        RB.NoForce(),
        :Prismatic;
    )
    PSJoint1 = RB.proto_joint_apparatus(
        nclusters_per_seg + ncables + 4,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb], 1, 1),
            RB.Anchor(rbs[1], 12, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :PlanarSpherical;
    )
    PSJoint2 = RB.proto_joint_apparatus(
        nclusters_per_seg + ncables + 5,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-1], 1, 1),
            RB.Anchor(rbs[1], 11, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :PlanarSpherical;
    )
    apparatuses = TypeSortedCollection(vcat([cst1, PrJoint1, PrJoint2, PSJoint1, PSJoint2], cables, clusters))
    cnt = RB.Connectivity(rigidbodies, apparatuses;) #？

    structure = RB.Structure(rigidbodies, apparatuses, cnt)
    gauges = Int[] # Gauges/sensors used to measure errors for cost functions, optimization, or feedback control; currently unused.
    actuators = reduce(vcat, [ # Actuators attached to one or more bodies or apparatuses; their types dispatch to execute! for different driving forces
        # Try whether ExternalForceActuator is sufficient; if not, create a new actuator
        # [RB.ExternalForceActuator(
        #     i,
        #     cluster,
        #     # Operator, 1 control input
        #     RB.NaiveOperator(1),
        #     # External force, may not be needed
        #     0.0,
        #     # Control input values
        #     [0.0],
        # )
        #  for (i, cluster) in enumerate(clusters)],
        [RB.RegisterActuator(
                1,# n * nclusters_per_seg +
                PSJoint1, # Create apparatus collection
                # Operator, 1 control input
                RB.FunctionOperator(
                    (values, u) -> values .+ u # Define operator behavior
                ),
                # Control input values, Reg
                (values=[0.0],)
            ) # Use this to drive the revolute joint
            ], 
        [RB.RegisterActuator(
                2,#n * nclusters_per_seg + 
                PSJoint2, # Create apparatus collection
                # Operator, 1 control input
                RB.FunctionOperator(
                    (values, u) -> values .+ u # Define operator behavior
                ),
                # Control input values, Reg
                (values=[0.0],)
            ) # Use this to drive the revolute joint
            ]
        ]
    )
    hub = RB.ControlHub(
        structure,
        gauges,
        actuators,
        RB.Coalition(structure, gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)

end

function build_2_seg_Prismatic(n, restlenc, rot_ang)# delta is the initial overlap distance between two units: 69.613
    ks1 = 1.195
    restlens1 = 25.0
    cs1 = 1.0
    ks2 = 0.037
    restlens2 = 55.0
    cs2 = 1.0
    @show kc = 10.0
    @show c_c = 0.0
    pre = 3.0
    μ = 0.00001
    deltas = restlenc .- 76.0
    prelenth = 60.0

    num_Y_modul = 4
    num_Y_modul2 = 6
    row_outer_spring = 3
    row_inter_spring = 3
    nb = num_Y_modul + n + num_Y_modul2 +2n
    # Cable-driven module: CoM-to-bottom distance; CoM-to-top distance
    len1_b = -30.54
    len1_t = -49.96
    # r1 outer spring, r2 cable hole, r4 inner spring
    r1 = 37.3
    r2 = 36.5
    r4 = 31
    # One Y-shaped part: CoM-to-top distance; CoM-to-bottom distance; bottom spring hole; r42 inner spring
    len2_t = -15.57
    len2_b = -60.64
    r3 = 5.0
    r42 = 30

    # CoM coordinates
    r̄g0 = [0.0, 0.0, 0.0]
    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]

    # Mass properties
    m0 = 0.0001
    m1 = 0.12
    m2 = 0.027

    # Inertia tensor
    i0 = [1.68e6 0 0; 0 1.68e6 0; 0 0 1.68e6]#
    i1 = [123.78 0 0; 0 113.959 0; 0 0 82.522]## base first section
    i2 = [19.333 0 0; 0 19.343 0; 0 0 15.492] #.* 0.1## single module
    # Key connection points
    r̄l0 = [
        [0.0, 0., 0.],
        [0., r2 * sqrt(2) / 2, 0.],
        [0., -r2 * sqrt(2) / 2, 0.]
    ]
    ## Cable-driven module
    r̄l1 = [
        # Top outer spring holes
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len1_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len1_t],
        [-r1, 0.0, len1_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len1_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len1_t],
        # Top inner spring holes
        [r4 * cos(pi / 3), r4 * sin(pi / 3), len1_t],
        [r4 * cos(pi / 3), -r4 * sin(pi / 3), len1_t],
        [-r4, 0.0, len1_t],
        # Rotation-axis points
        [0.0, 0.0, -len1_b],
        [0.0, 0.0, len1_t - deltas[2]],
        # Drive
        [r2 * sqrt(2) / 2, 0., -len1_b + prelenth],
        [-r2 * sqrt(2) / 2, 0., -len1_b + prelenth]
    ]
    # Y-shaped part
    r̄l2 = [ # Top outer springs
        [r1 * cos(pi / 3), r1 * sin(pi / 3), len2_t],
        [r1 * cos(pi / 3), -r1 * sin(pi / 3), len2_t],
        [-r1, 0.0, len2_t],
        # Top cable holes
        [r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        [r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, -r2 * sqrt(2) / 2, len2_t],
        [-r2 * sqrt(2) / 2, r2 * sqrt(2) / 2, len2_t],
        # Top inner spring holes
        [r42 * cos(pi / 3), r42 * sin(pi / 3), len2_t],
        [r42 * cos(pi / 3), -r42 * sin(pi / 3), len2_t],
        [-r42, 0.0, len2_t],
        # Bottom inner springs
        [r3 * cos(pi / 3), r3 * sin(pi / 3), -len2_b],
        [r3 * cos(pi / 3), -r3 * sin(pi / 3), -len2_b],
        [-r3, 0.0, -len2_b],
        # Rotation-axis points
        [0.0, 0.0, len2_t],
        [0.0, 0.0, -len2_b]
    ]

    np_rl0 = length(r̄l0)
    np_rl1 = length(r̄l1)
    np_rl2 = length(r̄l2)

    mesh0 = load(RB.assetpath("octopus/CTC_endpoint.STL"))#
    mesh1 = load(RB.assetpath("octopus/actuator-v2.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi))
    mesh2 = load(RB.assetpath("octopus/Y-module-v1.STL")) |> RB.make_patch(; rot=RotY(0)) |> RB.make_patch(; rot=RotX(pi))

    rd0 = RigidData(r̄g0, m0, i0, r̄l0, mesh0)
    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)

    ## Define the full arm; current issue: it is not defined by segments and has no rotational DOF

    if n == 1
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul]])
    else
        rds = reduce(vcat, [[rd1], [rd2 for _ in 1:num_Y_modul], [rd1], [rd2 for _ in 1:num_Y_modul2], [rd0, rd0, rd0, rd0]])
    end

    function get_r̄s(rds)# Define each module position in the initial configuration
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][3] = len1_b

        for i in 1:nb-1-2n
            i <= num_Y_modul + 1 ? delta = deltas[1] : delta = deltas[2]
            if length(rds[i].r̄l) == np_rl2 && i == num_Y_modul + 1 # If the previous section is Y-shaped and is the end of a segment, then...
                r̄s[i+1][3] += r̄s[i][3] + len1_b + len2_t
            elseif length(rds[i].r̄l) == np_rl2 && i != num_Y_modul + 1# 13 is the number of key node coordinates in the Y-shaped module
                r̄s[i+1][3] += r̄s[i][3] + len2_t + len2_b - delta
            else
                length(rds[i].r̄l) == np_rl1 # 11 is the number of key node coordinates in the actuator module
                r̄s[i+1][3] += r̄s[i][3] + len2_b + len1_t - delta

            end
        end

        r̄s[end-1][2] = -r2 * sqrt(2) / 2
        r̄s[end-1][3] = r̄s[num_Y_modul+2][3] + prelenth - len1_b 
        r̄s[end-0][2] = r2 * sqrt(2) / 2
        r̄s[end-0][3] = r̄s[num_Y_modul+2][3] + prelenth - len1_b

        r̄s[end-3][2] = -r2 * sqrt(2) / 2
        r̄s[end-3][3] = prelenth
        r̄s[end-2][2] = r2 * sqrt(2) / 2
        r̄s[end-2][3] = prelenth

        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        movable = true
        ci = Vector{Int}()
        Φi = collect(1:6)

        ri = r̄si # Initial CoM x-axis coordinates of each module
        # @show ri
        aps = rdsi.r̄l # Connection point positions

        if i >= num_Y_modul + 2
            theta = deg2rad(rot_ang[2]) # Initial deflection of each section; set to zero here if needed
        else
            theta = deg2rad(rot_ang[1]) # Initial deflection of each section; set to zero here if needed
        end
        R = [
            cos(theta) sin(theta) 0;
            -sin(theta) cos(theta) 0;
            0 0 1
        ]

        r̄g = rdsi.r̄g # CoM coordinates
        m = rdsi.m # Mass
        I = rdsi.i # Inertia tensor
        nrp = length(aps) # Number of nodes
        ṙo = zeros(3)
        ω = zeros(3)
        prop = RB.RigidBodyProperty(i, movable, m, SMatrix{3,3}(I), 
            RB.Locus(SVector{3}(r̄g)), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
                [SVector{3}(aps[i]) for i in 1:nrp], 
                [SVector{3}([0.0, 0.0, 1.0]) for i in 1:nrp], #normals, 
                zeros(nrp), 
                zeros(nrp)
            )]
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω,)

        nmcs = RB.NCF.NC1P3V(ri, ro,)# Initial natural coordinates?
        coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)

    #---Springs----##########################################################
    ncables = (num_Y_modul + num_Y_modul2) * (row_outer_spring + row_inter_spring) # This refers to inner and outer springs
    ks1 = ks1
    restlens1 = restlens1 # inner springs
    ks2 = ks2
    restlens2 = restlens2 # outer springs
    spring_dampers = [(i <= (num_Y_modul + num_Y_modul2) * row_inter_spring) ? RB.DistanceSpringDamper3D(restlens1, ks1, cs1) : RB.DistanceSpringDamper3D(restlens2, ks2, cs2) for i in 1:ncables] # First half are inner springs, second half are outer springs

    matrix_cnt_raw = Vector{Matrix{Int}}()
    # Define inner/outer spring connections
    for j in 1:num_Y_modul # connect inner springs
        # for j in 4i-3:4i
        s = zeros(3, nb-2n)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
        # end
    end
    for j in num_Y_modul+2:nb-1-2n # connect inner springs
        s = zeros(3, nb-2n)
        s[1, j] = 8
        s[1, j+1] = -11 # The numbers refer to the connection-point order defined above
        s[2, j] = 9
        s[2, j+1] = -12
        s[3, j] = 10
        s[3, j+1] = -13
        push!(matrix_cnt_raw, s)
    end
    for j in 1:num_Y_modul # connect outer springs
        # for j in 4i-3:4i
        s = zeros(3, nb-2n)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        push!(matrix_cnt_raw, s)
        # end
    end
    for j in num_Y_modul+2:nb-1-2n
        s = zeros(3, nb-2n)
        s[1, j] = 1
        s[1, j+1] = -1
        s[2, j] = 2
        s[2, j+1] = -2
        s[3, j] = 3
        s[3, j+1] = -3
        push!(matrix_cnt_raw, s)
    end
    connecting_matrix = RB.change_connecting_format(rigidbodies, reduce(vcat, matrix_cnt_raw))
    # Sliding cable - begin #############################################

    restlencs = zeros(nb - n -n)
    for i in 1:num_Y_modul+1 #nb-1
        if i ==1
            restlencs[i] = prelenth - len1_b - len1_t - 0.5
        else
            restlencs[i] = restlenc[1]
        end
    end
    for i in num_Y_modul+2:nb-n-n
        if i == num_Y_modul + 2
            restlencs[i] = prelenth - len1_b - len1_t - 0.5
        else
            restlencs[i] = restlenc[2]
        end
    end
    # @show restlencs
    nclusters_per_seg = 4
    cluster_segs = vcat(
        [[RB.DistanceSpringDamperSegment(restlencs[i], kc, c=c_c, prestress=pre) for i in 1:num_Y_modul+1]
         for _ = 1:nclusters_per_seg],
        [[RB.DistanceSpringDamperSegment(restlencs[num_Y_modul+1+i], kc, c=c_c, prestress=pre) for i in 1:num_Y_modul2+1]
         for _ = 1:nclusters_per_seg]
    )

    cluster_sps = vcat(
        [[RB.SlidingPoint(μ) for _ in 1:num_Y_modul]
         for _ = 1:nclusters_per_seg],
        [[RB.SlidingPoint(μ) for _ in 1:num_Y_modul2]
         for _ = 1:nclusters_per_seg]
    )

    cluster_matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:nclusters_per_seg*n]
    s1 = zeros(Int, num_Y_modul+1, 4, nclusters_per_seg)
    for k in 1:nclusters_per_seg
        if k == 1
            s1[1, :, k] = [nb - 3, 2, 1, 3 + k] # first CTC segment
        elseif k == 2
            s1[1, :, k] = [nb - 3, 3, 1, 3 + k] # first CTC segment
        elseif k == 3
            s1[1, :, k] = [nb-2, 3, 1, 3 + k] # first CTC segment
        elseif k == 4
            s1[1, :, k] = [nb-2, 2, 1, 3 + k] # first CTC segment
        end
        for l = 2:num_Y_modul+1
            s1[l, :, k] = [l-1, 3 + k, l, 3 + k]
        end
    end
    for k in 1:nclusters_per_seg
        push!(cluster_matrix_cnt_raw[k], s1[:, :, k])
    end

    s2 = zeros(Int, num_Y_modul2+1, 4, nclusters_per_seg)
    for k in 1:nclusters_per_seg
        if k == 1
            s2[1, :, k] = [nb - 1, 2, num_Y_modul+2, 3 + k] # first CTC segment
        elseif k == 2
            s2[1, :, k] = [nb - 1, 3, num_Y_modul+2, 3 + k] # first CTC segment
        elseif k == 3
            s2[1, :, k] = [nb, 3, num_Y_modul + 2, 3 + k] # first CTC segment
        elseif k == 4
            s2[1, :, k] = [nb, 2, num_Y_modul + 2, 3 + k] # first CTC segment
        end
        for l = 2:num_Y_modul2+1
            s2[l, :, k] = [l + num_Y_modul, 3 + k, l + num_Y_modul + 1, 3 + k]
        end
    end
    for k in nclusters_per_seg+1:nclusters_per_seg*2
        push!(cluster_matrix_cnt_raw[k], s2[:, :, k-nclusters_per_seg])
    end

    # Sliding cable - end ###############################################

    connecting_cluster_matrix = [reduce(vcat, cluster_matrix_cnt_raw[k]) for k in 1:nclusters_per_seg*n]
    cables, clusters = RB.connect_spring_and_frictionlessclusters(
        rigidbodies, spring_dampers,
        cluster_sps, cluster_segs,
        connecting_matrix, connecting_cluster_matrix,
        istart=nclusters_per_seg * 2
    )
    cst1 = RB.FixedBodyApparatus(nclusters_per_seg * 2 + ncables + 1, rbs[1])

    rj = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 2,
        RB.Hen2Egg(
            RB.Anchor(rbs[num_Y_modul+1], 14, 1),
            RB.Anchor(rbs[num_Y_modul+2], 11, 1)
        ),
        RB.NoForce(),
        :Revolute;
    )
    rjc = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 3,
        RB.Hen2Egg(
            RB.Anchor(rbs[num_Y_modul+1], 14, 1),
            RB.Anchor(rbs[num_Y_modul+2], 11, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :FloatingUniversal;
    )
    PrJoint1 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 4,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-2], 1, 1),
            RB.Anchor(rbs[1], 14, 1)
        ),
        RB.NoForce(),
        :Prismatic;
    )
    PrJoint2 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 5,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-3], 1, 1),
            RB.Anchor(rbs[1], 13, 1)
        ),
        RB.NoForce(),
        :Prismatic;
    )
    PSJoint1 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 6,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-2], 1, 1),
            RB.Anchor(rbs[1], 14, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :PlanarSpherical;
    )
    PSJoint2 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 7,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-3], 1, 1),
            RB.Anchor(rbs[1], 13, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :PlanarSpherical;
    )
    PrJoint3 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 8,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb], 1, 1),
            RB.Anchor(rbs[num_Y_modul+2], 14, 1)
        ),
        RB.NoForce(),
        :Prismatic;
    )
    PrJoint4 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 9,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-1], 1, 1),
            RB.Anchor(rbs[num_Y_modul+2], 13, 1)
        ),
        RB.NoForce(),
        :Prismatic;
    )
    PSJoint3 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 10,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb], 1, 1),
            RB.Anchor(rbs[num_Y_modul+2], 14, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :PlanarSpherical;
    )
    PSJoint4 = RB.proto_joint_apparatus(
        nclusters_per_seg * 2 + ncables + 11,
        RB.Hen2Egg(
            RB.Anchor(rbs[nb-1], 1, 1),
            RB.Anchor(rbs[num_Y_modul+2], 13, 1)
        ),
        RB.RheonomicJointForce(
        ),
        :PlanarSpherical;
    )
    # rj = RB.FixedJoint(nclusters_per_seg * 2 + ncables + 2, RB.Hen2Egg(RB.Anchor(rbs[num_Y_modul+1], 14, 1), RB.Anchor(rbs[num_Y_modul+2], 11, 1)))
    # rj = RB.FixedJoint(nclusters_per_seg * 2 + ncables + 2, RB.Hen2Egg(RB.Anchor(rbs[num_Y_modul+1], 14, 1), RB.Anchor(rbs[num_Y_modul+2], 11, 1)))
    # rj1 = RB.SphericalJoint(nclusters_per_seg * 2 + ncables + 3, RB.Hen2Egg(RB.Anchor(rbs[num_Y_modul+2], 12, 1), RB.Anchor(rbs[num_Y_modul+3], 15, 1)))
    apparatuses = TypeSortedCollection(vcat([cst1, rj, rjc, PrJoint1, PrJoint2, PSJoint1, PSJoint2, PrJoint3, PrJoint4, PSJoint3, PSJoint4], cables, clusters)) # middle revolute joint
    # apparatuses = TypeSortedCollection(vcat([cst1, rj], cables, clusters)) # use this when rigidly connecting the middle joint

    cnt = RB.Connectivity(rigidbodies, apparatuses;)
    structure = RB.Structure(rigidbodies, apparatuses, cnt)
    gauges = Int[] # Gauges/sensors used to measure errors for cost functions, optimization, or feedback control; currently unused.
    actuators = reduce(vcat, [ # Actuators attached to one or more bodies or apparatuses; their types dispatch to execute! for different driving forces
        # Try whether ExternalForceActuator is sufficient; if not, create a new actuator

        # [
        #     RB.ExternalForceActuator(
        #         i,
        #         cluster,
        #         # Operator, 1 control input
        #         RB.NaiveOperator(1),
        #         # External force, may not be needed
        #         0.0,
        #         # Control input values
        #         [0.0],
        #     )
        #     for (i, cluster) in enumerate(clusters)
        # ],
        RB.RegisterActuator(
            1,#n * nclusters_per_seg + 
            rjc, # Create apparatus collection
            # Operator, 1 control input
            RB.FunctionOperator(
                (values, u) -> values .+ u # Define operator behavior
            ),
            # Control input values, Reg
            (values=[0.0],)
        ), # Use this to drive the revolute joint
        # [RB.ExternalForceActuator(
        #     n*nclusters_per_seg+1,
        #     RB.Anchor(rbs[num_Y_modul+2], 11),
        #     # Operator, 1 control input
        #     RB.NaiveOperator(1),
        #     # External force, may not be needed
        #     [1.0,1.0,0.0],
        #     # Control input values
        #     [0.0],
        # )]
        [RB.RegisterActuator(
            2,# n * nclusters_per_seg +
            PSJoint1, # Create apparatus collection
            # Operator, 1 control input
            RB.FunctionOperator(
                (values, u) -> values .+ u # Define operator behavior
            ),
            # Control input values, Reg
            (values=[0.0],)
        ) # Use this to drive the revolute joint
        ],
        [RB.RegisterActuator(
            3,#n * nclusters_per_seg + 
            PSJoint2, # Create apparatus collection
            # Operator, 1 control input
            RB.FunctionOperator(
                (values, u) -> values .+ u # Define operator behavior
            ),
            # Control input values, Reg
            (values=[0.0],)
        ) # Use this to drive the revolute joint
        ],
        [RB.RegisterActuator(
            4,# n * nclusters_per_seg +
            PSJoint3, # Create apparatus collection
            # Operator, 1 control input
            RB.FunctionOperator(
                (values, u) -> values .+ u # Define operator behavior
            ),
            # Control input values, Reg
            (values=[0.0],)
        ) # Use this to drive the revolute joint
        ],
        [RB.RegisterActuator(
            5,#n * nclusters_per_seg + 
            PSJoint4, # Create apparatus collection
            # Operator, 1 control input
            RB.FunctionOperator(
                (values, u) -> values .+ u # Define operator behavior
            ),
            # Control input values, Reg
            (values=[0.0],)
        ) # Use this to drive the revolute joint
        ]

    ])

    hub = RB.ControlHub(
        structure,
        gauges,
        actuators,
        RB.Coalition(structure, gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)

end
