function uni(c=100.0;
        μ = 0.9,
        e = 0.0,
        z0 = 0.2,
        ωz = 5.0,
        mbar = 0.05,
        Rbar = RotZ(0.0),
        isbody = false,
    )
    lₛ = 0.104
    DE = 0.415
    θ = 2π/3
    p = SVector{3}.(
        [   
            [0,0,0.0],
            lₛ .*[cos(0θ), sin(0θ),0],
            lₛ .*[cos(1θ), sin(1θ),0],
            lₛ .*[cos(2θ), sin(2θ),0],
            [0,0,-DE/2],
            [0,0, DE/2]
        ]
    )
    # fig,ax,plt = scatter(p)
    # # ax.aspect = DataAspect()
    # fig
    ṙo_ref = SVector(0.0,0,0)
    function rigidbase(i,p)
        contactable = true
        visible = true
        ci = collect(1:12)
        cstr_idx = Int[]
        ri = p[1]
        ro = ri
        R = SMatrix{3,3}(Matrix(1.0I,3,3))
        ṙo = ṙo_ref
        ω = [1.0,0,ωz]
        # ω = zeros(3)
        mass_locus = SVector{3}(0.0,0.0,0.0)
        # @show p
        loci = p[[2,3,4,1]]
        axes = [SVector(0,0,1.0) for i = 1:4]
        m = inertia = 0.1
        Ī = SMatrix{3,3}(
            [
                inertia 0       0;
                0       inertia 0;
                0       0       inertia
            ]
        )
        prop = RB.RigidBodyProperty(
            i,
            contactable,
            m,
            Ī,
            mass_locus,
            loci,
            axes;
            visible = visible,
        )
        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        state = RB.RigidBodyState(prop, nmcs, ri, R, ṙo, ω, ci, cstr_idx)
        trimesh = load("BASE.STL") |> make_patch(;
            scale=1/1000,
            trans = [0,0,0.005],
            rot = RotZ(π/2),
        )
        # long = 25

        # b1 = endpoints2mesh(loci[1],loci[2];radius = norm(loci[2]-loci[1])/long)
        # b2 = endpoints2mesh(loci[2],loci[3];radius = norm(loci[3]-loci[2])/long)
        # b3 = endpoints2mesh(loci[3],loci[1];radius = norm(loci[3]-loci[1])/long)
        # trimesh = GB.merge([b1,b2,b3])
        RB.RigidBody(prop, state, trimesh)
    end

    rb1 = rigidbar(
        1, 
        Rbar*p[end-1],
        Rbar*p[end];
        ṙi = ṙo_ref,
        ṙj = ṙo_ref,
        m = mbar,
        contactable = true,
        visible = true,
        ci = Int[],
        cstr_idx = collect(1:6),
        # cstr_idx = [1],
        isbody,
    )
    rb2 = rigidbase(2, p)
    rbs = [
        rb1,
        rb2
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = RB.number(rigdibodies)
    indexedcoords = RB.index(rigdibodies)
    #
    ncables = 6

    original_restlens = zeros(ncables)
    original_restlens[1:3] .= 80e-3
    original_restlens[4:6] .= 80e-3

    ks = zeros(ncables)
    ks[1:3] .= 900.0
    ks[4:6] .= 1000.0

    cables = [
    RB.DistanceSpringDamper3D(i, original_restlens[i], ks[i], c;slack=false) for i = 1:ncables
    ]
    #
    force_elements = (cables = cables,)
    acs = [
    RB.ManualActuator(
        i,
        collect(1:6), [original_restlens[6*(i-1)+j] for j = 1:6],
    ) for i = [1]
    ]
    hub = (actuators = acs,)
    # #
    matrix_cnt = zeros(Int,ncables,2)
    for j = 1:ncables
    if j < 4
        matrix_cnt[j, 1:2] = [1, -j]
    else
        matrix_cnt[j, 1:2] = [2, -(j-3)]
    end
    end
    # display(matrix_cnt)
    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)

    # uj = RB.PinJoint(1,
    #     RB.Hen2Egg(
    #         1,
    #         RB.ID(rb2,4,4),
    #         RB.ID(rb1,3,3),
    #     )
    # )

    # uj = RB.UniversalJoint(1,
    #     RB.Hen2Egg(
    #         1,
    #         RB.ID(rb2,4,4),
    #         RB.ID(rb1,3,3),
    #     )
    # )

    uj = RB.UniversalPrismaticJoint(1,
        RB.Hen2Egg(
            1,
            RB.ID(rb2,4,4),
            RB.ID(rb1,3,3),
        )
    )

    jointed = RB.join([uj],indexedcoords)

    cnt = RB.Connectivity(numberedpoints, indexedcoords, tensioned, jointed)

    st = RB.Structure(rigdibodies, force_elements, cnt,)
    bot = RB.Robot(st, hub)
end