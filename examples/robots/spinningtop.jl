
function make_top(origin_position = [0.0,0.0,0.0],
        R = one(RotMatrix{3}),
        origin_velocity = [0.0,0.0,0.0],
        Ω = [0.0,0.0,5.0],
        cT = RB.QCF.QC;
        μ = 0.5,
        e = 0.9,
        color = :slategray,
        visible = true,
        loadmesh=false,
    )
    ω = R*Ω
    contactable = true
    pres_idx = Int[]

    m =  0.58387070
    mass_locus = @SVector zeros(3)
    # Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
    Ī = SMatrix{3,3}(
        Diagonal(SA[
            0.00022129,
            0.00022129,
            0.00030207
            ])
        )

    # h = 0.02292582
    radius = 0.044/√2
    h = 2*0.01897941
    loci = [radius.*[1,1,0] for i = 1:4]
    push!(loci,[0,0,-h])
    axes = [SVector(1.0,0,0) for i = 1:5]
    friction_coefficients = [μ for i = 1:5]
    restitution_coefficients = [e for i = 1:5]
    if loadmesh
        topmesh = load(
            RB.assetpath("Toupise2.STL")
        ) |> RB.make_patch(;
            scale=1/1000,
            color,
        )
    else
        pts = Point3.(loci)
        fcs = GB.TriangleFace.([
            [5,1,2],
            [5,4,3],
            [5,3,1],
            [5,2,4],
            [1,4,2],
            [4,1,3],
            [3,2,1],
            [2,3,4]
        ])
        nls = GB.normals(pts,fcs)
        topmesh = GB.Mesh(GB.meta(pts,normals=nls),fcs)
    end
    prop = RB.RigidBodyProperty(
        1,contactable,m,Ī,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients;
        visible
    )
    ri = origin_position+R*loci[5]
    @myshow ri
    if cT == RB.QCF.QC
        nmcs = RB.QCF.QC(m,Ī)
    else
        nmcs = RB.NCF.NC1P3V(ri,origin_position,R)
    end
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)
    coords = RB.NonminimalCoordinates(nmcs,pres_idx)
    rb1 = RB.RigidBody(prop,state,coords,topmesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
    force_elements = (cables = ss,)
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)
    cnt = RB.Connectivity(numberedpoints,indexedcoords,tensioned)
    st = RB.Structure(rbs,force_elements,cnt,)
    bot = RB.Robot(st)
end