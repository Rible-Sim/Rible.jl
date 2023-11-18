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
    ci = Int[],
    cstr_idx = collect(1:6),
    radius = 0.005,
    height = 1e-2,
    color = :darkorchid4,
    loadmesh = false,
    m = 3,
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
if id == 7
    mass_center = [0,0,0.1562442983-0.1*√2]
else
    mass_center = [0,0,0.1562442983-0.1*√2-0.1*√2/2]
end
# @show m,diag(Īg),mass_center

axes = [
    SVector(1.0,0,0)
    for i = eachindex(loci_positions)
]
prop = RB.RigidBodyProperty(
    id,
    contactable,
    mass,Īg,
    mass_center,
    Ref(inv(R)).*loci_positions,
    axes;
    visible=visible
)
ṙo = zero(ro)
ω = zero(ro)
u = R*(loci_positions[2] - loci_positions[1])
v = R*(loci_positions[3] - loci_positions[1])
w = R*(loci_positions[4] - loci_positions[1])
if rj isa Nothing
    nmcs = RB.NCF.NC1P3V(ri,ro,R)
elseif rk isa Nothing
    nmcs = RB.NCF.NC2P2V(ri,rj,ro,R)
elseif rl isa Nothing
    nmcs = RB.NCF.NC3P1V(ri,rj,rk,ro,R)
else
    nmcs = RB.NCF.NC4P(ri,rj,rk,rl,ro,R)
end
pretty_table(
    SortedDict(
        [
            ("id", id),
            ("density", density),
            ("mass", mass),
            ("inertia", Īg),
            ("mass center", mass_center),
            ("loci_positions[1]", loci_positions[1]),
            ("loci_positions[2]", loci_positions[2]),
            ("loci_positions[3]", loci_positions[3]),
            ("loci_positions[4]", loci_positions[4])
        ]
    )
)
# cf = RB.NCF.CoordinateFunctions(nmcs,q0,ci,free_idx)
# @show typeof(nmcs)
# radius = norm(loci_positions[2]-loci_positions[1])/32
# @show radius

if loadmesh
    trimesh = load(
        RB.assetpath("球铰/零件1.stl")
    ) |> make_patch(;
    scale = 1/200, 
    color = :tomato1,
    trans = [0,0,-0.39]
    )
else
    trimesh = GB.merge(
        [
            endpoints2mesh(loci_positions[i],loci_positions[j];
            radius,color)
            for (i,j) in [
                [1,2],[1,3],[1,4],
                [2,3],[3,4],[4,2]
            ]
        ]
    ) |> make_patch(;color = :darkslategrey)
    #     platemesh = endpoints2mesh(
    #         SVector(0.0,0.0,-height/2),
    #         SVector(0.0,0.0, height/2);
    #         radius=norm(loci[2]),n1=m,n2=2)
    #     trimesh = GB.merge(
    #         [
    #             trimesh,
    #             platemesh
    #         ]
    #     )
end
state = RB.RigidBodyState(prop,ro,R,ṙo,ω)
coords = RB.NonminimalCoordinates(nmcs,ci,cstr_idx)
body = RB.RigidBody(prop,state,coords,trimesh)
end
