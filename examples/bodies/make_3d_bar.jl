function make_3d_bar(
    id,
    ri,rj;
    ci = Int[],
    m = 0.080,
    radius_ratio = 1/30,
    mat_name = nothing, #"Teak",
    loadmesh = false,
    loadmesh2 = false,
)
# @show id,ri,rj
movable = true
if ci == Int[]
    constrained = false
else
    constrained = true
end
u = rj - ri
bar_length = norm(u)	
radius = bar_length*radius_ratio
û = u./bar_length
v̂,ŵ = RB.NCF.HouseholderOrthogonalization(û)
R = SMatrix{3,3}(hcat(û,v̂,ŵ))

mass_locus  = SVector{3}([ bar_length/2,0,0])
r̄p1 = SVector{3}([          0.0,0,0])
r̄p2 = SVector{3}([ bar_length,  0,0])
loci = [r̄p1,r̄p2]
if mat_name isa Nothing
    mass = m
else
    mat = filter(
        row -> row.name == mat_name, 
        material_properties
    )[1]
    density = mat.density |> ustrip
    sec_area = π*radius^2
    density_per_unit_len = density*sec_area
    mass = density_per_unit_len*bar_length
end
bar_inertia_exp = :(1/12*$mass*$bar_length^2)
Īg = SMatrix{3,3}(
    [	
        eval(bar_inertia_exp) 0.0 0.0;
        0.0 0.0 0.0;
        0.0 0.0 0.0
    ]
)
pretty_table(
    SortedDict(
        [
            ("id", id),
            ("radius", radius),
            # ("density", density),
            ("bar length", bar_length),
            ("mass", mass),
            ("inertia", Īg),
            ("bar_inertia_exp", bar_inertia_exp)
        ]
    )
)
axes = [SVector(1.0,0.0,0.0) for _ in eachindex(loci)]
prop = RB.RigidBodyProperty(
    id,
    movable,
    mass,
    Īg,
    mass_locus,
    loci,
    axes,;
    constrained=constrained
)
# @show prop.inertia
ṙo = zero(ri)
ω = zero(ri)
nmcs = RB.NCF.NC3D1P1V(ri,û,ri,R)
qcs = RB.QCF.QC(mass,Īg[1]*I(3))
@myshow id,û
# @show ri,rj,q0
# cf = RB.NCF.CoordinateFunctions(nmcs,q0,ci,free_idx)
# @show typeof(nmcs)
state = RB.RigidBodyState(prop,ri,R,ṙo,ω)
# coords = RB.NonminimalCoordinates(nmcs,ci)
coords = RB.NonminimalCoordinates(qcs,ci)
if loadmesh
    barmesh = load(RB.assetpath("装配体3.STL")) |> RB.make_patch(;
        # trans=[0,0,0.025],
        scale=1/500,
        color=:palegreen3,
    )
elseif loadmesh2
    barmesh = load(RB.assetpath("球铰/零件2.STL")) |> RB.make_patch(;
        scale=1/200,
        color=:palegreen3,
    )
else
    barmesh = endpoints2mesh(r̄p1,r̄p2;radius,)
end
body = RB.RigidBody(prop,state,coords,barmesh)
end