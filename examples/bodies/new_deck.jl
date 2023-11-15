function new_deck(id,loci,ro,R,ri,box;
    movable = true,
    constrained = false,
    ci = Int[],
    cstr_idx = collect(1:6),
)
mat = filter(
    row -> row.name == "Teak", 
    material_properties
)[1]
density = mat.density |> ustrip
box_volume = Meshes.volume(box)
mass = density*box_volume
# bar_inertia_exp = :(1/12*$mass*$bar_length^2)
wid,len,hei = Meshes.sides(box)
# m = 0.2999233976
Īg = SMatrix{3,3}(
    Matrix(
        Diagonal(
            [
                mass*(hei^2+len^2)/12, #Ix
                mass*(hei^2+wid^2)/12, #Iy
                mass*(len^2+wid^2)/12, #Iz
            ]
        )
    )
)
mass_locus = [0,0,0.0]
pretty_table(
    SortedDict(
        [
            ("id", id),
            ("mass", mass),
            ("density", density),
            ("inertia", Īg),
            ("wid,len,hei", (wid,len,hei)),
            ("box_volume", box_volume)
        ]
    )
)

prop = RB.RigidBodyProperty(
    id,
    movable,
    mass,
    Īg,
    mass_locus,
    loci;
    constrained=constrained
)
ṙo = zero(ro)
ω = zero(ro)
nmcs = RB.NCF.NC1P3V(ri,ro,R)
# cf = RB.NCF.CoordinateFunctions(nmcs,q0,ci,free_idx)
# @show typeof(nmcs)
boxmesh = Meshes.boundary(box) |> simple2mesh
state = RB.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci,cstr_idx)
RB.RigidBody(prop,state,boxmesh)
end