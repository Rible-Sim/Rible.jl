
function make_3d_plate(
    id,
    nodes_positions,
    ro,
    R,
    ri,
    rj=nothing,
    rk=nothing,
    rl=nothing;
    movable = true,
    m = 3,
    height=1e-2,
    radius=1e-3,
    ci = Int[],
    cstr_idx = collect(1:6),
    loadmesh = true,
    meshvisible = true,
)
# free_idx = collect(1:6)	
constrained = ci != Int[]
mat = filter(
    row -> row.name == "Teak", 
    material_properties
)[1]
density = mat.density |> ustrip
if meshvisible
    if loadmesh
        platemesh = load("柚木板3.STL") |> make_patch(;rot=RotZ(π/4))
    else
        platemesh = endpoints2mesh(
            SVector(0.0,0.0,-height/2),
            SVector(0.0,0.0,height/2);
            radius,n1=m,n2=2)
    end
else
    platemesh = nothing
end
# mass = density*GB.volume(platemesh)
if m == 4
    mass =  0.11008275 #kg
    Īg = SMatrix{3,3}(
        Matrix(Diagonal([
            0.00025403,
            0.00025403,
            0.00050623
        ])
        )
    )
else
    mass = 0.23569851 #kg
    Īg = SMatrix{3,3}(
        Matrix(Diagonal([
            0.00085048,
            0.00085048,
            0.00169704				
        ])
        )
    )
end
mass_center = SVector(0.0,0.0,0.0)
pretty_table(
    SortedDict(
        [
            ("id", id),
            ("shape", "plate"),
            ("No. vertices",m),
            ("radius", radius),
            ("height", height),
            ("density", density),
            # ("bar length", bar_length),
            ("mass", mass),
            ("inertia", Īg),
            # ("bar_inertia_exp", bar_inertia_exp)
        ]
    )
)

axes = [SVector(1.0,0.0,0.0) for _ in eachindex(nodes_positions)]
prop = RB.RigidBodyProperty(
    id,
    movable,mass,Īg,
    mass_center,
    nodes_positions,
    axes,
    ;
    constrained=constrained
)
ṙo = zero(ro)
ω = zero(ro)
# u = R*(nodes_positions[2] - nodes_positions[1])
# v = R*(nodes_positions[3] - nodes_positions[1])
# w = R*(nodes_positions[4] - nodes_positions[1])
nmcs = RB.NCF.NC1P3V(ri,ro,R)
state = RB.RigidBodyState(prop,ro,R,ṙo,ω)
coords = RB.NonminimalCoordinates(nmcs,ci,cstr_idx)
body = RB.RigidBody(prop,state,coords,platemesh)
end