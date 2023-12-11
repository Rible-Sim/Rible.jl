#--  painleve's Bar
function make_bar(;
        μ = 0.1,
        e = 0.0
    )
    contactable = true
    visible = true
    m = 0.402
    l = 1.0
    Ixx = 1/12*m*l^2
    Ī = SMatrix{3,3}(Diagonal([Ixx,0.0,0.0]))
    mass_locus = SVector{3}(0.0,0.0,0.0)
    loci = SVector{3}.([
        [ -l/2, 0.0, 0.0],
        [  l/2, 0.0, 0.0]
    ])
    θ = π/4
    zoffset = 1e-6
    xoffset = 1e-6
    ri = [     0.0-xoffset,0.0,sin(θ)*l-zoffset]
    rj = [cos(θ)*l-xoffset,0.0,     0.0-zoffset]
    origin_position = (ri + rj)./2
    origin_velocity = zero(origin_position)
    ω = zero(origin_position)
    R = Matrix(RotY(θ))

    # b = l/2*ω[2]*sin(θ)-9.81
    # p⁺ = 1 + 3*cos(θ)^2-μ*cos(θ)*sin(θ)
    # @show b,p⁺

    prop = RB.RigidBodyProperty(1,contactable,m,Ī,mass_locus,loci;visible)

    nmcs = RB.NCF.NC3D2P(ri,rj,origin_position,R,origin_velocity,ω)
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω)

    p1 = Meshes.Point(loci[1])
    p2 = Meshes.Point(loci[2])
    s = Meshes.Segment(p1,p2)
    cyl_bar = Meshes.Cylinder(Meshes.length(s)/40,s)
    cylsurf_bar = Meshes.boundary(cyl_bar)
    cyl_bar_simple = Meshes.triangulate(cylsurf_bar)
    cyl_bar_mesh = cyl_bar_simple |> simple2mesh
    rb1 = RB.RigidBody(prop,state,cyl_bar_mesh)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    ss = Int[]
    force_elements = (cables = ss,)
    connections = RB.connect(rbs,zeros(Int,0,0))
    cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
    st = RB.Structure(rbs,force_elements,cnt,)
    bot = RB.Robot(st)
end
