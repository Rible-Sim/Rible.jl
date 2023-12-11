function new_pointmass(;
        origin_position = [0.0,0.0,1.0],
        origin_velocity = zero(origin_position),
        m = 1.0,
        μ = 0.3,
        e = 0.9
    )
    contactable = false
    visible = true
    Ia = SMatrix{3,3}(Matrix(m*I,3,3))
    mass_locus  = SVector{3}([ 0.0, 0.0, 0.0])
    r̄p1 = SVector{3}([ 0.0, 0.0, 0.0])
    loci = [r̄p1]
    axes = [SVector{3}([ 1.0, 0.0, 0.0])]
    friction_coefficients = [μ]
    restitution_coefficients = [e]
    prop = RB.RigidBodyProperty(
        1,contactable,m,Ia,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients,
        ;visible=visible
    )
    ω = zero(origin_position)
    R = RotX(0.0)
    loci = Ref(origin_position) .+ Ref(R).*loci
    nmcs = RB.NCF.NC3D1P(loci[1],)
    ci = Int[]
    cstr_idx = Int[]
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)
    coords = RB.NonminimalCoordinates(nmcs,ci,cstr_idx)
    rb1 = RB.RigidBody(prop,state,coords)

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
