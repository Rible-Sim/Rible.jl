

function build_2d_ground(id)
    movable = false
    constrained = true
    h1 = 0.1
    b1 = 0.1
    a = 0.2
    m = 0.1
    Ia = SMatrix{2,2}(Matrix(m*a^2*I,2,2))
    r̄g_x = 0.0
    r̄g_y = 0.0
    mass_locus  = SVector{2}([ 0.0, 0.0])
    r̄p1 = SVector{2}([ 0.0, 0.0])
    r̄p2 = SVector{2}([ -b1, 0.0])
    r̄p3 = SVector{2}([ 2b1, 0.0])
    loci = [r̄p1,r̄p2,r̄p3]
    prop = RB.RigidBodyProperty(id,movable,m,Ia,
                mass_locus,loci;constrained=constrained
                )
    α = 0.0; ω = 0.0
    R = RB.rotation_matrix(α)
    ro = zeros(2)
    ṙo = zero(ro)
    loci_states = Ref(ro) .+ Ref(R).*loci
    nmcs = RB.NCF.NCMP(loci_states,ro,R)
    # cf = RB.NCF.CoordinateFunctions(nmcs,q0,ci,free_idx)
    # @show typeof(nmcs)
    nq = length(q0)
    ci = collect(1:nq)
    free_idx = Int[]
    cstr_idx = Int[]
    state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,cstr_idx)
    body = RB.RigidBody(prop,state)
end