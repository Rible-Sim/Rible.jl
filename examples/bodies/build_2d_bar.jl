function build_2d_bar(id,ri,rj;α = 0.0, ci = Int[])
    contactable = true
    if ci == Int[]
        visible = true
    else
        visible = true
    end
    u = rj - ri
    b = norm(u)
    α = atan(u[2],u[1])
    mass_locus  = SVector{2}([ b/2,0])
    r̄p1 = SVector{2}([ 0.0,0])
    r̄p2 = SVector{2}([   b,0])
    loci = [r̄p1,r̄p2]
    m = 2.6934977798E-02
    Īg = SMatrix{2,2}(
        [
            1/12*m*b^2 0.0;
                   0.0 0.0;
        ]
    )
    prop = RB.RigidBodyProperty(id,contactable,m,Īg,
                mass_locus,loci;visible=visible
                )
    ω = 0.0
    ro = ri
    ṙo = zero(ro)
    nmcs = RB.NCF.NC2D2P(ri,rj,ro,α)
    state = RB.RigidBodyState(prop,ro,α,ṙo,ω,ci)
    body = RB.RigidBody(prop,state)
end