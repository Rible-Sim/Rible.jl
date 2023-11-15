function dualtri(num_of_dof,;onedir=[1.0,0.0],θ=0.0,k=400.0,c=0.0,restlen=0.16)
    num_of_bodies = num_of_dof + 1
    nbp = 2nbodies - num_of_dof
    n_lower = count(isodd,1:num_of_bodies)
    n_upper = count(iseven,1:num_of_bodies)
    lower_index = 1:2:num_of_bodies
    upper_index = 2:2:num_of_bodies
    a = zeros(num_of_bodies)
    m = zeros(num_of_bodies)
    Ia = zeros(num_of_bodies)
    for (i,j) in enumerate(lower_index)
        a[j] = 12252e-5
        m[j] = 200.40071778e-3
        Ia[j] = 3738.8e-7
    end
    for (i,k) in enumerate(upper_index)
        a[k] = 12252e-5
        m[k] = 200.40071778e-3
        Ia[k] = 3738.8e-7
    end
    R = [
        cos(θ) -sin(θ);
        sin(θ) cos(θ)
    ]
    A = zeros(2,nbp)
    A[:,2] .= A[:,1] .+ a[1]*onedir
    for i in 3:nbp
        A[:,i] .= A[:,i-1] .+ a[i-1]*R^(i-2)*onedir
    end

    function rigidbody(i,m,a,Ia,ri,α)
        movable = true
        mass_center_position = SVector{2}([0.0,0.0])

        r̄p1 = SVector{2}([-a/2,0.0])
        r̄p2 = SVector{2}([ a/2,0.0])
        if isodd(i)
            r̄p3 = SVector{2}([0.0,-√3/2*a])
            r̄p4 = r̄p3 .+ SVector{2}([ a/2,tan(deg2rad(45))])
            r̄p5 = r̄p3 .+ SVector{2}([ a/2,tan(deg2rad(45))])
        else
            r̄p3 = SVector{2}([0.0,√3/2*a])
            r̄p4 = r̄p3 .- SVector{2}([ a/2,tan(deg2rad(45))])
            r̄p5 = r̄p3 .- SVector{2}([ a/2,tan(deg2rad(45))])
        end

        nodes_positions = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]

        axes = [SVector{2}([1.0,0.0]),]
        Ī = SMatrix{2, 2}([
            0.99Ia 0
            0 0.01Ia
        ])
        
        if i == 1
            constrained = true
            ci = collect(1:6)
            cstr_idx = Int[]
        else
            constrained = false
            ci = Int[]
            cstr_idx = collect(1:3)
        end
        
        ω = 0.0
        ro = ri
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
                    i,movable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    constrained=constrained
                    )

        nmcs = RB.NCF.NC1P2V(SVector{2}(ri), ro, α)

        state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)

        body = RB.RigidBody(prop,state)
    end
    rbs = [
        begin
            d = A[:,i+1]-A[:,i]
            rigidbody(
                i,m[i],a[i],
                Ia[i],
                (A[:,i]+A[:,i+1])./2,
                atan(d[2],d[1])
            )
        end
        for i = 1:num_of_bodies
    ]

    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)

    indexed = RB.index(rigdibodies)

    ncables = 2(num_of_bodies-1)
    upstringlen = restlen
    lostringlen = restlen
    original_restlens = zeros(ncables)
    restlens = zeros(ncables)
    actuallengths = zeros(ncables)
    ks = zeros(ncables)
    cs = zeros(ncables)
    cs .= c
    for i = 1:ncables
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],k,k)
    end
    ss = [RB.Cable2D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
    tensiles = (cables=ss,)

    matrix_cnt = zeros(Int,2(num_of_bodies-1),num_of_bodies)
    for i = 1:num_of_bodies-1
        ul = [1,-3]
        dl = [3,-2]
        if iseven(i)
            ul,dl = dl,ul
        end
        matrix_cnt[2(i-1)+1,i:i+1] = ul
        matrix_cnt[2(i-1)+2,i:i+1] = dl
    end

    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)

    function ganged_act(actid,id1,id2,original_restlens)
        ids = [id1,id2]
        original_values = original_restlens[[id1,id2]]
        RB.ManualActuator(actid,ids,original_values,RB.Ganged())
    end
    acs = [ifelse(isodd(i),ganged_act(i,2(i-1)+1,2i,original_restlens),
                           ganged_act(i,2i,2(i-1)+1,original_restlens)) for i = 1:num_of_bodies-1]
    hub = (actuators=acs,)
    pjs = [
        RB.PinJoint(i,RB.Hen2Egg(i,RB.ID(rbs[i],2,1),RB.ID(rbs[i+1],1)))
        for i = 1:num_of_bodies-1
    ]
    jointed = RB.join(pjs,indexed)

    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,hub)
end
