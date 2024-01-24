function man_nd(n;θ=0.0,k=0.0,c=0.0)
    nbodies = 2n - 1
    nbd = 2n
    A = zeros(2,nbd)
    unit_l = 0.2
    a = ones(nbodies)*unit_l
    A[1,2:end] = collect(1:2n-1)*unit_l
    m = ones(nbodies)*0.1
    Ia = m.*a.^2

    function rigidbody(i,m,a,Ia,ri,rj)
		contactable = true
        if i == 1
            visible = true
            pres_idx = [1,2,4]
        else
            visible = true
            pres_idx = Int[]
        end

        if isodd(i)
            r̄g_x = a/2
            r̄g_y = (sqrt(3)/2*(-sqrt(3)/6)+0.2*0.6)/(sqrt(3)/2 + 0.6)
        else
            r̄g_x = a/2
            r̄g_y = 0
        end
        mass_locus = SVector{2}([r̄g_x,r̄g_y])


        if isodd(i)
        ap1 = SVector{2}([0.0,0.0])
        ap2 = SVector{2}([a,0.0])
        ap3_x = a/2
        ap3_y = 0.6a
        ap4_x = a/2
        ap4_y = -√3/2*a
        ap3 = SVector{2}([ap3_x,ap3_y])
        ap4 = SVector{2}([ap4_x,ap4_y])
        aps = [ap1,ap2,ap3,ap4]
        else
            ap1 = SVector{2}([0.0,0.0])
            ap2 = SVector{2}([a,0.0])
            ap3 = SVector{2}([a/2,0.0])
            aps = [ap1,ap2,ap3]
        end


        prop = RB.RigidBodyProperty(i,contactable,m,Ia,
                    mass_locus,aps;visible=visible
                    )
        state = RB.RigidBodyState(prop,ri,rj,pres_idx)
        body = RB.RigidBody(prop,state)
    end

    rbs = [rigidbody(i,m[i],a[i],
            Ia[i],A[:,i],A[:,i+1]) for i = 1:nbodies]

    ncables = 4(n-1)
    upstringlen = 0.9norm(rbs[3].state.p[1] - rbs[1].state.p[3])
    lostringlen = 0.9norm(rbs[3].state.p[1] - rbs[1].state.p[4])
    original_restlens = zeros(ncables)
    restlens = zeros(ncables)
    actuallengths = zeros(ncables)
    ks = zeros(ncables)
    cs = zeros(ncables)
    cs .= c


    for i = 1:ncables
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,2],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,2],k,k)
    end
    ss = [RB.SString2D(original_restlens[i],ks[i],cs[i]) for i = 1:ncables]
    acs = [
        ifelse(isodd(i),RB.Actuator(SVector{2}(ss[2(i-1)+1:2i])),
                        RB.Actuator(SVector{2}(ss[2i:-1:2(i-1)+1])))
        for i = 1:2(n-1)
    ]

    rb2p = [
        [i,i+1] for i = 1:length(rbs)
    ]
    bodyid2q_raw = [
        [2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p
    ]
    bodyid2q = RB.filter_bodyid2q(bodyid2q_raw,rbs)
    string2ap_raw = [
        [zeros(Int,2),zeros(Int,2)] for i = 1:length(ss)
    ]
    string2ap = Vector{Tuple{RB.Signifier,RB.Signifier}}()
    for i = 1:n-1
        j = 2i - 1
        push!(string2ap,(RB.Signifier(j,3),RB.Signifier(j+2,1)))
        push!(string2ap,(RB.Signifier(j,2),RB.Signifier(j+2,3)))
        push!(string2ap,(RB.Signifier(j,4),RB.Signifier(j+2,1)))
        push!(string2ap,(RB.Signifier(j,2),RB.Signifier(j+2,4)))
    end
    cnt = RB.Connectivity(bodyid2q,string2ap)
    st = RB.Structure(rbs,ss,acs,cnt)
    RB.update_cables_apply_forces!(st)
    st
end
