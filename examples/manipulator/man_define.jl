
function man_ndof(ndof,θ=0.0)
    nbodies = ndof + 1
    nbp = 2nbodies - ndof
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    a_lower = fill(0.1,n_lower)
    a_upper = fill(0.08,n_upper)
    m_lower = fill(0.02,n_lower)
    m_upper = fill(0.013,n_upper)
    Ic_lower = fill(11.3,n_lower)
    Ic_upper = fill(5.4,n_upper)
    lower_index = 1:2:nbodies
    upper_index = 2:2:nbodies
    a = zeros(nbodies)
    m = zeros(nbodies)
    Ia = zeros(nbodies)
    for (i,j) in enumerate(lower_index)
        a[j] = a_lower[i]
        m[j] = m_lower[i]
        Ia[j] = Ic_lower[i] + m[j]*1/3*a[j]^2
    end
    for (i,k) in enumerate(upper_index)
        a[k] = a_upper[i]
        m[k] = m_lower[i]
        Ia[k] = Ic_upper[i] + m[k]*1/3*a[k]^2
    end
    A = zeros(2,nbp)
    A[:,2] .= A[:,1] .+ a[1]*[1.0,0.0]
    for i in 3:nbp
        A[:,i] .= A[:,i-1] .+ a[i-1]*[cos(θ*(i-2)),sin(θ*(i-2))]
    end

    function rigidbody(i,m,a,Ia,ri,rj)
        if i == 1
            movable = true
            constrained = true
            constrained_index = [1,2,4]
        else
            movable = true
            constrained = false
            constrained_index = Vector{Int}()
        end
        CoM_x = a/2
        CoM_y = √3/6*a
        if isodd(i)
            CoM_y = -CoM_y
        end
        CoM = SVector{2}([CoM_x,CoM_y])

        nap = 3
        ap1 = SVector{2}([0.0,0.0])
        ap2 = SVector{2}([a,0.0])
        ap3_x = a/2
        ap3_y = √3/2*a
        if isodd(i)
            ap3_y = -ap3_y
        end
        ap3 = SVector{2}([ap3_x,ap3_y])
        aps = [ap1,ap2,ap3]

        prop = TR.RigidBodyProperty(i,movable,m,Ia,
                    CoM,aps;constrained=constrained
                    )
        state = TR.RigidBodyState(prop,ri,rj,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,m[i],a[i],
            Ia[i],A[:,i],A[:,i+1]) for i = 1:nbodies]

    nstrings = 2(nbodies-1)
    upstringlen = 0.9norm(rbs[2].state.p[3] - rbs[1].state.p[1])
    lostringlen = 0.9norm(rbs[2].state.p[2] - rbs[1].state.p[3])
    original_restlens = zeros(nstrings)
    restlens = zeros(nstrings)
    actuallengths = zeros(nstrings)
    ks = zeros(nstrings)
    cs = zeros(nstrings)
    cs .= 0.0
    for i = 1:nstrings
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],1.0e2,1.0e2)
    end
    ss = [TR.SString2D(original_restlens[i],ks[i],cs[i]) for i = 1:nstrings]
    acs = [
        ifelse(isodd(i),TR.Actuator(SVector{2}(ss[2(i-1)+1:2i])),
                        TR.Actuator(SVector{2}(ss[2i:-1:2(i-1)+1])))
        for i = 1:nbodies-1
    ]

    rb2p = [
        [i,i+1] for i = 1:length(rbs)
    ]
    body2q_raw = [
        [2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p
    ]
    body2q = TR.filter_body2q(body2q_raw,rbs)
    string2ap_raw = [
        [zeros(Int,2),zeros(Int,2)] for i = 1:length(ss)
    ]
    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:length(rbs)-1
        push!(string2ap,(TR.ID(i,1),TR.ID(i+1,3)))
        push!(string2ap,(TR.ID(i,3),TR.ID(i+1,2)))
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
