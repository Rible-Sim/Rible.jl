
function man_ndof(ndof;θ=0.0,k=0.0,c=0.0,unit="mks")
    nbodies = ndof + 1
    nbp = 2nbodies - ndof
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    lower_index = 1:2:nbodies
    upper_index = 2:2:nbodies
    a = zeros(nbodies)
    m = zeros(nbodies)
    Ia = zeros(nbodies)
    if unit == "cgs"
        unit_L = 1e2
        unit_M = 1e3
        unit_I = unit_M*unit_L^2
    else
        unit_L = 1
        unit_M = 1
        unit_I = 1
    end
    for (i,j) in enumerate(lower_index)
        a[j] = 20.0e-2unit_L
        m[j] = 835.90254985e-3unit_M
        # Ia[j] = Ic_lower[i] + m[j]*1/3*a[j]^2
        Ia[j] = 28130.53053840*2e-7unit_I
        @show a[j],m[j],Ia[j]
    end
    for (i,k) in enumerate(upper_index)
        a[k] = 16.0e-2unit_L
        m[k] = 666.25659673e-3unit_M
        # Ia[k] = Ic_upper[i] + m[k]*1/3*a[k]^2
        Ia[k] = 14490.67513310*2e-7unit_I
        @show a[k],m[k],Ia[k]
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
    cs .= c
    for i = 1:nstrings
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],k,k)
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



function build_Y(tgstruct)
    @unpack nstrings,strings,actuators = tgstruct
    nact = length(actuators)
    ret = spzeros(Int,nstrings,nact)
    for (i,act) in enumerate(actuators)
        s1,s2 = act.strings
        is1 = findfirst((x)->x==s1,strings)
        is2 = findfirst((x)->x==s2,strings)
        ret[is1,i] = 1
        ret[is2,i] = -1
    end
    ret
end
