function make_tail(n)
    nver = n
    nhor = n + 1
    nb = nver + nhor
    n_revolute = 2 + 2(n-1)
    ver_lengths = zeros(nver)
    hor_lengths = zeros(nhor)
    ver_lengths .= 0.04
    hor_lengths .= 0.04
    O = zeros(2,nhor)
    for i in 2:nhor
        O[2,i] = O[2,i-1] - ver_lengths[i-1]
    end
    P = zeros(2,nhor)
    for i in 1:nhor
        P[:,i] = O[:,i] .+ [hor_lengths[i]/2,0.0]
    end

    m = fill(1.3e-2,nb)
    L = zeros(nb)

    ver_index = 2:2:nb
    hor_index = 1:2:nb
    for (i,j) in enumerate(ver_index)
        L[j] = ver_lengths[i]
    end
    for (i,j) in enumerate(hor_index)
        L[j] = hor_lengths[i]
    end
    r̄g = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        r̄g[j] .= [ver_lengths[i]/2,0.0]
    end
    inertia = zeros(nb)
    for (i,j) in enumerate(ver_index)
        inertia[j] = m[j]*L[j]^2/3
    end
    for (i,j) in enumerate(hor_index)
        inertia[j] = m[j]*L[j]^2/12
    end
    ri = [zeros(2) for i = 1:nb]
    rj = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        ri[j] .= O[:,i]
        rj[j] .= O[:,i+1]
    end
    for (i,j) in enumerate(hor_index)
        ri[j] .= O[:,i]
        rj[j] .= P[:,i]
    end
    p1 = [zeros(2) for i = 1:nb]
    p2 = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        p1[j] .= [0.0,0.0]
        p2[j] .= [ver_lengths[i],0.0]
    end
    for (i,j) in enumerate(hor_index)
        p1[j] .= [-hor_lengths[i]/2,0.0]
        p2[j] .= [ hor_lengths[i]/2,0.0]
    end
    function rigidbody(i,r̄g,m,inertia,ri,rj,aps)
        if i == 1
            movable = true
            constrained = true
            constrained_index = [1,2,4]
        else
            movable = true
            constrained = false
            constrained_index = Vector{Int}()
        end
        nap = length(aps)
        aps = [SVector{2}(aps[i]) for i = 1:nap]
        prop = TR.RigidBodyProperty(i,movable,
                    m,inertia,
                    SVector{2}(r̄g),
                    aps;constrained=constrained
                    )
        state = TR.RigidBodyState(prop,ri,rj,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,r̄g[i],m[i],inertia[i],ri[i],rj[i],[p1[i],p2[i]]) for i = 1:nb]

    ncables = 4n
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for i = 1:ncables
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],0.04,√3*0.02)
        ks[i] = ifelse(j∈[1,0],50.0,100.0)
    end
    ss = [TR.SString2D(i,0.5original_restlens[i],ks[i],0.0) for i = 1:ncables]

    tensiles = (cables=ss,)
	acs = [TR.ManualActuation(TR.SimpleActuator(4(i-1)+j,5original_restlens[4(i-1)+j])) for i = 1:n  for j = 1:4]
	hub = (actuators=acs,)

    rb2p = [
        ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
    ]
    body2q_raw = [[2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p]
    body2q = TR.filter_body2q(body2q_raw,rbs)
    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n
        push!(string2ap,(TR.ID(2i+1,1),TR.ID(2i-1,1)))
        push!(string2ap,(TR.ID(2i+1,1),TR.ID(2i  ,1)))
        push!(string2ap,(TR.ID(2i+1,2),TR.ID(2i  ,1)))
        push!(string2ap,(TR.ID(2i+1,2),TR.ID(2i-1,2)))
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    TR.update_cables_apply_forces!(tg)
	TR.jac_singularity_check(tg)
    tr = TR.TensegrityRobot(tg,hub)
end

function on_curve(n,d,f)
    x = zeros(n)
    y = zeros(n)
    y[1] = f(x[1])
    for i = 2:n
        function func!(func,unk)
            xi,yi = unk
            dy = yi - y[i-1]
            dx = xi - x[i-1]
            func[1] = dy^2 + dx^2 - d^2
            func[2] = yi-f(xi)
        end
        soli = nlsolve(func!,[x[i-1]+d,f(x[i-1]+d)])
        x[i],y[i] = soli.zero
    end
    x,y
end

function rot2dvec(v,β)
    x1,y1 = v
    x2 = cos(β)*x1 - sin(β)*y1
    y2 = sin(β)*x1 - cos(β)*y1
    [x2,y2]
end

function make_curve_tail(n)
    nver = n
    nhor = n + 1
    nb = nver + nhor
    n_revolute = 2 + 2(n-1)
    ver_lengths = zeros(nver)
    hor_lengths = zeros(nhor)
    ver_lengths .= 0.04
    hor_lengths .= 0.04
    O = zeros(2,nhor)
    O[2,1] = -0.02
    for i in 2:nhor
        O[2,i] = O[2,i-1] - ver_lengths[i-1]
    end
    amp = 0.05
    period = 2π/(0.04*(n+1))
    Ox,Oy = on_curve(n+1,ver_lengths[1],(x)->0.0*amp*sin(period*x))
    Vy = similar(Oy)
    for i in eachindex(Oy)
        # Vy[i] = -2amp*cos(period*Ox[i])
        Vy[i] = -2amp*cos(period*O[2,i])
    end
    # return Vy
    # O[1,:] = Oy
    # O[2,:] = -Ox
    P = zeros(2,nhor)
    for i in 1:nhor
        if i == nhor
            tangent_vec = O[:,i] - O[:,i-1]
        else
            tangent_vec = O[:,i+1] - O[:,i]
        end
        tangent_unit = tangent_vec/norm(tangent_vec)
        P[:,i] = O[:,i] .+ hor_lengths[i]*rot2dvec(tangent_unit,-π/2)
    end

    m = fill(1.3e-2,nb)
    L = zeros(nb)

    ver_index = 2:2:nb
    hor_index = 1:2:nb
    for (i,j) in enumerate(ver_index)
        L[j] = ver_lengths[i]
    end
    for (i,j) in enumerate(hor_index)
        L[j] = hor_lengths[i]
    end
    r̄g = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        r̄g[j] .= [ver_lengths[i]/2,0.0]
    end
    inertia = zeros(nb)
    for (i,j) in enumerate(ver_index)
        inertia[j] = m[j]*L[j]^2/3
    end
    for (i,j) in enumerate(hor_index)
        inertia[j] = m[j]*L[j]^2/12
    end
    ri = [zeros(2) for i = 1:nb]
    rj = [zeros(2) for i = 1:nb]
    vi = [zeros(2) for i = 1:nb]
    vj = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        ri[j] .= O[:,i]
        rj[j] .= O[:,i+1]
        vi[j] .= [Vy[i],0.0]
        vj[j] .= [Vy[i+1],0.0]
    end
    for (i,j) in enumerate(hor_index)
        ri[j] .= O[:,i]
        rj[j] .= P[:,i]
        vi[j] .= [Vy[i],0.0]
        vj[j] .= [Vy[i],0.0]
    end
    p1 = [zeros(2) for i = 1:nb]
    p2 = [zeros(2) for i = 1:nb]
    for (i,j) in enumerate(ver_index)
        p1[j] .= [0.0,0.0]
        p2[j] .= [ver_lengths[i],0.0]
    end
    for (i,j) in enumerate(hor_index)
        p1[j] .= [-hor_lengths[i]/2,0.0]
        p2[j] .= [ hor_lengths[i]/2,0.0]
    end
    function rigidbody(i,r̄g,m,inertia,ri,rj,vi,vj,aps)
        if i == 0
            movable = true
            constrained = true
            constrained_index = [1,2]
        elseif i in [7,21]
            movable = true
            constrained = true
            constrained_index = [1]
        else
            movable = true
            constrained = false
            constrained_index = Vector{Int}()
        end
        nap = length(aps)
        aps = [SVector{2}(aps[i]) for i = 1:nap]
        prop = TR.RigidBodyProperty(i,movable,
                    m,inertia,
                    SVector{2}(r̄g),
                    aps;constrained=constrained
                    )
        state = TR.RigidBodyState(prop,ri,rj,vi,vj,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,r̄g[i],m[i],inertia[i],
                    ri[i],rj[i],vi[i],vj[i],
                    [p1[i],p2[i]]) for i = 1:nb]

    ncables = 4n
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for i = 1:ncables
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],0.04,√3*0.02)
        ks[i] = ifelse(j∈[1,0],50.0,100.0)
    end
    ss = [TR.SString2D(0.5original_restlens[i],ks[i],0.0) for i = 1:ncables]
    # @code_warntype   TR.DString(k[i],original_restlen[i],
    #         restlen[i],actuallength[i],zeros(MVector{4}))

    acs = [TR.Actuator(SVector{4}(ss[4(i-1)+1:4i])) for i = 1:n]

    rb2p = [
        ifelse(isodd(i),[i,i+1],[i-1,i+1]) for i = 1:length(rbs)
    ]
    body2q_raw = [[2pid[1]-1,2pid[1],2pid[2]-1,2pid[2]] for pid in rb2p]
    body2q = TR.filter_body2q(body2q_raw,rbs)
    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n
        push!(string2ap,(TR.ID(2i+1,1),TR.ID(2i-1,1)))
        push!(string2ap,(TR.ID(2i+1,1),TR.ID(2i  ,1)))
        push!(string2ap,(TR.ID(2i+1,2),TR.ID(2i  ,1)))
        push!(string2ap,(TR.ID(2i+1,2),TR.ID(2i-1,2)))
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_cables_apply_forces!(tg)
    tg
end
