function move(trans,rot,o,X)
    x = @view X[1:3,1]
    y = @view X[1:3,2]
    z = @view X[1:3,3]
    f = trans ∘ rot
    f(o),hcat(rot(x),rot(y),rot(z))
end

function bend(n,dis,rm;
			o1 = SVector{3}(zeros(3)),
		    X1 = Matrix(1.0I,3,3),
			dir = [0,0,1]
			)
    rot = LinearMap(rm)
    trans = Translation(dis.*dir)

    os = [o1]
    Xs = [X1]
    for i = 2:n
        oi,Xi = move(trans,rot,os[end],Xs[end])
        push!(os,oi)
        push!(Xs,Xi)
    end
    os,Xs
end

function spine(n,dis,rm,heating_law;k=1000.0,c=0.0)
    nbodies = n
    nbp = 4*n

	mass = 1.6730422e-01 #kg
	inertia = Matrix(Diagonal([7.7035410e-01,7.7035412e-01,1.1875254]))*1.e-3
    r̄g = [0,0,-1.1037596].*1e-1

	pp = Vector{Vector{Float64}}(undef,13)
	pp[1] =  [ 16.7,    0,5]
	pp[2] =  [ -8.4, 14.5,5]
	pp[3] =  [ -8.4,-14.5,5]
	pp[4] =  [ 96.2,    0,108.7]
	pp[5] =  [-48.1, 83.3,108.7]
	pp[6] =  [-48.1,-83.3,108.7]
	pp[7] =  [ 96.2,    0,139.7]
	pp[8] =  [-48.1, 83.3,139.7]
	pp[9] =  [-48.1,-83.3,139.7]
	pp[10] = [ 65.2,    0,124.2]
	pp[11] = [-32.6, 56.5,124.2]
	pp[12] = [-32.6,-56.5,124.2]
	pp[13] = zeros(3)
    aps = SVector{3}.(map((x)->RotX(π)*x,pp)).*1e-3
	# scatter(aps)

    movable = ones(Bool,n)
	# movable[1] = false
    constrained = zeros(Bool,n)
    constrained[1] = true

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(r̄g),
                aps;constrained=constrained[i]) for i = 1:n]

    rs_raw,Rs = bend(n,dis,rm)
    ṙs = [zeros(3) for i = 1:n]
    ωs = [zeros(3) for i = 1:n]
	rs = map((r)->r+=[0,0,124.2/1000],rs_raw)
    function rigidbody(i,prop,aps,ro,R,ṙo,ω)
        ri,rj,rk,rl = [ro+R*ap for ap in [aps[13],aps[10],aps[11],aps[12]]]
        # lncs,q,q̇ = TR.NCF.NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
		# lncs,q,q̇ = TR.NCF.NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
		# lncs,q,q̇ = TR.NCF.NC2P2V(ri,rk,ro,R,ṙo,ω)
		lncs,q,q̇ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
		if i == 1
			pres_idx = TR.find_full_pres_indices(lncs,q)
		else
			pres_idx = Int[]
		end
        state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,q,q̇,pres_idx)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,props[i],aps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]

    nstrings = 6*(n-1)
	restlen = 3e-2
    stringlenH = 0.036
    stringlenR = 0.030 #m
    stringlens = repeat(vcat(fill(stringlenH,3),fill(stringlenR,3)),n-1)
    kH = k #N/m
    kR = k #N/m
    ks = repeat(vcat(fill(kH,3),fill(kR,3)),n-1)
    # c = 0.0
    cs = repeat(fill(c,6),n-1)
	SMA_strings_indices = collect([j for i = 1:n-1 for j in (6*(i-1)+1):(6*(i-1)+3)])
	strings_indices = collect([j for i = 1:n-1 for j in (6*(i-1)+4):(6*(i-1)+6)])
    SMA_strings = [TR.SMAString3D(i,stringlens[i],TR.LinearLaw(ks[i],0.0),cs[i]) for i in SMA_strings_indices]
    strings = [TR.SString3D(i,stringlens[i],ks[i],cs[i]) for i in strings_indices]
	tensiles = (strings = strings, SMA_strings = SMA_strings)
	acs = [TR.ManualActuation(TR.SimpleActuator(i,stringlens[i])) for i = strings_indices]
	heaters = [TR.ManualSerialHeating(TR.SimpleActuators(i.+[6j for j=0:n-2],fill(40.0,n-1)),fill(heating_law,n-1)) for i in 1:3]
	hub = (actuators=acs, heaters=heaters)
    bodynq = TR.get_nbodycoords(rbs[1])
    body2q_raw = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]
	# body2q = TR.filter_body2q(body2q_raw,rbs)
	body2q = body2q_raw

    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n-1
        for j = 4:6
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j+3)))
        end
        for j = 1:3
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j+9)))
        end
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	TR.jac_singularity_check(tg)
	tr = TR.TensegrityRobot(tg,hub)
	TR.heat!(tr,fill(40.0,length(heaters)))
    TR.update_strings_apply_forces!(tg)
    TR.update_SMA_strings_apply_forces!(tg)
	tr
end

function build_Y(tgstruct)
    (;nstrings,strings,actuators) = tgstruct
    nact = length(actuators)
    ret = spzeros(Int,nstrings,nact)
    for (i,act) in enumerate(actuators)
        s1 = act.strings[1]
        is1 = findfirst((x)->x==s1,strings)
        ret[is1,i] = 1
    end
    ret
end

function calculate_markers(p)
	pxy = [p[1],p[2],0]
	p̂xy = pxy/norm(pxy)
	θ = 0.0
	d = 0.010
	rot = RotZ(θ)
    p + d*rot*p̂xy
end

function spine_true(n,dis,rm;
					o1 = SVector{3}(zeros(3)),
					X1 = Matrix(1.0I,3,3),
					dir = [0,0,1])
    nbodies = n
    nbp = 4*n

	mass = 1.6730422e-01 #kg
	inertia = Matrix(Diagonal([7.7035410e-01,7.7035412e-01,1.1875254]))*1.e-3
    r̄g = [0,0,-1.1037596].*1e-1

	pp = Vector{Vector{Float64}}(undef,4)
	# pp[1] =  [ 16.7,    0,  5]
	# pp[2] =  [ -8.4, 14.5,  5]
	# pp[3] =  [ -8.4,-14.5,  5]
	# pp[4] =  [ 96.2,    0,108.7]
	# pp[5] =  [-48.1, 83.3,108.7]
	# pp[6] =  [-48.1,-83.3,108.7]
	pp[2] =  [ 96.2,    0,136.7] #7
	pp[3] =  [-48.1, 83.3,136.7] #8
	pp[4] =  [-48.1,-83.3,136.7] #9
	# pp[10] = [ 65.2,    0,124.2]
	# pp[11] = [-32.6, 56.5,124.2]
	# pp[12] = [-32.6,-56.5,124.2]
	# pp[13] = [    0,    0,-9.0]
	# pp[14] = [106.6,    0,127.2]
	# pp[15] = [-53.3, 92.3,127.2]
	# pp[16] = [-53.3,-92.3,127.2]
	pp[1] = [  0.0,  0.0,  5.0] #17
    aps = SVector{3}.(map((x)->RotX(π)*x.+[0.0,0.0,136.7],pp)).*1e-3
	# scatter(aps)
	# markers = calculate_markers.(aps[16:-1:14])
	# append!(aps,SVector{3}.(markers))
    movable = ones(Bool,n)
	# movable[1] = false
    constrained = zeros(Bool,n)
    constrained[1] = true

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(r̄g),
                aps;constrained=constrained[i]) for i = 1:n]

    rs_raw,Rs = bend(n,dis,rm;o1,X1,dir)
    ṙs = [zeros(3) for i = 1:n]
    ωs = [zeros(3) for i = 1:n]
	# rs = map((r)->r+=[0,0,124.2/1000],rs_raw)
	rs = rs_raw
    function rigidbody(i,prop,aps,ro,R,ṙo,ω)
        # ri,rj,rk,rl = [ro+R*ap for ap in [aps[1],aps[10],aps[11],aps[12]]]
		ri = copy(ro)
		# @show ri
        # lncs,q,q̇ = TR.NCF.NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
		# lncs,q,q̇ = TR.NCF.NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
		# lncs,q,q̇ = TR.NCF.NC2P2V(ri,rk,ro,R,ṙo,ω)
		lncs,q,q̇ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
		if i == 1
			ci = TR.find_full_pres_indices(lncs,q)
			Φi = collect(1:6)
		else
			ci = Int[]
			Φi = collect(1:6)
		end
        state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,ci,Φi)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,props[i],aps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]
	rigdibodies = TypeSortedCollection(rbs)
    numbered = TR.number(rigdibodies)
    indexed = TR.index(rigdibodies,)

    ncables = 6*(n-1)
	restlen = 3e-2
    cablelenH = 30.84760917385826e-3
    cablelenR = 22.75599611083790e-3 #m
    cablelens = repeat(vcat(fill(cablelenH,3),fill(cablelenR,3)),n-1)
    kH = 0.14290272007535298e3 #N/m
    kR = 1.132626507875379e3 #N/m
    ks = repeat(vcat(fill(kH,3),fill(kR,3)),n-1)
    c = 100.0
    cs = repeat(fill(c,6),n-1)
	cables_indices = collect([j for i = 1:n-1 for j in (6*(i-1)+1):(6*(i-1)+6)])
    cables = [TR.Cable3D(i,cablelens[i],ks[i],cs[i]) for i in cables_indices]
	tensiles = (cables = cables,)

	matrix_cnt = zeros(Int,6(n-1),n)
    for i = 1:n-1
        # for j = 4:6
        #     push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j+3)))
        # end
        # for j = 1:3
        #     push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j+9)))
        # end
        for j = 2:4
			matrix_cnt[6(i-1)+j-1,i:i+1] = [j,-j]
            # push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j)))
		end
		for j = 2:4
			matrix_cnt[6(i-1)+j+2,i:i+1] = [1,-j]
            # push!(string2ap,(TR.ID(i,1),TR.ID(i+1,j)))
        end
    end
	connected = TR.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	acs = [
		TR.ManualActuator(i,i,cablelens[i])
		for i = cables_indices
	]
	hub = (actuators=acs,)

	cnt = TR.Connectivity(numbered,indexed,tensioned,)
    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	bot = TR.TensegrityRobot(tg,hub)
end
