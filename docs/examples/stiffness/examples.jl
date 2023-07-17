function dualtri(ndof,;onedir=[1.0,0.0],θ=0.0,k=400.0,c=0.0,restlen=0.16)
    nbodies = ndof + 1
    nbp = 2nbodies - ndof
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    lower_index = 1:2:nbodies
    upper_index = 2:2:nbodies
    a = zeros(nbodies)
    m = zeros(nbodies)
    Ia = zeros(nbodies)
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
        r̄g = SVector{2}([0.0,0.0])

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

        r̄ps = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]

        ās = [SVector{2}([1.0,0.0]),]
		Ī = SMatrix{2, 2}([
			0.99Ia 0
			0 0.01Ia
		])
		
        if i == 1
			constrained = true
			ci = collect(1:6)
			Φi = Int[]
		else
			constrained = false
			ci = Int[]
			Φi = collect(1:3)
		end
		
		ω = 0.0
		ro = ri
		ṙo = zero(ro)

		prop = TR.RigidBodyProperty(
					i,movable,m,
					Ī,
                    r̄g,
					r̄ps,
                    ās;
					constrained=constrained
                    )

		lncs, q, _ = TR.NCF.NC1P2V(SVector{2}(ri), ro, α, ṙo, ω)

		state = TR.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)

        rb = TR.RigidBody(prop,state)
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
		for i = 1:nbodies
	]

	rigdibodies = TypeSortedCollection(rbs)
    numbered = TR.number(rigdibodies)

    indexed = TR.index(rigdibodies)

    ncables = 2(nbodies-1)
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
    ss = [TR.Cable2D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
	tensiles = (cables=ss,)

	matrix_cnt = zeros(Int,2(nbodies-1),nbodies)
    for i = 1:nbodies-1
		ul = [1,-3]
		dl = [3,-2]
		if iseven(i)
			ul,dl = dl,ul
		end
		matrix_cnt[2(i-1)+1,i:i+1] = ul
		matrix_cnt[2(i-1)+2,i:i+1] = dl
    end

	connected = TR.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	function ganged_act(actid,id1,id2,original_restlens)
		ids = [id1,id2]
		original_values = original_restlens[[id1,id2]]
		TR.ManualActuator(actid,ids,original_values,TR.Ganged())
	end
    acs = [ifelse(isodd(i),ganged_act(i,2(i-1)+1,2i,original_restlens),
						   ganged_act(i,2i,2(i-1)+1,original_restlens)) for i = 1:nbodies-1]
	hub = (actuators=acs,)
	pjs = [
		TR.PinJoint(i,TR.End2End(i,TR.ID(rbs[i],2,1),TR.ID(rbs[i+1],1)))
		for i = 1:nbodies-1
	]
	jointed = TR.join(pjs,indexed)

    cnt = TR.Connectivity(numbered,indexed,tensioned,jointed)
    tg = TR.TensegrityStructure(rigdibodies,tensiles,cnt)
    TR.TensegrityRobot(tg,hub)
end

function Tbars()
	SVo3 = SVector{3}([0.0,0.0,0.0])
	a = 1.0
	b = 1.0
	function make_base(i)
		movable = true
		constrained = true
        r̄g = SVo3
        r̄p1 = SVector{3}([ -a, b,0.0])
        r̄p2 = SVector{3}([ -a,-b,0.0])
        r̄p3 = SVector{3}([0.0, b,0.0])
        r̄p4 = SVector{3}([0.0,-b,0.0])
        r̄p5 = SVo3
        r̄ps = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5,]
        ās = [
			SVector{3}([1.0,0.0,0.0]),
			SVector{3}([0.0,1.0,0.0]),
		]
		m = 1.0
		Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
		ci = collect(1:12)
		Φi = Int[]	
		R = RotZ(0.0)	
		ω = SVo3
		ri = SVo3
		ro = ri
		ṙo = zero(ro)

		prop = TR.RigidBodyProperty(
					i,movable,m,
					Ī,
                    r̄g,
					r̄ps,
                    ās;
					constrained
                    )

		lncs, _, _ = TR.NCF.NC1P3V(ri, ro, R, ṙo, ω)

		state = TR.RigidBodyState(prop, lncs, ri, R, ṙo, ω, ci, Φi)

        TR.RigidBody(prop,state)
    end
	function make_slider(i;ri=SVo3,R=RotZ(0.0),)
		movable = true
		constrained = false
        r̄g = SVo3
        r̄p1 = SVo3
        r̄ps = [r̄p1,]
        ās = [SVector{3}([1.0,0.0,0.0]),]
		Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
		ω = SVo3
		ro = ri
		ṙo = zero(ro)
		m = 1.0

		prop = TR.RigidBodyProperty(
					i,movable,m,
					Ī,
                    r̄g,
					r̄ps,
                    ās;
					constrained
                    )

		lncs, _, _ = TR.NCF.NC1P3V(ri, ro, R, ṙo, ω)

		state = TR.RigidBodyState(prop, lncs, ri, R, ṙo, ω)

        TR.RigidBody(prop,state)
	end
	base = make_base(1)
	slider1 = make_slider(2;ri = SVector(-a,0.0,0.0))
	slider2 = make_slider(3;ri = SVo3,R=RotZ(π/2))
	bar = make_3d_bar(
		4,
		SVector(-a,0.0,0.0),SVo3;
	)
	rbs = [base,slider1,slider2,bar]
	rigdibodies = TypeSortedCollection(rbs)
    numbered = TR.number(rigdibodies)
	sm = [
		0 1 0 1;
		0 2 0 2;
		0 3 0 3;
		0 0 1 4;
		0 0 2 5;
		0 0 3 6;
	]
    indexed = TR.index(rigdibodies,sm)

    ncables = 5
	original_restlens = zeros(ncables)
	ks = zeros(ncables)
	cs = zeros(ncables)
    ss = [TR.Cable3D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
	tensiles = (cables=ss,)

	cm = [
		1 -1  0 0;
		2 -1  0 0;
		3  0 -1 0;
		4  0 -1 0;
		5 -1  0 0;
	]

	connected = TR.connect(rigdibodies, cm)
	tensioned = @eponymtuple(connected,)

	pj1 = TR.PrismaticJoint(1,TR.End2End(1,TR.ID(base,5,1),TR.ID(slider1,1,1)))
	pj2 = TR.PrismaticJoint(2,TR.End2End(2,TR.ID(base,5,2),TR.ID(slider2,1,1)))

	pjs = [
		pj1,pj2
	]
	jointed = TR.join(pjs,indexed)
	cnt = TR.Connectivity(numbered,indexed,tensioned,jointed)
    tg = TR.TensegrityStructure(rigdibodies,tensiles,cnt)
    TR.TensegrityRobot(tg,)
end