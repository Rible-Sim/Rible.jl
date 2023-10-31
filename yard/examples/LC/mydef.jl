function man_nd1(k=0.0,c=0.0;ratio=0.8)
    n = 3

	function make_rb1(ri)
		movable = false
		constrained = true
		ci = collect(1:6)
		cstr_idx = Int[]
		h1 = 0.1
		b1 = 0.1
		a = 0.2
		m = 0.1
		Ia = SMatrix{2,2}(Matrix(m*a^2*I,2,2))
		r̄g_x = 0.0
		r̄g_y = 0.0
		mass_locus  = SVector{2}([r̄g_x,r̄g_y])
		ap1 = SVector{2}([ 0.0, 0.0])
		ap2 = SVector{2}([-b1 , -h1])
		ap3 = SVector{2}([ b1 , -h1])
	    aps = [ap1,ap2,ap3]
	    prop = RB.RigidBodyProperty(1,movable,m,Ia,
	                mass_locus,aps;constrained=constrained
	                )
		α = 0.0; ω = 0.0
		ro = ri
		ṙo = zero(ro)
		nmcs = RB.NCF.NC1P2V(ri,ro,α)
		state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,cstr_idx)
	    body = RB.RigidBody(prop,state)
	end

	function make_rb2(ri,rj,α)
		movable = true
		constrained = true
		pres_idx = [1,2]
		# h2 = 0.1
		b2 = 0.2
		a = 1.0
		m = 0.1
		Ia = SMatrix{2,2}(Matrix(m*a^2*I,2,2))
		r̄g_x = 0.0
		r̄g_y = 0.0
		mass_locus = SVector{2}([r̄g_x,r̄g_y])
		ap1 = SVector{2}([0.0  ,0.0])
		ap2 = SVector{2}([ b2  ,0.0])
		ap3 = SVector{2}([ b2/2,0.0])
	    aps = [ap1,ap2,ap3]
	    prop = RB.RigidBodyProperty(2,movable,m,Ia,
	                mass_locus,aps;constrained=constrained
	                )
		ω = 0.0
		ro = ri
		ṙo = zero(ro)
		nmcs = RB.NCF.NC2P1V(ri,rj,ro,α)
		state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,pres_idx)
	    body = RB.RigidBody(prop,state)
	end

	function make_bar2(ri,rj,α)
		movable = true
		constrained = true
		pres_idx = [1,2]
		# h2 = 0.1
		b2 = 0.2
		a = 1.0
		m = 0.1
		Ia = SMatrix{2,2}(Matrix(m*a^2*I,2,2))
		r̄g_x = 0.0
		r̄g_y = 0.0
		mass_locus = SVector{2}([r̄g_x,r̄g_y])
		ap1 = SVector{2}([0.0  ,0.0])
		ap2 = SVector{2}([ b2  ,0.0])
	    aps = [ap1,ap2]
	    prop = RB.RigidBodyProperty(2,movable,m,Ia,
	                mass_locus,aps;constrained=constrained
	                )
		ω = 0.0
		ro = ri
		ṙo = zero(ro)
		nmcs = RB.NCF.NC2D2P(ri,rj,ro,α)
		state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,pres_idx)
	    body = RB.RigidBody(prop,state)
	end

	function make_rb3(ri,α)
		movable = true
		constrained = false
		h3 = 0.1
		b3 = 0.1
		a = 0.2
		m = 0.1
		Ia = SMatrix{2,2}(Matrix(m*a^2*I,2,2))
		r̄g_x = 0.0
		r̄g_y = 0.0
		mass_locus = SVector{2}([r̄g_x,r̄g_y])
		ap1 = SVector{2}([0.0,0.0])
		ap2 = SVector{2}([-b3,-h3])
		ap3 = SVector{2}([ b3,-h3])
	    aps = [ap1,ap2,ap3]
	    prop = RB.RigidBodyProperty(3,movable,m,Ia,
	                mass_locus,aps;constrained=constrained
	                )
		ω = 0.0
		ro = ri
		ṙo = zero(ro)
		nmcs = RB.NCF.NC1P2V(ri,ro,α)
		state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω)
	    body = RB.RigidBody(prop,state)
	end

	rI = SVector( 0.0,0.0)
	α1 = π/2
	rJ = RB.NCF.rotation_matrix(α1)*SVector( 0.2,0.0)
	rb1 = make_rb1(rI)
	rb2 = make_rb2(rI,rJ,α1)
	# rb2 = make_bar2(rI,rJ,α1)
	α2 = 0.0
	rb3 = make_rb3(rJ,α2)
    rigidbodies = TypeSortedCollection([rb1,rb2,rb3])
	numberedpoints = RB.number(rigidbodies)
	matrix_sharing = [
		1 1 0;
		2 2 0;
		0 3 1;
		0 4 2;
	]
	indexedcoords = RB.index(rigidbodies,matrix_sharing)
	#

    ncables = 4
	# ratio = 0.85
    α = ratio
    β = ratio
	innercableslen = α*norm(rb1.state.loci_states[3] - rb2.state.loci_states[2])
    outercableslen = β*norm(rb1.state.loci_states[2] - rb3.state.loci_states[2])
    original_restlens = [innercableslen,innercableslen,outercableslen,outercableslen]
    # restlens = zeros(ncables)
    # actuallengths = zeros(ncables)
    ks = fill(100.0,ncables)
    cs = zeros(ncables)
    # @show [(1-α)/α*lostringlen;(1-β)/β*upstringlen;(1-β)/β*upstringlen;(1-α)/α*lostringlen].*ks
    # original_restlens = [lostringlen;upstringlen;upstringlen;lostringlen]
    # @show  original_restlens
	#
    cables = [RB.Cable2D(i,original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
    acs = [RB.ManualActuator(i,i,original_restlens[i]) for i = 1:ncables]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)

	cnt_matrix_cables = [
		2 -2  0;
	    3 -2  0;
		2  0 -2;
		3  0 -3;
	]
	connected = RB.connect(rigidbodies,cnt_matrix_cables)
	tensioned = @eponymtuple(connected,)
	cnt = RB.Connectivity(numberedpoints,indexedcoords,tensioned)

	st = RB.Structure(rigidbodies,tensiles,cnt)
    bot = RB.Robot(st,hub)
end
