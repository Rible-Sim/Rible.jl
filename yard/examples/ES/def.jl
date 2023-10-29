function build_2d_bar(id,ri,rj;α = 0.0, ci = Int[])
	movable = true
	if ci == Int[]
		constrained = false
	else
		constrained = true
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
	prop = RB.RigidBodyProperty(id,movable,m,Īg,
				mass_locus,loci;constrained=constrained
				)
	ω = 0.0
	ro = ri
	ṙo = zero(ro)
	nmcs = RB.NCF.NC2D2P(ri,rj,ro,α)
	state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci)
	body = RB.RigidBody(prop,state)
end

function build_2d_tri(id,ri,rj=nothing,rk=nothing;
		α = 0.0, 
		ci = Int[], 
		constraints_indices = collect(1:3),		
		h = 0.05,
		b = 0.1,
		b1 = 0.05,
		m = 0.1271425597,
	)
	movable = true
	if isempty(ci)
		constrained = false
	else
		constrained = true
	end
	# unconstrained_indices = collect(1:6)
	b2 = b - b1
	a1 = 1.0b; b1 = √(2b^2-a1^2)
	a2 = 1.2b; b2 = √(2b^2-a2^2)
	mass_locus  = SVector{2}([ (b1+b)/3, h/3])
	r̄p1 = SVector{2}([      0.0, 0.0])
	r̄p2 = SVector{2}([       -a1, b2])
	r̄p3 = SVector{2}([        a1, b2])
	r̄p4 = SVector{2}([      -3a2, b2])
	r̄p5 = SVector{2}([       3a2, b2])
	loci = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]
	@myshow loci
	axes = [ SVector{2}([      1.0, 0.0])]
	Īgx = (b*h^3)/36
	Īgy = h*b*(b^2 - b*b1 + b1^2)/36
	A = b*h/2
	ρ = m/A
	Īgxy = -(b1^2*h^2)/24 + (b2^2*h^2)/24
	dx = -2/3*(b/2-b1)
	dy = -h/3
	Īgxy = Īgxy + A*dx*dy
	Īg = SMatrix{2,2}(
		ρ.*[
			Īgy Īgxy;
			Īgxy Īgx
		]
	)
	@show norm(mass_locus),m,Īg,tr(Īg)
	@show mass_locus,atan(mass_locus[2],mass_locus[1])
	prop = RB.RigidBodyProperty(id,movable,m,Īg,
				mass_locus,loci,;constrained=constrained
				)
	ω = 0.0
	ro = ri
	# ṙo = [0.0,0.0001]
	ṙo = zero(ro)
	if id == 7
		ṙo = SVector(0.01,0.0)
	# elseif id == 2
		ṙo = SVector(0.0,0.0001)
	end
	if rj isa Nothing
		nmcs = RB.NCF.NC1P2V(ri,ro,α)
	elseif rk isa Nothing
		nmcs = RB.NCF.NC2P1V(ri,rj,ro,α)
	else
		nmcs = RB.NCF.NC3P(ri,rj,rk,ro,α)
	end
	trimesh = begin
		if id == 1
			load("零件1 - 副本.STL") |> make_patch(;scale=1/1000,color=:mediumpurple4)
		else
			load("零件1.STL") |> make_patch(;scale=1/1000,color=:slategray4)
		end
	end
	state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,constraints_indices)
	body = RB.RigidBody(prop,state,trimesh)
end

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
	# cf = RB.NCF.CoordinateFunctions(nmcs,q0,ci,unconstrained_indices)
	# @show typeof(nmcs)
	nq = length(q0)
	ci = collect(1:nq)
	unconstrained_indices = Int[]
	constraints_indices = Int[]
	state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,constraints_indices)
	body = RB.RigidBody(prop,state)
end

function two_tri(;k=100.0,c=0.0,ratio=0.8)
	
	ro = SVector(0.0,0.0)
	α_tri1 = -π
	α_tri2 =  0
	rb1 = build_2d_tri(1,ro;
		α=α_tri1,
		ci=collect(1:6),
		constraints_indices=Int[]
	)
	rb2 = build_2d_tri(2,ro;
		α=α_tri2,
		ci=collect(1:2)
	)

	rbs = TypeSortedCollection((rb1,rb2))
	numberedpoints = RB.number(rbs)
	matrix_sharing = [
		1 1;
		2 2;
	]
	indexedcoords = RB.index(rbs,matrix_sharing)
	#
	ncables = 4
	restlen1 = 0.05
    restlens = fill(restlen1,ncables)
    ks = fill(k,ncables)
    cs = fill(c,ncables)
    cables = [RB.Cable2D(i,restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    acs = [RB.ManualActuator(1,1:ncables,restlens,RB.Uncoupled())]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)
	cnt_matrix_cables = [
		3 -2 ;
	   -2  3 ;
	    5 -2 ;
	   -4  3 ;
		# 5 -2 ;
	    # 4 -3 ;
	]
	connected = RB.connect(rbs,cnt_matrix_cables)
	#
	# cst1 = RB.PinJoint(
	# 	1,
	# 	RB.End2End(
	# 		1,
	# 		RB.ID(rb1,1),
	# 		RB.ID(rb2,1)
	# 	)
	# )

	# cst2 = RB.LinearJoint(
	# 	2,
	# 	2,
	# 	zeros(2),
	# 	begin
	# 		A = zeros(2,12)
	# 		A[1:2,1:2] = Matrix(1I,2,2)
	# 		A
	# 	end
	# )

	# jointedmembers = RB.join((cst1,),indexedcoords)

    cnt = RB.Connectivity(
		numberedpoints,
		indexedcoords,
		@eponymtuple(connected,),
		# jointedmembers
	)

    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,hub)
end

function one_bar(k=0.0,c=0.0;ratio=0.8)
    n = 1

	bps = SVector{2}.(
		[
		[0.0, 0.0],
		[0.1, 0.0]
		]
	)

	α1 = 0.0
	rb1 = build_2d_bar(1,bps[1],bps[2],α1)
    rbs = TypeSortedCollection((rb1,))
	numberedpoints = RB.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = RB.index(rbs,matrix_sharing)
	#
	ss = Int[]
    tensiles = (cables = ss,)
    hub = nothing
	#
	connections = RB.connect(rbs,zeros(Int,0,0))

	jointedmembers = RB.unjoin()

    cnt = RB.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers)

    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,hub)
end

function one_bar_one_tri()
    n = 2
	b1 = 0.05*√2
	bps = SVector{2}.(
		[
		   [  0.0,   0.0],
	       [b1/√2,-b1/√2]
		]
	)

	α1 = -π/4
	rb1 = build_2d_bar(1,bps[1],bps[2];α=α1,ci=collect(1:2))
	α2 =  0.0
	rb2 = build_2d_tri(2,bps[2];α=α2)
    rbs = TypeSortedCollection([rb1,rb2])
	numberedpoints = RB.number(rbs)
	matrix_sharing = [
		3 1;
		4 2;
	]
	indexedcoords = RB.index(rbs,matrix_sharing)
	# cables
	cables = Int[]
    tensiles = (cables = cables,)
    hub = nothing
	connected = RB.connect(rbs,zeros(Int,0,0))

    cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,hub)
end

function tower2d(;k=100.0,c=0.0,ratio=0.8,ratio1=ratio,slack=true)
    n = 8

	bps_raw = [zeros(2) for i = 1:11]
	for i = 1:4
		bps_raw[2i-1] .= [0.0,0.1*(i-1)]
		bps_raw[2i  ] .= [0.1,0.1*(i-1)]
	end
	for i = 1:3
		bps_raw[8+i] .= [0.05,0.15+0.1*(i-1)]
	end
	bps = SVector{2}.(bps_raw)
	# display(bps)
	α_bar = π/2
	α_tri = 0.0
	rb1 = build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	rb2 = build_2d_bar(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
	# rb1 = build_2d_tri(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# rb2 = build_2d_tri(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
	# rb3 = build_2d_tri(3,bps[3],bps[4],bps[9];α=α_tri)
	rb3 = build_2d_tri(3,bps[3],;α=α_tri)
	rb4 = build_2d_bar(4,bps[9],bps[10];α=α_bar)
	rb5 = build_2d_tri(5,bps[5],bps[6];α=α_tri)
	rb6 = build_2d_bar(6,bps[10],bps[11];α=α_bar)
	rb7 = build_2d_tri(7,bps[7],bps[8],bps[11];α=α_tri)
	rb8 = build_2d_ground(8)
    rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7,rb8))
	# rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7))
	numberedpoints = RB.number(rbs)
	matrix_sharing = [
		3 0 1 0 0 0 0;
		4 0 2 0 0 0 0;
		# 0 3 3 0 0 0 0;
		# 0 4 4 0 0 0 0;
		0 0 0 0 0 3 5;
		0 0 0 0 0 4 6;
	]
	indexedcoords = RB.index(rbs,matrix_sharing)
	#
	restlen4 = ratio1*0.1*√5
	restlen1 = ratio*0.1*√2
	restlen2 = ratio*0.1
	restlen3 = ratio*0.05*√2
    restlens = [
		restlen4,restlen4,
		restlen1,restlen1,
		restlen2,restlen2,
		restlen3,restlen3,
		restlen2,restlen2,
		restlen3,restlen3,
	]
	ncables = length(restlens)
	naux = 2
	ness = 10
    ks = vcat(fill(k,naux),fill(k,ness))
    cs = fill(c,ncables)
    cables = [RB.Cable2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
    acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)
	cnt_matrix_cables = [
		0  0 0 0 -1  0  0  2;
		0  0 0 0 -2  0  0  3;
		1 -2 0 0  0  0  0  0;
	   -2  1 0 0  0  0  0  0;
	    0  0 1 0 -1  0  0  0;
		0  0 2 0 -2  0  0  0;
		0  0 3 0 -1  0  0  0;
		0  0 3 0 -2  0  0  0;
 	    0  0 0 0  1  0 -1  0;
 		0  0 0 0  2  0 -2  0;
 		0  0 0 0  3  0 -1  0;
 		0  0 0 0  3  0 -2  0;
	]
	connected = RB.connect(rbs,cnt_matrix_cables)
	#

	cst1 = RB.PinJoint(RB.End2End(1,RB.ID(rb2,2),RB.ID(rb3,2)))
	cst2 = RB.PinJoint(RB.End2End(2,RB.ID(rb3,3),RB.ID(rb4,1)))
	cst3 = RB.PinJoint(RB.End2End(3,RB.ID(rb4,2),RB.ID(rb5,3)))
	cst4 = RB.PinJoint(RB.End2End(4,RB.ID(rb5,3),RB.ID(rb6,1)))
	jointedmembers = RB.join((cst1,cst2,cst3,cst4),indexedcoords)
	# jointedmembers = RB.join((cst1,cst2,cst3),indexedcoords)

    cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),jointedmembers)
	# cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,hub)
end

function plot_tower2d!(ax,
					st::RB.Structure;
					cablecolor = :dodgerblue,
					barcolor = :darkred,
					markercolor = :darkblue,
					refcolor = :darkgrey,
					isref = false,
					markit = false)
	if isref
		cablecolor = refcolor
		barcolor = refcolor
		markercolor = refcolor
	end
    (;nbodies) = st
    (;tensioned) = st.connectivity
	(;cables) = st.tensiles
	ncables = length(cables)
    ndim = RB.get_num_of_dims(st)
    T = RB.get_numbertype(st)

	rbs = RB.get_bodies(st)
	linesegs_bars = Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
	linesegs_tris =  Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
	ploys_tris = Vector{Vector{Point{ndim,T}}}()
	for rb in rbs
		if body.state.cache.funcs.nmcs isa RB.NCF.LNC2D2P
			push!(
				linesegs_bars,
				(Point(body.state.loci_states[1]),Point(body.state.loci_states[2])),
			)
		elseif body.state.cache.funcs.nmcs isa RB.NCF.LNCMP
		else
			append!(
				linesegs_tris,
				[
					(Point(body.state.loci_states[1]), Point(body.state.loci_states[2])),
					(Point(body.state.loci_states[2]), Point(body.state.loci_states[3])),
					(Point(body.state.loci_states[3]), Point(body.state.loci_states[1])),
				]
			)
			push!(
				ploys_tris,
				[Point(body.state.loci_states[i]) for i = 1:3]
			)
		end
	end
	linesegs_rbs = vcat(linesegs_bars,linesegs_tris)
	linesegments!(ax, linesegs_bars, color = barcolor, linewidth = barwidth)
	poly!(ax,ploys_tris, color = barcolor, transparency = false)
	if markit
		scatter!(ax,[rbs[3].state.rg], color = markercolor, markersize = (2 |> pt2px))
	end

	linesegs_noslack_cables = get_linesegs_cables(st;noslackonly=true)
	linesegs_slack_cables = get_linesegs_cables(st;slackonly=true)
	linesegments!(ax, linesegs_noslack_cables, color = cablecolor, linewidth = cablewidth)
	linesegments!(ax, linesegs_slack_cables, color = cablecolor, linewidth = cablewidth, linestyle = :dash)

end

function plot_one_bar_one_tri!(ax,
		tgob,
		sgi;
		cablecolor = :dodgerblue,
		barcolor = :darkred,
		markercolor = :darkblue,
		refcolor = :darkgrey,
		isref = false,
		markit = false)
	if isref
		cablecolor = refcolor
		barcolor = refcolor
		markercolor = refcolor
	end
	st = tgob[]
	(;nbodies) = st
	ndim = RB.get_num_of_dims(st)
	T = RB.get_numbertype(st)

	linesegs_bars = @lift begin
		ndim = RB.get_num_of_dims($tgob)
		T = RB.get_numbertype($tgob)
		ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
		foreach($tgob.bodies) do body
			if body.state.cache.funcs.nmcs isa RB.NCF.LNC2D2P
				push!(ret,
					Point(body.state.loci_states[1])=>
					Point(body.state.loci_states[2])
				)
			end
		end
		ret
	end

	ploys_tris = @lift begin
		ndim = RB.get_num_of_dims($tgob)
		T = RB.get_numbertype($tgob)
		ret = Vector{Vector{Point{ndim,T}}}()
		foreach($tgob.bodies) do body
			if body.state.cache.funcs.nmcs isa RB.NCF.LNC2D2P
				nothing
			else
				push!(ret,
					[Point(body.state.loci_states[i]) for i = 1:3]
				)
			end
		end
		ret
	end

	# linesegs_rbs = vcat(linesegs_bars,linesegs_tris)
	linesegments!(ax, linesegs_bars, color = barcolor, linewidth = barwidth)
	poly!(ax,ploys_tris, color = barcolor, transparency = false)
	# if markit
	# scatter!(ax,[rbs[3].state.rg], color = markercolor, markersize = (2 |> pt2px))
	# end

	# linesegs_noslack_cables = get_linesegs_cables(st;noslackonly=true)
	# linesegs_slack_cables = get_linesegs_cables(st;slackonly=true)
	# linesegments!(ax, linesegs_noslack_cables, color = cablecolor, linewidth = cablewidth)
	# linesegments!(ax, linesegs_slack_cables, color = cablecolor, linewidth = cablewidth, linestyle = :dash)
end

function scissor_tower(;
		θ = π/4,
		β = 0.0,#-π/4,
		k = nothing,
	)
	b = 1.0
	a = sqrt(2)/2*b
	c = cos(θ)*b
	d = sin(θ)*b
	n = 3
	bps_raw = [zeros(2) for i = 1:2n*3]
	for j = 1:n
		for i = 1:3
			is = 3(j-1)
			if j == n
				bps_raw[is+1] .= [-a/2,d*(j-1)]
				bps_raw[is+2] .= [ a/2,d*(j-1)]
				bps_raw[is+3] .= [ 0.0,d*(j-1)+a/2]
			else
				bps_raw[is+1] .= [ a/2-c,d*(j-1)]
				bps_raw[is+2] .= [-a/2+c,d*(j-1)]
				bps_raw[is+3] .= [ 0.0,d*(j-1)+a/2]
			end
		end
	end
	iss = n*3
	R = RB.rotation_matrix(β)
	x = SVector(0.0,a/2)
	for j = 1:n
		for i = 1:3
			is = 3(j-1)
			if j |> isodd
				bps_raw[iss+is+1] .= bps_raw[iss].+R*(x+[-a/2,d*(j-1)])
				bps_raw[iss+is+2] .= bps_raw[iss].+R*(x+[ a/2,d*(j-1)])
				bps_raw[iss+is+3] .= bps_raw[iss].+R*(x+[ 0.0,d*(j-1)+a/2])
			else
				bps_raw[iss+is+1] .= bps_raw[iss].+R*(x+[ a/2-c,d*(j-1)])
				bps_raw[iss+is+2] .= bps_raw[iss].+R*(x+[-a/2+c,d*(j-1)])
				bps_raw[iss+is+3] .= bps_raw[iss].+R*(x+[ 0.0,d*(j-1)+a/2])
			end
		end
	end
	bps = SVector{2}.(bps_raw)
	display(bps)
	# fig,_,_=scatter(bps)
	# return fig

	# α_tri = 0.0
	# bars = [
	# 	build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# 	build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# 	for i = 1:2
	# ]
	
	bar1 = build_2d_bar(1,bps[1],bps[5];α=θ,ci=collect(1:2))
	bar2 = build_2d_bar(2,bps[2],bps[4];α=π-θ,ci=collect(1:2))
	bar3 = build_2d_bar(3,bps[4],bps[8];α=θ,)
	bar4 = build_2d_bar(4,bps[5],bps[7];α=π-θ,)
	tri1 = build_2d_tri(5,bps[7],bps[8],bps[9];
		α = 0.0,
		b = b/sqrt(2),
		h = b/sqrt(2)/2,
		b1 = b/sqrt(2)/2
	)
	tri2 = build_2d_tri(6,bps[11],bps[10],bps[9];
		α = π+β,
		b = b/sqrt(2),
		h = b/sqrt(2)/2,
		b1 = b/sqrt(2)/2
	)

	# rb2 = build_2d_bar(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
	# rb1 = build_2d_tri(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# rb2 = build_2d_tri(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
	# rb3 = build_2d_tri(3,bps[3],bps[4],bps[9];α=α_tri)
	# rb3 = build_2d_tri(3,bps[3],;α=α_tri)
	# rb4 = build_2d_bar(4,bps[9],bps[10];α=α_bar)
	# rb5 = build_2d_tri(5,bps[5],bps[6];α=α_tri)
	# rb6 = build_2d_bar(6,bps[10],bps[11];α=α_bar)
	# rb7 = build_2d_tri(7,bps[7],bps[8],bps[11];α=α_tri)
	# rb8 = build_2d_ground(8)
    # rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7,rb8))
	rbs = TypeSortedCollection((bar1,bar2,bar3,bar4,tri1,tri2))
	numberedpoints = RB.number(rbs)
	sharing = [
		3 0 0 1 0 0;
		4 0 0 2 0 0;
		0 3 1 0 0 0;
		0 4 2 0 0 0;
		0 0 0 3 1 0;
		0 0 0 4 2 0;
		0 0 3 0 3 0;
		0 0 4 0 4 0;
	]
	indexedcoords = RB.index(rbs,sharing)
	# #
	# restlen4 = ratio1*0.1*√5
	# restlen1 = ratio*0.1*√2
	# restlen2 = ratio*0.1
	# restlen3 = ratio*0.05*√2
    # restlens = [
	# 	restlen4,restlen4,
	# 	restlen1,restlen1,
	# 	restlen2,restlen2,
	# 	restlen3,restlen3,
	# 	restlen2,restlen2,
	# 	restlen3,restlen3,
	# ]
	# ncables = length(restlens)
	# naux = 2
	# ness = 10
    # ks = vcat(fill(k,naux),fill(k,ness))
    # cs = fill(c,ncables)
    # cables = [RB.Cable2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
    # acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
    # tensiles = (cables = cables,)
    # hub = (actuators = acs,)
	cnt_matrix = [
		 1 -2  0  0 0  0;
	   	-2  1  0  0 0  0;
		 0  0  1 -2 0  0;
		 0  0 -2  1 0  0;
		 2 -2  0  0 0  0;
		 0  0  0  0 1 -2;
		 0  0  0  0 2 -1;
	]
	connected = RB.connect(rbs,cnt_matrix)
	# #

	# cst1 = RB.PinJoint(RB.End2End(1,RB.ID(rb2,2),RB.ID(rb3,2)))
	# cst2 = RB.PinJoint(RB.End2End(2,RB.ID(rb3,3),RB.ID(rb4,1)))
	# cst3 = RB.PinJoint(RB.End2End(3,RB.ID(rb4,2),RB.ID(rb5,3)))
	# cst4 = RB.PinJoint(RB.End2End(4,RB.ID(rb5,3),RB.ID(rb6,1)))
	# jointedmembers = RB.join((cst1,cst2,cst3,cst4),indexedcoords)
	# # jointedmembers = RB.join((cst1,cst2,cst3),indexedcoords)

	ncables = size(cnt_matrix,1)
	hncables = ncables
	if k isa Nothing
		cables = [RB.Cable2D(i,0.0,1000.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [RB.Cable2D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	acs = [
		RB.ManualActuator(i,[i],zeros(1))
		for i = 1:hncables
	]
	hub = (actuators = acs,)
	tensiles = (cables = cables,)
	connected = RB.connect(rbs,cnt_matrix)


    cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)
	# # cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
    st = RB.Structure(rbs,tensiles,cnt)
    RB.Robot(st,hub)
end

function trapezoid(id,ri,ro=ri,
		α = 0.0, 
		ci = Int[],
		constraints_indices = collect(1:3)

	)
	movable = true
	if isempty(ci)
		constrained = false
	else
		constrained = true
	end
	# unconstrained_indices = collect(1:6)
	b2 = b - b1
	r̄p1 = SVector{2}([ -r, 0.0])
	r̄p2 = SVector{2}([  r, 0.0])
	r̄p3 = SVector{2}([ -b/2,  h])
	r̄p4 = SVector{2}([  b/2,  h])
	loci = [r̄p1,r̄p2,r̄p3,r̄p4]
	ρ = 1.0
	m1 = ρ*2r
	m2 = ρ*h
	m3 = ρ*b
	m = m1+m2+m3
	r̄gy = (m2*h/2+m3*h)/m
	mass_locus  = SVector{2}([  0, r̄gy])
	I1_c = m1*(2r)^2/12
	I2_c = m2*h^2/12
	I3_c = m3*b^2/12
	I1_g = I1_c + m1*r̄gy^2
	I2_g = I2_c + m2*(h/2-r̄gy)^2
	I3_g = I3_c + m3*(h-r̄gy)^2
	I_g = I1_g + I2_g + I3_g
	Īg = SMatrix{2,2}(
		[
			0.99I_g 0.0;
			0.0 0.01I_g
		]
	)

	prop = RB.RigidBodyProperty(id,movable,m,Īg,
				mass_locus,loci;constrained=constrained
				)
	ω = 0.0
	# ṙo = [0.0,0.0001]
	ṙo = zero(ro)
	nmcs = RB.NCF.NC1P2V(ri,ro,α)
	state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,constraints_indices)
	body = RB.RigidBody(prop,state)

end

function birdneck(
		θ = π/4,
		β = 0.0,#-π/4,
		k = nothing,
	)
	L0 = 0.08
	L = 0.8L0
	b = 0.8L0/2
	r = r0 = 0.04
	h = 0.2*sqrt(L0^2-(L0/2)^2)/2

	b = 1.0
	a = sqrt(2)/2*b
	c = cos(θ)*b
	d = sin(θ)*b
	n = 3
	bps_raw = [zeros(2) for i = 1:2n*3]
	for j = 1:n
		for i = 1:3
			is = 3(j-1)
			if j |> isodd
				bps_raw[is+1] .= [-a/2,d*(j-1)]
				bps_raw[is+2] .= [ a/2,d*(j-1)]
				bps_raw[is+3] .= [ 0.0,d*(j-1)+a/2]
			else
				bps_raw[is+1] .= [ a/2-c,d*(j-1)]
				bps_raw[is+2] .= [-a/2+c,d*(j-1)]
				bps_raw[is+3] .= [ 0.0,d*(j-1)+a/2]
			end
		end
	end
	iss = n*3
	R = RB.rotation_matrix(β)
	x = SVector(0.0,a/2)
	for j = 1:n
		for i = 1:3
			is = 3(j-1)
			if j |> isodd
				bps_raw[iss+is+1] .= bps_raw[iss].+R*(x+[-a/2,d*(j-1)])
				bps_raw[iss+is+2] .= bps_raw[iss].+R*(x+[ a/2,d*(j-1)])
				bps_raw[iss+is+3] .= bps_raw[iss].+R*(x+[ 0.0,d*(j-1)+a/2])
			else
				bps_raw[iss+is+1] .= bps_raw[iss].+R*(x+[ a/2-c,d*(j-1)])
				bps_raw[iss+is+2] .= bps_raw[iss].+R*(x+[-a/2+c,d*(j-1)])
				bps_raw[iss+is+3] .= bps_raw[iss].+R*(x+[ 0.0,d*(j-1)+a/2])
			end
		end
	end
	bps = SVector{2}.(bps_raw)
	display(bps)
	# fig,_,_=scatter(bps)
	# return fig

	# α_tri = 0.0
	# bars = [
	# 	build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# 	build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# 	for i = 1:2
	# ]

	bar1 = build_2d_bar(1,bps[1],bps[5];α=θ,ci=collect(1:2))
	bar2 = build_2d_bar(2,bps[2],bps[4];α=π-θ,ci=collect(1:2))
	bar3 = build_2d_bar(3,bps[4],bps[8];α=θ,)
	bar4 = build_2d_bar(4,bps[5],bps[7];α=π-θ,)
	tri1 = build_2d_tri(5,bps[7],bps[8],bps[9];
		α = 0.0,
		b = b/sqrt(2),
		h = b/sqrt(2)/2,
		b1 = b/sqrt(2)/2
	)
	tri2 = build_2d_tri(6,bps[11],bps[10],bps[9];
		α = π+β,
		b = b/sqrt(2),
		h = b/sqrt(2)/2,
		b1 = b/sqrt(2)/2
	)

	# rb2 = build_2d_bar(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
	# rb1 = build_2d_tri(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
	# rb2 = build_2d_tri(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
	# rb3 = build_2d_tri(3,bps[3],bps[4],bps[9];α=α_tri)
	# rb3 = build_2d_tri(3,bps[3],;α=α_tri)
	# rb4 = build_2d_bar(4,bps[9],bps[10];α=α_bar)
	# rb5 = build_2d_tri(5,bps[5],bps[6];α=α_tri)
	# rb6 = build_2d_bar(6,bps[10],bps[11];α=α_bar)
	# rb7 = build_2d_tri(7,bps[7],bps[8],bps[11];α=α_tri)
	# rb8 = build_2d_ground(8)
	# rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7,rb8))
	rbs = TypeSortedCollection((bar1,bar2,bar3,bar4,tri1,tri2))
	numberedpoints = RB.number(rbs)
	sharing = [
		3 0 0 1 0 0;
		4 0 0 2 0 0;
		0 3 1 0 0 0;
		0 4 2 0 0 0;
		0 0 0 3 1 0;
		0 0 0 4 2 0;
		0 0 3 0 3 0;
		0 0 4 0 4 0;
	]
	indexedcoords = RB.index(rbs,sharing)
	# #
	# restlen4 = ratio1*0.1*√5
	# restlen1 = ratio*0.1*√2
	# restlen2 = ratio*0.1
	# restlen3 = ratio*0.05*√2
	# restlens = [
	# 	restlen4,restlen4,
	# 	restlen1,restlen1,
	# 	restlen2,restlen2,
	# 	restlen3,restlen3,
	# 	restlen2,restlen2,
	# 	restlen3,restlen3,
	# ]
	# ncables = length(restlens)
	# naux = 2
	# ness = 10
	# ks = vcat(fill(k,naux),fill(k,ness))
	# cs = fill(c,ncables)
	# cables = [RB.Cable2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
	# acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
	# tensiles = (cables = cables,)
	# hub = (actuators = acs,)
	cnt_matrix = [
		1  -2  0  0 0  0;
		-2  1  0  0 0  0;
		0   0  1 -2 0  0;
		0   0 -2  1 0  0;
		2  -2  0  0 0  0;
		0   0  0  0 1 -2;
		0   0  0  0 2 -1;
	]
	connected = RB.connect(rbs,cnt_matrix)
	# #

	# cst1 = RB.PinJoint(RB.End2End(1,RB.ID(rb2,2),RB.ID(rb3,2)))
	# cst2 = RB.PinJoint(RB.End2End(2,RB.ID(rb3,3),RB.ID(rb4,1)))
	# cst3 = RB.PinJoint(RB.End2End(3,RB.ID(rb4,2),RB.ID(rb5,3)))
	# cst4 = RB.PinJoint(RB.End2End(4,RB.ID(rb5,3),RB.ID(rb6,1)))
	# jointedmembers = RB.join((cst1,cst2,cst3,cst4),indexedcoords)
	# # jointedmembers = RB.join((cst1,cst2,cst3),indexedcoords)

	ncables = size(cnt_matrix,1)
	hncables = ncables
	if k isa Nothing
		cables = [RB.Cable2D(i,0.0,1000.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [RB.Cable2D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	acs = [
		RB.ManualActuator(i,[i],zeros(1))
		for i = 1:hncables
	]
	hub = (actuators = acs,)
	tensiles = (cables = cables,)
	connected = RB.connect(rbs,cnt_matrix)


	cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)
	# # cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
	st = RB.Structure(rbs,tensiles,cnt)
	RB.Robot(st,hub)
end