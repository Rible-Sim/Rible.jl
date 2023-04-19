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
	r̄g  = SVector{2}([ b/2,0])
	r̄p1 = SVector{2}([ 0.0,0])
	r̄p2 = SVector{2}([   b,0])
	r̄ps = [r̄p1,r̄p2]
	m = 2.6934977798E-02
	Īg = SMatrix{2,2}(
		[
			1/12*m*b^2 0.0;
			       0.0 0.0;
		]
	)
	@show b,m,Īg[1]
	prop = TR.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
				)
	ω = 0.0
	ro = ri
	ṙo = zero(ro)
	lncs,q0,q̇0 = TR.NCF.NC2D2P(ri,rj,ro,α,ṙo,ω)
	state = TR.RigidBodyState(prop,lncs,ro,α,ṙo,ω,ci)
	rb = TR.RigidBody(prop,state)
end

function build_2d_tri(id,ri,rj=nothing,rk=nothing;
		α = 0.0, 
		ci = Int[], 
		Φi = collect(1:3),		
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
	# uci = collect(1:6)
	b2 = b - b1
	r̄g  = SVector{2}([ (b1+b)/3, h/3])
	r̄p1 = SVector{2}([      0.0, 0.0])
	r̄p2 = SVector{2}([        b, 0.0])
	r̄p3 = SVector{2}([        b1,  h])
	r̄ps = [r̄p1,r̄p2,r̄p3]
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
	@show norm(r̄g),m,Īg,tr(Īg)
	@show r̄g,atan(r̄g[2],r̄g[1])
	prop = TR.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
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
		lncs,q0,q̇0 = TR.NCF.NC1P2V(ri,ro,α,ṙo,ω)
	elseif rk isa Nothing
		lncs,q0,q̇0 = TR.NCF.NC2P1V(ri,rj,ro,α,ṙo,ω)
	else
		lncs,q0,q̇0 = TR.NCF.NC3P(ri,rj,rk,ro,α,ṙo,ω)
	end
	state = TR.RigidBodyState(prop,lncs,ro,α,ṙo,ω,ci,Φi)
	rb = TR.RigidBody(prop,state)
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
	r̄g  = SVector{2}([ 0.0, 0.0])
	r̄p1 = SVector{2}([ 0.0, 0.0])
	r̄p2 = SVector{2}([ -b1, 0.0])
	r̄p3 = SVector{2}([ 2b1, 0.0])
	r̄ps = [r̄p1,r̄p2,r̄p3]
	prop = TR.RigidBodyProperty(id,movable,m,Ia,
				r̄g,r̄ps;constrained=constrained
				)
	α = 0.0; ω = 0.0
	R = TR.rotation_matrix(α)
	ro = zeros(2)
	ṙo = zero(ro)
	rps = Ref(ro) .+ Ref(R).*r̄ps
	lncs,q0 = TR.NCF.NCMP(rps,ro,R,ṙo,ω)
	# cf = TR.NCF.CoordinateFunctions(lncs,q0,ci,uci)
	# @show typeof(lncs)
	nq = length(q0)
	ci = collect(1:nq)
	uci = Int[]
	Φi = Int[]
	state = TR.RigidBodyState(prop,lncs,ro,α,ṙo,ω,ci,Φi)
	rb = TR.RigidBody(prop,state)
end

function two_tri(;k=100.0,c=0.0,ratio=0.8)
    n = 8

	bps_raw = [zeros(2) for i = 1:11]
	for i = 1:2
		bps_raw[2i-1] .= [0.0,0.1*(i-1)]
		bps_raw[2i  ] .= [0.1,0.1*(i-1)]
	end
	for i = 1:3
		bps_raw[4+i] .= [0.05,0.05+0.1*(i-1)]
	end
	bps = SVector{2}.(bps_raw)
	display(bps)
	α_tri1 = -π/2
	α_tri2 =  π/2
	rb1 = build_2d_tri(1,bps[3],bps[1],bps[5];α=α_tri1,ci=collect(1:6),Φi=Int[])
	rb2 = build_2d_tri(2,bps[2],bps[4],bps[5];α=α_tri2,ci=collect(5:6))

	rbs = TypeSortedCollection((rb1,rb2))
	numberedpoints = TR.number(rbs)
	matrix_sharing = [
		5 5;
		6 6;
	]
	indexedcoords = TR.index(rbs,matrix_sharing)
	#
	ncables = 2
	restlen1 = 0.05
    restlens = [
		restlen1,restlen1
	]
    ks = fill(k,ncables)
    cs = fill(c,ncables)
    cables = [TR.Cable2D(i,restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    acs = [TR.ManualActuator(TR.SimpleRegistor(i,restlens[i])) for i = 1:ncables]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)
	cnt_matrix_cables = [
		1 -2 0 0  0  0  0  0;
	   -2  1 0 0  0  0  0  0;
	]
	connected = TR.connect(rbs,cnt_matrix_cables)
	#

	# jointedmembers = TR.join((cst1,),indexedcoords)

    cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
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
	numberedpoints = TR.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = TR.index(rbs,matrix_sharing)
	#
	ss = Int[]
    tensiles = (cables = ss,)
    hub = nothing
	#
	connections = TR.connect(rbs,zeros(Int,0,0))

	jointedmembers = TR.unjoin()

    cnt = TR.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers)

    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
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
	numberedpoints = TR.number(rbs)
	matrix_sharing = [
		3 1;
		4 2;
	]
	indexedcoords = TR.index(rbs,matrix_sharing)
	# cables
	cables = Int[]
    tensiles = (cables = cables,)
    hub = nothing
	connected = TR.connect(rbs,zeros(Int,0,0))

    cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
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
	numberedpoints = TR.number(rbs)
	matrix_sharing = [
		3 0 1 0 0 0 0;
		4 0 2 0 0 0 0;
		# 0 3 3 0 0 0 0;
		# 0 4 4 0 0 0 0;
		0 0 0 0 0 3 5;
		0 0 0 0 0 4 6;
	]
	indexedcoords = TR.index(rbs,matrix_sharing)
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
    cables = [TR.Cable2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
    acs = [TR.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
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
	connected = TR.connect(rbs,cnt_matrix_cables)
	#

	cst1 = TR.PinJoint(TR.End2End(1,TR.ID(rb2,2),TR.ID(rb3,2)))
	cst2 = TR.PinJoint(TR.End2End(2,TR.ID(rb3,3),TR.ID(rb4,1)))
	cst3 = TR.PinJoint(TR.End2End(3,TR.ID(rb4,2),TR.ID(rb5,3)))
	cst4 = TR.PinJoint(TR.End2End(4,TR.ID(rb5,3),TR.ID(rb6,1)))
	jointedmembers = TR.join((cst1,cst2,cst3,cst4),indexedcoords)
	# jointedmembers = TR.join((cst1,cst2,cst3),indexedcoords)

    cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),jointedmembers)
	# cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
end

function plot_tower2d!(ax,
					tg::TR.TensegrityStructure;
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
    (;nbodies) = tg
    (;tensioned) = tg.connectivity
	(;cables) = tg.tensiles
	ncables = length(cables)
    ndim = TR.get_ndim(tg)
    T = TR.get_numbertype(tg)

	rbs = TR.get_rigidbodies(tg)
	linesegs_bars = Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
	linesegs_tris =  Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
	ploys_tris = Vector{Vector{Point{ndim,T}}}()
	for rb in rbs
		if rb.state.cache.funcs.nmcs isa TR.NCF.LNC2D2P
			push!(
				linesegs_bars,
				(Point(rb.state.rps[1]),Point(rb.state.rps[2])),
			)
		elseif rb.state.cache.funcs.nmcs isa TR.NCF.LNCMP
		else
			append!(
				linesegs_tris,
				[
					(Point(rb.state.rps[1]), Point(rb.state.rps[2])),
					(Point(rb.state.rps[2]), Point(rb.state.rps[3])),
					(Point(rb.state.rps[3]), Point(rb.state.rps[1])),
				]
			)
			push!(
				ploys_tris,
				[Point(rb.state.rps[i]) for i = 1:3]
			)
		end
	end
	linesegs_rbs = vcat(linesegs_bars,linesegs_tris)
	linesegments!(ax, linesegs_bars, color = barcolor, linewidth = barwidth)
	poly!(ax,ploys_tris, color = barcolor, transparency = false)
	if markit
		scatter!(ax,[rbs[3].state.rg], color = markercolor, markersize = (2 |> pt2px))
	end

	linesegs_noslack_cables = get_linesegs_cables(tg;noslackonly=true)
	linesegs_slack_cables = get_linesegs_cables(tg;slackonly=true)
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
	tg = tgob[]
	(;nbodies) = tg
	ndim = TR.get_ndim(tg)
	T = TR.get_numbertype(tg)

	linesegs_bars = @lift begin
		ndim = TR.get_ndim($tgob)
		T = TR.get_numbertype($tgob)
		ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
		foreach($tgob.bodies) do rb
			if rb.state.cache.funcs.nmcs isa TR.NCF.LNC2D2P
				push!(ret,
					Point(rb.state.rps[1])=>
					Point(rb.state.rps[2])
				)
			end
		end
		ret
	end

	ploys_tris = @lift begin
		ndim = TR.get_ndim($tgob)
		T = TR.get_numbertype($tgob)
		ret = Vector{Vector{Point{ndim,T}}}()
		foreach($tgob.bodies) do rb
			if rb.state.cache.funcs.nmcs isa TR.NCF.LNC2D2P
				nothing
			else
				push!(ret,
					[Point(rb.state.rps[i]) for i = 1:3]
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

	# linesegs_noslack_cables = get_linesegs_cables(tg;noslackonly=true)
	# linesegs_slack_cables = get_linesegs_cables(tg;slackonly=true)
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
	R = TR.rotation_matrix(β)
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
	numberedpoints = TR.number(rbs)
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
	indexedcoords = TR.index(rbs,sharing)
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
    # cables = [TR.Cable2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
    # acs = [TR.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
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
	connected = TR.connect(rbs,cnt_matrix)
	# #

	# cst1 = TR.PinJoint(TR.End2End(1,TR.ID(rb2,2),TR.ID(rb3,2)))
	# cst2 = TR.PinJoint(TR.End2End(2,TR.ID(rb3,3),TR.ID(rb4,1)))
	# cst3 = TR.PinJoint(TR.End2End(3,TR.ID(rb4,2),TR.ID(rb5,3)))
	# cst4 = TR.PinJoint(TR.End2End(4,TR.ID(rb5,3),TR.ID(rb6,1)))
	# jointedmembers = TR.join((cst1,cst2,cst3,cst4),indexedcoords)
	# # jointedmembers = TR.join((cst1,cst2,cst3),indexedcoords)

	ncables = size(cnt_matrix,1)
	hncables = ncables
	if k isa Nothing
		cables = [TR.Cable2D(i,0.0,1000.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [TR.Cable2D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	acs = [
		TR.ManualActuator(i,[i],zeros(1))
		for i = 1:hncables
	]
	hub = (actuators = acs,)
	tensiles = (cables = cables,)
	connected = TR.connect(rbs,cnt_matrix)


    cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)
	# # cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    TR.TensegrityRobot(tg,hub)
end

function trapezoid(id,ri,ro=ri,
		α = 0.0, 
		ci = Int[],
		Φi = collect(1:3)

	)
	movable = true
	if isempty(ci)
		constrained = false
	else
		constrained = true
	end
	# uci = collect(1:6)
	b2 = b - b1
	r̄p1 = SVector{2}([ -r, 0.0])
	r̄p2 = SVector{2}([  r, 0.0])
	r̄p3 = SVector{2}([ -b/2,  h])
	r̄p4 = SVector{2}([  b/2,  h])
	r̄ps = [r̄p1,r̄p2,r̄p3,r̄p4]
	ρ = 1.0
	m1 = ρ*2r
	m2 = ρ*h
	m3 = ρ*b
	m = m1+m2+m3
	r̄gy = (m2*h/2+m3*h)/m
	r̄g  = SVector{2}([  0, r̄gy])
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

	prop = TR.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
				)
	ω = 0.0
	# ṙo = [0.0,0.0001]
	ṙo = zero(ro)
	lncs,q0,q̇0 = TR.NCF.NC1P2V(ri,ro,α,ṙo,ω)
	state = TR.RigidBodyState(prop,lncs,ro,α,ṙo,ω,ci,Φi)
	rb = TR.RigidBody(prop,state)

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
	R = TR.rotation_matrix(β)
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
	numberedpoints = TR.number(rbs)
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
	indexedcoords = TR.index(rbs,sharing)
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
	# cables = [TR.Cable2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
	# acs = [TR.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
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
	connected = TR.connect(rbs,cnt_matrix)
	# #

	# cst1 = TR.PinJoint(TR.End2End(1,TR.ID(rb2,2),TR.ID(rb3,2)))
	# cst2 = TR.PinJoint(TR.End2End(2,TR.ID(rb3,3),TR.ID(rb4,1)))
	# cst3 = TR.PinJoint(TR.End2End(3,TR.ID(rb4,2),TR.ID(rb5,3)))
	# cst4 = TR.PinJoint(TR.End2End(4,TR.ID(rb5,3),TR.ID(rb6,1)))
	# jointedmembers = TR.join((cst1,cst2,cst3,cst4),indexedcoords)
	# # jointedmembers = TR.join((cst1,cst2,cst3),indexedcoords)

	ncables = size(cnt_matrix,1)
	hncables = ncables
	if k isa Nothing
		cables = [TR.Cable2D(i,0.0,1000.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [TR.Cable2D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	acs = [
		TR.ManualActuator(i,[i],zeros(1))
		for i = 1:hncables
	]
	hub = (actuators = acs,)
	tensiles = (cables = cables,)
	connected = TR.connect(rbs,cnt_matrix)


	cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)
	# # cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	TR.TensegrityRobot(tg,hub)
end