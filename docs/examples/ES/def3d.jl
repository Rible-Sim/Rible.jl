function make_3d_bar(id,ri,rj;ci = Int[])
	# @show id,ri,rj
	movable = true
	if ci == Int[]
		constrained = false
	else
		constrained = true
	end
	u = rj - ri
	b = norm(u)
	û = u./b
	v̂,ŵ = TR.NCF.HouseholderOrthogonalization(û)
	R = SMatrix{3,3}(hcat(û,v̂,ŵ))
	m = 0.080#19074114753529
	Īg = SMatrix{3,3}([1/12*m*b^2 0.0 0.0;
					   		  0.0 0.0 0.0;
					          0.0 0.0 0.0])
	@show b,m,Īg[1]
	r̄g  = SVector{3}([ 0.0,0,0])
	r̄p1 = SVector{3}([-b/2,0,0])
	r̄p2 = SVector{3}([ b/2,0,0])
	r̄ps = [r̄p1,r̄p2]
	prop = TR.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
				)
	# @show prop.inertia
	ro = (ri+rj)./2
	ṙo = zero(ro)
	ω = zero(ro)
	nmcs,q0,q̇0 = TR.NCF.NC3D2P(ri,rj,ro,R,ṙo,ω)
	# @show ri,rj,q0
	# cf = TR.NCF.CoordinateFunctions(nmcs,q0,ci,uci)
	# @show typeof(nmcs)
	state = TR.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci)
	radius = norm(r̄p2-r̄p1)/60
	@show radius
    barmesh = endpoints2mesh(r̄p1,r̄p2;radius,color=:darkred)
	rb = TR.RigidBody(prop,state,barmesh)
	# @show rb.state.cache.M
	rb
end

function make_3d_tri(
		id,
		r̄ps,ro,
		R,ri,rj=nothing,rk=nothing,rl=nothing;
		movable = true,
		constrained = false,
		ci = Int[],
		Φi = collect(1:6),
	)
	# uci = collect(1:6)
	m = 0.2999233976
	Īg = SMatrix{3,3}(
		Matrix(Diagonal([
			7.6639282053E-04,
			7.6638139752E-04,
			1.2464720496E-03
		])
		)
	)
	if id == 7
		r̄g = [0,0,0.1562442983-0.1*√2]
	else
		r̄g = [0,0,0.1562442983-0.1*√2-0.1*√2/2]
	end
	@show m,diag(Īg),r̄g
	prop = TR.RigidBodyProperty(id,
		movable,m,Īg,
		r̄g,r̄ps;
		constrained=constrained
	)
	ṙo = zero(ro)
	ω = zero(ro)
	u = R*(r̄ps[2] - r̄ps[1])
	v = R*(r̄ps[3] - r̄ps[1])
	w = R*(r̄ps[4] - r̄ps[1])
	if rj isa Nothing
		nmcs,_ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
	elseif rk isa Nothing
		nmcs,_ = TR.NCF.NC2P2V(ri,rj,ro,R,ṙo,ω)
	elseif rl isa Nothing
		nmcs,_ = TR.NCF.NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
	else
		nmcs,_ = TR.NCF.NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
	end
	# cf = TR.NCF.CoordinateFunctions(nmcs,q0,ci,uci)
	# @show typeof(nmcs)
	radius = norm(r̄ps[2]-r̄ps[1])/32
	@show radius
	trimesh = GB.merge(
		[
			endpoints2mesh(r̄ps[i],r̄ps[j];
			radius=0.0035,color=:darkorchid4)
			for (i,j) in [
				[1,2],[1,3],[1,4],
				[2,3],[3,4],[4,2]
			]
		]
	)
	state = TR.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci,Φi)
	rb = TR.RigidBody(prop,state,trimesh)
end

function tower3d(;
		k=500.0,k1=1000.0,c=0.0,
		d = 0.1*√2/2, r2 = 0.11,		
		α = π/6,
		ijkl=1
	)
	r1 = 0.1
	b = 0.22
	h = 0.1*√2/2
	γ = acos((d^2+r1^2+r2^2-b^2)/(2r1*r2))
	θ =  γ - 2π/3
	@show rad2deg.([γ,θ])
	deg120 = deg2rad(120)
	r̄pss = [
		SVector{3}.([
		       [           0,           0,  h],
		    r.*[           1,           0,  0],
		    r.*[cos( deg120), sin( deg120), 0],
		    r.*[cos(-deg120), sin(-deg120), 0]
		])
		for r in [r1,r2,r1,r1,r1,r1]
	]
	r̄ps = r̄pss[1]
	for i = 4:6
		r̄pss[i] .-= Ref(r̄pss[i][1])
	end
	# tetra_points = [GP.Point3D(tetra_coord...) for tetra_coord in tetra_coords]
	# tetra_centroid = GP.centroid(GP.Primitive(tetra_points...))
	# r̄ps = SVector{3}.(
	# 	vcat(
	# 		tetra_coords,
	# 		# [[GP.getx(tetra_centroid),GP.gety(tetra_centroid),GP.getz(tetra_centroid)]]
	# 	)
	# )
	ro_by_rbid = [
		SVector(0.0,0.0, 0),
		SVector(0.0,0.0, d),
		SVector(0.0,0.0,2d),
		SVector(0.0,0.0,2d-0.5h),
		SVector(0.0,0.0,2d-0.5h+2h),
		SVector(0.0,-0.5h*sin(α),2d-0.5h+2h+0.5h*cos(α))
	]
	R_by_rbid = [
		SMatrix(RotX(0.0)),
		SMatrix(RotZ(θ)),
		SMatrix(RotZ(2θ)),
		SMatrix(RotZ(2θ)),
		SMatrix(RotZ(2θ)*RotX(π)*RotY(α)),
		SMatrix(RotZ(2θ)*RotX(π)*RotY(α)),
	]
	rirjrkrl_by_rbid = [
		Ref(ro_by_rbid[i]) .+ Ref(R_by_rbid[i]).*r̄pss[i] for i = 1:6
	]
	# @show rirjrkrl_by_rbid[1]
	cycle3 = [2,3,4,2]
	# @show rirjrkrl_by_rbid[1]
	rb1_to_3 = [
		  make_3d_bar(i,rirjrkrl_by_rbid[1][cycle3[i  ]],
					    rirjrkrl_by_rbid[2][cycle3[i+1]]; ci = [1,2,3]) for i = 1:3
	]

	# rb4 = make_3d_tri(4,r̄ps,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1])
	# rb4 = make_3d_tri(4,r̄ps,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1:2]...)
	# rb4 = make_3d_tri(4,r̄ps,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1:3]...)
	# rb4 = make_3d_tri(4,r̄ps[2:4],ro_by_rbid[2],R_by_rbid[2],rirjrkrl_by_rbid[2][2:4]...)
	rb4_to_6 = [
		make_3d_bar(i+3,rirjrkrl_by_rbid[2][cycle3[i  ]],
					    rirjrkrl_by_rbid[3][cycle3[i+1]]; ci = Int[]) for i = 1:3
	]

	rb7 = make_3d_tri( 7,r̄ps,ro_by_rbid[3],R_by_rbid[3],rirjrkrl_by_rbid[3][1:4]...)

	rb8 = make_3d_tri( 8,r̄pss[4],ro_by_rbid[4],R_by_rbid[4],rirjrkrl_by_rbid[4][1])

	rb9 = make_3d_tri( 9,r̄pss[5],ro_by_rbid[5],R_by_rbid[5],rirjrkrl_by_rbid[5][1])

	rb10 = make_3d_tri(10,r̄pss[6],ro_by_rbid[6],R_by_rbid[6],rirjrkrl_by_rbid[6][1:ijkl]...)

	rbs = TypeSortedCollection(vcat(rb1_to_3,rb4_to_6,[rb7,rb8,rb9,rb10]))
	numberedpoints = TR.number(rbs)
	matrix_sharing = [
		4 0 0 0 1 0 0 0 0 0;
		5 0 0 0 2 0 0 0 0 0;
		6 0 0 0 3 0 0 0 0 0;
		0 4 0 0 0 1 0 0 0 0;
		0 5 0 0 0 2 0 0 0 0;
		0 6 0 0 0 3 0 0 0 0;
		0 0 4 1 0 0 0 0 0 0;
		0 0 5 2 0 0 0 0 0 0;
		0 0 6 3 0 0 0 0 0 0;

		# # #
		0 0 0 0 0 4 4 0 0 0;
		0 0 0 0 0 5 5 0 0 0;
		0 0 0 0 0 6 6 0 0 0;
		0 0 0 4 0 0 7 0 0 0;
		0 0 0 5 0 0 8 0 0 0;
		0 0 0 6 0 0 9 0 0 0;
		0 0 0 0 4 0 10 0 0 0;
		0 0 0 0 5 0 11 0 0 0;
		0 0 0 0 6 0 12 0 0 0;
		# #
		0 0 0 0 0 0 0 1 1 0;
		0 0 0 0 0 0 0 2 2 0;
		0 0 0 0 0 0 0 3 3 0;
	]
	indexedcoords = TR.index(rbs,matrix_sharing)
	# indexedcoords = TR.index(rbs)
	# #
	ndcables = 9
	nocables = 6
	nvcables = 3
	nocables = 6
	ncables = ndcables + nocables + nvcables + nocables
	restlend = 0.05
	restleno = 0.01
	restlenv = 0.05
	restleno = 0.01
    restlens = vcat(
		fill(restlend,ndcables),
		fill(restleno,nocables),
		fill(restlenv,nvcables),
		fill(restleno,nocables),
	)
    ks = vcat(
		fill(k1,ndcables),
		fill(k,nocables),
		fill(k,nvcables),
		fill(k,nocables),
	)
    cs = fill(c,ncables)
    cables = [TR.Cable3D(i,restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    acs = [TR.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)
	cnt_matrix_cables = [
		# triplex 1
		1 0 0 -1  0  0  0 0 0 0;
		0 1 0  0 -1  0  0 0 0 0;
		0 0 1  0  0 -1  0 0 0 0;
		# triplex intersect
		0 0 0  1 -1  0 0 0 0 0;
		0 0 0  0  1 -1 0 0 0 0;
		0 0 0 -1  0  1 0 0 0 0;
		# triplex 2
		0 0 0  1  0  0 -2 0 0 0;
		0 0 0  0  1  0 -3 0 0 0;
		0 0 0  0  0  1 -4 0 0 0;
		# Inner
		0 0 0 0 0 0 1 -2 0 0;
		0 0 0 0 0 0 1 -3 0 0;
		0 0 0 0 0 0 1 -4 0 0;
		# Outer
		0 0 0 0 0 0 2 -2 0 0;
		0 0 0 0 0 0 3 -3 0 0;
		0 0 0 0 0 0 4 -4 0 0;
		#
		0 0 0 0 0 0 0 2 -2 0;
		0 0 0 0 0 0 0 3 -4 0;
		0 0 0 0 0 0 0 4 -3 0;
		# Outer
		0 0 0 0 0 0 0 0 2 -2;
		0 0 0 0 0 0 0 0 3 -3;
		0 0 0 0 0 0 0 0 4 -4;
		# Inner
		0 0 0 0 0 0 0 0 2 -1;
		0 0 0 0 0 0 0 0 3 -1;
		0 0 0 0 0 0 0 0 4 -1;
	]
	connected = TR.connect(rbs,cnt_matrix_cables)
	# #
	#
	# cst1 = TR.PinJoint(TR.End2End(1,TR.ID(rb1_to_3[1],2),TR.ID(rb4,1)))
	# cst2 = TR.PinJoint(TR.End2End(2,TR.ID(rb1_to_3[2],2),TR.ID(rb4,2)))
	# cst3 = TR.PinJoint(TR.End2End(3,TR.ID(rb1_to_3[3],2),TR.ID(rb4,3)))
	# jointedmembers = TR.join((cst1,cst2,cst3),indexedcoords)
	#
    cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
end

function plot_decompose_tower3d(bot0,bot1;
							cablecolor = :dodgerblue,
							barcolor = :darkred,)
	fig_width = 1columnwidth
	fig_height = 0.8columnwidth
	fig = Figure(;resolution=(fig_width,fig_height),figure_padding=(-50,0,0,0))
	rows = (1,1,2,2)
	cols = (1,2,1,2)
	axs = [
		Axis3(fig[row,col],aspect=:data,
				# xlabel=L"x~(\mathrm{m})",
				# ylabel=L"y~(\mathrm{m})",
				# zlabel=L"z~(\mathrm{m})"
			)
		for (row,col) in zip(rows,cols)
	]
	linesegs_noslack_cables, linesegs_slack_cables = get_linesegs_cables(bot0.tg)
	linesegs_cables = vcat(linesegs_noslack_cables, linesegs_slack_cables)
    linesegments!(axs[1], linesegs_cables[10:15], color = cablecolor, linewidth = cablewidth)
    linesegments!(axs[2], linesegs_cables[16:18], color = cablecolor, linewidth = cablewidth)
    linesegments!(axs[3], linesegs_cables[ 1: 9], color = cablecolor, linewidth = cablewidth)

	linesegs_bars = get_linesegs_bars(bot0.tg)
	linesegments!(axs[1], linesegs_bars[ 6+1: 6+12], color = barcolor, linewidth = barwidth)
	linesegments!(axs[2], linesegs_bars[12+1:12+12], color = barcolor, linewidth = barwidth)
	linesegments!(axs[3], linesegs_bars[   1:    9], color = barcolor, linewidth = barwidth)

	linesegments!(axs[4], reduce(vcat,get_linesegs_cables(bot1.tg))[ 1: 9], color = cablecolor, linewidth = cablewidth)
	linesegments!(axs[4],              get_linesegs_bars(bot1.tg)[   1: 9], color = barcolor, linewidth = barwidth)

	# zlims!(ax,0,0.5)
	# ax.xlabel = "x (m)"
	# ax.ylabel = "y (m)"
	# ax.zlabel = "z (m)"
	# ax.yticks = ax.xticks = [-0.1,0,0.1]
	# ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
	# ax.ylabeloffset = ax.xlabeloffset = 40
	# ax.alignmode =  Mixed(;left = 25, right = -25)

	for (iax,ax) in enumerate(axs)
		if iax in [1,2,3,4]
			hidedecorations!(ax)
			hidespines!(ax)
		end
		if iax in [6]
			ax.zticklabelsvisible = false
			ax.zticksvisible = false
			ax.zlabelvisible = false
		end
		xlims!(ax,-0.125,0.125)
		ylims!(ax,-0.125,0.125)
		if iax in [3,4]
			zlims!(ax,0,0.3)
		end
		ax.azimuth = 7.595530633326982
		ax.elevation = 0.12269908169872408
	end
	leftpads = [150,150,100,100]
	toppads = [20,20,-0,-0]
	subtitles = [
		"Class-1 Module",
		"Class-2 Module",
		"2-Stage Triplex (Folded)",
		"2-Stage Triplex (Unfolded)"
	]
	for (ilabel,label) in enumerate(alphabet[1:4])
		Label(fig.layout[rows[ilabel], cols[ilabel], Bottom()],
			("($label) "),
			font = "CMU Serif Bold",
			textsize = fontsize,
			padding = (leftpads[ilabel]-1.5fontsize, 0, 0, toppads[ilabel]),
			halign = :left,
			valign = :top
		)
	end
	for (ilabel,label) in enumerate(alphabet[1:4])
		Label(fig.layout[rows[ilabel], cols[ilabel], Bottom()],
			("$(subtitles[ilabel])"),
			textsize = fontsize,
			padding = (leftpads[ilabel], 0, 0, toppads[ilabel]),
			halign = :left,
			valign = :top
		)
	end
    rowsize!(fig.layout, 2, Relative(0.64))
	colgr̄p!(fig.layout,0)
	rowgr̄p!(fig.layout,0)

	fig
end

function plot_compose_tower3d(bot0,bot1;
							cablecolor = :dodgerblue,
							barcolor = :darkred,)
	fig_width = 0.7testwidth
	fig_height = 0.6testwidth
	fig = Figure(;resolution=(fig_width,fig_height),figure_padding=(fontsize,0,0,-3fontsize))
	axs = [
		Axis3(fig[1,i],aspect=:data,
				# xlabel=L"x~(\mathrm{m})",
				# ylabel=L"y~(\mathrm{m})",
				# zlabel=L"z~(\mathrm{m})"
			)
		for i = 1:2
	]
	plot_tower3d!(axs[1],bot0.tg)
	plot_tower3d!(axs[2],bot1.tg)
	# zlims!(ax,0,0.5)
	# ax.xlabel = "x (m)"
	# ax.ylabel = "y (m)"
	# ax.zlabel = "z (m)"
	# ax.yticks = ax.xticks = [-0.1,0,0.1]
	# ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
	# ax.ylabeloffset = ax.xlabeloffset = 40
	# ax.alignmode =  Mixed(;left = 25, right = -25)

	for (iax,ax) in enumerate(axs)
		# if iax in [2]
		# 	ax.zticklabelsvisible = false
		# 	ax.zticksvisible = false
		# 	ax.zlabelvisible = false
		# end
		xlims!(ax,-0.125,0.125)
		ylims!(ax,-0.125,0.125)
		ax.azimuth = 7.595530633326982
		ax.elevation = 0.12269908169872408
	end
	leftpad = 50
	leftpadplus = -25
	toppad = 0
	subtitles = [
		"3D Tensegrity Tower (Initial Configuration)",
		"3D Tensegrity Tower (Target Configuration)",
	]

	# ndim = TR.get_ndim(bot0.tg)
	# T = TR.get_numbertype(bot0.tg)
	# bot0_rcs_by_cables = Vector{MVector{ndim,T}}()
	# foreach(bot0.tg.connectivity.tensioned) do scnt
	# 	push!(bot0_rcs_by_cables,
	# 		(
	# 			scnt.end1.rbsig.state.rps[scnt.end1.pid].+
	# 			scnt.end2.rbsig.state.rps[scnt.end2.pid]
	# 		)./2
	# 	)
	# end
	# text!(axs[1],
	# 	["l$(i)" for (i,rc) in enumerate(bot0_rcs_by_cables)];
	# 	position = bot0_rcs_by_cables,
	# 	textsize = fontsize,
	# 	color = cablecolor,
	# 	align = (:left, :center),
	# 	# offset = (-5, -10)
	# )
	for (ilabel,label) in enumerate(alphabet[1:2])
		Label(fig.layout[1, ilabel, Bottom()],
			("($label) "),
			font = "CMU Serif Bold",
			textsize = fontsize,
			padding = (leftpad-1.8fontsize, 0, 0, toppad),
			halign = :left,
			valign = :bottom
		)
	end
	for (ilabel,label) in enumerate(alphabet[1:2])
		Label(fig.layout[1, ilabel, Bottom()],
			("$(subtitles[ilabel])"),
			textsize = fontsize,
			padding = (leftpad, 0, 0, toppad),
			halign = :left,
			valign = :bottom
		)
	end
	colgr̄p!(fig.layout,2fontsize)
	fig
end

function class1(;
		k=500.0,c=0.0,
		d = 0.1*√2/2, r2 = 0.11,
		ijkl=1,
		R1 = RotX(0),
		)
	r1 = 0.1
	b = 0.22
	h = 0.1*√2/2
	γ = acos((d^2+r1^2+r2^2-b^2)/(2r1*r2))
	θ =  γ - 2π/3
	@show rad2deg.([γ,θ])
	deg120 = deg2rad(120)
	r̄pss = [
		SVector{3}.([
		       [           0,           0, 1.5h],
		    r.*[           1,           0,  0],
		    r.*[cos( deg120), sin( deg120), 0],
		    r.*[cos(-deg120), sin(-deg120), 0]
		])
		for r in [r1,r1]
	]
	r̄ps = r̄pss[1]

	ro_by_rbid = [
		SVector(0.0,0.0,0.0),
		SVector(0.0,0.0,1.0h)
	]
	R_by_rbid = [
		SMatrix(RotZ(0.0)),
		SMatrix(R1),
	]
	rirjrkrl_by_rbid = [
		Ref(ro_by_rbid[i]) .+ Ref(R_by_rbid[i]).*r̄pss[i] for i = 1:2
	]
	
	cycle3 = [2,3,4,2]
	rb1 = make_3d_tri(1,r̄ps,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1];			
		constrained = true,
		ci = collect(1:12),
		Φi = Int[]
	)

	rb2 = make_3d_tri(2,r̄ps,ro_by_rbid[2],R_by_rbid[2],rirjrkrl_by_rbid[2][1:ijkl]...)

	rbs = TypeSortedCollection([rb1,rb2])
	numberedpoints = TR.number(rbs)
	# matrix_sharing = [
	# ]
	# indexedcoords = TR.index(rbs,matrix_sharing)
	indexedcoords = TR.index(rbs)
	# #
	ncables = 6
	restlen = 0.01
    restlens = fill(restlen,ncables)
    ks = fill(k,ncables)
    cs = fill(c,ncables)
    cables = [TR.Cable3D(i,restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    acs = [TR.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)
	cnt_matrix_cables = [
		# Outer
		2 -2;
		3 -3;
		4 -4;
		# Inner
		-1 2;
		-1 3;
		-1 4;
	]
	connected = TR.connect(rbs,cnt_matrix_cables)
	#
    cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

    tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    TR.TensegrityRobot(tg,hub)
end

function spine3d(n;
		c=0.0,
		RR=RotX(0.0)
	)
    nbp = 4*n

    a = 0.04 #m
    h = 0.04 #m
    θ = π/4

    mass = 0.1 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-2
    r̄g = [0.0, 0.0, 0.0]

    r̄p_x = cos(θ)*a
    r̄p_y = sin(θ)*a
    r̄p5 = SVector{3}([  0.0,   0.0,   0.0])
    r̄p4 = SVector{3}([  0.0,  r̄p_x,  r̄p_y])
    r̄p3 = SVector{3}([  0.0, -r̄p_y,  r̄p_y])
    r̄p2 = SVector{3}([ r̄p_x,   0.0, -r̄p_y])
    r̄p1 = SVector{3}([-r̄p_x,   0.0, -r̄p_y])
    r̄ps = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]

    movable = ones(Bool,n)
    movable[1] = false
    constrained = zeros(Bool,n)
    constrained[1] = true

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(r̄g),
                r̄ps;constrained=constrained[i]) for i = 1:n]

    rs = [[0.0,0.0,i*h] for i = 0:n-1]
    Rs = [Matrix(RR)^i for i = 0:n-1]
    ṙs = [zeros(3) for i = 0:n-1]
    ωs = [zeros(3) for i = 0:n-1]

    function rigidbody(i,prop,r̄ps,r,R,ṙ,ω)
        if i == 1
            ci = collect(1:12)
            Φi = Int[]
        else
            ci = Int[]
			Φi = collect(1:6)
        end

        # ri,rj,rk,rl = [r+R*r̄p for r̄p in r̄ps]
        # nmcs,q,q̇ = TR.NCF.NC4P(ri,rj,rk,rl,r,R,ṙ,ω)
        # nmcs,q,q̇ = TR.NCF.NC3P1V(ri,rk,rl,r,R,ṙ,ω)
        # nmcs,q,q̇ = TR.NCF.NC2P2V(rk,rl,r,R,ṙ,ω)
        nmcs,_ = TR.NCF.NC1P3V(r,r,R,ṙ,ω)
        state = TR.RigidBodyState(prop,nmcs,r,R,ṙ,ω,ci,Φi)
		radius = norm(r̄ps[1]-r̄ps[5])/25
		vertmesh = merge([
			endpoints2mesh(r̄ps[i],r̄ps[5];radius)
			for i in 1:4
		])
        rb = TR.RigidBody(prop,state,vertmesh)
    end
    rbs = [rigidbody(i,props[i],r̄ps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]
	rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigdibodies)
	indexedcoords = TR.index(rigdibodies)

    ncables = 8*(n-1)
    cablelenH = 0.6h
    cablelenR = 0.04
    cablelens = repeat(vcat(fill(cablelenH,4),fill(cablelenR,4)),n-1)
    kH = 400e1
    kR = 400e1
    ks = repeat(vcat(fill(kH,4),fill(kR,4)),n-1)
    # c = 0.0
    cs = repeat(fill(c,8),n-1)
    cables = [TR.Cable3D(i,cablelens[i],ks[i],cs[i]) for i = 1:ncables]
    tensiles = (cables=cables,)	
    acs = [TR.ManualActuator(1,collect(1:ncables),cablelens[1:ncables])]
    hub = (actuators=acs,)

	matrix_cnt_raw = Vector{Matrix{Int}}()
    for i = 2:n
        s = zeros(Int,8, n)
		for j = 1:4
        	s[j,   i-1] = j;  s[j, i] = -j
        end
		s[5, i-1] = 3;  s[5, i] = -1
		s[6, i-1] = 4;  s[6, i] = -1
		s[7, i-1] = 3;  s[7, i] = -2
		s[8, i-1] = 4;  s[8, i] = -2

		# s = zeros(Int,6, n)
		# for j = 1:2
        # 	s[j,   i-1] = j+2;  s[j, i] = -(j+2)
        # end
		# s[3, i-1] = 3;  s[3, i] = -1
		# s[4, i-1] = 4;  s[4, i] = -1
		# s[5, i-1] = 3;  s[5, i] = -2
		# s[6, i-1] = 4;  s[6, i] = -2

        push!(matrix_cnt_raw, s)
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)
	# display(matrix_cnt)
    connected = TR.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)
    cnt = TR.Connectivity(numberedpoints, indexedcoords, tensioned)
    tg = TR.TensegrityStructure(rigdibodies,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
end

function newspine3d(n;
		c=0.0,
		RR=RotX(0.0)
	)
    nbp = 4 * n

    a = 0.04 #m
    h = 0.04 #m
    θ = 35*π / 180

    mass = 0.1 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174, 45.174, 25.787])) * 1e-2
    r̄g = [0.0, 0.0, 0.0]

    r̄p_x = cos(θ) * a
    r̄p_y = sin(θ) * a
    r̄p5 = SVector(  0.0,   0.0,   0.0)
    r̄p4 = SVector(  0.0,  r̄p_x,  r̄p_y)
    r̄p3 = SVector(  0.0, -r̄p_x,  r̄p_y)
    r̄p2 = SVector( r̄p_x,   0.0, -r̄p_y)
    r̄p1 = SVector(-r̄p_x,   0.0, -r̄p_y)
    r̄ps = [r̄p1, r̄p2, r̄p3, r̄p4, r̄p5]

    movable = ones(Bool, n)
    movable[1] = false
    constrained = zeros(Bool, n)
    constrained[1] = true

    props = [
		TR.RigidBodyProperty(
			i, movable[i], mass,
			SMatrix{3,3}(inertia),
			SVector{3}(r̄g),
			r̄ps; 
			constrained=constrained[i]
		) for i = 1:n
	]

    rs = [[0.0, 0.0, i * h] for i = 0:n-1]
    Rs = [Matrix(RR)^i for i = 0:n-1]
    ṙs = [zeros(3) for i = 0:n-1]
    ωs = [zeros(3) for i = 0:n-1]

    function rigidbody(i, prop, r̄ps, r, R, ṙ, ω)
        if i == 1
            ci = collect(1:12)
            Φi = Int[]
        else
            ci = Int[]
            Φi = collect(1:6)
        end

        # ri,rj,rk,rl = [r+R*r̄p for r̄p in r̄ps]
        # nmcs,q,q̇ = TR.NCF.NC4P(ri,rj,rk,rl,r,R,ṙ,ω)
        # nmcs,q,q̇ = TR.NCF.NC3P1V(ri,rk,rl,r,R,ṙ,ω)
        # nmcs,q,q̇ = TR.NCF.NC2P2V(rk,rl,r,R,ṙ,ω)
        nmcs, _ = TR.NCF.NC1P3V(r, r, R, ṙ, ω)
        state = TR.RigidBodyState(prop, nmcs, r, R, ṙ, ω, ci, Φi)
        radius = norm(r̄ps[1] - r̄ps[5]) / 25
        vertmesh = merge([
            endpoints2mesh(r̄ps[i], r̄ps[j]; radius)
            for (i,j) in [
				(1,5),(2,5),(3,5),(4,5),(1,2),(3,4)
			]
        ])
        rb = TR.RigidBody(prop, state, vertmesh)
    end
    rbs = [rigidbody(i, props[i], r̄ps, rs[i], Rs[i], ṙs[i], ωs[i]) for i = 1:n]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigdibodies)
    indexedcoords = TR.index(rigdibodies)

    ncables = 8 * (n - 1)
    cablelenH = 0.6h
    cablelenR = 0.04
    cablelens = repeat(vcat(fill(cablelenH, 4), fill(cablelenR, 4)), n - 1)
    kH = 400e1
    kR = 400e1
    ks = repeat(vcat(fill(kH, 4), fill(kR, 4)), n - 1)
    # c = 0.0
    cs = repeat(fill(c, 8), n - 1)
    cables = [TR.Cable3D(i, cablelens[i], ks[i], cs[i]) for i = 1:ncables]
    tensiles = (cables=cables,)
    acs = [TR.ManualActuator(1, collect(1:ncables), cablelens[1:ncables])]
    hub = (actuators=acs,)

    matrix_cnt_raw = Vector{Matrix{Int}}()
    for i = 2:n
        s = zeros(Int, 8, n)
        for j = 1:4
            s[j, i-1] = j
            s[j, i] = -j
        end
        s[5, i-1] = 3
        s[5, i] = -1
        s[6, i-1] = 4
        s[6, i] = -1
        s[7, i-1] = 3
        s[7, i] = -2
        s[8, i-1] = 4
        s[8, i] = -2

        # s = zeros(Int,6, n)
        # for j = 1:2
        # 	s[j,   i-1] = j+2;  s[j, i] = -(j+2)
        # end
        # s[3, i-1] = 3;  s[3, i] = -1
        # s[4, i-1] = 4;  s[4, i] = -1
        # s[5, i-1] = 3;  s[5, i] = -2
        # s[6, i-1] = 4;  s[6, i] = -2

        push!(matrix_cnt_raw, s)
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)
    # display(matrix_cnt)
    connected = TR.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    cnt = TR.Connectivity(numberedpoints, indexedcoords, tensioned)
    tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt)
    bot = TR.TensegrityRobot(tg, hub)
end

function new_deck(id,r̄ps,ro,R,ri,box;
		movable = true,
		constrained = false,
		ci = Int[],
		Φi = collect(1:6),
	)
	# uci = collect(1:6)
	m = 0.2999233976
	Īg = SMatrix{3,3}(
		Matrix(Diagonal([7.6639282053E-04,7.6638139752E-04,1.2464720496E-03]))
	)
	r̄g = [0,0,0.0]
	@show m,diag(Īg),r̄g
	prop = TR.RigidBodyProperty(
		id,movable,m,Īg,
		r̄g,r̄ps;constrained=constrained
	)
	ṙo = zero(ro)
	ω = zero(ro)
	nmcs,_ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
	# cf = TR.NCF.CoordinateFunctions(nmcs,q0,ci,uci)
	# @show typeof(nmcs)
	boxmesh = Meshes.boundary(box) |> simple2mesh
	state = TR.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci,Φi)
	TR.RigidBody(prop,state,boxmesh)
end

function bridge3d(;
		k = nothing,
		n = 2,		
		d = 10.0,
		hw = 9.0,
		h = 4.5,
		b = 6.0,
		hww = 3.0,
		l =  (hw + hww )^2+(h+b)^2 |> sqrt
	)

	hw = sqrt(l^2 - (h+b)^2 ) - hww

	up = [
		[ hw*i,6+3*i+12*j, h+b]
		for i in [-1,1], j = 0:n-1
	] .|> SVector{3}
	reverse!(@view up[2,:])
	lo = [
		[-hww*i,6+3*i+12*j,0]
		for i in [-1,1], j = 0:n-1
	] .|> SVector{3}
	reverse!(@view lo[2,:])
	mid = begin
		ret = [
			[3.0*i,6*j,b]
			for i in [-1,1], j = 0:2n
		]
		reverse!(@view ret[2,:])
		@show ret
		ret
	end |> transpose |> vec .|> SVector{3}
	display(mid)
	lbars = [
		make_3d_bar(
		i,
		lo[1,i],
		up[1,i]; ci = [1,2,3]) 
		for i = 1:n
	]
	rbars = [
		make_3d_bar(
		n+i,
		lo[2,i],
		up[2,i]; ci = [1,2,3]) 
		for i = 1:n
	]
	deckcenter = SVector(0.0,n/2*12,b)
	r̄ps_deck = mid.-Ref(deckcenter)
	厚度=0.5
	box= Meshes.Box(
		Meshes.Point(r̄ps_deck[begin]), 
		Meshes.Point(r̄ps_deck[2n+2]+SVector(0,0,厚度))
	)
	deck = new_deck(
		2n+1,r̄ps_deck, # r̄ps
		deckcenter,#ro
		RotX(0.0),#R
		deckcenter,#ri
		box;
	)
	# ;k=500.0,k1=1000.0,c=0.0,
	# 			d = 0.1*√2/2, r2 = 0.11,
	# 			ijkl=1)

	rbs = TypeSortedCollection(vcat(lbars,rbars,[deck]))
	numberedpoints = TR.number(rbs)
	indexedcoords = TR.index(rbs)
	# # #
	cnt_matrix_elas = ElasticArray{Int}(undef, 2n+1, 0)
	# left
	for i = 1:n
		if i == 1
			js, je = 2i-1, 2i+2
		elseif i == n
			js, je = 2i-3, 2i+1
		else
			js, je = 2i-3, 2i+2
		end
		for j = js:je
			row = zeros(Int,2n+1)
			row[i] = 2
			row[2n+1] = -j
			append!(cnt_matrix_elas,row)
		end
		for j = js:je
			row = zeros(Int,2n+1)
			row[i] = 1
			row[2n+1] = -j
			append!(cnt_matrix_elas,row)
		end
		for j = js:je
			row = zeros(Int,2n+1)
			row[i] = 1
			row[2n+1] =-( 4n+3-j)
			append!(cnt_matrix_elas,row)
		end
		for j = 1:n
			row = zeros(Int,2n+1)
			row[i] = 1
			row[n+j] = -2
			append!(cnt_matrix_elas,row)
		end
		# if i != 1
		# 	row = zeros(Int,2n+1)
		# 	row[[i,i+n-1]] .= [2,-1]
		# 	append!(cnt_matrix_elas,row)
		# end
		# if i != n
		# 	row = zeros(Int,2n+1)
		# 	row[[i,i+1]] .= [2,-2]
		# 	append!(cnt_matrix_elas,row)
		# end
	end
	# Right
	for i = n+1:2n
		if i == n+1
			js, je = 2i, 2i+3
		elseif i == 2n
			js, je = 2i-2, 2i+2
		else
			js, je = 2i-1, 2i+3
		end
		for j = js:je
			row = zeros(Int,2n+1)
			row[i] = 2
			row[2n+1] = -j
			append!(cnt_matrix_elas,row)
		end
		
		for j = js:je
			row = zeros(Int,2n+1)
			row[i] = 1
			row[2n+1] = -j
			append!(cnt_matrix_elas,row)
		end
		for j = js:je
			row = zeros(Int,2n+1)
			row[i] = 1
			row[2n+1] = -( 4n+3 - j)
			append!(cnt_matrix_elas,row)
		end
		for j = 1:n
			row = zeros(Int,2n+1)
			row[i] = 1
			row[j] = -2
			append!(cnt_matrix_elas,row)
		end
	# 	if i != 1
	# 		row = zeros(Int,2n+1)
	# 		row[[i,i+n-1]] .= [2,-1]
	# 		append!(cnt_matrix_elas,row)
	# 	end
		# if i != 2n
		# 	row = zeros(Int,2n+1)
		# 	row[[i,i+1]] .= [2,-2]
		# 	append!(cnt_matrix_elas,row)
		# end
	end
	display(cnt_matrix_elas)
	cnt_matrix = Matrix(transpose(cnt_matrix_elas))
	ncables = size(cnt_matrix,1)
	hncables = div(ncables,2)
	if k isa Nothing
		cables = [TR.Cable3D(i,0.0,100.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [TR.Cable3D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	acs = [
		TR.ManualActuator(i,[i,i+hncables],zeros(2))
		for i = 1:hncables
	]
	tensiles = (cables = cables,)
	hub = (actuators = acs,)
	# 	# triplex 1
	# 	1 0 0 -1  0  0  0 0 0 0;
	# 	0 1 0  0 -1  0  0 0 0 0;
	# 	0 0 1  0  0 -1  0 0 0 0;
	# 	# triplex intersect
	# 	0 0 0  1 -1  0 0 0 0 0;
	# 	0 0 0  0  1 -1 0 0 0 0;
	# 	0 0 0 -1  0  1 0 0 0 0;
	# 	# triplex 2
	# 	0 0 0  1  0  0 -2 0 0 0;
	# 	0 0 0  0  1  0 -3 0 0 0;
	# 	0 0 0  0  0  1 -4 0 0 0;
	# 	# Inner
	# 	0 0 0 0 0 0 1 -2 0 0;
	# 	0 0 0 0 0 0 1 -3 0 0;
	# 	0 0 0 0 0 0 1 -4 0 0;
	# 	# Outer
	# 	0 0 0 0 0 0 2 -2 0 0;
	# 	0 0 0 0 0 0 3 -3 0 0;
	# 	0 0 0 0 0 0 4 -4 0 0;
	# 	#
	# 	0 0 0 0 0 0 0 2 -2 0;
	# 	0 0 0 0 0 0 0 3 -4 0;
	# 	0 0 0 0 0 0 0 4 -3 0;
	# 	# Outer
	# 	0 0 0 0 0 0 0 0 2 -2;
	# 	0 0 0 0 0 0 0 0 3 -3;
	# 	0 0 0 0 0 0 0 0 4 -4;
	# 	# Inner
	# 	0 0 0 0 0 0 0 0 2 -1;
	# 	0 0 0 0 0 0 0 0 3 -1;
	# 	0 0 0 0 0 0 0 0 4 -1;
	# 	]
	connected = TR.connect(rbs,cnt_matrix)
    # ss = Int[]
	# tensiles = (cables = ss,)
	# connected = TR.connect(rbs,zeros(Int,0,0))
	# #
	#
	# cst1 = TR.PinJoint(TR.End2End(1,TR.ID(rb1_to_3[1],2),TR.ID(rb4,1)))
	# cst2 = TR.PinJoint(TR.End2End(2,TR.ID(rb1_to_3[2],2),TR.ID(rb4,2)))
	# cst3 = TR.PinJoint(TR.End2End(3,TR.ID(rb1_to_3[3],2),TR.ID(rb4,3)))
	# jointedmembers = TR.join((cst1,cst2,cst3),indexedcoords)
	#
	cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	bot = TR.TensegrityRobot(tg,hub)
end

function lander(;k=nothing)
	a = 1.0 #tetrahedra base triangle length
	h = 1.5 #tegrahedra height
	l = 2.0 #legs lengths
	θ = 2π/3 #triangle angle
	d = 1.0
	α = π/3
	P = vcat(
		[[2a/sqrt(3)*cos(i*θ),2a/sqrt(3)*sin(i*θ),0.0] for i = -1:1],
		[[sin(α).*l*cos(i*θ+π/3),sin(α).*l*sin(i*θ+π/3),-d-cos(α)*l] for i = -1:1],
		[[0.0,0.0,-h],[0.0,0.0,-d]]
	) .|> SVector{3}
	
	# scatter(P)
	legs = [
		make_3d_bar(
			i,
			P[end],
			P[i+3],; ci = Int[]
		)
		for i = 1:3
	]

	base = make_3d_tri(
		4,
		P[[7,1,2,3]],
		zero(P[7]),
		RotX(0.0),P[7];
		movable = true,
		constrained = true,
		ci = collect(1:12),
		Φi = Int[],
	)
	
	# # #
	nb = 4
	sharing_elas = ElasticArray{Int}(undef, nb, 0)
	# left
	for i = 2:3
		for j = 1:3
			row = zeros(Int,nb)
			row[i] = j
			row[1] = j
			append!(sharing_elas,row)
		end
	end
	sharing = Matrix(sharing_elas')
	display(sharing)
	rbs = TypeSortedCollection(vcat(legs,[base]))
	numberedpoints = TR.number(rbs)
	indexedcoords = TR.index(rbs,sharing)
	
	# 
	cnt_matrix_elas = ElasticArray{Int}(undef, nb, 0)
	cnt_circular = CircularArray([2,3,4])
	# left
	for i = 1:3		
		row = zeros(Int,nb)
		row[i] = 2
		row[nb] = -1
		append!(cnt_matrix_elas,row)
		for j = 0:1		
			row = zeros(Int,nb)
			row[i] = 2
			row[nb] = -cnt_circular[i+j]
			append!(cnt_matrix_elas,row)
		end		
		row = zeros(Int,nb)
		row[i] = 2
		row[cnt_circular[i+1]-1] = -2
		append!(cnt_matrix_elas,row)
	end
	for i = 2:4		
		row = zeros(Int,nb)
		row[1] = 1
		row[nb] = -i
		append!(cnt_matrix_elas,row)
	end	
	row = zeros(Int,nb)
	row[1] = 1
	row[nb] = -1
	append!(cnt_matrix_elas,row)

	display(cnt_matrix_elas)
	cnt_matrix = Matrix(cnt_matrix_elas')
	ncables = size(cnt_matrix,1)
	hncables = ncables
	if k isa Nothing
		cables = [TR.Cable3D(i,0.0,100.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [TR.Cable3D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	acs = [
		TR.ManualActuator(i,[i],zeros(1))
		for i = 1:hncables
	]
	hub = (actuators = acs,)
	tensiles = (cables = cables,)
	connected = TR.connect(rbs,cnt_matrix)
	# #
	#
	# cst1 = TR.PinJoint(TR.End2End(1,TR.ID(rb1_to_3[1],2),TR.ID(rb4,1)))
	# cst2 = TR.PinJoint(TR.End2End(2,TR.ID(rb1_to_3[2],2),TR.ID(rb4,2)))
	# cst3 = TR.PinJoint(TR.End2End(3,TR.ID(rb1_to_3[3],2),TR.ID(rb4,3)))
	# jointedmembers = TR.join((cst1,cst2,cst3),indexedcoords)
	#

	cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	bot = TR.TensegrityRobot(tg,hub)
end