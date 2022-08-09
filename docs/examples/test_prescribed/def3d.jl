function make_3d_bar(id,ri,rj;ci = Vector{Int}())
	# @show id,ri,rj
	movable = true
	if ci == Vector{Int}()
		constrained = false
	else
		constrained = true
	end
	u = rj - ri
	b = norm(u)
	û = u./b
	v̂,ŵ = TR.NaturalCoordinates.HouseholderOrthogonalization(û)
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
	lncs,q0,q̇0 = TR.NaturalCoordinates.NC3D2P(ri,rj,ro,R,ṙo,ω)
	# @show ri,rj,q0
	# cf = TR.NaturalCoordinates.CoordinateFunctions(lncs,q0,ci,uci)
	# @show typeof(lncs)
	state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,ci)
	rb = TR.RigidBody(prop,state)
	# @show rb.state.cache.M
	rb
end

function make_3d_tri(id,r̄ps,ro,R,ri,rj=nothing,rk=nothing,rl=nothing;)
	movable = true
	constrained = false
	ci = Vector{Int}()
	# uci = collect(1:6)
	m = 0.2999233976
	Īg = SMatrix{3,3}(
		Matrix(Diagonal([7.6639282053E-04,7.6638139752E-04,1.2464720496E-03]))
	)
	r̄g = [0,0,0.1562442983-0.1*√2]
	@show m,diag(Īg),r̄g
	prop = TR.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
				)
	ṙo = zero(ro)
	ω = zero(ro)
	u = R*(r̄ps[2] - r̄ps[1])
	v = R*(r̄ps[3] - r̄ps[1])
	w = R*(r̄ps[4] - r̄ps[1])
	if rj isa Nothing
		lncs,_ = TR.NaturalCoordinates.NC1P3V(ri,ro,R,ṙo,ω)
	elseif rk isa Nothing
		lncs,_ = TR.NaturalCoordinates.NC2P2V(ri,rj,ro,R,ṙo,ω)
	elseif rl isa Nothing
		lncs,_ = TR.NaturalCoordinates.NC3P1V(ri,rj,rk,ro,R,ṙo,ω)
	else
		lncs,_ = TR.NaturalCoordinates.NC4P(ri,rj,rk,rl,ro,R,ṙo,ω)
	end
	# cf = TR.NaturalCoordinates.CoordinateFunctions(lncs,q0,ci,uci)
	# @show typeof(lncs)
	state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,ci)
	rb = TR.RigidBody(prop,state)
end

function tower3d(;k=500.0,k1=1000.0,c=0.0,
					d = 0.1*√2/2, r2 = 0.11,
					ijkl=1)
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
		SVector(0.0,0.0,2d+0.5h),
		SVector(0.0,0.0,2d+0.5h+2h),
		SVector(0.0,0.0,2d+0.5h+2h+0.5h)
	]
	R_by_rbid = [
		SMatrix(RotX(0.0)),
		SMatrix(RotZ(θ)),
		SMatrix(RotZ(2θ)),
		SMatrix(RotZ(2θ)),
		SMatrix(RotZ(2θ)*RotX(π)),
		SMatrix(RotZ(2θ)*RotX(π)),
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
					    rirjrkrl_by_rbid[3][cycle3[i+1]]; ci = Vector{Int}()) for i = 1:3
	]

	 rb7 = make_3d_tri( 7,r̄ps,ro_by_rbid[3],R_by_rbid[3],rirjrkrl_by_rbid[3][1:4]...)

	 rb8 = make_3d_tri( 8,r̄ps,ro_by_rbid[4],R_by_rbid[4],rirjrkrl_by_rbid[4][1])

	 rb9 = make_3d_tri( 9,r̄ps,ro_by_rbid[5],R_by_rbid[5],rirjrkrl_by_rbid[5][1])

	rb10 = make_3d_tri(10,r̄ps,ro_by_rbid[6],R_by_rbid[6],rirjrkrl_by_rbid[6][1:ijkl]...)

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

function get_linesegs_cables(tg)
	(;nrigids) = tg
	(;connected) = tg.connectivity.tensioned
	(;cables) = tg.tensiles
	ndim = TR.get_ndim(tg)
	T = TR.get_numbertype(tg)
	# cables
	# cables
	linesegs_noslack_cables = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
	linesegs_slack_cables = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
	foreach(connected) do scnt
		scable = cables[scnt.id]
		ret = Point(scnt.end1.rbsig.state.rps[scnt.end1.pid]) =>
				Point(scnt.end2.rbsig.state.rps[scnt.end2.pid])
		if scable.state.tension > 0
			push!(linesegs_noslack_cables,ret)
		else
			push!(linesegs_slack_cables,ret)
		end
	end
	linesegs_noslack_cables, linesegs_slack_cables
end

function get_linesegs_bars(tg)
	ndim = TR.get_ndim(tg)
	T = TR.get_numbertype(tg)
	rbs = TR.get_rigidbodies(tg)
	linesegs = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
	for rb in rbs
		if rb.state.cache.funcs.lncs isa TR.NaturalCoordinates.LNC3D2P
			push!(
				linesegs,
				Point(rb.state.rps[1]) => Point(rb.state.rps[2]),
			)
		elseif rb.state.cache.funcs.lncs isa TR.NaturalCoordinates.LNCMP
		else
			# @show typeof(rb.state.cache.funcs.lncs)
			append!(
				linesegs,
				[
					Point(rb.state.rps[2]) => Point(rb.state.rps[3]),
					Point(rb.state.rps[3]) => Point(rb.state.rps[4]),
					Point(rb.state.rps[4]) => Point(rb.state.rps[2]),
					Point(rb.state.rps[2]) => Point(rb.state.rps[1]),
					Point(rb.state.rps[3]) => Point(rb.state.rps[1]),
					Point(rb.state.rps[4]) => Point(rb.state.rps[1]),
				]
			)
		end
	end
	linesegs
end

function plot_tower3d!(ax,
					tg::TR.TensegrityStructure;
					cablecolor = :dodgerblue,
					barcolor = :darkred,
					markercolor = :darkblue,
					markit = false,
					textsize = 10)

	linesegs_noslack_cables, linesegs_slack_cables = get_linesegs_cables(tg)
    linesegments!(ax, linesegs_noslack_cables, color = cablecolor, linewidth = cablewidth)
    linesegments!(ax, linesegs_slack_cables, color = cablecolor, linewidth = cablewidth, linestyle = :dash)

	linesegs_bars = get_linesegs_bars(tg)
	linesegments!(ax, linesegs_bars, color = barcolor, linewidth = barwidth)

	rbs = TR.get_rigidbodies(tg)
	if markit
		meshscatter!(ax,[rbs[end].state.rg], color = markercolor, markersize = 0.008)
	end

	xlims!(ax,-0.125,0.125)
	ylims!(ax,-0.125,0.125)
	zlims!(ax,0,0.5)
	ax.xlabel = "x (m)"
	ax.ylabel = "y (m)"
	ax.zlabel = "z (m)"
	ax.yticks = ax.xticks = [-0.1,0,0.1]
	ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
	ax.ylabeloffset = ax.xlabeloffset = 40
	ax.alignmode =  Mixed(;left = 25, right = -25)
	ax.azimuth = 7.525530633326982
	ax.elevation = 0.12269908169872408
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
	colgap!(fig.layout,0)
	rowgap!(fig.layout,0)

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
	colgap!(fig.layout,2fontsize)
	fig
end
