function make_3d_bar(
		id,
		ri,rj;
		ci = Int[],
		m = 0.080,
		radius_ratio = 1/30,
		mat_name = "Teak",
		loadmesh = false,
	)
	# @show id,ri,rj
	movable = true
	if ci == Int[]
		constrained = false
	else
		constrained = true
	end
	u = rj - ri
	bar_length = norm(u)	
	radius = bar_length*radius_ratio
	û = u./bar_length
	v̂,ŵ = TR.NCF.HouseholderOrthogonalization(û)
	R = SMatrix{3,3}(hcat(û,v̂,ŵ))

	r̄g  = SVector{3}([ 0.0,0,0])
	r̄p1 = SVector{3}([ 0.0,0,0])
	r̄p2 = SVector{3}([ bar_length,0,0])
	r̄ps = [r̄p1,r̄p2]
	mat = filter(
		row -> row.name == mat_name, 
		material_properties
	)[1]
	density = mat.density |> ustrip
	sec_area = π*radius^2
	density_per_unit_len = density*sec_area
	mass = density_per_unit_len*bar_length
	bar_inertia_exp = :(1/12*$mass*$bar_length^2)
	Īg = SMatrix{3,3}(
		[	
			eval(bar_inertia_exp) 0.0 0.0;
					   		  0.0 0.0 0.0;
					          0.0 0.0 0.0
		]
	)
	pretty_table(
		SortedDict(
			[
				("id", id),
				("radius", radius),
				("density", density),
				("bar length", bar_length),
				("mass", mass),
				("inertia", Īg),
				("bar_inertia_exp", bar_inertia_exp)
			]
		)
	)
	ās = [SVector(1.0,0.0,0.0)]
	prop = TR.RigidBodyProperty(
		id,
		movable,
		mass,
		Īg,
		r̄g,
		r̄ps,
		ās,;
		constrained=constrained
	)
	# @show prop.inertia
	ṙo = zero(ri)
	ω = zero(ri)
	nmcs,q0,q̇0 = TR.NCF.NC3D1P1V(ri,û,ri,R,ṙo,ω)
	@myshow id,û
	# @show ri,rj,q0
	# cf = TR.NCF.CoordinateFunctions(nmcs,q0,ci,uci)
	# @show typeof(nmcs)
	state = TR.RigidBodyState(prop,nmcs,ri,R,ṙo,ω,ci)
	if loadmesh
		barmesh = load("装配体3.STL") |> make_patch(;
			trans=[0,0,0.025],
			scale=1/500,
			color=:palegreen3,
		)
	else
		barmesh = endpoints2mesh(r̄p1,r̄p2;radius,)
	end
	rb = TR.RigidBody(prop,state,barmesh)
	# @show rb.state.cache.M
	rb
end

function make_3d_tri(
		id,
		r̄ps,
		ro,
		R,
		ri,
		rj = nothing,
		rk = nothing,
		rl = nothing;
		movable = true,
		constrained = false,
		ci = Int[],
		Φi = collect(1:6),
		radius = 0.005,
		height = 1e-2,
		color = :darkorchid4,
		loadmesh = false,
		m = 3,
	)
	# uci = collect(1:6)
	mass = 0.2999233976
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
	# @show m,diag(Īg),r̄g

	ās = [
		SVector(1.0,0,0)
		for i = axes(r̄ps)
	]
	prop = TR.RigidBodyProperty(
		id,
		movable,
		mass,Īg,
		r̄g,r̄ps,ās;
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
	pretty_table(
		SortedDict(
			[
				("id", id),
				("density", density),
				("mass", mass),
				("inertia", Īg),
				("mass center", r̄g),
				("r̄ps[1]", r̄ps[1]),
				("r̄ps[2]", r̄ps[2]),
				("r̄ps[3]", r̄ps[3]),
				("r̄ps[4]", r̄ps[4])
			]
		)
	)
	# cf = TR.NCF.CoordinateFunctions(nmcs,q0,ci,uci)
	# @show typeof(nmcs)
	# radius = norm(r̄ps[2]-r̄ps[1])/32
	# @show radius
	trimesh = GB.merge(
		[
			endpoints2mesh(r̄ps[i],r̄ps[j];
			radius,color)
			for (i,j) in [
				[1,2],[1,3],[1,4],
				[2,3],[3,4],[4,2]
			]
		]
	) |> make_patch(;color = :darkslategrey)
	if loadmesh
		platemesh = endpoints2mesh(SVector(0.0,0.0,-height/2),SVector(0.0,0.0,height/2);radius=norm(r̄ps[2]),n1=m,n2=2)
		trimesh = GB.merge(
			[
				trimesh,
				platemesh
			]
		)
	end
	state = TR.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci,Φi)
	rb = TR.RigidBody(prop,state,trimesh)
end

function make_3d_plate(
		id,
		r̄ps,
		ro,
		R,
		ri,
		rj=nothing,
		rk=nothing,
		rl=nothing;
		movable = true,
		m = 3,
		height=1e-2,
		radius=1e-3,
		ci = Int[],
		Φi = collect(1:6),
		loadmesh = true,
		meshvisible = true,
	)
	# uci = collect(1:6)	
	constrained = ci != Int[]
	mat = filter(
		row -> row.name == "Teak", 
		material_properties
	)[1]
	density = mat.density |> ustrip
	if meshvisible
		if loadmesh
			platemesh = load("柚木板3.STL") |> make_patch(;rot=RotZ(π/4))
		else
			platemesh = endpoints2mesh(
				SVector(0.0,0.0,-height/2),
				SVector(0.0,0.0,height/2);
				radius,n1=m,n2=2)
		end
	else
		platemesh = nothing
	end
	# mass = density*GB.volume(platemesh)
	if m == 4
		mass =  0.11008275 #kg
		Īg = SMatrix{3,3}(
			Matrix(Diagonal([
				0.00025403,
				0.00025403,
				0.00050623
			])
			)
		)
	else
		mass = 0.23569851 #kg
		Īg = SMatrix{3,3}(
			Matrix(Diagonal([
				0.00085048,
				0.00085048,
				0.00169704				
			])
			)
		)
	end
	r̄g = SVector(0.0,0.0,0.0)
	pretty_table(
		SortedDict(
			[
				("id", id),
				("shape", "plate"),
				("No. vertices",m),
				("radius", radius),
				("height", height),
				("density", density),
				# ("bar length", bar_length),
				("mass", mass),
				("inertia", Īg),
				# ("bar_inertia_exp", bar_inertia_exp)
			]
		)
	)
	
	ās = [SVector(1.0,0.0,0.0)]
	prop = TR.RigidBodyProperty(
		id,
		movable,mass,Īg,
		r̄g,
		r̄ps,
		ās,
		;
		constrained=constrained
	)
	ṙo = zero(ro)
	ω = zero(ro)
	# u = R*(r̄ps[2] - r̄ps[1])
	# v = R*(r̄ps[3] - r̄ps[1])
	# w = R*(r̄ps[4] - r̄ps[1])
	nmcs,_ = TR.NCF.NC1P3V(ri,ro,R,ṙo,ω)
	state = TR.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci,Φi)
	rb = TR.RigidBody(prop,state,platemesh)
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
	for i = 4:6
		r̄pss[i] .-= Ref(r̄pss[i][1])
	end
	ro_by_rbid = [
		SVector(0.0,0.0, 0),
		SVector(0.0,0.0, d),
		SVector(0.0,0.0,2d),
		SVector(0.0,0,2d+1.5h),
		SVector(0.0,0,2d+1.5h),
		SVector(0.0,-0.5h*sin(α),2d+1.5h+0.5h*cos(α))
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

	rb7  = make_3d_tri( 7,r̄pss[3],ro_by_rbid[3],R_by_rbid[3],rirjrkrl_by_rbid[3][1:4]...)

	rb8  = make_3d_tri( 8,r̄pss[4],ro_by_rbid[4],R_by_rbid[4],rirjrkrl_by_rbid[4][1])

	rb9  = make_3d_tri( 9,r̄pss[5],ro_by_rbid[5],R_by_rbid[5],rirjrkrl_by_rbid[5][1])

	rb10 = make_3d_tri(10,r̄pss[6],ro_by_rbid[6],R_by_rbid[6],rirjrkrl_by_rbid[6][1:ijkl]...)

	rbs = TypeSortedCollection(vcat(rb1_to_3,rb4_to_6,[rb7,rb8,rb9,]))
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
	indexedcoords = TR.index(rbs,matrix_sharing[begin:end,begin:end-1])
	# indexedcoords = TR.index(rbs)
	# #
	ndcables = 9
	nocables = 6
	nvcables = 3
	nocables = 6
	ncables = ndcables + nocables + nvcables #+ nocables
	restlend = 0.05
	restleno = 0.01
	restlenv = 0.05
	restleno = 0.01
    restlens = vcat(
		fill(restlend,ndcables),
		fill(restleno,nocables),
		fill(restlenv,nvcables),
		# fill(restleno,nocables),
	)
    ks = vcat(
		fill(k1,ndcables),
		fill(k,nocables),
		fill(k,nvcables),
		# fill(k,nocables),
	)
    cs = fill(c,ncables)
	pretty_table(
		SortedDict(
			[
				("k1 for first $ndcables", k1),
				("k", k),
				("c", c),
			]
		)
	)
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
		# # Outer
		# 0 0 0 0 0 0 0 0 2 -2;
		# 0 0 0 0 0 0 0 0 3 -3;
		# 0 0 0 0 0 0 0 0 4 -4;
		# # Inner
		# 0 0 0 0 0 0 0 0 2 -1;
		# 0 0 0 0 0 0 0 0 3 -1;
		# 0 0 0 0 0 0 0 0 4 -1;
	]
	connected = TR.connect(rbs,cnt_matrix_cables[begin:end,begin:end-1])
	#
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

function twotre3d(;
		κ = nothing,
		r = 1.0,
		r1 = 1.0,
		b = 2.0,
		m = 4,
		α = 2π/m,
		θ = 1.25α, 
		n = 1,
		outer = false,
		loadmesh = true,
		isspine = false,
	)
	# hh = 0.5
	@assert θ>α
	c = sqrt(
		r^2 + r1^2 -2*r*r1*cos(θ)
	)
	h² = b^2-c^2
	@assert h² > 0
	h = sqrt(b^2-c^2)
	@show r,r1
	bps = [
		vcat(
			[
				[0,0,h]
			],
			[
				[3r*cos(i*α),3r*sin(i*α),0.0]
				for i = 1:m
			]
		) .|> SVector{3}
	]
	# fig = Figure()
	# ax = Axis3(fig[1,1])
	# scatter!(ax,reduce(vcat,bps))
	# scatter!(ax,reduce(vcat,midps))
	# fig
	plates = [
		begin
			if j == 1				
				ci = collect(1:12)
				Φi = Int[]
			else
				ci = Int[]
				Φi = collect(1:6)
			end
			id = j
			if isspine
				ro = SVector(0.0,0.0,0.9(j-1)*h)
				R = RotX(0.0)
			else
				ro = SVector(0.0,0.0,(j-1)*2h)
				R = RotX(((id-1)*π))
			end
			ri = ro
			r̄ps_raw = bps[1] |> Array
			# r̄ps_ext = [
			# 	SVector(3r̄p[1],3r̄p[2],0.0)
			# 	for r̄p in r̄ps_raw
			# ]
			r̄ps = vcat(r̄ps_raw,)
			make_3d_tri(
				id,
				r̄ps,ro,
				R,ri,;
				# movable = true,
				# constrained = false,
				# ci = Int[],
				# Φi = collect(1:6),
				color = :slategrey,
				loadmesh,
			)
		end
		for j = 1:n+1
	]
	nb =  length(plates)
	rbs = plates |> TypeSortedCollection
	numberedpoints = TR.number(rbs)
	# indexedcoords = TR.index(rbs,sharing)
	indexedcoords = TR.index(rbs,)

	connecting_elas = ElasticArray{Int}(undef, nb, 0)
	# outer
	if isspine
		ncables = 2m
		for i = 2:m+1
			row = zeros(Int,nb)
			row[1] =  i
			row[2] = -i
			append!(connecting_elas,row)
		end
		for i = 2:m+1
			row = zeros(Int,nb)
			row[1] =  1
			row[2] = -i
			append!(connecting_elas,row)
		end
	else
		ncables = m
		row = zeros(Int,nb)
		row[1] =   2
		row[2] = -(3)
		append!(connecting_elas,row)
		row = zeros(Int,nb)
		row[1] =   3
		row[2] = -(2)
		append!(connecting_elas,row)
		row = zeros(Int,nb)
		row[1] =   4
		row[2] = -(4)
		append!(connecting_elas,row)
	end
	connecting = Matrix(connecting_elas')
	display(connecting)
	connected = TR.connect(rbs,connecting)

	ncables = size(connecting,1)
	# @assert ncables == ncables_prism + ncables
	κ0 = 72e9*π*(3e-3)^2/1.0
	@show κ0

	cables = [TR.Cable3D(i,0.8*2h,1e3,0.0;slack=true) for i = 1:ncables]
	
	acs = [
		TR.ManualActuator(
			2,
			collect(1:ncables),
			zeros(ncables)
		),
	]
	tensiles = (cables = cables,)
	hub = (actuators = acs,)

	# jointedmembers = TR.join(csts,indexedcoords)
	cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)
	
	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	TR.TensegrityRobot(tg,hub)
end

function prisms(;
		κ = nothing,
		r = 1.0,
		r1 = 1.0,
		b = 2.0,
		m = 4,
		α = 2π/m,
		θ = 1.25α, 
		n = 1,
	)
	# hh = 0.5
	@assert θ>α
	c = sqrt(
		r^2 + r1^2 -2*r*r1*cos(θ)
	)
	h² = b^2-c^2
	@assert h² > 0
	h = sqrt(b^2-c^2)
	@show r,r1
	bps = [
		[
			[r*cos(i*α),r*sin(i*α),j*2h]
			for i = 1:m
		] .|> SVector{3} |> CircularArray
		for j = 0:n
	]
	midps = [
		[
			[r1*cos(i*α+θ),r1*sin(i*α+θ),j*2h-h]
			for i = 1:m
		] .|> SVector{3} |> CircularArray
		for j = 1:n
	]
	# fig = Figure()
	# ax = Axis3(fig[1,1])
	# scatter!(ax,reduce(vcat,bps))
	# scatter!(ax,reduce(vcat,midps))
	# fig	

	bars = [
		vcat(
			[
				make_3d_bar(
					2m*(j-1)+i,
					bps[j][i],
					midps[j][i]; 
					#ci = ifelse(j==1,[1,2,3],Int[]),
					radius_ratio = 1/60,
				) for i = 1:m
			],
			[
				make_3d_bar(
					2m*(j-1)+m+i,
					midps[j][i-1],
					bps[j+1][i-1];
					radius_ratio = 1/60,
				) for i = 1:m
			]
		)
		for j = 1:n
	]
	nbars = length.(bars) |> sum
	nb = nbars + 1
	
	ro = SVector(0.0,0.0,0.0)
	r̄ps_raw = bps[1] |> Array
	r̄ps = vcat(r̄ps_raw,)
	# r̄ps_ext = [
	# 	SVector(3r̄p[1],3r̄p[2],0.0)
	# 	for r̄p in r̄ps_raw
	# ]
	# r̄ps = vcat(r̄ps_raw,r̄ps_ext)
	plate = make_3d_plate(
		nb,
		r̄ps,
		ro,
		RotZ(0.0),
		ro;
		radius=r,
		movable = true,
		m=3,
		height=1e-2,
		ci = collect(1:12),
		Φi = Int[],
		loadmesh = false,
		meshvisible = true,
	)
	rbs = vcat(
		reduce(vcat,bars),
		plate
	 ) |> TypeSortedCollection
	numberedpoints = TR.number(rbs)
	sharing_elas = ElasticArray{Int}(undef, nb, 0)
	cm = CircularArray(collect(1:m))
	for j = 1:n
		is = 2m*(j-1)
		for i = 1:m
			for k = 1:3
				row = zeros(Int,nb)
				row[is+cm[i]  ]   = k+3
				row[is+m+cm[i+1]] = k
				append!(sharing_elas,row)		
			end
		end
	end	
	for j = 2:n
		is = 2m*(j-2)+m
		for i = 1:m
			for k = 1:3
				row = zeros(Int,nb)
				row[is+cm[i]  ]   = k+3
				row[is+m+cm[i+2]] = k
				append!(sharing_elas,row)		
			end
		end
	end
	sharing = Matrix(sharing_elas')
	# display(sharing)
	# indexedcoords = TR.index(rbs,sharing)
	indexedcoords = TR.index(rbs,)

	connecting_elas = ElasticArray{Int}(undef, nb, 0)
	for j = 1:n
		is = 2m*(j-1)
		# lower cross
		for i = 1:m
			row = zeros(Int,nb)
			row[is+cm[i  ]] =  1
			row[is+cm[i-1]] = -2
			append!(connecting_elas,row)
		end
		# additional lower cross
		# for i = 1:m
		# 	row = zeros(Int,nb)
		# 	row[is+cm[i  ]] =  1
		# 	row[is+cm[i-2]] = -2
		# 	append!(connecting_elas,row)
		# end
		# upper cross
		for i = 1:m
			row = zeros(Int,nb)
			row[is+cm[i+1]] = -2
			row[is+m+cm[i]] =  2
			append!(connecting_elas,row)
		end
		# addtional upper cross
		# for i = 1:m
		# 	row = zeros(Int,nb)
		# 	row[is+cm[i+1]] = -2
		# 	row[is+m+cm[i+1]] =  2
		# 	append!(connecting_elas,row)
		# end
		# mid
		for i = 1:m
			row = zeros(Int,nb)
			row[is+cm[i  ]] =  2
			row[is+cm[i+1]] = -2
			append!(connecting_elas,row)
		end
		# # upper
		# for i = 1:m
		# 	row = zeros(Int,nb)
		# 	row[is+m+cm[i  ]] =  2
		# 	row[is+m+cm[i+1]] = -2
		# 	append!(connecting_elas,row)
		# end
	end
	ncables_prism = size(connecting_elas,2)
	for j = 1:n
		is = 0
		for i = 1:m
			if j > 1
				row = zeros(Int,nb)
				row[is+cm[i]] =    1
				row[is+cm[i+1]] = -1
				append!(connecting_elas,row)
			end
			row = zeros(Int,nb)
			row[is+m+cm[i]] =    2
			row[is+m+cm[i+1]] = -2
			append!(connecting_elas,row)
		end
	end
	ncables_prism += (n)*m
	connecting = Matrix(connecting_elas')
	# display(connecting)
	connected = TR.connect(rbs,connecting)

	ncables = size(connecting,1)
	# @assert ncables == ncables_prism + ncables_outer

	mat_cable = filter(
		row->row.name == "Nylon 66",
		material_properties
	)[1]
	diameter = 1e-3Unitful.m
	cable_length = 0.1Unitful.m
	κ = (mat_cable.modulus_elas)*π*(diameter/2)^2/cable_length
	# @show κ
	@show uconvert(Unitful.N/Unitful.m,κ),ustrip(Unitful.N/Unitful.m,κ)
	cables_prism = [TR.Cable3D(i,0.0,   ustrip(Unitful.N/Unitful.m,κ),0.0;slack=true) for i = 1:ncables_prism]
	cables = cables_prism
	acs = [
		TR.ManualActuator(
			1,
			collect(1:ncables_prism),
			zeros(ncables_prism)
		),
	]
	tensiles = (cables = cables,)
	hub = (actuators = acs,)

	
	csts_bar2plate = [
		TR.PinJoint(
			i,
			TR.End2End(
				i,
				TR.ID(plate,i,1),
				TR.ID(bars[1][cm[i]],1,1),
			)
		)
		for i = 1:m
	]

	csts_bar2bar = [
		TR.PinJoint(
			m+i+m*(j-1),
			TR.End2End(
				m+i+m*(j-1),
				TR.ID(bars[j][m+cm[i+1]],1,1),
				TR.ID(bars[j][cm[i]]    ,2,1)
			)
		)
		for j = 1:n for i = 1:m
	]

	csts = vcat(
		csts_bar2plate,
		csts_bar2bar,
		
	)
	
	jointedmembers = TR.join(csts,indexedcoords)

	cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),jointedmembers)

	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	bot = TR.TensegrityRobot(tg,hub)
end
function embed3d(;
		κ = nothing,
		r = 1.0,
		r1 = 1.0,
		b = 2.0,
		m = 4,
		α = 2π/m,
		θ = 1.25α, 
		n = 1,
		outer = false,
		isprism = false,
	)
	# hh = 0.5
	@assert θ>α
	c = sqrt(
		r^2 + r1^2 -2*r*r1*cos(θ)
	)
	h² = b^2-c^2
	@assert h² > 0
	h = sqrt(b^2-c^2)
	@show r,r1
	bps = [
		[
			[r*cos(i*α),r*sin(i*α),j*2h]
			for i = 1:m
		] .|> SVector{3} |> CircularArray
		for j = 0:n
	]
	midps = [
		[
			[r1*cos(i*α+θ),r1*sin(i*α+θ),j*2h-h]
			for i = 1:m
		] .|> SVector{3} |> CircularArray
		for j = 1:n
	]
	# fig = Figure()
	# ax = Axis3(fig[1,1])
	# scatter!(ax,reduce(vcat,bps))
	# scatter!(ax,reduce(vcat,midps))
	# fig
	bars = [
		vcat(
			[
				make_3d_bar(
					2m*(j-1)+i,
					bps[j][i],
					midps[j][i]; ci = ifelse(j==1,[1,2,3],Int[])
				) for i = 1:m
			],
			[
				make_3d_bar(
					2m*(j-1)+m+i,
					midps[j][i-1],
					bps[j+1][i-1];
				) for i = 1:m
			]
		)
		for j = 1:n
	]
	plates = [
		begin
			if j == 1				
				ci = collect(1:12)
				Φi = Int[]
			else
				ci = Int[]
				Φi = collect(1:6)
			end
			id = (2n)*m+j
			ro = SVector(0.0,0.0,(j-1)*2h)
			ri = ro
			R = RotZ(0.0)
			r̄ps_raw = bps[1] |> Array
			r̄ps_ext = [
				SVector(3r̄p[1],3r̄p[2],0.0)
				for r̄p in r̄ps_raw
			]
			r̄ps = vcat(r̄ps_raw,r̄ps_ext)
			make_3d_plate(
				id,
				r̄ps,
				ro,
				R,
				ri;
				radius=3r,
				movable = true,
				m,
				height=1e-2,
				ci,
				Φi,
				loadmesh = false,
				meshvisible = !isprism,
			)
		end
		for j = 1:n+1
	]
	nb = sum(length.(bars)) + length(plates)
	rbs = vcat(reduce(vcat,bars),plates) |> TypeSortedCollection
	numberedpoints = TR.number(rbs)
	sharing_elas = ElasticArray{Int}(undef, nb, 0)
	cm = CircularArray(collect(1:m))
	for j = 1:n
		is = 2m*(j-1)
		for i = 1:m
			for k = 1:3
				row = zeros(Int,nb)
				row[is+cm[i]  ]   = k+3
				row[is+m+cm[i+1]] = k
				append!(sharing_elas,row)		
			end
		end
	end	
	for j = 2:n
		is = 2m*(j-2)+m
		for i = 1:m
			for k = 1:3
				row = zeros(Int,nb)
				row[is+cm[i]  ]   = k+3
				row[is+m+cm[i+2]] = k
				append!(sharing_elas,row)		
			end
		end
	end
	sharing = Matrix(sharing_elas')
	# display(sharing)
	indexedcoords = TR.index(rbs,sharing)
	# indexedcoords = TR.index(rbs,)

	connecting_elas = ElasticArray{Int}(undef, nb, 0)
	for j = 1:n
		is = 2m*(j-1)
		# lower cross
		for i = 1:m
			row = zeros(Int,nb)
			row[is+cm[i  ]] =  1
			row[is+cm[i-1]] = -2
			append!(connecting_elas,row)
		end
		# additional lower cross
		# for i = 1:m
		# 	row = zeros(Int,nb)
		# 	row[is+cm[i  ]] =  1
		# 	row[is+cm[i-2]] = -2
		# 	append!(connecting_elas,row)
		# end
		# upper cross
		for i = 1:m
			row = zeros(Int,nb)
			row[is+cm[i+1]] = -2
			row[is+m+cm[i]] =  2
			append!(connecting_elas,row)
		end
		# addtional upper cross
		# for i = 1:m
		# 	row = zeros(Int,nb)
		# 	row[is+cm[i+1]] = -2
		# 	row[is+m+cm[i+1]] =  2
		# 	append!(connecting_elas,row)
		# end
		# mid
		for i = 1:m
			row = zeros(Int,nb)
			row[is+cm[i  ]] =  2
			row[is+cm[i+1]] = -2
			append!(connecting_elas,row)
		end
		# # upper
		# for i = 1:m
		# 	row = zeros(Int,nb)
		# 	row[is+m+cm[i  ]] =  2
		# 	row[is+m+cm[i+1]] = -2
		# 	append!(connecting_elas,row)
		# end
	end
	ncables_prism = size(connecting_elas,2)
	# outer
	ncables_outer = 0
	for j = 1:n
		if outer
			is = (2n)*m
			for i = 1:m
				row = zeros(Int,nb)
				row[is+j] =  m+i
				row[is+j+1] = -i-m
				append!(connecting_elas,row)
			end
			ncables_outer = n*m
		end
	end
	if isprism
		for j = 1:n
			is = 0
			for i = 1:m
				row = zeros(Int,nb)
				row[is+cm[i]] =    1
				row[is+cm[i+1]] = -1
				append!(connecting_elas,row)
				row = zeros(Int,nb)
				row[is+m+cm[i]] =    2
				row[is+m+cm[i+1]] = -2
				append!(connecting_elas,row)
			end
		end
		ncables_prism += (n+1)*m
	end
	connecting = Matrix(connecting_elas')
	# display(connecting)
	connected = TR.connect(rbs,connecting)

	ncables = size(connecting,1)
	# @assert ncables == ncables_prism + ncables_outer

	mat_cable = filter(
		row->row.name == "Nylon 66",
		material_properties
	)[1]
	diameter = 1e-3Unitful.m
	cable_length = 0.1Unitful.m
	κ = (mat_cable.modulus_elas)*π*(diameter/2)^2/cable_length
	# @show κ
	@show uconvert(Unitful.N/Unitful.m,κ),ustrip(Unitful.N/Unitful.m,κ)
	cables_prism = [TR.Cable3D(i,0.0,   ustrip(Unitful.N/Unitful.m,κ),0.0;slack=true) for i = 1:ncables_prism]
	cables_outer = [
		TR.Cable3D(i,[0.95,0.85,0.75][((i-ncables_prism) % 3)+1]*0.2,ustrip(Unitful.N/Unitful.m,κ),0.0;slack=true) 
		for i = ncables_prism+1:ncables
	]
	@show 2h
	cables = vcat(
		cables_prism,
		cables_outer
	)
	acs = [
		TR.ManualActuator(
			1,
			collect(1:ncables_prism),
			zeros(ncables_prism)
		),
		TR.ManualActuator(
			2,
			collect(ncables_prism+1:ncables),
			zeros(ncables_outer)
		),
	]
	tensiles = (cables = cables,)
	hub = (actuators = acs,)
	
	csts = [
		TR.PinJoint(i+m*(j-1),TR.End2End(i,TR.ID(bars[j][m+cm[i]],2),TR.ID(plates[j+1],cm[i-1])))
		 for j = 1:n for i = 1:m
	]

	jointedmembers = TR.join(csts,indexedcoords)

	cnt = TR.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),jointedmembers)

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
	# 			scnt.hen.rbsig.state.rps[scnt.hen.pid].+
	# 			scnt.egg.rbsig.state.rps[scnt.egg.pid]
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
	mat = filter(
		row -> row.name == "Teak", 
		material_properties
	)[1]
	density = mat.density |> ustrip
	box_volume = Meshes.volume(box)
	mass = density*box_volume
	# bar_inertia_exp = :(1/12*$mass*$bar_length^2)
	wid,len,hei = Meshes.sides(box)
	# m = 0.2999233976
	Īg = SMatrix{3,3}(
		Matrix(
			Diagonal(
				[
					mass*(hei^2+len^2)/12, #Ix
					mass*(hei^2+wid^2)/12, #Iy
					mass*(len^2+wid^2)/12, #Iz
				]
			)
		)
	)
	r̄g = [0,0,0.0]
	pretty_table(
		SortedDict(
			[
				("id", id),
				("mass", mass),
				("density", density),
				("inertia", Īg),
				("wid,len,hei", (wid,len,hei)),
				("box_volume", box_volume)
			]
		)
	)

	prop = TR.RigidBodyProperty(
		id,
		movable,
		mass,
		Īg,
		r̄g,
		r̄ps;
		constrained=constrained
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
		# @show ret
		ret
	end |> transpose |> vec .|> SVector{3}
	# display(mid)
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
	r̄ps_deck = vcat(mid,[SVector(0.0,0.0,b)]).-Ref(deckcenter)
	hei = 0.25
	box = Meshes.Box(
		Meshes.Point(r̄ps_deck[begin]), 
		Meshes.Point(r̄ps_deck[2n+2]+SVector(0,0,hei))
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
	# cross
	for i = 1:n
		if  i == 1
			js, je = 2n+1-i, 2n+1-i
		else
			js, je = 2n+1-i, 2n+1-i+1
		end
		for j = js:je		
			row = zeros(Int,2n+1)
			row[i] = 2
			row[j] = -2
			append!(cnt_matrix_elas,row)
		end
	end
	# display(cnt_matrix_elas)
	cnt_matrix = Matrix(transpose(cnt_matrix_elas))
	ncables = size(cnt_matrix,1)
	hncables = div(ncables,2)
	mat_cable = filter(
		row->row.name == "Nylon 66",
		material_properties
	)[1]
	diameter = 1e-2
	κ = (mat_cable.modulus_elas |> ustrip)*1e9*π*(diameter/2)^2/10
	@show κ
	if k isa Nothing
		cables = [TR.Cable3D(i,0.0,κ,0.0;slack=false) for i = 1:ncables]
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

function tower(;k=nothing)
	a = 1.0 #tetrahedra base triangle length
	h = 0.6 #tegrahedra height
	l = 2.0 #legs lengths
	θ = 2π/3 #triangle angle
	d = 1.0
	α = π/3
	P = vcat(
		[
			[
				2a/sqrt(3)*cos(i*θ),
				2a/sqrt(3)*sin(i*θ),
				-0.5h
			] 
			for i = -1:1
		],
		[
			[
				2a/sqrt(3)*cos((i+0.5)*θ),
				2a/sqrt(3)*sin((i+0.5)*θ),
				h
			] 
			for i = -1:1
		],
		[
			[0.0,0.0,d],
			[0.0,0.0,0.0],
		]
	) .|> SVector{3}
	
	# scatter(P)
	
	base = make_3d_tri(
		1,
		P[[8,1,2,3]],
		zero(P[8]),
		RotX(0.0),
		P[8];
		movable = true,
		constrained = true,
		ci = collect(1:12),
		Φi = Int[],
		radius = 0.04,
	)

	bar = make_3d_bar(
		2,
		P[end],
		P[end-1],;		
		radius_ratio = 1/40,
		ci = Int[]
	)

	
	top = make_3d_tri(
		3,
		P[[7,4,5,6]],
		zero(P[7]),
		RotX(0.0),
		P[7];
		radius = 0.04,
		movable = true,
		constrained = false,
	)

	# # #
	nb = 3
	rbs = TypeSortedCollection([base,bar,top])
	numberedpoints = TR.number(rbs)
	indexedcoords = TR.index(rbs,)
	
	# 
	cnt_matrix_elas = ElasticArray{Int}(undef, nb, 0)
	cnt_circular = CircularArray([2,3,4])
	# left
	for i = 1:3	
		row = zeros(Int,nb)
		row[1] = -(i+1)
		row[2] = 2
		append!(cnt_matrix_elas,row)
	end
	for i = 1:3	
		row = zeros(Int,nb)
		row[2] = 1
		row[3] = -(i+1)
		append!(cnt_matrix_elas,row)
	end
	for i = 1:3	
		row = zeros(Int,nb)
		row[1] = cnt_circular[i]
		row[3] = -cnt_circular[i]
		append!(cnt_matrix_elas,row)
		row = zeros(Int,nb)
		row[1] = cnt_circular[i]
		row[3] = -cnt_circular[i-1]
		append!(cnt_matrix_elas,row)
	end	

	display(cnt_matrix_elas)
	cnt_matrix = Matrix(cnt_matrix_elas')
	ncables = size(cnt_matrix,1)
	if k isa Nothing
		cables = [TR.Cable3D(i,0.0,100.0,0.0;slack=false) for i = 1:ncables]
	else
		cables = [TR.Cable3D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
	end
	tensiles = (cables = cables,)
	connected = TR.connect(rbs,cnt_matrix)

	cst1 = TR.PinJoint(
		1,
		TR.End2End(
			1,
			TR.ID(base,1),
			TR.ID(bar,1),
		)
	)
	
	cst2 = TR.PinJoint(
		2,
		TR.End2End(
			2,
			TR.ID(bar,2),
			TR.ID(top,1)
		)
	)

	jointedmembers = TR.join((cst1,cst2,),indexedcoords)

	cnt = TR.Connectivity(
			numberedpoints,
			indexedcoords,
			@eponymtuple(connected,),
			jointedmembers
		)

	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
	bot = TR.TensegrityRobot(tg,)
end

function simple(;
		c=0.0,
		m = 4,
		α = 2π/m,
		k = 100.0,
		z0 = 0.2,
		ωz = 5.0,
		mbar = 0.05,
		free = false,
	)
	lₛ = 15e-3
	DE = 300e-3
	d = (150*√2)*1e-3
	r = d/2
	bps = [
		[r*cos(i*α),r*sin(i*α),0.0]
		for i = 1:m
	] .|> SVector{3}
	ro = SVector(0,0,z0)
	Rplate = RotX(π/4)
	Rbar = RotY(-π/12)
	p = Ref(Rbar) .* [
		[0,0, DE/2],
		[0,0,-DE/2]
	] .+ Ref(ro) .|> SVector{3}
	
	# fig,ax,plt = scatter(p)
	# # ax.aspect = DataAspect()
	# fig
	ṙo_ref = SVector(0.0,0,0)

	rbs = [
		make_3d_bar(
			1, 
			p[end-1],
			p[end];
			# ṙi = ṙo_ref, 
			# ṙj = ṙo_ref,
			m = mbar,
			mat_name = "Teak"
		),
		make_3d_plate(
			2, 
			bps,ro,
			Rplate,ro;
			m,
			radius = r,
			# constrained = !free,	
			ci = ifelse(free,Int[],collect(1:12)),
			Φi = ifelse(free,collect(1:6),Int[]),
		)
	]

	rigdibodies = TypeSortedCollection(rbs)
	numberedpoints = TR.number(rigdibodies)
	indexedcoords = TR.index(rigdibodies)
	#
	ncables = 2m

	original_restlens = zeros(ncables)
	original_restlens[  1: m] .= 0.05
	original_restlens[m+1:2m] .= 0.1

	ks = zeros(ncables)
	ks[  1: m] .= ks[m+1:2m] .= k

	cables = [
		TR.Cable3D(i, original_restlens[i], ks[i], c;slack=true) for i = 1:ncables
	]

	pretty_table(
		reduce(vcat,[
				[
					"No.$i" "Rest length" "$(original_restlens[i])" "Stiffness" "$(ks[i])" "Damping" "$c";
				] 
				for i = 1:2m
			]
		)
	)
	#
	tensiles = (cables = cables,)
	acs = [
	TR.ManualActuator(
		i,
		collect(1:2m), [original_restlens[2m*(i-1)+j] for j = 1:6],
	) for i = [1]
	]
	hub = (actuators = acs,)
	# #
	matrix_cnt = zeros(Int,ncables,2)
	for j = 1:ncables
	if j <= m
		matrix_cnt[j, 1:2] = [1, -j]
	else
		matrix_cnt[j, 1:2] = [2, -(j-m)]
	end
	end
	# display(matrix_cnt)
	connected = TR.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)
	#
	cnt = TR.Connectivity(numberedpoints, indexedcoords, tensioned)
	# #
	tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt, )
	bot = TR.TensegrityRobot(tg, hub)
end