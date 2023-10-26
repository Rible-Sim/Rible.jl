function make_strut(id,ri,rj;ci = Int[])

    movable = true
	if ci == Int[]
		constrained = false
	else
		constrained = true
	end
	u = rj - ri
	b = norm(u)
	û = u./b
	v̂,ŵ = RB.NCF.HouseholderOrthogonalization(û)
	R = SMatrix{3,3}(hcat(û,v̂,ŵ))
	d = 28/1000
	r = d/2
	Δ = 1.5/1000
	ca = π*(r^2-(r-Δ)^2)
	vol = ca*b
	ρ =  7850.0
	m = ρ*vol
	Īg = SMatrix{3,3}([1/12*m*b^2 0.0 0.0;
					   		  0.0 0.0 0.0;
					          0.0 0.0 0.0])
	# @show id,b#,m,Īg[1]
	r̄g  = SVector{3}([b/2,0,0])
	r̄p1 = SVector{3}([0.0,0,0])
	r̄p2 = SVector{3}([  b,0,0])
	r̄ps = [r̄p1,r̄p2]
	prop = RB.RigidBodyProperty(id,movable,m,Īg,
				r̄g,r̄ps;constrained=constrained
				)
	# @show prop.inertia
	ro = ri
	ṙo = zero(ro)
	ω = zero(ro)
	lncs,q0,q̇0 = RB.NCF.NC3D2P(ri,rj,ro,R,ṙo,ω)
	state = RB.RigidBodyState(prop,lncs,ro,R,ṙo,ω,ci)
	strutmesh = endpoints2mesh(r̄p1,r̄p2)
	rb = RB.RigidBody(prop,state,strutmesh)
end

function nbridge(n,m=2;θ=missing,r=missing,c=0.0,h=1.0,o2=[0,4.0,0],right=false)
	l = Int(right)
    α = 2π/n
    b = 1.35
	a = sqrt(b^2-h^2)
	if (θ isa Missing) ⊻ (r isa Missing)
		if θ isa Missing
			@assert abs(a/(2r)) < 1
			θ = 2asin(a/(2r)) - α
			@show rad2deg(θ)
		else
			r = a/2/sin((α+θ)/2)
			@show r
		end
	else
		error("Either one of θ or r should be specified.")
	end

    as = [
        SVector(
			r*cos(i*α+j*θ+π/2),j*h,
			r*sin(i*α+j*θ+π/2)
		)
        for i = 0:n-1, j = 0:m
    ]

	# a = sqrt((2*sin(α)*r)^2 + h^2)
	# @show a
	a = norm(as[3,2]-as[1,1])
	s = sqrt(b^2 - (a/2)^2)
	r1 = s - norm((as[3,2].+as[1,1])./2-[0,h/2,0])
	@show r, r1
	@assert r1 > 0
	bs = [
		SVector(
			r1*cos((i+1//2)*α+(j+1//2)*θ+π/2),(j+1/2)*h,
		    r1*sin((i+1//2)*α+(j+1//2)*θ+π/2)
		)
		for i = 0:n-1, j = 0:m-1
	]

	rot180 = LinearMap(RotZ(π))
	translate2o2 = Translation(o2)
	mirror = LinearMap(
		SMatrix{3,3}([
		-1 0 0;
		 0 1 0;
		 0 0 1;
		])
	)
	transform2right = translate2o2 ∘ mirror ∘ rot180
	asr = transform2right.(as)
	bsr = transform2right.(bs)

	aslr = hcat(as,asr)
	bslr = hcat(bs,bsr)
	nhalf = 3*n*m
	rbs = reduce(vcat,[
		[
			make_strut(k*nhalf+(j-1)*3*n+(i-1)*3+1,CA(aslr[:,k*(m+1)+j  ])[i  ],
										  		   CA(aslr[:,k*(m+1)+j+1])[i+1];
										  		   ci=ifelse(j==1,[1,2,3],Int[])),
			make_strut(k*nhalf+(j-1)*3*n+(i-1)*3+2,CA(aslr[:,k*(m+1)+j])[i  ],
										  	   	   CA(bslr[:,k*m+j])[i-2];
										  		   ci=ifelse(j==1,[1,2,3],Int[])),
			make_strut(k*nhalf+(j-1)*3*n+(i-1)*3+3,CA(bslr[:,k*m+j  ])[i-2],
										  	       CA(aslr[:,k*(m+1)+j+1])[i+2];)
		]
		for k = 0:l
		for j = 1:m
		for i = 1:n
	])

	nb = length(rbs)
	rigidbodies = TypeSortedCollection(rbs)
	numbered = RB.number(rigidbodies)
	sharing_rows = [
		begin
			sharing_row = zeros(Int,ifelse(j==m,3*3,4*3),nb)
			ca = CA((j-1)*3n+1:    j*3n)
			cb = CA(    j*3n+1:(j+1)*3n)
			sharing_row[1:3,k*nhalf+ca[(i-1)*3+1]] .= 1:3
			sharing_row[1:3,k*nhalf+ca[(i-1)*3+2]] .= 1:3
			sharing_row[4:6,k*nhalf+ca[(i-1)*3+2]] .= 4:6
			sharing_row[4:6,k*nhalf+ca[(i-1)*3+3]] .= 1:3
			sharing_row[7:9,k*nhalf+ca[(i-1)*3+3]] .= 4:6
			sharing_row[7:9,k*nhalf+ca[(i-1)*3+4]] .= 4:6
			if j !== m
				sharing_row[10:12,k*nhalf+ca[(i-1)*3+1 ]] .= 4:6
				sharing_row[10:12,k*nhalf+cb[(i-1)*3+19]] .= 1:3
			end
			sharing_row
		end
		for k = 0:l
		for j = 1:m
		for i = 1:n
	]
	sharing = reduce(vcat,sharing_rows)
	# browse(matrix_sharing)
	indexed = RB.index(rigidbodies,sharing)

	nhcables = 5n*m
	ncables = (l+1)*nhcables
	active_cable_indices = [
		k*nhcables + CA((j-1)*5*n+1:j*5n)[(j-1)*5+(i-1)*5+v]
		for k = 0:l for j = 1:m for i = 1:n for v = 2:3
	]
	stale_cable_indices = [
		k*nhcables + CA((j-1)*5*n+1:j*5n)[-(j-1)*5+(i-1)*5+v]
		for k = 0:l for j = 1:m for i = 1:n for v = [1,19]
	]

	restlens = fill(0.0,ncables)
	ks = fill(2000.0,ncables)
	ks[stale_cable_indices] .= 3000.0
    cs = fill(c,ncables)
    cables = [RB.Cable3D(i,restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    tensiles = (cables = cables,)
	cnt_rows = [
		begin
			cnt_row = zeros(Int,5,nb)
			# cnt_row[1,CA(3n,(i-1)*3:1)] )[ 1; cnt_row[1,CA(3n,      i:3+1))[ = -1
			ca = CA((j-1)*3n+1:j*3n)
			cnt_row[1,k*nhalf + ca[  (i-1)*3+1]] =  1
			cnt_row[1,k*nhalf + ca[n+(i-1)*3+1]] = -1
			cnt_row[2,k*nhalf + ca[  (i-1)*3+1]] =  1
			cnt_row[2,k*nhalf + ca[n+(i-1)*3+4]] = -1
			cnt_row[3,k*nhalf + ca[  (i-1)*3+1]] =  2
			cnt_row[3,k*nhalf + ca[n+(i-1)*3+4]] = -1
			cnt_row[4,k*nhalf + ca[  (i-1)*3+1]] =  2
			cnt_row[4,k*nhalf + ca[n+(i-1)*3+7]] = -1

			cnt_row[5,k*nhalf + ca[  (i-1)*3+1]] =  2
			cnt_row[5,k*nhalf + ca[      i*3+1]] = -2

			cnt_row
		end
		for k in 0:l
		for j in 1:m
		for i in 1:n
	]
	cnt_cables = reduce(vcat,cnt_rows)
	# display(cnt_cables)
	connected = RB.connect(rigidbodies,cnt_cables)
	tensioned = @eponymtuple(connected,)
    cnt = RB.Connectivity(numbered,indexed,tensioned)
    st = RB.Structure(rigidbodies,tensiles,cnt)

	actuators_active =
		RB.ManualActuator(
			1,
			         active_cable_indices,
			restlens[active_cable_indices],
			RB.Serial()
		)
	actuators_stale =
		RB.ManualActuator(
			2,
			stale_cable_indices,
			restlens[stale_cable_indices],
			RB.Serial()
		)
	actuators_rings = [
		RB.ManualActuator(
			2+j,
				     [k*nhcables + CA((j-1)*5*n+1:j*5n)[i*5] for i = 1:n],
			restlens[[k*nhcables + CA((j-1)*5*n+1:j*5n)[i*5] for i = 1:n]],
			RB.Serial()
		)
		for k = 0:l
		for j = 1:m
	]

	hub = (actuators = vcat(actuators_active,actuators_stale,actuators_rings),)

	# hub = (actuators = vcat(actuators_helixes,actuators_rings),)

    bot = RB.Robot(st,hub)

end
# 参数：
# n：边数
# m：段数
# θ：段中内节旋转角度
# r：段端外节半径
# h：外节间距离
# c：绳索阻尼
# right：有否右对接桥
# o2：右对接桥原点
function nrailbridge(n,m=2;θ=missing,r=missing,c=0.0,h=1.0,o2=[0,4.0,0],right=false)
	l = Int(right)
    α = 2π/n
	# 杆长
    b = 1.35
	a = sqrt(b^2-h^2)
	# θ 和 r 只能给定其中一个
	if (θ isa Missing) ⊻ (r isa Missing)
		if θ isa Missing
			@assert abs(a/(2r)) < 1
			θ = 2asin(a/(2r)) - α
			@show rad2deg(θ)
		else
			r = a/2/sin((α+θ)/2)
			@show r
		end
	else
		error("Either θ or r should be specified.")
	end
	# 外节杆端点
    as = [
        SVector(
			r*cos(i*α+j*θ+π/2),j*h,
			r*sin(i*α+j*θ+π/2)
		)
        for i = 0:n-1, j = 0:m
    ]

	aa = norm(as[3,2]-as[1,1])
	s = sqrt(b^2 - (aa/2)^2)
	# 内节杆半径
	r1 = s - norm((as[3,2].+as[1,1])./2-[0,h/2,0])
	@show r, r1
	@assert r1 > 0
	# 内节杆端点
	bs = [
		SVector(
			r1*cos((i+1//2)*α+(j+1//2)*θ+π/2),(j+1/2)*h,
		    r1*sin((i+1//2)*α+(j+1//2)*θ+π/2)
		)
		for i = 0:n-1, j = 0:m-1
	]

	# 右对接桥
	rot180 = LinearMap(RotZ(π))
	translate2o2 = Translation(o2)
	mirror = LinearMap(
		SMatrix{3,3}([
		-1 0 0;
		 0 1 0;
		 0 0 1;
		])
	)
	transform2right = translate2o2 ∘ mirror ∘ rot180
	asr = transform2right.(as)
	bsr = transform2right.(bs)

	# 左右对接
	aslr = hcat(as,asr)
	bslr = hcat(bs,bsr)
	# 每边端点数
	nhalf = 3*n*m
	# 杆件
	rbs = reduce(vcat,[
		[
			make_strut(k*nhalf+(j-1)*3*n+(i-1)*3+1,CA(aslr[:,k*(m+1)+j  ])[i  ],
										  		   CA(aslr[:,k*(m+1)+j+1])[i+1];
										  		   ci=Int[]),
			make_strut(k*nhalf+(j-1)*3*n+(i-1)*3+2,CA(aslr[:,k*(m+1)+j])[i  ],
										  	   	   CA(bslr[:,k*m+j])[i-2];
										  		   ci=Int[]),
			make_strut(k*nhalf+(j-1)*3*n+(i-1)*3+3,CA(bslr[:,k*m+j  ])[i-2],
										  	       CA(aslr[:,k*(m+1)+j+1])[i+2];)
		]
		for k = 0:l
		for j = 1:m
		for i = 1:n
	])
	# 杆数
	nb = length(rbs)
	rigidbodies = TypeSortedCollection(rbs)
	numbered = RB.number(rigidbodies)
	# 共享坐标
	sharing_rows = [
		begin
			sharing_row = zeros(Int,ifelse(j==m,3*3,4*3),nb)
			ca = CA((j-1)*3n+1:    j*3n)
			cb = CA(    j*3n+1:(j+1)*3n)
			sharing_row[1:3,k*nhalf+ca[(i-1)*3+1]] .= 1:3
			sharing_row[1:3,k*nhalf+ca[(i-1)*3+2]] .= 1:3
			sharing_row[4:6,k*nhalf+ca[(i-1)*3+2]] .= 4:6
			sharing_row[4:6,k*nhalf+ca[(i-1)*3+3]] .= 1:3
			sharing_row[7:9,k*nhalf+ca[(i-1)*3+3]] .= 4:6
			sharing_row[7:9,k*nhalf+ca[(i-1)*3+4]] .= 4:6
			if j !== m
				sharing_row[10:12,k*nhalf+ca[(i-1)*3+1 ]] .= 4:6
				sharing_row[10:12,k*nhalf+cb[(i-1)*3+19]] .= 1:3
			end
			sharing_row
		end
		for k = 0:l
		for j = 1:m
		for i = 1:n
	]
	sharing = reduce(vcat,sharing_rows)
	# browse(matrix_sharing)
	indexed = RB.index(rigidbodies,sharing)

	# 单边绳索数
	nhcables = 5n*m
	# 总绳索数
	ncables = (l+1)*nhcables
	# 主动绳索编号
	active_cable_indices = [
		k*nhcables + CA((j-1)*5*n+1:j*5n)[(j-1)*5+(i-1)*5+v]
		for k = 0:l for j = 1:m for i = 1:n for v = 2:3
	]
	# 被动绳索编号
	stale_cable_indices = [
		k*nhcables + CA((j-1)*5*n+1:j*5n)[-(j-1)*5+(i-1)*5+v]
		for k = 0:l for j = 1:m for i = 1:n for v = [1,19]
	]

	restlens = fill(0.0,ncables)
	ks = fill(2000.0,ncables)
	ks[stale_cable_indices] .= 3000.0
    cs = fill(c,ncables)
    cables = [RB.Cable3D(i,restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    tensiles = (cables = cables,)
	# 绳索连接
	cnt_rows = [
		begin
			cnt_row = zeros(Int,5,nb)
			# cnt_row[1,CA(3n,(i-1)*3:1)] )[ 1; cnt_row[1,CA(3n,      i:3+1))[ = -1
			ca = CA((j-1)*3n+1:j*3n)
			cnt_row[1,k*nhalf + ca[  (i-1)*3+1]] =  1
			cnt_row[1,k*nhalf + ca[n+(i-1)*3+1]] = -1
			cnt_row[2,k*nhalf + ca[  (i-1)*3+1]] =  1
			cnt_row[2,k*nhalf + ca[n+(i-1)*3+4]] = -1
			cnt_row[3,k*nhalf + ca[  (i-1)*3+1]] =  2
			cnt_row[3,k*nhalf + ca[n+(i-1)*3+4]] = -1
			cnt_row[4,k*nhalf + ca[  (i-1)*3+1]] =  2
			cnt_row[4,k*nhalf + ca[n+(i-1)*3+7]] = -1

			cnt_row[5,k*nhalf + ca[  (i-1)*3+1]] =  2
			cnt_row[5,k*nhalf + ca[      i*3+1]] = -2

			cnt_row
		end
		for k in 0:l
		for j in 1:m
		for i in 1:n
	]
	cnt_cables = reduce(vcat,cnt_rows)
	# display(cnt_cables)
	connected = RB.connect(rigidbodies,cnt_cables)
	tensioned = @eponymtuple(connected,)

	# 导轨约束
	A_cst = zeros(2n,indexed.nfull)
	for i = 1:n
	    j = 3(i-1)+1
	    inds = indexed.mem2sysfull[j]
	    indx, indy, indz = inds[1:3]
		rb = rbs[j]
		q, _ = RB.NCF.rigidstate2naturalcoords(
						rb.state.cache.funcs.nmcs,
						rb.state.ro,
						rb.state.R,
						rb.state.ṙo,
						rb.state.ω)
	    A_cst[2(i-1)+1,indx], A_cst[2(i-1)+1,indz] = q[3], -q[1]
		A_cst[2(i-1)+2,indy] = 1
	end
	railcst = RB.LinearJoint(A_cst,zeros(2n))
	jointed = RB.join([railcst],indexed)

	# 连接信息
    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)

	# 张拉整体结构
    st = RB.Structure(rigidbodies,tensiles,cnt)

	# 作动器
	actuators_active =
		RB.ManualActuator(
			1,
			         active_cable_indices,
			restlens[active_cable_indices],
			RB.Serial()
		)
	actuators_stale =
		RB.ManualActuator(
			2,
			stale_cable_indices,
			restlens[stale_cable_indices],
			RB.Serial()
		)
	actuators_rings = [
		RB.ManualActuator(
			2+j,
				     [k*nhcables + CA((j-1)*5*n+1:j*5n)[i*5] for i = 1:n],
			restlens[[k*nhcables + CA((j-1)*5*n+1:j*5n)[i*5] for i = 1:n]],
			RB.Serial()
		)
		for k = 0:l
		for j = 1:m
	]

	hub = (actuators = vcat(actuators_active,actuators_stale,actuators_rings),)

	# hub = (actuators = vcat(actuators_helixes,actuators_rings),)

    bot = RB.Robot(st,hub)

end
