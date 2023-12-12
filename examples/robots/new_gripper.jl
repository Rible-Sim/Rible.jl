

function new_gripper(;c=0.0,x1 = 41.02, ϕ1 = 0.7652, ϕ2 = 0.0)
	function rigidbody(i,m,Ia,ri,rj,aps)
		contactable = true
		# end
		u = rj - ri
		d = norm(u)
		mass_locus = SVector{2}([d/2,0])
	
		Ī = SMatrix{2, 2}([
			0.99Ia 0
			0 0.01Ia
		])
	
		α = atan(u[2],u[1])
		ω = 0.0
		ro = ri
		ṙo = zero(ro)
	
	
		if i in [1,2]
			visible = true
			if i == 1
				ci = collect(1:6)
				cstr_idx = Int[]
			else
				ci = [2]
				cstr_idx = collect(1:3)
			end
		else
			visible = true
			ci = Int[]
			cstr_idx = collect(1:3)
		end
		prop = RB.RigidBodyProperty(
					i,contactable,m,
					Ī,
					mass_locus,
					aps;
					visible
					)
		nmcs = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α)
	
		state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)
	
		body = RB.RigidBody(prop,state)
	end
	
    barlengths = [140,50,86.29165124598852]
	ext2 = 15
	rs = SVector{2}.([
		[0.0           ,0.0],
		[barlengths[1] ,0.0],
		[x1            ,0.0],
		[x1+(barlengths[2]+ext2)*cos(-ϕ1),
		(barlengths[2]+ext2)*sin(-ϕ1)],
		[x1+(barlengths[2]+ext2)*cos(-ϕ1)+barlengths[3]*cos(-ϕ2),
		(barlengths[2]+ext2)*sin(-ϕ1)+barlengths[3]*sin(-ϕ2)],
	])


	m = 100.0
	Ia = 1/12*m*barlengths[2]^2

	ap1 = SVector{2}([          0.0,0.0])
	ap2 = SVector{2}([barlengths[1],0.0])
	aps = [ap1,ap2]
    rb1 = rigidbody(1,m,Ia,rs[1],rs[2],aps)
	ap1 = SVector{2}([          0.0,0.0])
	ap2 = SVector{2}([barlengths[2],0.0])
	aps = [ap1,ap2]
    rb2 = rigidbody(2,m,Ia,rs[3],rs[4],aps)
	ap1 = SVector{2}([          -20,0.0])
	ap2 = SVector{2}([           15,0.0])
	ap3 = SVector{2}([barlengths[3],0.0])
	aps = [ap1,ap2,ap3]
    rb3 = rigidbody(3,m,Ia,rs[4],rs[5],aps)
	rbs = [rb1,rb2,rb3]

	rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)

    matrix_sharing = [
		0 3 1;
		0 4 2;
	]
    indexed = RB.index(rigdibodies, matrix_sharing)
    # indexed = RB.index(rigdibodies)


	# restlengths = [59,59,64.2,64.2,36.7,36.7]
    # ks = [0.482,0.482,0.568,0.568,0.256,0.256]
	restlengths = [59,64.2,36.7,36.7,18.0,50.0]
    ks = [0.482,0.568,0.256/2,0.256/2,0.03,0.01]
	ncables = length(ks)
    cs = zeros(ncables)
    ss = [RB.DistanceSpringDamper2D( restlengths[i],ks[i],cs[i];slack=false) for i = 1:ncables]
	apparatuses = (cables=ss,)

	# matrix_cnt = [
	# 	1 -2  0
	# 	1  0 -2
	# 	2 -2  0
	# 	2  0 -2
	# 	1 -1  0
	# 	2 -1  0
	# ]
	matrix_cnt = [
		1 -2  0;
		2 -2  0;
		1 -1  0;
		2 -1  0;
		0  2 -1;
		2  0 -2
	]
	connected = RB.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	hub = nothing

    cnt = RB.Connectivity(numbered,indexed,tensioned)
    st = RB.Structure(rigdibodies,apparatuses,cnt)
    bot = RB.Robot(st,hub)
end