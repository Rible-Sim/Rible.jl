
function rigidbody(i,m,Ia,ri,rj,aps)
	movable = true
	# end
	u = rj - ri
	d = norm(u)
	r̄g = SVector{2}([d/2,0])

	Ī = SMatrix{2, 2}([
		0.99Ia 0
		0 0.01Ia
	])

	α = atan(u[2],u[1])
	ω = 0.0
	ro = ri
	ṙo = zero(ro)


	if i in [1,2]
		constrained = true
		if i == 1
			ci = collect(1:6)
			Φi = Int[]
		else
			ci = [2]
			Φi = collect(1:3)
		end
	else
		constrained = false
		ci = Int[]
		Φi = collect(1:3)
	end
	prop = TR.RigidBodyProperty(
				i,movable,m,
				Ī,
				r̄g,
				aps;
				constrained
				)
	lncs, q, _ = TR.NaturalCoordinates.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α, ṙo, ω)

	state = TR.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)

	rb = TR.RigidBody(prop,state)
end

function new_gripper(;c=0.0,x1 = 41.02, ϕ1 = 0.7652, ϕ2 = 0.0)
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
    numbered = TR.number(rigdibodies)

    matrix_sharing = [
		0 3 1;
		0 4 2;
	]
    indexed = TR.index(rigdibodies, matrix_sharing)
    # indexed = TR.index(rigdibodies)


	# restlengths = [59,59,64.2,64.2,36.7,36.7]
    # ks = [0.482,0.482,0.568,0.568,0.256,0.256]
	restlengths = [59,64.2,36.7,36.7,18.0,50.0]
    ks = [0.482,0.568,0.256/2,0.256/2,0.03,0.01]
	ncables = length(ks)
    cs = zeros(ncables)
    ss = [TR.Cable2D(i, restlengths[i],ks[i],cs[i];slack=false) for i = 1:ncables]
	tensiles = (cables=ss,)

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
	connected = TR.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	hub = nothing

    cnt = TR.Connectivity(numbered,indexed,tensioned)
    tg = TR.TensegrityStructure(rigdibodies,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
end


function skew(a)
	[-a[2],a[1]]
end

function build_tgg_nullspace(bot,q̌)
	(;tg) = bot
	(;ndof) = tg
	(;indexed) = tg.connectivity
	(;nfree,mem2sysfree) = indexed
	q̌2 = @view q̌[mem2sysfree[2]]
	q̌3 = @view q̌[mem2sysfree[3]]
	x2 = q̌2[1]
	r2i = [x2,0.0]
	r2j = q̌2[2:3]
	u2 = r2j - r2i
	v2 = q̌2[4:5]
	r3j = q̌3[2:3]
	u3 = r3j - r2i
	v3 = q̌3[4:5]
	ret = zeros(eltype(q̌),nfree,ndof)
	ret .= [
		 1            0  0;
		[1,0] -skew(u2)  [0,0];
		[0,0] -skew(v2)  [0,0];
	    [1,0]     [0,0] -skew(u3);
	    [0,0]     [0,0] -skew(v3);
	]
end

function build_halftggriper_nullspace(bot,q̌)
	(;tg) = bot
	(;ndof) = tg
	(;indexed) = tg.connectivity
	(;nfree,mem2sysfree) = indexed
	q̌2 = @view q̌[mem2sysfree[2]]
	q̌3 = @view q̌[mem2sysfree[3]]
	x2 = q̌2[1]
	r2i = [x2,0.0]
	r2j = q̌2[2:3]
	u2 = r2j - r2i
	v2 = q̌2[4:5]
	r3j = q̌3[3:4]
	u3 = r3j - r2j
	v3 = q̌3[5:6]
	ret = zeros(eltype(q̌),nfree,ndof)
	ret .= [
		 1            0         0;
		[1,0] -skew(u2)     [0,0];
		[0,0] -skew(v2)     [0,0];
	    [1,0] -skew(u2) -skew(u3);
	    [0,0]     [0,0] -skew(v3);
	]
end

function sup!(ax,tgob)
	revert_y = [1  0;
				0 -1]
	bars = @lift begin
		rbs = TR.get_rigidbodies($tgob)
		ndim = TR.get_ndim($tgob)
        T = TR.get_numbertype($tgob)
        ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
		push!(ret,
			Point2($tgob.state.rigids[1].q[1:2])=>
			Point2($tgob.state.rigids[1].q[3:4])
		)
		push!(ret,
			Point2($tgob.state.rigids[2].q[1:2])=>
			Point2($tgob.state.rigids[2].q[3:4])
		)
		push!(ret,
			Point2($tgob.state.rigids[2].q[1:2])=>
			Point2(revert_y*$tgob.state.rigids[2].q[3:4]
			)
		)
		push!(ret,
			Point(rbs[3].state.rps[1])=>
			Point(rbs[3].state.rps[3])
		)
		push!(ret,
			Point2(revert_y*rbs[3].state.rps[1])=>
			Point2(revert_y*rbs[3].state.rps[3])
		)

		ret
	end
    linesegments!(ax, bars, color = :slategrey, linewidth = 10)
end
