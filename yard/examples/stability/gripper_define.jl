
function rigidbody(i,m,Ia,ri,rj,aps)
	movable = true
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
	prop = RB.RigidBodyProperty(
				i,movable,m,
				Ī,
				mass_locus,
				aps;
				constrained
				)
	lncs, q, _ = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α, ṙo, ω)

	state = RB.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)

	body = RB.RigidBody(prop,state)
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
    ss = [RB.Cable2D(i, restlengths[i],ks[i],cs[i];slack=false) for i = 1:ncables]
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
	connected = RB.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	hub = nothing

    cnt = RB.Connectivity(numbered,indexed,tensioned)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    bot = RB.Robot(st,hub)
end

function skew(a)
	[-a[2],a[1]]
end

function build_tgg_nullspace(bot,q̌)
	(;st) = bot
	(;ndof) = st
	(;indexed) = st.connectivity
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

function make_halftggriper_nullspace(bot,q̌)
	(;st) = bot
	(;ndof) = st
	(;indexed) = st.connectivity
	(;nfree,mem2sysfree) = indexed
	function inner_Ň(q̌,c)
		ret = zeros(eltype(q̌),nfree,ndof)
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
		ret .= [
			1            0         0;
			[1,0] -skew(u2)     [0,0];
			[0,0] -skew(v2)     [0,0];
			[1,0] -skew(u2) -skew(u3);
			[0,0]     [0,0] -skew(v3);
		]
	end
end

function sup_ggriper!(ax,tgob,sgi;
		color = :slategrey,
		linewidth = 10
	)
	revert_y = [1  0;
				0 -1]
	bars = @lift begin
		rbs = RB.get_bodies($tgob)
		ndim = RB.get_ndim($tgob)
        T = RB.get_numbertype($tgob)
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
			Point(rbs[3].state.loci_states[1])=>
			Point(rbs[3].state.loci_states[3])
		)
		push!(ret,
			Point2(revert_y*rbs[3].state.loci_states[1])=>
			Point2(revert_y*rbs[3].state.loci_states[3])
		)

		ret
	end
	points = @lift [Point2f($tgob.state.rigids[2].q[1:2]) ]
	# @show points[]
    linesegments!(ax, bars; color, linewidth)
	if sgi in [2]
		arrows!(ax, points, [Vec2f(-1,0)]; 
			color = :red, 
			align = :tailend, 
			arrowsize = 50,
			lengthscale = sgi*20,
			linewidth = 8
		)
	end
end

function make_spine(n,θ=0.0)
    a = 0.1
    b = 0.1*√2
    α = 3/4*π
    d = 0.15

    function vert(i,r,θ,a,b,α)
        if i == 1
            movable = false
            constrained = true
            ci = collect(1:6)
			Φi = Int[]
        else
            movable = true
            constrained = false
            ci = Int[]
			Φi = collect(1:3)
	        # ap1 = b*[cos( α),sin( α)]
	        # ap2 = b*[cos( α),sin( α)]
	        # ap3 = b*[cos(-α),sin(-α)]
	        # ap4 = b*[cos(-α),sin(-α)]
        end		
		ap1 = [a,0.0]
		ap2 = b*[cos( α),sin( α)]
		ap3 = b*[cos(-α),sin(-α)]
		# ap4 = [a,0.0]
        mass_locus = zeros(2)
        aps = [ap1,ap2,ap3]
        m = 0.495 #kg
        inertia = SMatrix{2,2}(Diagonal(ones(2)))
        prop = RB.RigidBodyProperty(
			i,
			movable,m,inertia,mass_locus,aps;
			constrained
		)

        ri = SVector{2}(r)
        ro = SVector{2}(r)
        ω = 0.0
        ṙo = @SVector zeros(2)
        lncs,_,_ = RB.NCF.NC1P2V(ri,ro,θ,ṙo,ω)
        state = RB.RigidBodyState(prop,lncs,ro,θ,ṙo,ω, ci, Φi)
        RB.RigidBody(prop,state)
    end

    rs = [[i*d,0.0] for i = 0:n-1]
    θs = [i*θ for i = 0:n-1]
    rbs = [vert(i,rs[i],θs[i],a,b,α) for i = 1:n]

	rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
    indexed = RB.index(rigdibodies)

    ncables = 4*(n-1)
    k = 840.0 #N/m
    c = 100.0
    original_restlens = repeat([0.15, 0.11, 0.11, 0.15]./2,n-1)
    ss = [RB.Cable2D(i,original_restlens[i],k,c) for i = 1:ncables]
    tensiles = (cables = ss,)
    hub = nothing
	

	matrix_cnt = zeros(Int,ncables,n)
	for i = 1:n-1
		matrix_cnt[4(i-1)+1,[i,i+1]] = [2,-2]
		matrix_cnt[4(i-1)+2,[i,i+1]] = [3,-3]
		matrix_cnt[4(i-1)+3,[i,i+1]] = [1,-2]
		matrix_cnt[4(i-1)+4,[i,i+1]] = [1,-3]
	end
	matrix_cnt |> display
	connected = RB.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

    cnt = RB.Connectivity(numbered,indexed,tensioned)
	# rbs
    st = RB.Structure(rigdibodies,tensiles,cnt)
    bot = RB.Robot(st,hub)
end

function sup_spine2d!(ax,tgob,sgi)
	bars = @lift begin
		ndim = RB.get_ndim($tgob)
        T = RB.get_numbertype($tgob)
        ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
		foreach($tgob.rigidbodies) do body
			
			push!(ret,
				Point(body.state.ro)=>
				Point(body.state.loci_states[1])
			)
			push!(ret,
				Point(body.state.ro)=>
				Point(body.state.loci_states[2])
			)
			push!(ret,
				Point(body.state.ro)=>
				Point(body.state.loci_states[3])
			)

		end
		ret
	end
    linesegments!(ax, bars, linewidth = 10)
end

function dualtri(ndof,onedir=[1.0,0.0];θ=0.0,k=400.0,c=0.0,restlen=0.16)
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
        mass_locus = SVector{2}([0.0,0.0])

        ap1 = SVector{2}([-a/2,0.0])
        ap2 = SVector{2}([ a/2,0.0])
		if isodd(i)
			ap3 = SVector{2}([0.0,-√3/2*a])
		else
        	ap3 = SVector{2}([0.0,√3/2*a])
		end

        aps = [ap1,ap2,ap3]

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

		prop = RB.RigidBodyProperty(
					i,movable,m,
					Ī,
                    mass_locus,
					aps;
					constrained=constrained
                    )

		lncs, q, _ = RB.NCF.NC1P2V(SVector{2}(ri), ro, α, ṙo, ω)

		state = RB.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)

        body = RB.RigidBody(prop,state)
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
    numbered = RB.number(rigdibodies)

   
    indexed = RB.index(rigdibodies)

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
    ss = [RB.Cable2D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
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

	connected = RB.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	function ganged_act(actid,id1,id2,original_restlens)
		ids = [id1,id2]
		original_values = original_restlens[[id1,id2]]
		RB.ManualActuator(actid,ids,original_values,RB.Ganged())
	end
    acs = [ifelse(isodd(i),ganged_act(i,2(i-1)+1,2i,original_restlens),
						   ganged_act(i,2i,2(i-1)+1,original_restlens)) for i = 1:nbodies-1]
	hub = (actuators=acs,)
	pjs = [
		RB.PinJoint(RB.End2End(i,RB.ID(rbs[i],2),RB.ID(rbs[i+1],1)))
		for i = 1:nbodies-1
	]
	jointed = RB.join(pjs,indexed)

    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,hub)
end


function sup_dualtri!(ax,tgob;
		color = :black,
		linewidth = 10
	)
	bars = @lift begin
		ndim = RB.get_ndim($tgob)
        T = RB.get_numbertype($tgob)
        ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
		foreach($tgob.rigidbodies) do body
			
			push!(ret,
				Point(body.state.loci_states[1])=>
				Point(body.state.loci_states[2])
			)
			push!(ret,
				Point(body.state.loci_states[2])=>
				Point(body.state.loci_states[3])
			)
			push!(ret,
				Point(body.state.loci_states[3])=>
				Point(body.state.loci_states[1])
			)

		end
		ret
	end
    linesegments!(ax, bars; color, linewidth)
end
function get_angle(v,w)
	atan(w[2]*v[1]-w[1]*v[2],w[1]*v[1]+w[2]*v[2])
end

function get_angles(st)
    rbs = RB.get_bodies(st)
    angles = zeros(st.nrigids-1)
    for (bodyid,rb) in enumerate(rbs)
        if bodyid > 1
			state0 = rbs[bodyid-1].state
            v = state0.loci_states[2]-state0.loci_states[1]
            state1 = rbs[bodyid].state
            w = state1.loci_states[2]-state1.loci_states[1]
            angles[bodyid-1] = get_angle(v,w)
        end
    end
    rad2deg.(angles)
end

