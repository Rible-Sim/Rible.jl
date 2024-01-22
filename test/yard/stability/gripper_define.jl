function skew(a)
	[-a[2],a[1]]
end

function build_tgg_nullspace(bot,q̌)
	(;st) = bot
	(;num_of_dof) = st
	(;indexed) = st.connectivity
	(;num_of_free_coords,bodyid2sys_free_coords) = indexed
	q̌2 = @view q̌[bodyid2sys_free_coords[2]]
	q̌3 = @view q̌[bodyid2sys_free_coords[3]]
	x2 = q̌2[1]
	r2i = [x2,0.0]
	r2j = q̌2[2:3]
	u2 = r2j - r2i
	v2 = q̌2[4:5]
	r3j = q̌3[2:3]
	u3 = r3j - r2i
	v3 = q̌3[4:5]
	ret = zeros(eltype(q̌),num_of_free_coords,num_of_dof)
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
	(;num_of_dof) = st
	(;indexed) = st.connectivity
	(;num_of_free_coords,bodyid2sys_free_coords) = indexed
	function inner_Ň(q̌,c)
		ret = zeros(eltype(q̌),num_of_free_coords,num_of_dof)
		q̌2 = @view q̌[bodyid2sys_free_coords[2]]
		q̌3 = @view q̌[bodyid2sys_free_coords[3]]
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
		ndim = RB.get_num_of_dims($tgob)
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

function sup_spine2d!(ax,tgob,sgi)
	bars = @lift begin
		ndim = RB.get_num_of_dims($tgob)
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

function dualtri(num_of_dof,onedir=[1.0,0.0];θ=0.0,k=400.0,c=0.0,restlen=0.16)
    nbodies = num_of_dof + 1
    nbp = 2nbodies - num_of_dof
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
		contactable = true
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
			visible = true
			ci = collect(1:6)
			cstr_idx = Int[]
		else
			visible = true
			ci = Int[]
			cstr_idx = collect(1:3)
		end
		
		ω = 0.0
		ro = ri
		ṙo = zero(ro)

		prop = RB.RigidBodyProperty(
					i,contactable,m,
					Ī,
                    mass_locus,
					aps;
					visible=visible
                    )

		nmcs = RB.NCF.NC1P2V(SVector{2}(ri), ro, α)

		state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)

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
    ss = [RB.DistanceSpringDamper2D( original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
	apparatuses = (cables=ss,)

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
		RB.PinJoint(RB.Hen2Egg(RB.Signifier(rbs[i],2),RB.Signifier(rbs[i+1],1)))
		for i = 1:nbodies-1
	]
	jointed = RB.join(pjs,indexed)

    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,apparatuses,cnt)
    RB.Robot(st,hub)
end


function sup_dualtri!(ax,tgob;
		color = :black,
		linewidth = 10
	)
	bars = @lift begin
		ndim = RB.get_num_of_dims($tgob)
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

