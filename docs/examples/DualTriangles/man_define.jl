
function man_ndof(ndof,onedir=[0.0,-1.0];θ=0.0,k=0.0,c=0.0,unit="mks",restlen=0.16,isvirtual=true)
    nbodies = ndof + 1
    nbp = nbodies + 1
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    ends_indices = [1,nbodies]
    inner_indices = 2:nbodies-1
    a = zeros(nbodies)
    m = zeros(nbodies)
    Ia = zeros(nbodies)
    if unit == "cgs"
        unit_L = 1e2
        unit_M = 1e3
        unit_I = unit_M*unit_L^2
    else # "mks"
        unit_L = 1
        unit_M = 1
        unit_I = 1
    end
	for j in ends_indices
        a[j] = 12252e-5unit_L
		if isvirtual
			m[j] = 94.873e-3unit_M
		else
        	m[j] = 198.553e-3unit_M
		end
        Ia[j] = 3068.8e-7unit_I  #??2e-7??
    end
    for k in inner_indices
        a[k] = 12252e-5unit_L
		if isvirtual
        	m[k] = 2*94.873e-3unit_M
		else
			m[k] = 294.873e-3unit_M
		end
        Ia[k] = 2*3738.8e-7unit_I
    end
    R = [
		cos(θ) -sin(θ);
		sin(θ) cos(θ)
	]
    A = zeros(2,nbp)
    A[:,2] .= A[:,1] .+ √3/2*a[1]*onedir
    for i in 3:nbp
        A[:,i] .= A[:,i-1] .+ √3*a[i-1]*R^(i-2)*onedir
    end
	B = zeros(2,nbodies)
	B[:,1] = A[:,1]
	for i in 2:nbodies
        B[:,i] .= (A[:,i] .+ A[:,i+1])./2
    end

    function rigidbody(i,m,a,Ī,ri,rj)

		movable = true
		sti_l = 0.01832
		rao_l = 0.00475
		ap3 = SVector{2}([0.0, a/2])
		ap4 = SVector{2}([0.0,-a/2])
		if i == 1
			CoM_x = a/(2*√3)
	        CoM_y = 0.0
			ap1 = SVector{2}([   0.0,0.0])
	        ap2 = SVector{2}([√3/2*a,0.0])
			ap4_y = -a/2
			ap5_x = -rao_l*sin(pi/6)
			ap5_y = a/2+rao_l*cos(pi/6)
			ap6_x = -rao_l*sin(pi/6)
			ap6_y = -a/2-rao_l*cos(pi/6)
			ap7_x = sti_l*sin(pi/6)
			ap7_y = a/2-sti_l*cos(pi/6)
			ap8_x = sti_l*sin(pi/6)
			ap8_y = -a/2+sti_l*cos(pi/6)
		elseif i == nbodies
			if isvirtual
				CoM_x = 0.05704695*sin(pi/3) + 0.02546222*cos(pi/3) - √3/2*a
				CoM_y = 0.05704695*cos(pi/3) - 0.02546222*sin(pi/3)
			else
				CoM_x = 0.05704695*sin(pi/3) + 0.02546222*cos(pi/3) - √3/2*a
				CoM_y = 0.05704695*cos(pi/3) - 0.02546222*sin(pi/3)
			end
	        ap1 = SVector{2}([-√3/2*a,0.0])
			ap2 = SVector{2}([    0.0,0.0])
			ap5_x = rao_l*sin(pi/6)
			ap5_y = a/2+rao_l*cos(pi/6)
			ap6_x = rao_l*sin(pi/6)
			ap6_y = -a/2-rao_l*cos(pi/6)
			ap7_x = -sti_l*sin(pi/6)
			ap7_y = a/2-sti_l*cos(pi/6)
			ap8_x = -sti_l*sin(pi/6)
			ap8_y = -a/2+sti_l*cos(pi/6)
		else
	        CoM_x = 0.10215556 - √3/2*a
	        CoM_y =  0.0020302
			ap1 = SVector{2}([-√3/2*a,0.0])
	        ap2 = SVector{2}([ √3/2*a,0.0])
			ap5_x = 0.0
			ap5_y = a/2+rao_l
			ap6_x = 0.0
			ap6_y = -a/2-rao_l
			ap7_x = 0.0
			ap7_y = a/2-sti_l
			ap8_x = 0.0
			ap8_y = -a/2+sti_l
		end

        r̄g = SVector{2}([CoM_x,CoM_y])
		ap5 = SVector{2}([ap5_x,ap5_y])
		ap6 = SVector{2}([ap6_x,ap6_y])

		ap7 = SVector{2}([ap7_x,ap7_y])
		ap8 = SVector{2}([ap8_x,ap8_y])
        aps = [ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8]

		# only to get lncs
        prop = TR.RigidBodyProperty(
					i,movable,m,
					Ī,
                    r̄g,
					aps;
					constrained=ifelse(i==1,true,false)
                )
        ro = copy(ri)
        ω = 0.0
        ṙo = @SVector zeros(2)
		α = get_angle([1.0,0.0],rj-ri)
        lncs, q0, _ = TR.NaturalCoordinates.NC1P2V(ri,ro,α,ṙo,ω)
		if i == 1
			ci = collect(1:6)
			Φi = Int[]
		else
			ci = Int[]
			Φi = collect(1:3)
		end
        state = TR.RigidBodyState(prop,lncs,ro,α,ṙo,ω,ci,Φi)

        rb = TR.RigidBody(prop,state)
    end
    rbs = [
		rigidbody(i,m[i],a[i],
		SMatrix{2,2}(Diagonal([Ia[i]/2,Ia[i]]/2)),
		SVector{2}(B[:,i]),SVector{2}(A[:,i+1])) for i = 1:nbodies
	]
	rigdibodies = TypeSortedCollection(rbs)
    numbered = TR.number(rigdibodies)
    indexed = TR.index(rigdibodies,)

    nstrings = 2(nbodies-1)
    # upstringlen = 0.182-5.2024/292.14
    # lostringlen = 0.182-5.2024/292.14
    upstringlen = restlen
    lostringlen = restlen
    original_restlens = zeros(nstrings)
    restlens = zeros(nstrings)
    actuallengths = zeros(nstrings)
    ks = zeros(nstrings)
    cs = zeros(nstrings)
    cs .= c
	stringlen = [0.1898-14.5/k,0.1769-14.5/k]
	snum = 1
	for i = 1:nstrings
		if isvirtual
	        j = i % 4
	        original_restlens[i] = ifelse(j∈[1,0],upstringlen,lostringlen)
	        ks[i] = ifelse(j∈[1,0],k,k)
		else
			original_restlens[i] = stringlen[snum]
			if i%2 == 0
				snum = snum + 1
			end
			ks[i] = k
		end
    end

	# ks[1:2] .= 1e10
	# ks[end-1:end] .= 1e-8
    ss = [TR.Cable2D(i, original_restlens[i],ks[i],cs[i]) for i = 1:nstrings]
	tensiles = (cables=ss,)

	matrix_cnt = zeros(Int,2(nbodies-1),nbodies)
	for i = 1:nbodies-1
		if isvirtual
			matrix_cnt[2(i-1)+1,i:i+1] = [3,-3]
			matrix_cnt[2(i-1)+2,i:i+1] = [4,-4]
		else
			matrix_cnt[2(i-1)+1,i:i+1] = [5,-7]
			matrix_cnt[2(i-1)+2,i:i+1] = [6,-8]
		end
	end

	connected = TR.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

	function ganged_act(actid,id1,id2,original_restlens)
		ids = [id1,id2]
		original_values = original_restlens[[id1,id2]]
		TR.ManualActuator(actid,ids,original_values,TR.Ganged())
	end
    acs = [ifelse(isodd(i),ganged_act(i,2(i-1)+1,2i,original_restlens),
						   ganged_act(i,2i,2(i-1)+1,original_restlens)) for i = 1:nbodies-1]
	hub = (actuators=acs,)

	pinjoints = [
		TR.PinJoint(
			TR.End2End(
				i-1,TR.ID(rbs[i-1],2),TR.ID(rbs[i],1)
			)
		)
		for i in 2:nbodies
	]

	jointed = TR.join(pinjoints,indexed)

    cnt = TR.Connectivity(numbered,indexed,tensioned,jointed)
    tg = TR.TensegrityStructure(rigdibodies,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,hub)
end



function get_angle(v,w)
	atan(w[2]*v[1]-w[1]*v[2],w[1]*v[1]+w[2]*v[2])
end

function get_angles(bot)
    (;tg) = bot
    rbs = TR.get_rigidbodies(tg)
    angles = zeros(tg.nrigids-1)
    for (rbid,rb) in enumerate(rbs)
        if rbid > 1
			state0 = rbs[rbid-1].state
            v = state0.rps[2]-state0.rps[1]
            state1 = rbs[rbid].state
            w = state1.rps[2]-state1.rps[1]
            angles[rbid-1] = get_angle(v,w)
        end
    end
    rad2deg.(angles)
end
