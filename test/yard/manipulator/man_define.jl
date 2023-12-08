
function man_ndof(num_of_dof,onedir=[1.0,0.0];θ=0.0,k=0.0,c=0.0,unit="mks",restlen=0.16,isvirtual=true)
    nbodies = num_of_dof + 1
    nbp = 2nbodies - num_of_dof
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    lower_index = 1:2:nbodies
    upper_index = 2:2:nbodies
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
	for (i,j) in enumerate(lower_index)
        # a[j] = 20.0e-2unit_L
		# if j==nbodies
		# 	a[j] = 12252e-5unit_L
		# 	m[j] = 90.153058e-3unit_M
		# 	Ia[j] = 2580.2e-7unit_I
		# else
	        a[j] = 12252e-5unit_L
	        #m[j] = 835.90254985e-3unit_M
	        m[j] = 200.40071778e-3unit_M
	        # Ia[j] = Ic_lower[i] + m[j]*1/3*a[j]^2
	        #Ia[j] = 28130.53053840*2e-7unit_I
	        Ia[j] = 3738.8e-7unit_I  #??2e-7??
		# end
        # @show a[j],m[j],Ia[j]
    end
    for (i,k) in enumerate(upper_index)
		# if k==nbodies
		# 	a[k] = 12252e-5unit_L
		# 	m[k] = 90.153058e-3unit_M
		# 	Ia[k] = 2580.2e-7unit_I
		# else
	        # a[k] = 16.0e-2unit_L
	        a[k] = 12252e-5unit_L
	        #m[k] = 666.25659673e-3unit_M
	        m[k] = 200.40071778e-3unit_M
	        # Ia[k] = Ic_upper[i] + m[k]*1/3*a[k]^2
	        Ia[k] = 3738.8e-7unit_I
	        # @show a[k],m[k],Ia[k]
		# end
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

    function rigidbody(i,m,a,Ia,ri,rj)
		contactable = true
		# if i == nbodies
		# 	CoM_x = 0.05933624
		# 	CoM_y = 0.03109915
		# else
        CoM_x = 0.06535974
        CoM_y = 0.03465684
		# end

        if isodd(i)
            CoM_y = -CoM_y
        end
        mass_locus = SVector{2}([CoM_x,CoM_y])

        # nap = 3 #?
		sti_l = 0.01832
		rao_l = 0.00475
        ap1 = SVector{2}([0.0,0.0])
        ap2 = SVector{2}([a,0.0])
        ap3_x = a/2
        ap3_y = √3/2*a
		ap4_x = a-sti_l*cos(pi/6)
		ap4_y = sti_l*sin(pi/6)
		ap5_x = a/2
		ap5_y = √3/2*a-sti_l
		ap6_x = -rao_l*cos(pi/6)
		ap6_y = -rao_l*sin(pi/6)
		ap7_x = a/2
		ap7_y = √3/2*a+rao_l

        if isodd(i)
            ap3_y = -ap3_y
			ap4_y = -ap4_y
			ap5_y = -ap5_y
			ap6_y = -ap6_y
			ap7_y = -ap7_y

        end
        ap3 = SVector{2}([ap3_x,ap3_y])
		ap4 = SVector{2}([ap4_x,ap4_y])
		ap5 = SVector{2}([ap5_x,ap5_y])
		ap6 = SVector{2}([ap6_x,ap6_y])
		ap7 = SVector{2}([ap7_x,ap7_y])

        aps = [ap1,ap2,ap3,ap4,ap5,ap6,ap7]

		Ī = SMatrix{2, 2}([
			0.99Ia 0
			0 0.01Ia
		])
		# only to get nmcs
		visible = true
        ci = Int[]
        prop = RB.RigidBodyProperty(
					i,contactable,m,
					Ī,
                    mass_locus,
					aps;
					visible=visible
                    )
		α = 0.0
		ω = 0.0
		ro = ri
		ṙo = zero(ro)

		visible = true
		ci = Int[]
		cstr_idx =  collect(1:3)
		prop = RB.RigidBodyProperty(
					i,contactable,m,
					Ī,
					mass_locus,
					aps;
					visible=ifelse(i in [1,2],true,false)
					)
		nmcs = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α)

        if i in [1,2]
			visible = true
			if i == 1
            	ci = RB.find_full_pres_idx(nmcs, q)
				display(ci)
			else
				ci = [1]
			end
		end
		state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)

        body = RB.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,m[i],a[i],
            Ia[i],A[:,i],A[:,i+1]) for i = 1:nbodies]

	rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)

    matrix_sharing = zeros(Int,2nbodies-2,nbodies)
	for i = 1:nbodies-1
		matrix_sharing[2(i-1)+1,i:i+1] = [3, 1]
		matrix_sharing[2(i-1)+2,i:i+1] = [4, 2]
	end
	# display(matrix_sharing)
    indexed = RB.index(rigdibodies, matrix_sharing)

    ncables = 2(nbodies-1)
    # upstringlen = 0.182-5.2024/292.14
    # lostringlen = 0.182-5.2024/292.14
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
    ss = [RB.DistanceSpringDamper2D(i, original_restlens[i],ks[i],cs[i]) for i = 1:ncables]
	tensiles = (cables=ss,)

	matrix_cnt = zeros(Int,2(nbodies-1),nbodies)
    for i = 1:nbodies-1
		if isvirtual
			ul = [1,-3]
			dl = [3,-2]
		else
			ul = [6,-5]
			dl = [7,-4]
		end
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

    cnt = RB.Connectivity(numbered,indexed,tensioned)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    bot = RB.Robot(st,hub)
end



# function build_Y(tr)
# 	@unpack st, hub = tr
# 	@unpack actuators = hub
#     @unpack ncables,cables = st
#     nact = length(actuators)
#     ret = spzeros(Int,ncables,nact)
#     for (i,ctrller) in enumerate(actuators)
#         is1, is2 = ctrller.acts.id_cables
#         ret[is1,i] = 1
#         ret[is2,i] = -1
#     end
#     ret
# end

function get_angle(v,w)
	atan(w[2]*v[1]-w[1]*v[2],w[1]*v[1]+w[2]*v[2])
end

function get_angles(bot)
    (;st) = bot
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
