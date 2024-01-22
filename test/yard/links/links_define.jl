function move(trans,rot,o,X)
    x = @view X[1:3,1]
    y = @view X[1:3,2]
    z = @view X[1:3,3]
    f = trans ∘ rot
    f(o),hcat(rot(x),rot(y),rot(z))
end

function bend(n,di,rm)
    o1 = @SVector zeros(3)
    X1 = Matrix(1.0I,3,3)
    rot = LinearMap(rm)
    trans = Translation(0,0,di)

    os = [o1]
    Xs = [X1]
    for i = 2:n
        oi,Xi = move(trans,rot,os[end],Xs[end])
        push!(os,oi)
        push!(Xs,Xi)
    end
    os,Xs
end

function links(n,di,rm;k=3e1,c=0.0)
    nbodies = n
    nbp = 4*n

    a = 8.0e-2 #m
    h = 10.0e-2 #m
    θ = 2π/3
    l = √(a^2+h^2)
    b = √3a

    mass = 93.562831e-3 #kg
    # inertia = [2175.129	1.164e-05 6.58e-06;
	# 		   1.164e-05 2175.129 8.24e-06;
    #            6.58e-06	8.24e-06 2460.3283]
   	inertia = Matrix(Diagonal([2175.129,2175.129,2460.3283]))*1.e-7
    mass_locus = [0.0, 0.0, 7.5932384]*1.e-2


    ap1 = SVector{3}([0.0, 0.0, 0.0])
    ap2 = SVector{3}([a,   0.0,   h])
    ap3 = SVector{3}([a*cos(θ),  a*sin(θ), h])
    ap4 = SVector{3}([a*cos(θ), -a*sin(θ), h])
    aps = [ap1,ap2,ap3,ap4]

    contactable = ones(Bool,n)
	# contactable[1] = false
    visible = zeros(Bool,n)
    visible[1] = true

    props = [RB.RigidBodyProperty(i,contactable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(mass_locus),
                aps;visible=visible[i]) for i = 1:n]

    rs,Rs = bend(n,di,rm)
    ṙs = [zeros(3) for i = 1:n]
    ωs = [zeros(3) for i = 1:n]

    function rigidbody(i,prop,aps,ro,R,ṙo,ω)
        ri,rj,rk,rl = [ro+R*ap for ap in aps]
        # nmcs = RB.NCF.NC4P(ri,rj,rk,rl,ro,R)
        # nmcs = RB.NCF.NC3P1V(ri,rj,rk,ro,R)
		# nmcs = RB.NCF.NC2P2V(ri,rj,rk-ri,rl-ri,ro,R)
		nmcs = RB.NCF.NC1P3V(ri,ro,R)
		if i == 1
			pres_idx = RB.find_full_pres_idx(nmcs,q)
		else
			pres_idx = Int[]
		end
        state = RB.RigidBodyState(prop,ro,R,ṙo,ω,q,q̇,pres_idx)
        body = RB.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,props[i],aps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]

    nstrings = 6*(n-1)
	restlen = 3e-2
    stringlenH = restlen
    stringlenR = restlen #m
    stringlens = repeat(vcat(fill(stringlenH,3),fill(stringlenR,3)),n-1)
    kH = k #N/m
    kR = k #N/m
    ks = repeat(vcat(fill(kH,3),fill(kR,3)),n-1)
    # c = 0.0
    cs = repeat(fill(c,6),n-1)
    ss = [RB.SString3D(i,stringlens[i],ks[i],cs[i]) for i = 1:nstrings]
	apparatuses = (strings=ss,)
	acs = [RB.ManualActuator(RB.SimpleRegistor(6(i-1)+j,stringlens[6(i-1)+j])) for i = 1:n-1  for j = 1:6]
	hub = (actuators=acs,)
	bodynq = RB.get_num_of_coords(rbs[1])
    bodyid2q_raw = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]
	# bodyid2q = RB.filter_bodyid2q(bodyid2q_raw,rbs)
	bodyid2q = bodyid2q_raw

    string2ap = Vector{Tuple{RB.Signifier,RB.Signifier}}()
    for i = 1:n-1
        for j = 2:4
            push!(string2ap,(RB.Signifier(i,j),RB.Signifier(i+1,j)))
        end
        for j = 2:4
            push!(string2ap,(RB.Signifier(i,j),RB.Signifier(i+1,1)))
        end
    end
    cnt = RB.Connectivity(bodyid2q,string2ap)
    st = RB.Structure(rbs,apparatuses,cnt)
    RB.update_strings_apply_forces!(st)
	RB.jac_singularity_check(st)
    tr = RB.Robot(st,hub)
end

function build_Y(tr)
    @unpack st, hub = tr
	@unpack nstrings, strings = st
	actuators = hub.actuators
    nact = length(actuators)
    ret = spzeros(Int,nstrings,nact)
    for (i,actuator) in enumerate(actuators)
        ret[actuator.reg.id_string,i] = 1
    end
    ret
end


# function inverse2actuation(tgstruct,refstruct=deepcopy(tgstruct);gravity=false)
#     # Rx = RotY(π/18)
#     # Rx = RotY(0.0)
#     # reflinkn = links(n,h,Rx)
#     refq0,_,_ = RB.get_initial(refstruct)
#     refλ0,Δu,a = RB.inverse(tgstruct,refstruct,build_Y(tgstruct),gravity=false)
#
#     u0 = [s.original_restlen for s in tgstruct.strings]
#     rl = u0 + Δu
#     if any((x)->x<0,rl)
#         @error "Negative rest lengths"
#     end
#     ℓ = [s.state.length for s in tgstruct.strings]
#     if any((x)->x<0,ℓ-rl)
#         @error "Zero tension"
#     end
#     s = 1 ./ℓ
#     rl,a
# end
