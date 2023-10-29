function spine3d(n;c=0.0)
    nbodies = n
    nbp = 4*n

    a = 0.04 #m
    h = 0.04 #m
    θ = π/4

    mass = 0.1 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-2
    mass_locus = [0.0, 0.0, 0.0]

    ap_x = cos(θ)*a
    ap_y = sin(θ)*a
    ap4 = SVector{3}([  0.0,  ap_x,  ap_y])
    ap3 = SVector{3}([  0.0, -ap_y,  ap_y])
    ap2 = SVector{3}([ ap_x,   0.0, -ap_y])
    ap1 = SVector{3}([-ap_x,   0.0, -ap_y])
    aps = [ap1,ap2,ap3,ap4]

    movable = ones(Bool,n)
    movable[1] = false
    constrained = zeros(Bool,n)
    constrained[1] = true

    props = [RB.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(mass_locus),
                aps;constrained=constrained[i]) for i = 1:n]

    rs = [[0.0,0.0,i*h] for i = 0:n-1]
    Rs = [RotX(0.0) for i = 1:n]
    ṙs = [zeros(3) for i = 1:n]
    ωs = [zeros(3) for i = 1:n]

    function rigidbody(i,prop,aps,r,R,ṙ,ω)
        if i == 1
            ci = collect(1:12)
            constraints_indices = Int[]
        else
            ci = Int[]
			constraints_indices = collect(1:6)
        end

        # ri,rj,rk,rl = [r+R*ap for ap in aps]
        # nmcs = RB.NCF.NC4P(ri,rj,rk,rl,r,R,ṙ,ω)
        # nmcs = RB.NCF.NC3P1V(ri,rk,rl,r,R,ṙ,ω)
        # nmcs = RB.NCF.NC2P2V(rk,rl,r,R,ṙ,ω)
        nmcs = RB.NCF.NC1P3V(r,r,R,ṙ,ω)
        state = RB.RigidBodyState(prop,nmcs,r,R,ṙ,ω,ci,constraints_indices)
        body = RB.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,props[i],aps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]
	rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = RB.number(rigdibodies)
	indexedcoords = RB.index(rigdibodies)

    ncables = 8*(n-1)
    stringlenH = 0.6h
    stringlenR = 0.04
    stringlens = repeat(vcat(fill(stringlenH,4),fill(stringlenR,4)),n-1)
    kH = 1e1
    kR = 4e1
    ks = repeat(vcat(fill(kH,4),fill(kR,4)),n-1)
    # c = 0.0
    cs = repeat(fill(c,8),n-1)
    cables = [RB.Cable3D(i,stringlens[i],ks[i],cs[i]) for i = 1:ncables]
    tensiles = (cables=cables,)
    acs = [RB.ManualActuator(RB.SimpleRegistor(8(i-1)+j,stringlens[8(i-1)+j])) for i = 1:n-1  for j = 1:6]
    hub = (actuators=acs,)


    string2ap = Vector{Tuple{RB.ID,RB.ID}}()
	asslist = [(j=j,k=k) for j = 3:4 for k = 1:2]
    string2ap = reduce(vcat,
        [
        vcat(
            [
            (id=8(i-1)+j,hen=RB.ID(rbs[i],j),egg=RB.ID(rbs[i+1],j))
            for j = 1:4
            ],
            [
            (id=8(i-1)+4+a,hen=RB.ID(rbs[i],jk.j),egg=RB.ID(rbs[i+1],jk.k))
            for (a,jk) in enumerate(asslist)
            ]
        )
        for i = 1:n-1
        ]
    )
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
        push!(matrix_cnt_raw, s)
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)
	# display(matrix_cnt)
    connections = RB.connect(rigdibodies, matrix_cnt)
    cnt = RB.Connectivity(numberedpoints, indexedcoords, connections)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    bot = RB.Robot(st,hub)
end
