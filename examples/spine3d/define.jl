function spine3d(n)
    nbodies = n
    nbp = 4*n

    a = 0.04 #m
    h = 0.04 #m
    θ = π/4

    mass = 0.1 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-2
    CoM = [0.0, 0.0, 0.0]

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

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(CoM),
                aps;constrained=constrained[i]) for i = 1:n]

    rs = [[0.0,-i*h,0.0] for i = 0:n-1]
    Rs = [RotX(π/2) for i = 1:n]
    ṙs = [zeros(3) for i = 1:n]
    ωs = [zeros(3) for i = 1:n]

    function rigidbody(i,prop,aps,r,R,ṙ,ω)
        if i == 1
            constrained_index = collect(1:2:11)
            constrained_index = Vector{Int}()
        else
            constrained_index = Vector{Int}()
        end

        ri,rj,rk,rl = [r+R*ap for ap in aps]
        bps,q,q̇ = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r,R,ṙ,ω)
        state = TR.RigidBodyState(prop,bps,r,R,ṙ,ω,q,q̇,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,props[i],aps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]

    nstrings = 8*(n-1)
    stringlenH = 0.6h
    stringlenR = 0.04
    stringlens = repeat(vcat(fill(stringlenH,4),fill(stringlenR,4)),n-1)
    kH = 1e1
    kR = 4e1
    ks = repeat(vcat(fill(kH,4),fill(kR,4)),n-1)
    c = 0.0
    cs = repeat(fill(c,8),n-1)
    ss = [TR.SString3D(stringlens[i],ks[i],cs[i]) for i = 1:nstrings]
    acs = [TR.Actuator(SVector{1}(ss[8(i-1)+j])) for i = 1:n-1  for j = 1:6]
    bodynq = TR.get_nbodycoords(rbs[1])
    body2q_raw = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]
    body2q = TR.filter_body2q(body2q_raw,rbs)
    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n-1
        for j = 1:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j)))
        end
        for j = 3:4
            for k = 1:2
                push!(string2ap,(TR.ID(i,j),TR.ID(i+1,k)))
            end
        end
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
