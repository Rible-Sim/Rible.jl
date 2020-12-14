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
    CoM = [0.0, 0.0, 7.5932384]*1.e-2


    ap1 = SVector{3}([0.0, 0.0, 0.0])
    ap2 = SVector{3}([a,   0.0,   h])
    ap3 = SVector{3}([a*cos(θ),  a*sin(θ), h])
    ap4 = SVector{3}([a*cos(θ), -a*sin(θ), h])
    aps = [ap1,ap2,ap3,ap4]

    movable = ones(Bool,n)
	# movable[1] = false
    constrained = zeros(Bool,n)
    constrained[1] = true

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(CoM),
                aps;constrained=constrained[i]) for i = 1:n]

    rs,Rs = bend(n,di,rm)
    ṙs = [zeros(3) for i = 1:n]
    ωs = [zeros(3) for i = 1:n]

    function rigidbody(i,prop,aps,r,R,ṙ,ω)
        if i == 1
            constrained_index = collect(2:2:12)
        else
            constrained_index = Vector{Int}()
        end

        ri,rj,rk,rl = [r+R*ap for ap in aps]
        bps,q,q̇ = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r,R,ṙ,ω)
        state = TR.RigidBodyState(prop,bps,r,R,ṙ,ω,q,q̇,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,props[i],aps,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]

    nstrings = 6*(n-1)
    stringlenH = 0.6h
    stringlenR = 6e-2 #m
    stringlens = repeat(vcat(fill(stringlenH,3),fill(stringlenR,3)),n-1)
    kH = k #N/m
    kR = k #N/m
    ks = repeat(vcat(fill(kH,3),fill(kR,3)),n-1)
    # c = 0.0
    cs = repeat(fill(c,6),n-1)
    ss = [TR.SString3D(stringlens[i],ks[i],cs[i]) for i = 1:nstrings]
    acs = [TR.Actuator(SVector{1}(ss[6(i-1)+j])) for i = 1:n-1  for j = 1:6]
    bodynq = TR.get_nbodycoords(rbs[1])
    body2q_raw = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]
	# body2q = TR.filter_body2q(body2q_raw,rbs)
	body2q = body2q_raw

    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n-1
        for j = 2:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j)))
        end
        for j = 2:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,1)))
        end
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end

function build_Y(tgstruct)
    @unpack nstrings,strings,actuators = tgstruct
    nact = length(actuators)
    ret = spzeros(Int,nstrings,nact)
    for (i,act) in enumerate(actuators)
        s1 = act.strings[1]
        is1 = findfirst((x)->x==s1,strings)
        ret[is1,i] = 1
    end
    ret
end


# function inverse2actuation(tgstruct,refstruct=deepcopy(tgstruct);gravity=false)
#     # Rx = RotY(π/18)
#     # Rx = RotY(0.0)
#     # reflinkn = links(n,h,Rx)
#     refq0,_,_ = TR.get_initial(refstruct)
#     refλ0,Δu,a = TR.inverse(tgstruct,refstruct,build_Y(tgstruct),gravity=false)
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
