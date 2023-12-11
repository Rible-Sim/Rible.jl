
function man_ndof(num_of_dof,onedir=[0.0,-1.0];θ=0.0,k=0.0,c=0.0,unit="mks",restlen=0.16,isvirtual=true)
    nbodies = num_of_dof + 1
    nbp = nbodies + 1
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    ends_idx = [1,nbodies]
    inner_idx = 2:nbodies-1
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
    for j in ends_idx
        a[j] = 12252e-5unit_L
        if isvirtual
            m[j] = 94.873e-3unit_M
        else
            m[j] = 198.553e-3unit_M
        end
        Ia[j] = 3068.8e-7unit_I  #??2e-7??
    end
    for k in inner_idx
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

        contactable = true
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

        mass_locus = SVector{2}([CoM_x,CoM_y])
        ap5 = SVector{2}([ap5_x,ap5_y])
        ap6 = SVector{2}([ap6_x,ap6_y])

        ap7 = SVector{2}([ap7_x,ap7_y])
        ap8 = SVector{2}([ap8_x,ap8_y])
        aps = [ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8]

        # only to get nmcs
        prop = RB.RigidBodyProperty(
                    i,contactable,m,
                    Ī,
                    mass_locus,
                    aps;
                    visible=ifelse(i==1,true,false)
                )
        ro = copy(ri)
        ω = 0.0
        ṙo = @SVector zeros(2)
        α = get_angle([1.0,0.0],rj-ri)
        nmcs = RB.NCF.NC1P2V(ri,ro,α)
        if i == 1
            ci = collect(1:6)
            cstr_idx = Int[]
        else
            ci = Int[]
            cstr_idx = collect(1:3)
        end
        state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,cstr_idx)

        body = RB.RigidBody(prop,state)
    end
    rbs = [
        rigidbody(i,m[i],a[i],
        SMatrix{2,2}(Diagonal([Ia[i]/2,Ia[i]]/2)),
        SVector{2}(B[:,i]),SVector{2}(A[:,i+1])) for i = 1:nbodies
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
    indexed = RB.index(rigdibodies,)

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
    ss = [RB.DistanceSpringDamper2D(i, original_restlens[i],ks[i],cs[i]) for i = 1:nstrings]
    force_elements = (cables=ss,)

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

    pinjoints = [
        RB.PinJoint(
            RB.Hen2Egg(
                i-1,RB.ID(rbs[i-1],2),RB.ID(rbs[i],1)
            )
        )
        for i in 2:nbodies
    ]

    jointed = RB.join(pinjoints,indexed)

    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,force_elements,cnt)
    bot = RB.Robot(st,hub)
end


function man_ndof_2022(num_of_dof,onedir=[1.0,0.0];θ=0.0,k=1250.0,c=0.0,unit="mks")

    edi = 2 # 1 为HUI版，2 为ZHAO版

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
        if j==nbodies
            a[j] = 12250e-5unit_L#12252
            m[j] = 99.55e-3unit_M#99.26,149.263,198.55
            Ia[j] =  2977e-7unit_I#2977,3029,3136
        else
            a[j] = 12250e-5unit_L
            #m[j] = 835.90254985e-3unit_M
            # m[j] = 200.40071778e-3unit_M
            m[j] = 273.11e-3unit_M
            # Ia[j] = Ic_lower[i] + m[j]*1/3*a[j]^2
            #Ia[j] = 28130.53053840*2e-7unit_I
            Ia[j] =  9907e-7unit_I
            # Ia[j] = 3738.8e-7unit_I  #??2e-7??
        end
        # @show a[j],m[j],Ia[j]
    end
    for (i,k) in enumerate(upper_index)
        if k==nbodies
            a[k] = 12250e-5unit_L
            m[k] = 99.263e-3unit_M
            Ia[k] = 2977e-7unit_I
        # elseif k==nbodies-1
        #     a[k] = 12252e-5unit_L
        #     m[k] = 211.285e-3unit_M
        #     Ia[k] = 400700.4159e-7unit_I
        else
            # a[k] = 16.0e-2unit_L
            a[k] = 12250e-5unit_L
            #m[k] = 666.25659673e-3unit_M
            # m[k] = 200.40071778e-3unit_M
            m[k] = 273.11e-3unit_M
            # Ia[k] = Ic_upper[i] + m[k]*1/3*a[k]^2
            Ia[k] =  9907e-7unit_I#8270.5
            # Ia[k] = 3738.8e-7unit_I
            # @show a[k],m[k],Ia[k]
        end
    end
    R = [
        cos(θ) -sin(θ);
        sin(θ) cos(θ)
    ]
    A = zeros(2,nbp)
    # A[:,2] .= A[:,1] .+ a[1]*[1.0,0.0]#HUI版
    # A[:,2] .= A[:,1] .+ √3/2*a[1]*onedir#zhao版
    # edi == 1 ? A[:,2] .= A[:,1] .+ a[1]*[1.0,0.0] : A[:,2] .= A[:,1] .+ √3/2*a[1]*onedir
    edi == 1 ? A[:,2] .= A[:,1] .+ a[1]*[0.0,-1.0] : A[:,2] .= A[:,1] .+ √3/2*a[1]*onedir#竖直
    A[:,3] .= A[:,2] .+ √3*a[2]*R*onedir
    θ=0.0
    R = [
        cos(θ) -sin(θ);
        sin(θ) cos(θ)
    ]
    A[:,4] .= A[:,3] .+ √3*a[3]*R*onedir

    # for i in 3:nbp
    #     # A[:,i] .= A[:,i-1] .+ a[i-1]*[sin(θ*(i-2)),-cos(θ*(i-2))]
    #     edi == 1 ? A[:,i] .= A[:,i-1] .+ a[i-1]*[sin(θ*(i-2)),-cos(θ*(i-2))] : A[:,i] .= A[:,i-1] .+ √3*a[i-1]*R^(i-2)*onedir#原*R^(i-2)*onedir
    #
    # end

    function rigidbody(i,m,a,Ia,ri,rj)
        contactable = true
        if i == 1
            visible = true
            # edi == 1 ? ci = [1,3,4] : ci = [1,2,4]
            ci = collect(1:6)
            cstr_idx = Int[]
        elseif i == 2
            visible = true
            ci = collect(1:2)
            cstr_idx = collect(1:3)
        else
            visible = true
            ci = Int[]
            cstr_idx = collect(1:3)
        end

        sti_l = 0.01832 #弹簧偏移顶点的量
        rao_l = 0.00475 #绕线滑轮
        #---HUI版
        if edi == 1
            if i == nbodies
                CoM_x = 0.06383
                CoM_y = 0.03174
            else
                CoM_x = 0.06383
                CoM_y = 0.03174
            end

            if isodd(i)
                CoM_x =  0.06389
                CoM_y = -0.03188
            end
            # nap = 3 #?
            ap1 = SVector{2}([0.0,0.0])
            ap2 = SVector{2}([a,0.0])
            ap3 = SVector{2}([a/2,0.0])
            ap4_x = a/2
            ap4_y = √3/2*a
            ap5_x = a-sti_l*cos(pi/6)
            ap5_y = sti_l*sin(pi/6)
            ap6_x = a/2
            ap6_y = √3/2*a-sti_l
            ap7_x = -rao_l*sin(pi/6)
            ap7_y = -rao_l*cos(pi/6)
            ap8_x = a/2+rao_l*sin(pi/6)
            ap8_y = √3/2*a+rao_l*cos(pi/6)

            if isodd(i)
                ap4_y = -ap4_y
                ap5_y = -ap5_y
                ap6_y = -ap6_y
                ap7_y = -ap7_y
                ap8_y = -ap8_y

            end
        # zhao版
        else
            ap1 = SVector{2}([0.0,0.0])
            if i == 1
                CoM_x = 0.06049371*sin(pi/3)-0.0242439*cos(pi/3)
                CoM_y = a/2-0.06049371*cos(pi/3)-0.0242439*sin(pi/3)
                # CoM_x = 0.06535974*sin(pi/3)-0.03465684*cos(pi/3)
                # CoM_y = a/2-0.06535974*cos(pi/3)-0.03465684*sin(pi/3)
                ap2 = SVector{2}([√3/2*a,0.0])
                ap3_x =  0.0
                ap3_y =  a/2
                ap4_x =  0.0
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
                # CoM_x = 0.05797798*sin(pi/3) + 0.03024948*cos(pi/3)
                # CoM_y = 0.05797798*cos(pi/3) - 0.03024948*sin(pi/3)
                CoM_x = 0.06522388#0.07629847，50g0.063166
                CoM_y = 0.0025259#0.00336616,50g0.0051033
                ap2 = SVector{2}([√3/2*a,0.0])
                ap3_x =  √3/2*a
                ap3_y =  a/2
                ap4_x =  √3/2*a
                ap4_y = -a/2
                ap5_x = √3/2*a+rao_l*sin(pi/6)
                ap5_y = a/2+rao_l*cos(pi/6)
                ap6_x = √3/2*a+rao_l*sin(pi/6)
                ap6_y = -a/2-rao_l*cos(pi/6)
                ap7_x = √3/2*a-sti_l*sin(pi/6)
                ap7_y = a/2-sti_l*cos(pi/6)
                ap8_x = √3/2*a-sti_l*sin(pi/6)
                ap8_y = -a/2+sti_l*cos(pi/6)

            else
                CoM_x = 0.1015888
                CoM_y =  0.0003247
                ap2 = SVector{2}([√3*a,0.0])
                ap3_x =  √3/2*a
                ap3_y =  a/2
                ap4_x =  √3/2*a
                ap4_y = -a/2
                ap5_x = √3/2*a
                ap5_y = a/2+rao_l
                ap6_x = √3/2*a
                ap6_y = -a/2-rao_l
                ap7_x = √3/2*a
                ap7_y = a/2-sti_l
                ap8_x = √3/2*a
                ap8_y = -a/2+sti_l
            end
            ap3 = SVector{2}([ap3_x,ap3_y])
        end
        ##--------------------

        mass_locus  = SVector{2}([CoM_x,CoM_y])
        ap4 = SVector{2}([ap4_x,ap4_y])
        ap5 = SVector{2}([ap5_x,ap5_y])
        ap6 = SVector{2}([ap6_x,ap6_y])
        ap7 = SVector{2}([ap7_x,ap7_y])
        ap8 = SVector{2}([ap8_x,ap8_y])
        aps = [ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8]

        prop = RB.RigidBodyProperty(i,contactable,m,Ia,
                    mass_locus,aps;visible=visible
                    )

        α = get_angle([1.0,0.0],rj-ri)
        ω = 0.0
        ro = ri
        ṙo = zero(ro)

        nmcs = RB.NCF.NC2P1V(ri,rj,ro,α)

        state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,cstr_idx)

        body = RB.RigidBody(prop,state)

    end
    rbs = [rigidbody(i,m[i],a[i],
            SMatrix{2,2}(Diagonal([Ia[i]*0.01,Ia[i]]*0.99)),
            A[:,i],A[:,i+1]) for i = 1:nbodies]

    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)

    matrix_sharing = zeros(Int,2nbodies-2,nbodies)
    for i = 1:nbodies-1
        matrix_sharing[2i-1,i:i+1] = [3, 1]
        matrix_sharing[2i  ,i:i+1] = [4, 2]
    end
    # display(matrix_sharing)
    indexed = RB.index(rigdibodies, matrix_sharing)

    nstrings = 2(nbodies-1)

    upstringlen = 0.178-14.5/1275
    lostringlen = upstringlen
    #此处定义每个关节的初始静止长度
    # stringlen = [0.1733,0.1714,0.1738,0.1725,0.1698,0.171,0.1698,0.171].-14.5/k
    stringlen = [0.1918,0.1875,0.18,0.18,0.18,0.18].-14.5/k
    original_restlens = zeros(nstrings)
    restlens = zeros(nstrings)
    actuallengths = zeros(nstrings)
    ks = zeros(nstrings)
    cs = zeros(nstrings)
    kz = [1250,1250,310,310]
    cs .= [3,3,1.,1.]#c
    if edi == 1
            #这里定义每关节的绳长不同
        snum = 1
        for i = 1:nstrings
            # j = i % 4
            original_restlens[i] = stringlen[snum]#ifelse(j∈[1,0],stringlen[snum],stringlen[snum])
            ks[i] = k #ifelse(j∈[1,0],k,k)
            if i%2 == 0
                snum = snum + 1
            end
        end
        # for i = 1:nstrings
        #     j = i % 4
        #     original_restlens[i] = ifelse(j∈[1,0],upstringlen,lostringlen)
        #     ks[i] = ifelse(j∈[1,0],k,k)
        # end
    else
        # stringlen = [0.1898-14.5/k,0.1769-14.5/k]
        snum = 1
        for i = 1:nstrings
            # j = i % 4
            original_restlens[i] = stringlen[snum]#ifelse(j∈[1,0],stringlen[snum],stringlen[snum])
            ks[i] = kz[i] #ifelse(j∈[1,0],k,k)
            if i%2 == 0
                snum = snum + 1
            end
        end
    end

    ss = [RB.DistanceSpringDamper2D(i, original_restlens[i],ks[i],cs[i]) for i = 1:nstrings]
    force_elements = (cables=ss,)

    matrix_cnt = zeros(Int,2(nbodies-1),nbodies)

    for i = 1:length(rbs)-1
        if edi == 1
            # push!(string2ap,(RB.ID(i,7),RB.ID(i+1,6)))
            # push!(string2ap,(RB.ID(i,8),RB.ID(i+1,5)))
            matrix_cnt[2(i-1)+1,i:i+1] = [7,-6]
            matrix_cnt[2(i-1)+2,i:i+1] = [8,-5]
        else
            # push!(string2ap,(RB.ID(i,5),RB.ID(i+1,7)))
            # push!(string2ap,(RB.ID(i,6),RB.ID(i+1,8)))
            matrix_cnt[2(i-1)+1,i:i+1] = [5,-7]
            matrix_cnt[2(i-1)+2,i:i+1] = [6,-8]
        end

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
    st = RB.Structure(rigdibodies,force_elements,cnt)
    bot = RB.Robot(st,hub)
end


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
