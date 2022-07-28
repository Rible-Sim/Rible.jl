function make_new_tail(n)

    nver = n
    nhor = n + 1
    nb = nver + nhor
    ver_lengths = fill(0.2,nver)
    hor_lengths = fill(0.1,nhor)
    O = zeros(2, 2, nhor)
    for i = 2:nhor
        O[:, 1, i] = O[:, 1, i-1] + [0.0,-ver_lengths[i-1]]
    end
    for i = 1:nhor
        O[:, 2, i] = O[:, 1, i] + [hor_lengths[i], 0.0]
    end

    function rigidbody(i,O)
        lev,pos = divrem(i,2)
        if i == 1
            movable = false
            constrained = true
            ci = collect(1:6)
            Φi = Vector{Int}()
        else
            movable = true
            if i == 2
                constrained = true
                ci = collect(1:2)
                Φi = [1]
            else
                constrained = false
                ci = Vector{Int}()
                if isodd(pos)
                    Φi = collect(1:3)
                else
                    Φi = [1]
                end
            end
        end

        if isodd(pos)
            α = 0.0
            ri = SVector{2}(O[:,1, lev+1])
            rj = SVector{2}(O[:,2, lev+1])
            b = norm(rj-ri)
            r̄g = SVector{2}(0.0,0.0)
            r̄p1 = SVector{2}(-b,0.0)
            r̄p2 = SVector{2}( b,0.0)
            m = inertia = 0.1
            Ī = SMatrix{2,2}([
                inertia 0
                0 inertia
            ])
        else
            α = -π / 2
            ri = SVector{2}(O[:,1, lev])
            rj = SVector{2}(O[:,1, lev+1])
            b = norm(rj-ri)
            r̄g = SVector{2}(b/2,0.0)
            r̄p1 = SVector{2}(0.0,0.0)
            r̄p2 = SVector{2}(  b,0.0)
            m = inertia = 0.1
            Ī = SMatrix{2,2}([
                inertia 0
                0 0
            ])
        end
        r̄ps = [r̄p1,r̄p2]
        nr̄p = length(r̄ps)
        ṙo = zeros(2); ω = 0.0
        prop = TR.RigidBodyProperty(
            i,
            movable,
            m,
            Ī,
            SVector{2}(r̄g),
            r̄ps;
            constrained = constrained,
        )
        ro = ri
        if isodd(pos)
            lncs, _ = TR.NaturalCoordinates.NC2P1V(ri, rj, ro, α, ṙo, ω)
        else
            lncs, _ = TR.NaturalCoordinates.NC2D2P(ri, rj, ro, α, ṙo, ω)
        end
        state = TR.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)
        rb = TR.RigidBody(prop, state)
    end
    rbs = [
        rigidbody(i, O) for i = 1:nb
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigdibodies)
    matrix_sharing_raw = Vector{Matrix{Int}}()
    for i = 1:2n
        s = zeros(2, nb)
        if isodd(i)
            s[1:2, i] = 1:2
            s[1:2, i+1] = 1:2
        else
            s[1:2, i] = 3:4
            s[1:2, i+1] = 1:2
        end
        push!(matrix_sharing_raw, s)
    end
    matrix_sharing = reduce(vcat, matrix_sharing_raw)
    display(matrix_sharing)
    indexedcoords = TR.index(rigdibodies, matrix_sharing)

    ncables = 4n
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for i = 1:ncables
        lev,j = divrem(i-1,4)
        original_restlens[i] = ifelse(j ∈ [1, 0], ver_lengths[lev+1], ver_lengths[lev+1])
        ks[i] = ifelse(j ∈ [1, 0], 100, 100.0)
    end
    cables =
        [TR.Cable2D(i, 0.5original_restlens[i], ks[i], 0.0) for i = 1:ncables]  #
    tensiles = (cables = cables,)
    acs = [
        TR.ManualActuator(1,
            [4(i-1)+j for i = 1:n for j = 1:4],
            5original_restlens[[4(i-1)+j for i = 1:n for j = 1:4]],
        )
    ]
    hub = (actuators = acs,)

    matrix_cnt_raw = Vector{Matrix{Int}}()
    for i = 1:n
        s = zeros(4, nb)
        s[1, 2i-1] = 1
        s[1, 2i+1] = -1
        s[2, 2i] = 1
        s[2, 2i+1] = -1
        s[3, 2i] = 1
        s[3, 2i+1] = -2
        s[4, 2i-1] = 2
        s[4, 2i+1] = -2
        push!(matrix_cnt_raw, s)
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)
    display(matrix_cnt)
    connections = TR.connect(rigdibodies, matrix_cnt)

    cnt = TR.Connectivity(numberedpoints, indexedcoords, connections)

    tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt)
    bot = TR.TensegrityRobot(tg, hub)
end
