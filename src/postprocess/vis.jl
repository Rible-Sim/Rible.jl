MakieCore.@recipe(Viz, structure) do scene
    # theme_pub,
    MakieCore.Attributes(        
        isref=false,
        showlabels=false,
        show_cable_labels = false,
        show_mass_center_labels = false,
        show_node_labels = false,
        showpoints=false,
        show_mass_centers = false,
        show_nodes = false,
        showarrows=false,
        showmesh=true,
        showwire=false,
        showcables=true,
        pointcolor=:black,
        slack_linestyle = :dash,
        cablecolor=:deepskyblue,
        cablelabelcolor=:darkgreen,
        rigidlabelcolor=:darkblue,
        refcolor=:lightgrey,
        cablewidth=0.6,
        meshcolor=:slategrey,
        fontsize = 12,
    )
end

function MakieCore.plot!(viz::Viz{Tuple{S}};
    ) where S <:AbstractBody
    body_ob = viz[:structure]
    # body decorations
    mass_center_ob = lift(body_ob) do body_ob
        body_ob.state.mass_locus_state.frame.position |> GB.Point
    end
    nodes_ob = lift(body_ob) do body_ob
        [
            locus_state.frame.position |> GB.Point
            for locus_state in body_ob.state.loci_states
        ]
    end

    # arrows_ob = lift(body_ob) do body_ob
    #     body_ob.state.as .|> GB.Point
    # end
    id  = lift(body_ob) do body_ob
        body_ob.prop.id
    end
    if viz.show_mass_centers[]
        MakieCore.scatter!(viz,mass_center_ob;)
        if viz.show_mass_center_labels[]
            MakieCore.text!(viz,
                "r$(id[])g" ,
                position = mass_center_ob,
                color = viz.rigidlabelcolor[],
                align = (:left, :top),
                offset = (-5, -10)
            )
        end
    end
    if viz.show_nodes[]
        MakieCore.scatter!(viz,nodes_ob; color=viz.pointcolor[])
        if viz.show_node_labels[]
            MakieCore.text!(viz,
                ["r$(id[])p$pid" for (pid,rp) in enumerate(nodes_ob[])],
                position = nodes_ob,
                color = :darkred,
                align = (:left, :top),
                offset = (0,2viz.fontsize[]*(rand()-0.5))
            )
        end
    end
    # body mesh
    if viz.showmesh[] || viz.showwire[]
        meshes_ob = lift(body_ob) do body_ob
            build_mesh(body_ob,color=viz.meshcolor[])
        end
        if viz.showwire[]
            strokewidth = linewidth
        else
            strokewidth = 0
        end
        if viz.isref[]
            MakieCore.mesh!(viz, meshes_ob; shading = MakieCore.Automatic())
        else
            MakieCore.poly!(viz, meshes_ob; shading = MakieCore.Automatic(),
                # strokewidth
            )
        end
    end
    viz
end

function MakieCore.plot!(viz::Viz{Tuple{S}};
    ) where S <: Apparatus{<:CableJoint}
    cable_appar_ob = viz[:structure]

    seg_ob = lift(cable_appar_ob) do cab
        point_start = cab.force.state.start |> GB.Point
        point_stop = cab.force.state.stop |> GB.Point
        [(point_start,point_stop)]
    end

    # slackonly=false,
    # noslackonly=true
    linestyle_ob = lift(cable_appar_ob) do cab
        if cab.force.state.length > cab.force.state.restlen
            return :solid
        else
            return viz.slack_linestyle[]
        end
    end

    MakieCore.linesegments!(
        viz, seg_ob, 
        color = viz.cablecolor[], 
        linewidth = viz.cablewidth[], 
        linestyle = linestyle_ob
    )

    point_mid_ob = lift(cable_appar_ob) do cab
        point_start = cab.force.state.start
        point_stop = cab.force.state.stop 
        (point_start .+ point_stop)./2 |> GB.Point
    end


    id_ob = lift(cable_appar_ob) do cab
        "c$(cab.id)"
    end
  

    # show cable labels
    if viz.show_cable_labels[]
        MakieCore.text!(viz,
            id_ob,
            position = point_mid_ob,
            color = viz.cablelabelcolor[],
            align = (:left, :top),
            offset = (-5, -10)
        )
    end
    viz
end

function MakieCore.plot!(viz::Viz{Tuple{S}};
    ) where S <:AbstractStructure
    if viz.isref[]
        showlabels = false
        cablecolor=
        cablelabelcolor=viz.refcolor[]
        meshcolor=viz.refcolor[]
        viz.showcables[] = false
    else
        meshcolor=viz.meshcolor[]
        cablecolor = viz.cablecolor[]
        cablelabelcolor = viz.cablelabelcolor[]
    end
    if viz.showlabels[]
        show_cable_labels = 
        show_mass_center_labels = 
        show_node_labels = true
    else
        show_cable_labels = viz.show_cable_labels[]
        show_mass_center_labels = viz.show_mass_center_labels[]
        show_node_labels = viz.show_node_labels[]
    end
    if viz.showpoints[]
        show_mass_centers = show_nodes = true
    else
        show_mass_centers = viz.show_mass_centers[]
        show_nodes = viz.show_nodes[]
    end
    tgob = viz[:structure]
    (;connectivity) = tgob[]
    (;num_of_bodies,) = connectivity.indexed
    cable_apparatuses_ob = lift(tgob) do tgob
        filter(sort(tgob.apparatuses)) do appar
            appar.joint isa CableJoint
        end
    end

    ncables = length(cable_apparatuses_ob[])

    cables_array_ob = [
        begin
            lift(cable_apparatuses_ob) do cable_apparatuses_ob
                cable_apparatuses_ob[i]
            end
        end
        for i = 1:ncables
    ]

    if ncables > 0 && viz.showcables[]
        for cable_ob in cables_array_ob
            viz!(viz,cable_ob;
                cablecolor,
                cablewidth = viz.cablewidth[],
                cablelabelcolor = viz.cablelabelcolor[],
                slack_linestyle = viz.slack_linestyle[],
                show_cable_labels,
            )
        end
    end
    if num_of_bodies > 0
        bodies_ob = lift(tgob) do tgob
            get_bodies(tgob)
        end
        bodies_array_ob = [
            begin
                lift(bodies_ob) do bodies_ob
                    bodies_ob[i]
                end
            end
            for i = 1:num_of_bodies
        ]
        for body_ob in bodies_array_ob
            viz!(viz,body_ob;
                show_mass_centers,
                show_mass_center_labels,
                show_nodes,
                show_node_labels,
                meshcolor,
                pointcolor = viz.pointcolor[],
                showmesh = viz.showmesh[],
                # showwire = viz.showwire[],
            )
        end
    end
    viz
end

function build_mesh(body::AbstractRigidBody;update=true,color=nothing)
    (;mesh) = body
    @assert !(mesh isa Nothing)
    if update
        frame = to_3D(body.state.origin_frame)
        origin_position = frame.position
        R = frame.axes.X
    else
        origin_position = SVector(0,0,0)
        R = Matrix(1I,3,3)
    end
    trans = Translation(origin_position)
    rot = LinearMap(R)
    ct = trans ∘ rot
    updated_pos = GB.Point3f.(ct.(mesh.position))
    fac = GB.faces(mesh)
    nls = GB.normals(updated_pos,fac)
    
    
    if !(color isa Nothing)
        parsedcolor = parse(CT.RGB{Float32},color)
        colors = fill(parsedcolor,length(updated_pos))
    elseif hasproperty(mesh,:color)
        colors = mesh.color
    else
        parsedcolor = parse(CT.RGB{Float32},:slategrey) 
        colors = fill(parsedcolor,length(updated_pos))
    end
    GB.Mesh(GB.meta(updated_pos,normals=nls,color=colors),fac)
    # GB.Mesh(GB.meta(coloredpoints,normals=nls),fac)
end

function make_patch(;trans=[0.0,0,0],rot=RotX(0.0),scale=1,color=:slategrey)
    parsedcolor = parse(CT.RGBA{Float32},color)
    function patch(mesh)
        ct = Translation(trans) ∘ LinearMap(rot)
        updated_pos = ct.(mesh.position.*scale)
        fac = GB.faces(mesh)
        nls = GB.normals(updated_pos,fac)
        colors = fill(parsedcolor,length(updated_pos))
        GB.Mesh(GB.meta(updated_pos,normals=nls,color=colors),fac)
    end
end

function get_groundmesh(f::Function,rect)
    GeometryBasics.Mesh(f, rect, NaiveSurfaceNets()) |> make_patch(;color = :snow)
end

function get_groundmesh(plane::Plane,rect)
    GB.Mesh(rect, Meshing.MarchingCubes()) do v
        signed_distance(v,plane)
    end  |> make_patch(;color = :snow)
end

function get_groundmesh(static_env::StaticContactSurfaces,rect)
    map(static_env.surfaces) do surface
        GB.Mesh(rect, Meshing.MarchingCubes()) do v
            signed_distance(v,surface)
        end 
    end |> GB.merge |> make_patch(;color = :snow)
end

function hidex(ax)
    ax.xticklabelsvisible = false
    ax.xlabelvisible = false
end

function hidey(ax)
    ax.yticklabelsvisible = false
    ax.ylabelvisible = false
end

function hidez(ax)
    ax.zticklabelsvisible = false
    ax.zlabelvisible = false
end

function hidexyz(ax)
    hidex(ax)
    hidey(ax)
    hidez(ax)
end

function time2step(at,t)
    atstep = findfirst((x)->x>=at, t)
    @assert !isa(atstep,Nothing)
    atstep
end

function get_groundmesh(::Nothing,rect)
    plane_n = [0,0,1.0]
    plane_r = zeros(3)
    plane = Plane(plane_n,plane_r)
    get_groundmesh(plane,rect)
end

function get_linesegs_cables(structure;slackonly=false,noslackonly=false)
    (;connected) = structure.connectivity.tensioned
    (;cables) = structure.apparatuses
    ndim = get_num_of_dims(structure)
    T = get_numbertype(structure)
    linesegs_cables = Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
    foreach(connected) do scnt
        scable = cables[scnt.id]
        ret = (Point(scnt.hen.body.state.loci_states[scnt.hen.pid].position),
                Point(scnt.egg.body.state.loci_states[scnt.egg.pid].position))
        slacking = scable.state.tension <= 0
        if (slackonly && slacking) ||
           (noslackonly && !slacking) ||
           (!slackonly && !noslackonly)
            push!(linesegs_cables,ret)
        end
    end
    linesegs_cables
end

function endpoints2mesh(
        p1,p2;
        radius=norm(p2-p1)/40,
        n1=10,n2=2,
        color=:slategrey
    )
    cyl_bar = Meshes.Cylinder(
        Meshes.Point(p1),
        Meshes.Point(p2),
        radius
    )
    cylsurf_bar = Meshes.boundary(cyl_bar)
    # Meshes.sample(cylsurf_bar,Meshes.RegularSampling(10,3))
    cyl_bar_simple = Meshes.discretize(cylsurf_bar,Meshes.RegularDiscretization(n1,n2))
    simple2mesh(cyl_bar_simple,color)
end

function spbasis(n)
    a = abs.(n)
    if (a[1]≥0 && a[2]≥0) || (a[1]≤0 && a[2]≤0) 
        v = SVector(n[1]+1,n[2]-1,n[3])
    else
        v = SVector(n[1]-1,n[2]-1,n[3])
    end
    t = n×v |> normalize
    b = n×t |> normalize
    t,b
end

function build_mesh(fb::FlexibleBody,nsegs=100;color=:slategrey)
    (;state,coords,cache) = fb
    (;e) = cache
    ancs = coords.nmcs
    (;L,radius) = ancs
    T = typeof(L)
    V = T <: AbstractFloat ? T : Float64
    _r = ANCF.make_r(ancs,e)
    _rₓ = ANCF.make_rₓ(ancs,e)

    sz = (10,nsegs)
    φmin, φmax = V(0), V(2π)
    xmin, xmax = V(0), V(L)
    δφ = (φmax - φmin) / sz[1]
    φrange = range(φmin, stop=φmax-δφ, length=sz[1])
    xrange = range(xmin, stop=xmax,    length=sz[2])

    function point(φ, x)
        o = _r(x) |> Meshes.Point
        n = _rₓ(x) |> normalize
        u,v = spbasis(n)
        R = [u v n]
        R*Meshes.Vec(radius*cos(φ), radius*sin(φ), 0.0) + o
    end

    points = IterTools.ivec(point(φ, x) for φ in φrange, x in xrange) |> collect
     # connect regular samples with quadrangles
    nx, ny = sz
    topo   = Meshes.GridTopology((nx-1, ny-1))
    middle = collect(Meshes.elements(topo))
    for j in 1:ny-1
        u = (j  )*nx
        v = (j-1)*nx + 1
        w = (j  )*nx + 1
        z = (j+1)*nx
        quad = Meshes.connect((u, v, w, z))
        push!(middle, quad)
    end

    connec = middle
    Meshes.SimpleMesh(points, connec) |> (x)->simple2mesh(x,color)
end

function simple2mesh(sp,color=:slategrey)
    dim   = Meshes.embeddim(sp)
    nvert = Meshes.nvertices(sp)
    nelem = Meshes.nelements(sp)
    verts = Meshes.vertices(sp)
    topo  = Meshes.topology(sp)
    elems = Meshes.elements(topo)

    # coords of vertices
    coords = Meshes.coordinates.(verts)
    # fan triangulation (assume convexity)
    tris4elem = map(elems) do elem
      I = Meshes.indices(elem)
      [[I[1], I[i], I[i+1]] for i in 2:length(I)-1]
    end

    # flatten vector of triangles
    tris = [tri for tris in tris4elem for tri in tris]
    points  = GB.Point.(coords)
    faces  = GB.TriangleFace{UInt64}.(tris)
    nls = GB.normals(points,faces)
    parsedcolor = parse(CT.RGB{Float32},color)
    colors = fill(parsedcolor,length(points))
    GB.Mesh(GB.meta(points,normals=nls,color=colors),faces)
end

function get_err_avg(bots;bid=1,pid=1,di=1,field=:traj)
    nbots = length(bots)
    lastbot = bots[nbots]
    if field == :traj
        get! = get_trajectory!
    elseif field == :vel
        get! = get_velocity!
    elseif field == :midvel
        get! = get_mid_velocity!
    elseif field == :ang
        get!(bot,bid,pid) = get_orientation!(bot,bid)
    end
    lastbot_traj = get!(lastbot,bid,pid)
    if field == :midvel
        lastbot_t = get_mid_times(lastbot)
    else
        lastbot_t = lastbot.traj.t
    end
    lastbot_dt = lastbot.traj.t[begin+1]-lastbot.traj.t[begin]
    lastbot_traj_itp = begin 
        itp = interpolate(lastbot_traj[di,:], BSpline(Linear()))
        scale(
            itp,
            range(lastbot_t[begin],lastbot_t[end];step=lastbot_dt)
        )
    end
    dts = [
        bot.traj.t[begin+1]-bot.traj.t[begin]
        for bot in bots[begin:end-1]
    ]
    err_avg = [
        begin
            if field == :midvel
                t = get_mid_times(bot)
            else
                (;t) = bot.traj
            end
            bot_traj_di = get!(
                bots[i],bid,pid
            )[di,begin:end-1]
            ref_traj_di = lastbot_traj_itp(
                t[begin:end-1]
            ) 
            bot_traj_di .- ref_traj_di .|> abs |> mean
        end
        for (i,bot) in enumerate(bots[1:nbots-1])
    ]
    dts,err_avg
end
