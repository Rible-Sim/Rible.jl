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
        cablewidth=2,
        meshcolor=nothing,
        fontsize = 12,
    )
end

function MakieCore.plot!(viz::Viz{Tuple{S}};
    ) where S <:AbstractBody
    body_ob = viz[:structure]
    # body decorations
    mass_center_ob = lift(body_ob) do body_ob
        body_ob.state.mass_locus_state.position |> GB.Point
    end
    nodes_ob = lift(body_ob) do body_ob
        [
            locus_state.position |> GB.Point
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
        MakieCore.scatter!(viz,nodes_ob;color=viz.pointcolor[])
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

function MakieCore.plot!(viz::Viz{Tuple{Vector{S}}};
    ) where S <: Cable
    cables_ob = viz[:structure]
    point_mid_ob = lift(cables_ob) do cables_ob
        [
            begin 
                point_start = cab.state.start
                point_stop = cab.state.stop 
                (point_start .+ point_stop)./2 |> GB.Point
            end
            for cab in cables_ob
        ]
    end

    slackseg_ob = lift(cables_ob) do cables_ob
        [
            begin 
                point_start = cab.state.start |> GB.Point
                point_stop = cab.state.stop |> GB.Point
                (point_start,point_stop)
            end
            for cab in cables_ob if cab.state.length <= cab.state.restlen
        ]
    end

    noslackseg_ob = lift(cables_ob) do cables_ob
        [
            begin 
                point_start = cab.state.start |> GB.Point
                point_stop = cab.state.stop |> GB.Point
                (point_start,point_stop)
            end
            for cab in cables_ob if cab.state.length > cab.state.restlen
        ]
    end

    ids_ob = lift(cables_ob) do cables_ob
        [
            "c$(cab.id)"
            for cab in cables_ob
        ]
    end
    # slackonly=false,
    # noslackonly=true
    MakieCore.linesegments!(
        viz, noslackseg_ob, 
        color = viz.cablecolor[], 
        linewidth = viz.cablewidth[], 
        linestyle = :solid
    )
    MakieCore.linesegments!(
        viz, slackseg_ob, 
        color = viz.cablecolor[], 
        linewidth = viz.cablewidth[], 
        linestyle = viz.slack_linestyle[]
    )
    # show cable labels
    if viz.show_cable_labels[]
        MakieCore.text!(viz,
            ids_ob,
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
    (;tensiles,num_of_bodies) = tgob[]
    ncables = length(tensiles.cables)
    if ncables > 0 && viz.showcables[]
        cables_ob = lift(tgob) do tgob
            tgob.tensiles.cables
        end
        viz!(viz,cables_ob;
            cablecolor,
            cablewidth = viz.cablewidth[],
            cablelabelcolor = viz.cablelabelcolor[],
            slack_linestyle = viz.slack_linestyle[],
            show_cable_labels,
        )
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
        origin_position,R,_ = get3Dstate(body)
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

theme_pub = MakieCore.Attributes(;
    size = (1280,720),
    palette = (
        vlinecolor = [:slategrey],
        linestyle = [
            :solid,
            :dash, 
            :dot,
            :dashdot,
            :dashdotdot,
        ],
        marker = [
            :xcross,:cross,
            :utriangle,:dtriangle,
            :ltriangle,:rtriangle,
            :diamond,:hexagon,
            :star8,:star5
        ],
        markercolor = [:red, :blue],
    ),
    Axis = (
        # titlefont = "CMU Serif Bold",
        # titlesize = fontsize,
        titlegap = 0,
    ),
    Axis3 = (
        # titlefont = "CMU Serif Bold",
        # titlesize = fontsize,
        titlegap = 0,
    ),
    Label = (
        # fontsize = fontsize,
        # font = "CMU Serif Bold",
        halign = :left,
        padding = (0, 0, 0, 0),
    ),
    VLines = (
        cycle = [:color => :vlinecolor],
    ),
    Mesh = (
        color = :slategrey,
        transparency = false
    ),
    Poly = (
        color = :slategrey,
        transparency = false,
    )
)

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
