@unit pt "pt" Point (100//7227)*u"inch" false
Unitful.register(@__MODULE__)
@unit px "px" Pixel (3//4)*u"pt" false
Unitful.register(@__MODULE__)
to_resolution(dpi,len) = uconvert(Unitful.NoUnits,dpi*len)

pt2px(x, ppi = 300*u"px"/1u"inch") = ustrip(u"px",x*u"pt"*ppi)

fontsize::Float64 = 8 |> pt2px
markersize::Float64 = 0.5fontsize
linewidth::Float64 = 0.5 |> pt2px
# cablewidth = 0.75 |> pt2px
# barwidth = 1.5 |> pt2px
cw::Float64 = 455 |> pt2px
tw::Float64 = 455 |> pt2px
th::Float64 = 688.5 |> pt2px
455.24411
mks_cyc::Base.Iterators.Cycle{Vector{String}} = Iterators.cycle(["o","v","^","s","P","X","d","<",">","h"])
lss_cyc::Base.Iterators.Cycle{Vector{String}} = Iterators.cycle(["-","--","-.",":"])
alphabet::String = join('a':'z')
tenmarkers::Vector{Symbol} = [
    :xcross,:cross,
    :utriangle,:dtriangle,
    :ltriangle,:rtriangle,
    :diamond,:hexagon,
    :star8,:star5
]
tlabel::LaTeXString = L"t~(\mathrm{s})"

function hidex(ax)
    ax.xticklabelsvisible = false
    ax.xlabelvisible = false
end

function hidey(ax)
    ax.yticklabelsvisible = false
    ax.ylabelvisible = false
end

theme_try = Theme(;
    font = "CMU Serif",
    fontsize,
    markersize,
    linewidth,
    figure_padding = (2fontsize,fontsize,fontsize,fontsize),
    palette = (
        vlinecolor = [:slategrey],
    ),
    Axis3 = (
        titlefont = "CMU Serif Bold",
        titlesize = fontsize,
        titlegap = 0
    ),
    Mesh = (
        color = :slategrey,
        transparency = true
    )
)

set_theme!(theme_try)
theme_pub = Theme(;
    font = "CMU Serif",
    fontsize,
    markersize,
    linewidth,
    figure_padding = (fontsize,fontsize,fontsize,fontsize),
    resolution=(0.9tw,0.8tw),
    palette = (
        vlinecolor = [:slategrey],
    ),
    Axis = (
        titlefont = "CMU Serif Bold",
        titlesize = fontsize,
        titlegap = 0
    ),
    Axis3 = (
        titlefont = "CMU Serif Bold",
        titlesize = fontsize,
        titlegap = 0
    ),
    Label = (
        textsize = fontsize,
        font = "CMU Serif Bold",
        halign = :left,
        padding = (0, 0, 0, 0),
    ),
    VLines = (
        cycle = [:color => :vlinecolor],
    ),
    Mesh = (
        color = :slategrey,
        transparency = true
    )
)

function plot_rigid(rb::TR.AbstractRigidBody;showmesh=true,showupdatemesh=false)
    (;r̄g,r̄ps,nr̄ps) = rb.prop
    (;mesh) = rb
    fig = Figure(resolution=(1920,1080))
    ndim = TR.get_ndim(rb)
    if ndim == 2 && !showmesh
        ax = Axis(fig[1,1])
        ax.aspect = DataAspect()
    else
        ax = LScene(fig[1,1], scenekw = (clear=true,))
    end
    scatter!(ax,r̄ps,markersize=5)
    scatter!(ax,r̄g,markersize=5)
    text!(ax,
        ["$i $(string(r̄p))" for (i,r̄p) in enumerate(r̄ps)] ,
        position = r̄ps,
        color = :darkblue,
        align = (:left, :baseline),
        offset = (-5, 10)
    )
    text!(ax,
        "r̄g $(string(r̄g))",
        position = r̄g,
        color = :darkred,
        align = (:left, :baseline),
        offset = (-5, 10)
    )
    if !(mesh isa Nothing)
        if showmesh
            mesh!(ax,mesh;transparency=true)
        end
        if showupdatemesh
            mesh!(ax,update_rigidmesh(rb))
        end
    end
    fig
end

function plot_tg!(tg::TR.TensegrityStructure)
    ndim = TR.get_ndim(tg)
    fig = Figure(resolution=(1280,720))
    if ndim == 2
        ax = Axis(fig[1,1], aspect=:data)
    else
        ax = Axis3(fig[1,1], aspect=:data)
    end
    tgob = Observable(deepcopy(tg))
    init_plot!(ax,tgob)
    fig
end

function match_figsize(figsize)
    @match figsize begin
        :UHD => (3840,2160)
        :FHD => (1920,1080)
        :HD  => (1280,720)
        fs::Tuple => fs
    end
end

function time2step(at,t)
    atstep = findfirst((x)->x>=at, t)
    @assert !isa(atstep,Nothing)
    atstep
end

function get_groundmesh(f::Function,rect)
    GeometryBasics.Mesh(f, rect, NaiveSurfaceNets()) |> make_patch()
end

function get_groundmesh(plane::TR.Plane,rect)
    GeometryBasics.Mesh(rect, NaiveSurfaceNets()) do v
        TR.signed_distance(v,plane)
    end |> make_patch()
end

function get_groundmesh(::Nothing,rect)
    plane_n = [0,0,1.0]
    plane_r = zeros(3)
    plane = TR.Plane(plane_n,plane_r)
    get_groundmesh(plane,rect)
end

function fullupdate!(tg)
    TR.update!(tg)
end

function plot_traj!(bot::TR.TensegrityRobot;
            AxisType=LScene,
            figsize=:FHD,
            fig = Figure(resolution=match_figsize(figsize)),
            gridsize=(1,1),
            attimes=nothing,
            atsteps=nothing,
            doslide=true,
            showinit=false,
            showlabels=true,
            showpoints=true,
            showmesh=true,
            showwire=false,
            showinfo=true,
            xlims=(-1.0,1.0),
            ylims=(-1.0,1.0),
            zlims=(-1e-3,1.0),
            showground=true,
            ground=nothing,
            fontsize=20,
            actuate=false,
            figname=nothing,
            colorbar=nothing,
            showcb=false,
            rowgap=2fontsize,
            titleformatfunc = (sgi,tt)-> begin
                @sprintf "(%s) t = %.10G (s)" alphabet[sgi] tt
            end,
            sup! = (ax,tgob,sgi)->nothing,
            kargs...)
    (;tg,traj) = bot
    tg.state.system.q .= traj.q[begin]
    fullupdate!(tg)
    tgobini = Observable(deepcopy(tg))
    showmesh = mapreduce(&, TR.get_rigidbodies(tg)) do rb
            rb.mesh isa GeometryBasics.Mesh
    end
    ndim = TR.get_ndim(tg)
    xmin,xmax = xlims
    ymin,ymax = ylims
    zmin,zmax = zlims
    xwid = xmax - xmin
    ywid = ymax - ymin
    zwid = zmax - zmin
    grid1 = fig[1,1] = GridLayout()
    subgrids = [
        grid1[i,j] = GridLayout()
        for i in 1:gridsize[1], j = 1:gridsize[2]
    ] |> permutedims
    parsed_steps = @match (attimes,atsteps) begin
        (a::Nothing,b::Nothing) => [
                1 for sgi in eachindex(subgrids)
            ]
    	(a,b::Nothing) => [
                time2step(a[sgi],traj.t)
                for sgi in eachindex(subgrids)
            ]
        (a::Nothing,b) => b
        (a,b) => begin
            @warn "Ignoring `attimes`"
            b
        end
    end

    # @warn "Overwriting `atsteps`"
    # @warn "Ignoring `attimes`"
    rowgap!(grid1,rowgap)
	# colgap!(grid1,10fontsize)
    for sgi in eachindex(parsed_steps)
        if sgi > length(subgrids)
            sg = subgrids[end]
        else
            sg = subgrids[sgi]
        end
        this_step = parsed_steps[sgi]
        this_time = Observable(traj.t[this_step])
        tg.state.system.c .= traj.c[this_step]
        tg.state.system.q .= traj.q[this_step]
        fullupdate!(tg)
        tgob = Observable(deepcopy(tg))
        axtitle = lift(this_time) do tt
            titleformatfunc(sgi,tt)
        end
        if ndim == 2 && !showmesh
            # showinfo = false
            showground = false
            ax = Axis(sg[1,1],title=axtitle)
            ax.aspect = DataAspect()
            xlims!(ax,xmin,xmax)
            ylims!(ax,ymin,ymax)
        elseif AxisType <: Axis3
            ax = Axis3(sg[1,1],
                title=axtitle,
                aspect=:data,
            )
            xlims!(ax,xmin,xmax)
            ylims!(ax,ymin,ymax)
            zlims!(ax,zmin,zmax)
        elseif AxisType <: LScene
            showinfo = false
            # ax = LScene(fig[1,1], show_axis=false, scenekw = (clear=true,))
            ax = LScene(fig[1,1], ) #scenekw = (clear=true,)
            # cam = Makie.camera(ax.scene)
            # cam = cam3d!(ax.scene, projectiontype=Makie.Orthographic)
            # update_cam!(ax.scene, cam)
        else
            error("Unknown AxisType")
        end
        if showground
            rect = Rect3f((xmin,ymin,zmin),(xwid,ywid,zwid))
            groundmesh = get_groundmesh(ground,rect)
            mesh!(ax,groundmesh;color = :snow)
        end
        if showinit
            # init_plot!(ax,deepcopy(tgob);showmesh,showwire,isref=true,showpoints,kargs...)
            sup!(ax,tgobini,sgi)
            init_plot!(ax,tgobini;showmesh,showwire,isref=true,showpoints,kargs...)
        end
        if colorbar isa Nothing
            rigidcolor = :slategray4
        else
            rigidcolor = colorbar.colormap[sgi]
            if showcb
                Colorbar(sg[1, 2]; 
                        size = fontsize,
                        colorbar...
                    )
            end
        end
        sup!(ax,tgob,sgi)
        init_plot!(ax,tgob;
            showlabels,showmesh,showwire,fontsize,showpoints,
            rigidcolor,
            kargs...
        )

        if doslide
            if showinfo
                grid2 = sg[:,2] = GridLayout(;tellheight=false)
                grid_info = grid2[1,1] = GridLayout(;tellheight=false)
                dict_info = [
                    "fig. height" => lift(string,ax.height),
                    "fig. width" => lift(string,ax.width)
                ]
                if ndim == 3 && AxisType == Axis3
                    cam_info = [
                        "azimuth" => lift(string,ax.azimuth),
                        "elevation" => lift(string,ax.elevation)
                    ]
                    append!(dict_info,cam_info)
                end
                for (i,(infoname,infovalue)) in enumerate(dict_info)
                    Label(grid_info[i,1],
                        lift(infovalue) do iv
                            "$infoname = $iv"
                        end,
                        justification = :left,
                    )
                end
                grid_button = grid2[2,1] = GridLayout(;tellheight=false)
                # camera = cameracontrols(ax.scene)
                # lookat = camera.lookat
                # eyepos = camera.eyeposition
                botton_printinfo = Button(grid_button[1,1], label = "Print Infos")
                on(botton_printinfo.clicks) do _
                    for (infoname,infovalue) in dict_info
                        println("$infoname = $(infovalue[])")
                    end
                end
            end
            grid3 = sg[2,:] = GridLayout()
            slidergrid = SliderGrid(
                grid3[1,1],
                (label = "Step", range = 1:length(traj), startvalue = this_step),
            )
            on(slidergrid.sliders[1].value) do this_step
                if actuate
                    TR.actuate!(bot,[traj.t[this_step]])
                end
                tg.state.system.q .= traj.q[this_step]
                tg.state.system.c .= traj.c[this_step]
                fullupdate!(tg)
                # @show TR.mechanical_energy(tg)
                #
                this_time[] = traj.t[this_step]
                TR.analyse_slack(tg,true)
                tgob[] = tg
            end
        end
    end
    if fig isa Figure
        savefig(fig,figname)
        DataInspector(fig)
    end
    fig
end

function savefig(fig,figname=nothing)
    if !isa(figname,Nothing)
        if isdefined(Main,:figdir)
            figpath = joinpath(figdir,figname)
        else
            figpath = figname
        end
        if Makie.current_backend() == CM
            @info "Saving to $figpath.pdf"
            CM.save(figpath*".pdf",fig)
        else
            @info "Saving to $figpath.png"
            GM.save(figpath*".png",fig)
        end
    end
    fig
end

function init_plot!(ax,tgob;isref=false,
        showlabels=true,
        showpoints=true,
        showmesh=true,
        showwire=false,
        fontsize=10,
        cablecolor=:deepskyblue,
        cablelabelcolor=:darkgreen,
        rigidcolor=:slategray4,
        rigidlabelcolor=:darkblue
    )
    (;tensiles,nrigids) = tgob[]
    ncables = length(tensiles.cables)
    ndim = TR.get_ndim(tgob[])
    if isref
        showlabels = false
        cablecolor =
        cablelabelcolor =
        rigidcolor =
        rigidlabelcolor = :grey
    end
    # cables
    if ncables !== 0
    # if false
        linesegs_cables = @lift begin
            (;tensioned) = $tgob.connectivity
            ndim = TR.get_ndim($tgob)
            T = TR.get_numbertype($tgob)
            ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
            mapreduce(
                (scnt)->
                Point(scnt.end1.rbsig.state.rps[scnt.end1.pid]) =>
                Point(scnt.end2.rbsig.state.rps[scnt.end2.pid]),
                vcat,
                tensioned.connected
                ;init=ret
            )
        end
        linesegments!(ax, linesegs_cables, color = cablecolor, linewidth = 2)
        rcs_by_cables = @lift begin
            (;tensioned) = $tgob.connectivity
            ndim = TR.get_ndim($tgob)
            T = TR.get_numbertype($tgob)
            ret = Vector{MVector{ndim,T}}()
            mapreduce(
                (scnt)->
                [(
                    scnt.end1.rbsig.state.rps[scnt.end1.pid].+
                    scnt.end2.rbsig.state.rps[scnt.end2.pid]
                )./2],
                vcat,
                tensioned.connected
                ;init=ret
            )
        end
        # @show rcs_by_cables
        if showlabels
            text!(ax,
                ["c$(i)" for (i,rc) in enumerate(rcs_by_cables[])] ,
                position = rcs_by_cables,
                textsize = fontsize,
                color = cablelabelcolor,
                align = (:left, :top),
                offset = (-5, -10)
            )
        end
    end

    rbsob = @lift begin
        TR.get_rigidbodies($tgob)
    end
    nb = length(rbsob[])
    if showmesh || showwire
        for rbid = 1:nb
            rigid_mesh = @lift begin
                rb = ($rbsob)[rbid]
                update_rigidmesh(rb)
            end
            if showwire
                strokewidth = linewidth
            else
                strokewidth = 0
            end
            poly!(ax, rigid_mesh;
                    color = rigidcolor,
                    shading = true,
                    strokewidth
            )
            # mesh!(ax, rigid_mesh;
            #         color = rigidcolor,
            #         shading = true
            # )
        end
    end
    # mass centers
    rg_by_rbs = @lift begin
        [rb.state.rg for rb in $rbsob]
    end
    # points on rigidbodies
    rps_by_rbs = [
        @lift begin
            rbs = TR.get_rigidbodies($tgob)
            rbs[rbid].state.rps
        end
        for rbid = 1:nrigids
    ]
    if showpoints
        scatter!(ax,rg_by_rbs; color = rigidlabelcolor)
        for (rbid,rps) in enumerate(rps_by_rbs)
            scatter!(ax,rps)
        end
    end
    if showlabels
        text!(ax,
            ["r$(i)g" for (i,rg) in enumerate(rg_by_rbs[])] ,
            position = rg_by_rbs,
            textsize = fontsize,
            color = rigidlabelcolor,
            align = (:left, :top),
            offset = (-5, -10)
        )
        for (rbid,rps) in enumerate(rps_by_rbs)
            text!(ax,
                # ["r$(rbid)p$pid $(string(rp))" for (pid,rp) in enumerate(rps[])],
                ["r$(rbid)p$pid" for (pid,rp) in enumerate(rps[])],
                position = rps,
                textsize = fontsize,
                color = :darkred,
                align = (:left, :top),
                offset = (20(rand()-0.5), 20(rand()-0.5))
            )
        end
    end
end

function get3Dstate(rb)
    (;state) = rb
    (;ro,R,ṙo,ω) = state
    ndim = TR.get_ndim(rb)
    T = TR.get_numbertype(rb)
    o = zero(T)
    i = one(T)
    if ndim == 3
        return ro, R, ṙo, ω
    else
        ro3 = MVector{3}(ro[1],ro[2],o)
        ṙo3 = MVector{3}(ṙo[1],ṙo[2],o)
        R3 = MMatrix{3,3}(
            [
                R[1,1] -R[2,1] o;
                -R[1,2] R[2,2] o;
                o      o      i;
            ]
        )
        ω3 = MVector{3}(o,o,ω[1])
        return ro3, R3, ṙo3, ω3
    end
end

function update_rigidmesh(rb)
    (;mesh) = rb
    @assert !(mesh isa Nothing)
    ro,R,ṙo,ω = get3Dstate(rb)
    trans = Translation(ro)
    rot = LinearMap(R)
    ct = trans ∘ rot
    updated_pos = ct.(mesh.position)
    fac = GeometryBasics.faces(mesh)
    nls = normals(updated_pos,fac)
    GeometryBasics.Mesh(meta(updated_pos,normals=nls),fac)
end


function make_patch(;trans=[0.0,0,0],rot=RotX(0.0))
    function patch(mesh)
        ct = Translation(trans) ∘ LinearMap(rot)
        updated_pos = ct.(mesh.position)
        fac = GeometryBasics.faces(mesh)
        nls = normals(updated_pos,fac)
        GeometryBasics.Mesh(meta(updated_pos,normals=nls),fac)
    end
end

function endpoints2mesh(p1,p2;long=400)
	s = Meshes.Segment(Meshes.Point(p1),Meshes.Point(p2))
	cyl_bar = Meshes.Cylinder(Meshes.length(s)/long,s)
    cylsurf_bar = Meshes.boundary(cyl_bar)
    # Meshes.sample(cylsurf_bar,Meshes.RegularSampling(10,3))
	cyl_bar_simple = Meshes.discretize(cylsurf_bar,Meshes.RegularDiscretization(10,2))
    cyl_bar_simple |> simple2mesh
end


function parse_Adams_dynamic(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    entities = findall("//adams:Entity[@entType='ADAMS_Variable']", adams_xml.root, ns)
    names = Symbol.(vcat("time",[replace(ele["entity"],"."=>"_") for ele in entities]))
    step_initialConditions = findfirst("//adams:Step[@type='initialConditions']",adams_xml.root, ns)
    initialConditions = (;zip(names,parse.(Float64,split(step_initialConditions.content,"\n")[2:end-1]))...)
    ret = StructArray([initialConditions])
    step_dynamic = findall("//adams:Step[@type='dynamic']",adams_xml.root, ns)
    for step in step_dynamic
        step_values = parse.(Float64,split(step.content,"\n")[2:end-1])
        # @assert length(variable_names) == length(step_data)
        push!(ret,(;zip(names,step_values)...))
    end
    ret
end

function parse_Adams_static(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    entities = findall("//adams:Entity[@entType='ADAMS_Variable']", adams_xml.root, ns)
    names = Symbol.(vcat("time",[replace(ele["entity"],"."=>"_") for ele in entities]))
    step_initialConditions = findfirst("//adams:Step[@type='initialConditions']",adams_xml.root, ns)
    step_static = findfirst("//adams:Step[@type='static']",adams_xml.root, ns)
    step_values = parse.(Float64,split(step_static.content,"\n")[2:end-1])
    (;zip(names,step_values)...)
end

function parse_Adams_frequency(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    frequency_node = findfirst("//adams:ModalEnergyHeader[@type='undamped_natural_freqs']",adams_xml.root, ns)
    frequency_string = split(frequency_node.content,"\n")[2]
    frequency = parse.(Float64,split(frequency_string," "))
end

function simple2mesh(sp)
	dim   = Meshes.embeddim(sp)
	nvert = Meshes.nvertices(sp)
	nelem = Meshes.nelements(sp)
	verts = Meshes.vertices(sp)
	topo  = Meshes.topology(sp)
	elems = Meshes.elements(topo)

	# coordinates of vertices
	coords = Meshes.coordinates.(verts)
	# fan triangulation (assume convexity)
	tris4elem = map(elems) do elem
	  I = Meshes.indices(elem)
	  [[I[1], I[i], I[i+1]] for i in 2:length(I)-1]
	end

	# flatten vector of triangles
	tris = [tri for tris in tris4elem for tri in tris]

	points  = Point.(coords)
	faces  = TriangleFace{UInt64}.(tris)
	# tcolors  = colorant

	# convert connectivities to matrix format
	# tmatrix = reduce(hcat, tconnec) |> transpose |> Matrix
    nls = normals(points,faces)
    GeometryBasics.Mesh(meta(points,normals=nls),faces)
end

function plotsave_energy(bot,figname=nothing)
	(;traj) = bot
	(;t) = traj
	tmin,tmax = t[begin], t[end]
	titles = [
		"机械能",
		"动能",
		"势能"
	]
	with_theme(theme_pub; resolution = (0.9tw,0.6tw)) do
		ME = TR.mechanical_energy!(bot;gravity=true)
		fig = Figure()
		axs = [
			begin
				ax = Axis(
					fig[i,1],
					xlabel=L"t~\mathrm{(s)}",
					ylabel="Energy (J)",
					title=titles[i]
				)
				xlims!(ax,tmin,tmax)
				# ylims!(ax,,)
				if i !== 3
					ax.xlabelvisible = false
					ax.xticklabelsvisible = false
				end
				Label(fig[i,1,TopLeft()], "($(alphabet[i]))")
				ax
			end
			for i = 1:3
		]

		lines!(axs[1], t, ME.E)
		lines!(axs[2], t, ME.T)
		lines!(axs[3], t, ME.V)


		savefig(fig,figname)

		fig
	end
end


function plotsave_friction_direction(bots,x,xs,figname=nothing;
        resolution = (0.9tw,0.4tw),
        mo = 8,
        mo_α = mo,
        vtol=1e-5,
    )
	with_theme(theme_pub;
        resolution
	) do
		fig = Figure()
		for (botid,bot) in enumerate(bots)
			ax1 = Axis(fig[botid,1],
						xlabel=tlabel,
						ylabel=L"\delta\alpha~(\mathrm{Rad})",
						title=latexstring("$x=$(xs[botid])")
				)
			ax2 = Axis(fig[botid,2],
						xlabel=tlabel,
						ylabel=L"\delta\theta~(\mathrm{Rad})",
						title=latexstring("$x=$(xs[botid])")
				)
			Label(fig[botid,1,TopLeft()],"($(alphabet[2botid-1]))")
			Label(fig[botid,2,TopLeft()],"($(alphabet[2botid]))")
			(;t) = bot.traj
			contacts_traj_voa = VectorOfArray(bot.contacts_traj)
            if eltype(xs) <: Int
			    c1s = contacts_traj_voa[xs[botid],:]
            else
                c1s = contacts_traj_voa[1,:]
            end
			idx_sli = findall(c1s) do c
						issliding(c;vtol)
					end
			idx_imp = findall(isimpact,c1s) ∩ idx_sli
			δα_imp = map(c1s[idx_imp]) do c
						get_contact_angle(c)
					end
			idx_per = findall(doespersist,c1s) ∩ idx_sli
            # @show idx_per[begin]
			δα_per = map(c1s[idx_per]) do c
						get_contact_angle(c)
					end
			scaling = 10.0^(-mo_α)
			Label(fig[botid,1,Top()],latexstring("\\times 10^{-$(mo_α)}"))
			markersize = fontsize
			scatter!(ax1,idx_per,δα_per./scaling;
						marker=:diamond, markersize)
			scatter!(ax1,idx_imp,δα_imp./scaling;
						marker=:xcross, markersize)
			ylims!(ax1,-1.0,1.0)
			# xlims!(ax1,extrema(t)...)

			θ_imp = map(c1s[idx_imp]) do c
					get_friction_direction(c)
				end
			θ_per = map(c1s[idx_per]) do c
					get_friction_direction(c)
				end
			scaling = 10.0^(-mo)
			Label(fig[botid,2,Top()],latexstring("\\times 10^{-$(mo)}"))
			scatter!(ax2,t[idx_per],(θ_per.-π/4)./scaling, label="Persistent";
						marker=:diamond, markersize)
			scatter!(ax2,t[idx_imp],(θ_imp.-π/4)./scaling, label="Impact";
						marker=:xcross, markersize)
			xlims!(ax2,extrema(t)...)
			ylims!(ax2,-1.0,1.0)
			if botid !== length(bots)
				hidex(ax1)
				hidex(ax2)
			else
				Legend(fig[length(bots)+1,:],ax2;
					orientation=:horizontal
				)
			end
		end
		rowgap!(fig.layout,fontsize/2)
		savefig(fig,figname)
		fig
	end
end


function plotsave_error(mbots,dts,figname=nothing;
        resolution = (0.8tw,0.3tw),
    )
	with_theme(theme_pub;
            resolution,
            figure_padding = (0,fontsize,0,0),
		) do
		fig = Figure()
		for j in axes(mbots,2)
			bots = mbots[:,j]
			nbots = length(bots)
			itp = interpolate(get_trajectory!(bots[nbots],1,0)[1,:], BSpline(Linear()))
			ref_traj = scale(itp,0:dts[nbots]:bots[nbots].traj.t[end])
			# ref_traj(0:0.1:2)
			err_avg = [
				get_trajectory!(bots[i],1,0)[1,begin:end-1].-ref_traj(bots[i].traj.t[begin:end-1]) .|> abs |> mean
				for (i,dt) in enumerate(dts[1:nbots-1])
			]

			ax = Axis(fig[1,j];
				xlabel = L"h",
				ylabel = "Avg. Err.",
				yscale = Makie.log10,yminorticksvisible = true, yminorgridvisible = true,yminorticks = IntervalsBetween(8),
				xscale = Makie.log10,xminorticksvisible = true, xminorgridvisible = true,xminorticks = IntervalsBetween(8),
				)
			Label(fig[1,j,TopLeft()],"($(alphabet[j]))")
			scatterlines!(ax,dts[1:nbots-1],err_avg;marker=:rect,color=:red,label="Error of NMSI")
			lines!(ax,dts[1:nbots-1],err_avg[1] .*(dts[1:nbots-1]./dts[1]).^2,label="2nd-Order")
			lines!(ax,dts[1:nbots-1],err_avg[1] .*(dts[1:nbots-1]./dts[1]),label="1st-Order")
			# ylims!(ax,1e-11,1e-2)
			if j == 1
				axislegend(ax;position=:rb)
			else
				hidey(ax)
			end
		end
		savefig(fig,figname)
		fig
	end
end