Unitful.register(@__MODULE__)
@unit pt "pt" Point (100//7227)*u"inch" false
@unit px "px" Pixel (3//4)*u"pt" false
to_resolution(dpi,len) = uconvert(Unitful.NoUnits,dpi*len)

pt2px(x, ppi = 300*u"px"/1u"inch") = ustrip(u"px",x*u"pt"*ppi)
px2pt(x, ppi = 300*u"px"/1u"inch") = ustrip(u"pt",x*u"px"/ppi)

fontsize::Float64 = 8 |> RB.pt2px
markersize::Float64 = 0.5fontsize
linewidth::Float64 = 0.5 |> RB.pt2px
# cablewidth = 0.75 |> RB.pt2px
# barwidth = 1.5 |> RB.pt2px
cw::Float64 = 455 |> RB.pt2px
tw::Float64 = 455 |> RB.pt2px
th::Float64 = 688.5 |> RB.pt2px
# 455.24411
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

macro myshow(exs...)
    blk = Expr(:block)
    for ex in exs
        push!(blk.args, :(print($(sprint(Base.show_unquoted,ex)*" = "))))
        push!(blk.args, :(show(stdout, "text/plain", begin value=$(esc(ex)) end)))
        push!(blk.args, :(println()))
    end
    isempty(exs) || push!(blk.args, :value)
    return blk
end

set_theme!(RB.theme_pub)

# fonts = (; 
#     regular = "CMU Serif", 
#     bold = "CMU Serif Bold",
#     italic = "CMU Serif Italic",
#     math = "NewComputerModern 10 Italic"
# )
# fontsize,
# markersize,
# linewidth,
# figure_padding = (fontsize,fontsize,fontsize,fontsize),
    
function RB.get_groundmesh(f::Function,rect)
    GB.Mesh(f, rect, NaiveSurfaceNets()) |> make_patch(;color = :snow)
end

function RB.get_groundmesh(plane::RB.Plane,rect)
    GB.Mesh(rect, MarchingCubes()) do v
        RB.signed_distance(v,plane)
    end |> RB.make_patch(;color = :snow)
end

function plot_traj!(bot::RB.Robot;
        AxisType=LScene,
        figsize=:FHD,
        fig = Figure(),
        gridsize=(1,1),
        attimes=nothing,
        atsteps=nothing,
        doslide=true,
        dorecord=false,
        speedup=1,
        showlabels=true,
        showpoints=true,
        showmesh=true,
        showwire=false,
        showtitle=true,
        showcables=true,
        showinfo=true,
        xlims=(-1.0,1.0),
        ylims=(-1.0,1.0),
        zlims=(-1e-3,1.0),
        showground=true,
        ground=nothing,
        showback=false,
        actuate=false,
        figname=nothing,
        showinit=false,
        auto=false,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font=:bold),
                (@sprintf "t = %.10G (s)" tt)
            )
        end,
        sup! = (ax,tgob,sgi)->nothing,
        kargs...
    )
    (;structure,traj) = bot
    ndim = RB.get_num_of_dims(structure)
    structure.state.system.q .= traj.q[begin]
    RB.update!(structure)
    tgobini = Observable(deepcopy(structure))
    xmin, xmax = xlims
    ymin, ymax = ylims
    zmin, zmax = zlims
    xwid = xmax - xmin
    ywid = ymax - ymin
    zwid = zmax - zmin
    grid1 = fig[1,1] = GridLayout(;tellheight=false)
    subgrids = [
        grid1[i,j] = GridLayout(;tellheight=false)
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
    for sgi in eachindex(parsed_steps)
        if sgi > length(subgrids)
            sg = subgrids[end]
        else
            sg = subgrids[sgi]
        end
        this_step = parsed_steps[sgi]
        this_time = Observable(traj.t[this_step])
        RB.goto_step!(bot,this_step;actuate)
        tgob = Observable(deepcopy(structure))
        axtitle = map(this_time) do tt
            titleformatfunc(sgi,tt)
        end
        if ndim == 2 && !showmesh
            # showinfo = false
            showground = false
            ax = Axis(sg[1,1],)
            ax.aspect = DataAspect()
            xlims!(ax,xmin,xmax)
            ylims!(ax,ymin,ymax)
        elseif AxisType <: Axis3
            ax = Axis3(sg[1,1],
                # title=axtitle,
                aspect=:data,
            )
            xlims!(ax,xmin,xmax)
            ylims!(ax,ymin,ymax)
            zlims!(ax,zmin,zmax)
        elseif AxisType <: LScene
            showinfo = false
            # ax = LScene(fig[1,1], show_axis=false, scenekw = (clear=true,))
            ax = LScene(fig[1,1],) #
            # cam = Makie.camera(ax.scene)
            # cam = cam3d!(ax.scene, projectiontype=Makie.Orthographic)
            # update_cam!(ax.scene, cam)
        else
            error("Unknown AxisType")
        end
        if showground
            rect = Rect3f((xmin,ymin,zmin),(xwid,ywid,zwid))
            groundmesh = RB.get_groundmesh(ground,rect)
            mesh!(ax,groundmesh;color = :snow)
        end
        if showwire || showmesh || showcables || showlabels || showpoints
            if showinit
                viz!(ax,tgobini;
                    showmesh,
                    showwire,
                    isref=true,
                    showcables=false,
                    showpoints,
                    kargs...
                )
            end
            viz!(ax,tgob;
                showmesh,
                showwire,
                showlabels,
                showpoints,
                showcables,
                kargs...
            )
        end
        sup!(ax,tgob,sgi)
        if showtitle    
            Label(
                sg[1, 1, Top()],
                axtitle,
                padding = (0,0,0,0),
                justification = :right,
                lineheight = 1.0,
                halign = :center,
                valign = :bottom,
            )
        end
        # add background
        if showback
            subWinWidth,subWinHeight = match_figsize(figsize)
            subWindow0 = Scene(ax.scene, 
                px_area=Rect(0, 0, subWinWidth, subWinHeight), 
                clear=true, 
                backgroundcolor=:white
            )
            campixel!(subWindow0;farclip=1)
            bg = load(RB.assetpath("stars.jpg")) |> rotr90
            image!(subWindow0,
                LinRange(0,subWinWidth,size(bg,1)),
                LinRange(0,subWinHeight,size(bg,2)), 
                bg;
                inspectable=false
            )
        end
        if doslide || dorecord
            if showinfo
                grid2 = sg[:,2] = GridLayout(;tellheight=false)
                grid_info = grid2[1,1] = GridLayout(;tellheight=false)
                dict_info = [
                    "fig. height" => map(string,ax.height),
                    "fig. width" => map(string,ax.width)
                ]
                if ndim == 3 && AxisType == Axis3
                    cam_info = [
                        "azimuth" => map(string,ax.azimuth),
                        "elevation" => map(string,ax.elevation)
                    ]
                    append!(dict_info,cam_info)
                end
                for (i,(infoname,infovalue)) in enumerate(dict_info)
                    Label(grid_info[i,1],
                        map(infovalue) do iv
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
            if dorecord
                dt = traj.t[begin+1] - traj.t[begin]
                framerate = 30 
                skipstep = round(Int,1/framerate/dt*speedup)
                recordsteps = 1:skipstep:length(traj.t)
                record(fig, figname, recordsteps;
                    framerate) do this_step
                    if actuate
                        RB.actuate!(bot,[traj.t[this_step]])
                    end
                    structure.state.system.q .= traj.q[this_step]
                    structure.state.system.c .= traj.c[this_step]
                    RB.update!(structure)
                    # @show RB.mechanical_energy(structure)
                    #
                    this_time[] = traj.t[this_step]
                    RB.analyse_slack(structure,true)
                    tgob[] = structure
                end
            else
                grid3 = sg[2,:] = GridLayout()
                slidergrid = SliderGrid(
                    grid3[1,1],
                    (label = "Step", range = 1:length(traj), startvalue = this_step),
                )
                on(slidergrid.sliders[1].value) do this_step
                    if actuate
                        RB.actuate!(bot,[traj.t[this_step]])
                    end
                    structure.state.system.q .= traj.q[this_step]
                    structure.state.system.c .= traj.c[this_step]
                    RB.update!(structure)
                    # @show RB.mechanical_energy(structure)
                    #
                    this_time[] = traj.t[this_step]
                    RB.analyse_slack(structure,true)
                    tgob[] = structure
                    if auto
                        if AxisType <: LScene
                            center!(ax.scene)
                        elseif AxisType <: Axis3
                            autolimits!(ax)
                        end
                    end
                end
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

function get_linesegs_cables(structure;slackonly=false,noslackonly=false)
    (;connected) = structure.connectivity.tensioned
    (;cables) = structure.tensiles
    ndim = RB.get_num_of_dims(structure)
    T = RB.get_numbertype(structure)
    linesegs_cables = Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
    foreach(connected) do scnt
        scable = cables[scnt.id]
        ret = (Point(scnt.hen.rbsig.state.loci_states[scnt.hen.pid].position),
                Point(scnt.egg.rbsig.state.loci_states[scnt.egg.pid].position))
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

function RB.build_mesh(fb::RB.FlexibleBody,nsegs=100;color=:slategrey)
    (;state) = fb
    (;cache) = state
    (;funcs,e) = cache
    (;ancs) = funcs
    (;L,radius) = ancs
    T = typeof(L)
    V = T <: AbstractFloat ? T : Float64
    _r = RB.ANCF.make_r(ancs,e)
    _rₓ = RB.ANCF.make_rₓ(ancs,e)

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
    parsedcolor = parse(Makie.ColorTypes.RGB{Float32},color)
    colors = fill(parsedcolor,length(points))
    GB.Mesh(GB.meta(points,normals=nls,color=colors),faces)
end

function plotsave_friction_direction(bots,x,xs,figname=nothing;
        resolution = (0.9tw,0.4tw),
        mo = 8,
        mo_α = mo,
        vtol=1e-5,
        cid = 1,
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
                c1s = contacts_traj_voa[cid,:]
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
            scatter!(ax1,t[idx_per],δα_per./scaling;
                        marker=:diamond, markersize)
            scatter!(ax1,t[idx_imp],δα_imp./scaling;
                        marker=:xcross, markersize)
            ylims!(ax1,-1.0,1.0)
            xlims!(ax1,extrema(t)...)

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

function get_err_avg(bots;bid=1,pid=1,di=1,field=:traj)
    nbots = length(bots)
    lastbot = bots[nbots]
    if field == :traj
        get! = RB.get_trajectory!
    elseif field == :vel
        get! = RB.get_velocity!
    end
    lastbot_traj = get!(lastbot,bid,pid)
    lastbot_t = lastbot.traj.t
    lastbot_dt = lastbot_t[begin+1]-lastbot_t[begin]
    lastbot_traj_itp = begin 
        itp = interpolate(lastbot_traj[di,:], BSpline(Linear()))
        scale(
            itp,
            lastbot_t[begin]:lastbot_dt:lastbot_t[end]
        )
    end
    dts = [
        bot.traj.t[begin+1]-bot.traj.t[begin]
        for bot in bots[begin:end-1]
    ]
    err_avg = [
        begin
            (;t) = bot.traj
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

function plot_convergence_order!(ax,dts,err_avg;show_orders=true)

    scatterlines!(ax,
        dts,
        err_avg;
        marker=:rect,
        color=:red,
        label="NMSI"
    )
    

    if show_orders
        o2 = err_avg[1] .*(dts./dts[1]).^2
        lines!(ax,dts,o2,label=L"\mathcal{O}(h^2)")
        o1 = err_avg[1] .*(dts./dts[1])
        lines!(ax,dts,o1,label=L"\mathcal{O}(h)")
    end

    ax.xlabel = L"h~(\mathrm{s})"
    ax.yscale = Makie.log10
    ax.yminorticksvisible = true 
    ax.yminorgridvisible = true 
    ax.yminorticks = IntervalsBetween(8)
    ax.xscale = Makie.log10
    ax.xminorticksvisible = true 
    ax.xminorgridvisible = true 
    ax.xminorticks = IntervalsBetween(8)

    ax
end

function plotsave_contact_persistent(bot,figname=nothing;
    cid=1,
    tol = 1e-7 # rule out false active
    )
    with_theme(theme_pub;) do
        contacts_traj_voa = VectorOfArray(bot.contacts_traj)
        c1_traj = contacts_traj_voa[cid,:]
        # steps = 1:length(c1_traj)
        active_nonpersist = findall(c1_traj) do c
            c.state.active && !c.state.persistent && norm(c.state.Λ) > tol
        end
        active_persist = findall(c1_traj) do c
            c.state.active && c.state.persistent && norm(c.state.Λ) > tol
        end
        # @myshow active_persist
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = "Step", ylabel = "Contact Type")
        scatter!(ax,active_nonpersist,one.(active_nonpersist),label="Impact")
        scatter!(ax,active_persist,zero.(active_persist),label="Persistent")
        axislegend(ax)
        savefig(fig,figname)
        fig
    end
end

function plot_self_stress_states(
        botinput,
        S;
        rtol = 1e-14
        # Ň = build_nullspace_on_free(bot.structure)
    )
    bot = deepcopy(botinput)
    @myshow S
    maxS = maximum(abs.(S))
    Sbool = S.> maxS*rtol
    S[.!Sbool] .= 0.0
    ns = size(S,2)
    with_theme(theme_pub;
        ) do 
        plot_traj!(
            bot;
            AxisType=Axis3,
            gridsize=(1,ns), 
            doslide=false,
            showlabels=false,
            showpoints=false,
            showcables = false,
            xlims = (-4e-1,4e-1),
            ylims = (-4e-1,4e-1),
            zlims = (-4e-1,4e-1),
            showground = false,
            sup! = (ax,tgob,sgi)-> begin
            # cables
                hidexyz(ax)
                @myshow Sbool[:,sgi]
                linesegs_cables = @lift begin
                    get_linesegs_cables($tgob;)[Sbool[:,sgi]]
                end
                linesegments!(ax, 
                    linesegs_cables, 
                    color = :red, 
                    # linewidth = cablewidth
                    )
                rcs_by_cables = @lift begin
                    (;tensioned) = $tgob.connectivity
                    ndim = RB.get_num_of_dims($tgob)
                    T = RB.get_numbertype($tgob)
                    ret = Vector{MVector{ndim,T}}()
                    mapreduce(
                        (scnt)->
                        [(
                            scnt.hen.rbsig.state.loci_states[scnt.hen.pid].+
                            scnt.egg.rbsig.state.loci_states[scnt.egg.pid]
                        )./2],
                        vcat,
                        tensioned.connected
                        ;init=ret
                    )
                end
                # @show rcs_by_cables
                Stext = [
                        @sprintf "%4.2f"  S[i,sgi] 
                        for i in axes(S,1)
                        if Sbool[i,sgi]
                    ]
                @myshow Stext
                text!(
                    ax,
                    Stext,
                    position = rcs_by_cables[][Sbool[:,sgi]],
                    fontsize = fontsize,
                    color = :red,
                    align = (:right, :top),
                    offset = (-fontsize/2, 0)
                )
            end
        )
    end
end

function plot_kinematic_indeterminacy(
        botinput,
        D,
        Ň,
        # Ň = build_nullspace_on_free(bot.structure)
    )
    bot = deepcopy(botinput)
    nk = size(D,2)
    δq̌ = [Ň*D[:,i] for i in axes(D,2)]
    scaling=0.3
    for i = 1:nk
        push!(bot.traj,deepcopy(bot.traj[end]))
        bot.traj.t[end] = i
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        bot.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    plot_traj!(
        bot;
        AxisType=Axis3,
        gridsize=(1,nk),        
        atsteps=2:nk+1,
        doslide=false,
        showlabels=false,
        showpoints=false,
        # showcables = false,
        showground = false,
        # xlims = (-4e-2,4e-2),
        # ylims = (-4e-2,4e-2),
        # zlims = (-4e-2,4e-2),        
        slack_linestyle = :solid,
        showinit = true,
    )
end

function plotsave_contactpoints(bot,figname=nothing)
    contacts_traj_voa = VectorOfArray(bot.contacts_traj)
    (;t) = bot.traj
    with_theme(theme_pub;
            figure_padding = (0,fontsize,0,0),
            resolution = (0.9tw,0.3tw),
        ) do 
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = tlabel, ylabel = "Contact No.")
        markersize = fontsize
        nc = size(contacts_traj_voa,1)
        for ic in eachindex(contacts_traj_voa[1])
            c_traj = contacts_traj_voa[ic,:]
            idx_imp = findall((x)->isimpact(x;Λtol=0), c_traj)
            idx_per = findall((x)->doespersist(x;Λtol=0),c_traj)
            scatter!(ax,t[idx_per],fill(ic,length(idx_per)); marker=:diamond, markersize)
            scatter!(ax,t[idx_imp],fill(ic,length(idx_imp)); marker=:xcross, markersize)
        end
        xlims!(ax,t[begin],t[end])
        ylims!(ax,0.5,nc+0.5)
        ax.yticks = 1:nc
        elem_1 = [MarkerElement(;color = :blue, marker = :xcross, markersize)]
        elem_2 = [MarkerElement(;color = :blue, marker = :diamond, markersize)]
        Legend(fig[1, 1],
            [elem_1, elem_2],
            ["Impact", "Persistent"],
            tellheight = false,
            tellwidth = false,
            orientation = :horizontal, 
            halign = :right,
            valign = :top,
            margin = (0,fontsize,fontsize,2fontsize)
        )
        savefig(fig,figname)
        fig
    end
end