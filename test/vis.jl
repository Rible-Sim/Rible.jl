Unitful.register(@__MODULE__)
@unit pt "pt" Point (100//7227)*u"inch" false
@unit px "px" Pixel (3//4)*u"pt" false
to_resolution(dpi,len) = uconvert(Unitful.NoUnits,dpi*len)

pt2px(x, ppi = 300*u"px"/1u"inch") = ustrip(u"px",x*u"pt"*ppi)
px2pt(x, ppi = 300*u"px"/1u"inch") = ustrip(u"pt",x*u"px"/ppi)

fontsize::Float64 = 8 #|> pt2px
markersize::Float64 = 2.0#fontsize
linewidth::Float64 = 0.5 #|> pt2px
tw::Float64 = 455.24411 #|> pt2px
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

theme_pub = Makie.Attributes(;
    fonts = (; 
        regular = "CMU Serif", 
        bold = "CMU Serif Bold",
        italic = "CMU Serif Italic",
        math = "NewComputerModern 10 Italic"
    ),
    size = (1280,720),
    fontsize,#pt
    linewidth,
    markersize = 3.2,
    rowgap = fontsize,
    colgap = fontsize,
    figure_padding = (0,fontsize,0,0),
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
        spinewidth=0.2,
        xticksize = 0.4fontsize, yticksize = 0.4fontsize, 
        xticklabelpad = 0.05fontsize, yticklabelpad = 0.1fontsize, 
        xgridwidth=0.2,ygridwidth=0.2,
        xminorgridwidth=0.2,yminorgridwidth=0.2,
        xtickwidth = 0.2, ytickwidth = 0.2,
        xminorticksize = 2.0, yminorticksize = 2.0,
        xminortickwidth = 0.1, yminortickwidth = 0.1,
        xlabelpadding = 0.1fontsize, ylabelpadding = 0.1fontsize,
    ),
    Axis3 = (
        # titlefont = "CMU Serif Bold",
        # titlesize = fontsize,
        titlegap = 0,
        xticksize = 0.4fontsize, yticksize = 0.4fontsize, zticksize = 0.4fontsize,
        xspinewidth=0.2,yspinewidth=0.2,zspinewidth=0.2,
        xgridwidth=0.2,ygridwidth=0.2,zgridwidth=0.2,
        xminorgridwidth=0.2,yminorgridwidth=0.2,zminorgridwidth=0.2,
        xtickwidth = 0.2, ytickwidth = 0.2,ztickwidth=0.2,
        xminortickwidth = 0.2, yminortickwidth = 0.2, zminortickwidth=0.2,
        xlabeloffset = 1.5fontsize, ylabeloffset = 1.5fontsize, zlabeloffset = 2fontsize,
        xticklabelpad = 0.05fontsize, yticklabelpad = 0.05fontsize, zticklabelpad = 0.05fontsize, 
        protrusions = (2fontsize,2fontsize,2fontsize,0)
    ),
    Label = (
        # fontsize = fontsize,
        # font = "CMU Serif Bold",
        halign = :left,
        padding = (0, 0, 0, 0),
        justification = :right,
        lineheight = 1.0,
        valign = :bottom,
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
    ),
    Legend = (
        # alignmode, 
        # backgroundcolor, bgcolor, 
        # framecolor, framevisible, 
        framewidth = 0.2, 
        # gridshalign, gridsvalign, 
        groupgap = 0, 
        # tellheight, tellwidth, 
        # halign, height, 
        # valign, width,
        # margin, 
        # nbanks, orientation, 
        padding = (0.2fontsize, 0.2fontsize, 0.1fontsize, 0.1fontsize,), 
        # label, labelcolor, labelfont, labelhalign, labelvalign, labeljustification, 
        labelsize = 0.8fontsize, 
        # linecolor, linepoints, linestyle, 
        # linewidth = 0.2, 
        # marker, markercolor, markerpoints, markerstrokecolor, 
        # markersize = 1,
        # markerstrokewidth = fontsize, 
        # patchcolor, 
        patchlabelgap = 0.1fontsize, 
        patchsize = (0.6fontsize, 0.6fontsize), 
        # patchstrokecolor, patchstrokewidth, 
        # polycolor, polypoints,  polystrokecolor, 
        # polystrokewidth = 1, 
        rowgap = 0, 
        colgap = 0.2fontsize, 
        # titlecolor, titlefont, titlegap, titlehalign, titleposition, titlesize, titlevalign, titlevisible, 

    ),
    # Scatterlines = (
    #     markersize = 4,
    # )
)

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
        titleformatfunc = (subgrid_idx,tt)-> begin
            rich(
                rich("($(alphabet[subgrid_idx])) ", font=:bold),
                (@sprintf "t = %.10G (s)" tt)
            )
        end,
        sup! = (ax,tgob,subgrid_idx)->nothing,
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
                1 for subgrid_idx in eachindex(subgrids)
            ]
        (a,b::Nothing) => [
                time2step(a[subgrid_idx],traj.t)
                for subgrid_idx in eachindex(subgrids)
            ]
        (a::Nothing,b) => b
        (a,b) => begin
            @warn "Ignoring `attimes`"
            b
        end
    end
    # @warn "Overwriting `atsteps`"
    # @warn "Ignoring `attimes`"
    for subgrid_idx in eachindex(parsed_steps)
        if subgrid_idx > length(subgrids)
            subgrid = subgrids[end]
        else
            subgrid = subgrids[subgrid_idx]
        end
        this_step = parsed_steps[subgrid_idx]
        this_time = Observable(traj.t[this_step])
        RB.goto_step!(bot,this_step;actuate)
        tgob = Observable(deepcopy(structure))
        axtitle = map(this_time) do tt
            titleformatfunc(subgrid_idx,tt)
        end
        if ndim == 2 && !showmesh
            # showinfo = false
            showground = false
            ax = Axis(subgrid[1,1],)
            ax.aspect = DataAspect()
            xlims!(ax,xmin,xmax)
            ylims!(ax,ymin,ymax)
        elseif AxisType <: Axis3
            ax = Axis3(subgrid[1,1],
                # title=axtitle,
                aspect=:data,
                # tellheight=false,
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
                RB.viz!(ax,tgobini;
                    showmesh,
                    showwire,
                    isref=true,
                    showcables=false,
                    showpoints,
                    kargs...
                )
            end
            RB.viz!(ax,tgob;
                showmesh,
                showwire,
                showlabels,
                showpoints,
                showcables,
                kargs...
            )
        end
        sup!(ax,tgob,subgrid_idx)
        if showtitle    
            Label(
                subgrid[1, 1, Top()],
                axtitle,
                padding = (0,0,0,0),
                justification = :right,
                lineheight = 1.0,
                halign = :center,
                valign = :bottom,
            )
            # Label(subgrid[1,1,TopLeft()],"($(alphabet[1]))",font=:bold)
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
                grid2 = subgrid[:,2] = GridLayout(;tellheight=false)
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
                grid3 = subgrid[2,:] = GridLayout()
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
                    # RB.analyse_slack(structure,true)
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

function savefig(fig,figname=nothing;backend=GM)
    if !isa(figname,Nothing)
        if isdefined(Main,:figdir)
            figpath = joinpath(figdir,figname)
        else
            figpath = figname
        end
        if backend == CM
            @info "Saving to $figpath.pdf"
            save(figpath*".pdf",fig;backend,pt_per_unit=1)
        else
            px_per_unit = ustrip(u"inch",tw*u"pt")*600/tw # 300dpi
            @info "Saving to $figpath.png"
            save(figpath*".png",fig;backend,px_per_unit)
        end
    end
    fig
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
            sup! = (ax,tgob,subgrid_idx)-> begin
            # cables
                hidexyz(ax)
                @myshow Sbool[:,subgrid_idx]
                linesegs_cables = @lift begin
                    get_linesegs_cables($tgob;)[Sbool[:,subgrid_idx]]
                end
                linesegments!(ax, 
                    linesegs_cables, 
                    color = :red, 
                    # linewidth = cablewidth
                    )
                rcs_by_cables = @lift begin
                    (;tensioned) = $tgob.connectivity
                    ndim = get_num_of_dims($tgob)
                    T = get_numbertype($tgob)
                    ret = Vector{MVector{ndim,T}}()
                    mapreduce(
                        (scnt)->
                        [(
                            scnt.hen.bodysig.state.loci_states[scnt.hen.pid].+
                            scnt.egg.bodysig.state.loci_states[scnt.egg.pid]
                        )./2],
                        vcat,
                        tensioned.connected
                        ;init=ret
                    )
                end
                # @show rcs_by_cables
                Stext = [
                        @sprintf "%4.2f"  S[i,subgrid_idx] 
                        for i in axes(S,1)
                        if Sbool[i,subgrid_idx]
                    ]
                @myshow Stext
                text!(
                    ax,
                    Stext,
                    position = rcs_by_cables[][Sbool[:,subgrid_idx]],
                    fontsize = fontsize,
                    color = :red,
                    align = (:right, :top),
                    offset = (-fontsize/2, 0)
                )
            end
        )
    end
end

function plot_convergence_order!(ax,dts,err_avg;
    show_orders=true,
    orders = [1,2],
    marker=:rect,
    color=:red,
    label="NMSI"
)

    scatterlines!(ax,dts,err_avg;marker,color,label)

    if show_orders
        for order = orders
            o = err_avg[1] .*(dts./dts[1]).^order
            label = ifelse(order==1,L"\mathcal{O}(h)",latexstring("\\mathcal{O}(h^$order)"))
            lines!(ax,dts,o;label)
        end
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

function plotsave_friction_direction(bots,x,xs,figname=nothing;
        size = (0.9tw,0.4tw),
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

function plotsave_contactpoints(bot,figname=nothing)
    contacts_traj_voa = VectorOfArray(bot.contacts_traj)
    (;t) = bot.traj
    with_theme(theme_pub;
            figure_padding = (0,fontsize,0,0),
            size = (0.9tw,0.3tw),
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

