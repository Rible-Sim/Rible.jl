const alphabet::String = join('a':'z')

function time2step(at, t)
    atstep = findfirst((x) -> x >= at, t)
    @assert !isa(atstep, Nothing)
    atstep
end

function plot_traj!(sim_result::SimulationResult;kwargs...)
    (;simulator,solver_cache) = sim_result
    (;prob,controller,tspan,restart,totalstep,solver_history) = simulator
    (;bot,env, policy) = prob
    (;geometry, field,) = env
    plot_traj!(bot;ground = geometry,policy,kwargs...)
end

"""
    VisConfig

Visualization configuration struct.
"""
@kwdef struct VisConfig
    "Whether to show labels."
    show_labels::Bool = true
    "Whether to show points."
    show_loci::Bool = true
    "Whether to show the mesh."
    show_mesh::Bool = true
    "Whether to show the wireframe."
    show_wire::Bool = false
    "Whether to show the title."
    show_title::Bool = true
    "Whether to show cables."
    show_cables::Bool = true
    "Whether to show the info panel."
    show_info::Bool = true
    "Whether to show axes."
    show_axis::Bool = true
    "Whether to show the ground."
    show_ground::Bool = false
    "Whether to show the background."
    show_background::Bool = false
    "Whether to show the initial state."
    show_init::Bool = false
    "Title formatting function."
    title_format_func::Function = (subgrid_idx,tt)-> begin
        Makie.rich(
            Makie.rich("($(alphabet[subgrid_idx])) ", font=:bold),
            (@sprintf "t = %.10G (s)" tt)
        )
    end
end

"""
    AxisConfig

Axis configuration struct.
"""
@kwdef struct AxisConfig
    "Axis type, e.g. Makie.LScene or Makie.Axis3."
    AxisType::Type = Makie.LScene
    "Figure size."
    figsize::Symbol = :FHD
    "Grid size."
    gridsize::Tuple{Int,Int} = (1,1)
    "X-axis limits."
    xlims = (-1.0,1.0)
    "Y-axis limits."
    ylims = (-1.0,1.0)
    "Z-axis limits."
    zlims = (-1e-3,1.0)
    "Row gap."
    rowgap::Real = 0
    "Column gap."
    colgap::Real = 0
    "Whether to auto-adjust."
    auto::Bool = false
end

"""
    PlayConfig

Playback/recording configuration struct.
"""
@kwdef struct PlayConfig
    "Whether to enable the slider."
    do_slide::Bool = true
    "Whether to record video."
    do_record::Bool = false
    "Speedup factor."
    speedup::Real = 1
    "File name."
    figname::Union{String,Nothing} = nothing
    "Specific times."
    at_times::Union{AbstractVector,Nothing} = nothing
    "Specific steps."
    at_steps::Union{AbstractVector,Nothing} = nothing
end

function Base.show(io::IO, c::VisConfig)
    println(io, "VisConfig(")
    for n in fieldnames(VisConfig)
        println(io, "    $n = $(getfield(c, n)),")
    end
    print(io, ")")
end

function Base.show(io::IO, c::AxisConfig)
    println(io, "AxisConfig(")
    for n in fieldnames(AxisConfig)
        println(io, "    $n = $(getfield(c, n)),")
    end
    print(io, ")")
end

function Base.show(io::IO, c::PlayConfig)
    println(io, "PlayConfig(")
    for n in fieldnames(PlayConfig)
        println(io, "    $n = $(getfield(c, n)),")
    end
    print(io, ")")
end

"""
    setup_axis!(subgrid, config::AxisConfig, ndim::Int, show_mesh::Bool)

Set up the axis.

# Arguments
- `subgrid`: subgrid
- `config`: axis configuration
- `ndim`: dimension
- `show_mesh`: whether to show the mesh
"""
function setup_axis!(subgrid, config::AxisConfig, ndim::Int, show_mesh::Bool)
    (;AxisType, xlims, ylims, zlims, auto) = config
    xmin, xmax = xlims
    ymin, ymax = ylims
    zmin, zmax = zlims

    @debug "Setting up axis" AxisType ndim show_mesh

    if ndim == 2 && !show_mesh
        @debug "Creating 2D Axis"
        # Use Makie.Axis for 2D plots without a mesh.
        ax = Makie.Axis(subgrid[1,1])
        ax.aspect = Makie.DataAspect()
        Makie.xlims!(ax,xmin,xmax)
        Makie.ylims!(ax,ymin,ymax)
        return ax
    elseif AxisType <: Makie.Axis3
        @debug "Creating Axis3"
        # Use Makie.Axis3.
        ax = Makie.Axis3(subgrid[1,1], aspect=:data)
        Makie.xlims!(ax,xmin,xmax)
        Makie.ylims!(ax,ymin,ymax)
        Makie.zlims!(ax,zmin,zmax)
        if auto
            @debug "Auto-limiting Axis3"
            Makie.autolimits!(ax)
        end
        return ax
    elseif AxisType <: Makie.LScene
        @debug "Creating LScene"
        # Use Makie.LScene.
        ax = Makie.LScene(subgrid[1,1], show_axis=true, scenekw = (show_axis = true,))
        if auto
            @debug "Centering LScene"
            Makie.center!(ax.scene)
        end
        return ax
    else
        error("Unknown AxisType: $AxisType")
    end
end

"""
    draw_ground!(ax, ground, config::AxisConfig)

Draw the ground.
"""
function draw_ground!(ax, ground, config::AxisConfig)
    (;xlims, ylims, zlims) = config
    xmin, xmax = xlims
    ymin, ymax = ylims
    zmin, zmax = zlims
    xwid = xmax - xmin
    ywid = ymax - ymin
    zwid = zmax - zmin
    
    # Create the ground rectangle.
    rect = Makie.Rect3f((xmin,ymin,zmin),(xwid,ywid,zwid))
    groundmesh = get_groundmesh(ground,rect)
    Makie.mesh!(ax,groundmesh;color = :snow)
end

"""
Draw the background scene (e.g., starry sky).
$(TYPEDSIGNATURES)
"""
function draw_background!(ax, config::AxisConfig)
    (;figsize) = config
    subWinWidth,subWinHeight = match_figsize(figsize)
    # Create the background scene.
    subWindow0 = Makie.Scene(ax.scene, 
        px_area=Makie.Rect(0, 0, subWinWidth, subWinHeight), 
        clear=true, 
        backgroundcolor=:white
    )
    Makie.campixel!(subWindow0;farclip=1)
    # Load the star background.
    bg = Makie.load(assetpath("stars.jpg")) |> rotr90
    Makie.image!(subWindow0,
        LinRange(0,subWinWidth,size(bg,1)),
        LinRange(0,subWinHeight,size(bg,2)), 
        bg;
        inspectable=false
    )
end

"""
    setup_info_panel!(subgrid, ax, axis_config::AxisConfig)

Set up the info panel.
"""
function setup_info_panel!(subgrid, ax, axis_config::AxisConfig)
    (;AxisType) = axis_config
    
    @debug "Showing info panel"
    # Show the info panel.
    grid2 = subgrid[:,2] = Makie.GridLayout(;tellheight=false)
    grid_info = grid2[1,1] = Makie.GridLayout(;tellheight=false)
    dict_info = [
        "fig. height" => map(string,ax.height),
        "fig. width" => map(string,ax.width)
    ]
    if AxisType == Makie.Axis3
        @debug "Adding camera info for Axis3"
        cam_info = [
            "azimuth" => map(string,ax.azimuth),
            "elevation" => map(string,ax.elevation)
        ]
        append!(dict_info,cam_info)
    end
    for (i,(infoname,infovalue)) in enumerate(dict_info)
        Makie.Label(grid_info[i,1],
            map(infovalue) do iv
                "$infoname = $iv"
            end,
            justification = :left,
        )
    end
    grid_button = grid2[2,1] = Makie.GridLayout(;tellheight=false)
    botton_printinfo = Makie.Button(grid_button[1,1], label = "Print Infos")
    on(botton_printinfo.clicks) do _
        for (infoname,infovalue) in dict_info
            println("$infoname = $(infovalue[])")
        end
    end
end

"""
    setup_interaction!(subgrid, traj, bot, structure, viz_st, ax, this_time, this_step, policy, config::PlayConfig, axis_config::AxisConfig, vis_config::VisConfig)

Set up interaction.
"""
function setup_interaction!(subgrid, traj, bot, structure, viz_st, ax, this_time, this_step, policy, config::PlayConfig, axis_config::AxisConfig, vis_config::VisConfig)
    (;do_slide, do_record, ) = config
    (;AxisType) = axis_config
    (;show_info) = vis_config

    if show_info
        setup_info_panel!(subgrid, ax, axis_config)
    end

    if !do_record
        @debug "Setting up slider (not recording)"
        # Add the slider.
        grid3 = subgrid[2,:] = Makie.GridLayout()
        slidergrid = Makie.SliderGrid(
            grid3[1,1],
            (label = "Step", range = 1:length(traj), startvalue = this_step),
        )
        on(slidergrid.sliders[1].value) do step
            update_robot_vis!(viz_st, structure, bot, step, traj, this_time, policy, axis_config, ax)
        end
    end
end

"""
    record_trajectory!(fig, traj, bot, structure, viz_st, this_time, policy, config::PlayConfig)

Record the trajectory.
"""
function record_trajectory!(fig, traj, bot, structure, viz_st, this_time, policy, config::PlayConfig)
    (;figname, speedup) = config
    dt = traj.t[begin+1] - traj.t[begin]
    framerate = 30 
    skipstep = round(Int,1/framerate/dt*speedup)
    recordsteps = 1:skipstep:length(traj.t)
    
    @info "Recording video" figname framerate speedup
    
    Makie.record(fig, figname, recordsteps; px_per_unit = 2, framerate) do step
        @show step / length(traj.t)
        # We don't need to pass axis_config or ax here as we are not re-centering in record mode usually, 
        # or we can pass them if needed. For now let's assume no auto-centering during record for performance/stability or keep it simple.
        # Actually the original code didn't auto-center in record loop.
        
        # Manually inline update logic or call a simplified update function
        structure.state.system.q .= traj.q[step]
        structure.state.system.c .= traj.c[step]
        update!(structure)
        this_time[] = traj.t[step]
        foreach(structure.bodies) do body
            (;id) = body.prop
            viz_st.bodies[id][] = body
        end
        foreach(structure.apparatuses) do appar
            (;id) = appar
            viz_st.apparatuses[id][] = appar
        end
    end
end

"""
    update_robot_vis!(viz_st, structure, bot, step, traj, this_time, policy, config::AxisConfig, ax)

Update the robot visualization.
"""
function update_robot_vis!(viz_st, structure, bot, step, traj, this_time, policy, config::AxisConfig, ax)
    inst_state = traj[step]
    actuate!(bot,policy,inst_state)
    structure.state.system.q .= traj.q[step]
    structure.state.system.c .= traj.c[step]
    update!(structure, NoField())
    this_time[] = traj.t[step]
    foreach(structure.bodies) do body
        (;id) = body.prop
        viz_st.bodies[id][] = body
    end
    foreach(structure.apparatuses) do appar
        (;id) = appar
        viz_st.apparatuses[id][] = appar
    end
    
    if config.auto
        if config.AxisType <: Makie.LScene
            @debug "Auto-centering LScene in update"
            Makie.center!(ax.scene)
        elseif config.AxisType <: Makie.Axis3
            @debug "Auto-limiting Axis3 in update"
            Makie.autolimits!(ax)
        end
    end
end

function plot_traj!(bot::Robot;
        AxisType=Makie.LScene,
        figsize=:FHD,
        fig = Makie.Figure(),
        gridsize=(1,1),
        at_times=nothing,
        at_steps=nothing,
        do_slide=true,
        do_record=false,
        speedup=1,
        show_labels=true,
        show_loci=true,
        show_mesh=true,
        show_wire=false,
        show_title=true,
        show_cables=true,
        show_info=true,
        xlims=(-1.0,1.0),
        ylims=(-1.0,1.0),
        zlims=(-1e-3,1.0),
        show_axis = true,
        show_ground=false,
        ground=nothing,
        show_background=false,
        figname=nothing,
        rowgap=0,
        colgap=0,
        show_init=false,
        policy= NoPolicy(),
        auto=false,
        title_format_func = (subgrid_idx,tt)-> begin
            Makie.rich(
                Makie.rich("($(alphabet[subgrid_idx])) ", font=:bold),
                (@sprintf "t = %.10G (s)" tt)
            )
        end,
        sup! = (ax,viz_st,subgrid_idx)->nothing,
        kargs...
    )
    
    vis_config = VisConfig(;
        show_labels, show_loci, show_mesh, show_wire, show_title, 
        show_cables, show_info, show_axis, show_ground, show_background, show_init, title_format_func
    )
    
    axis_config = AxisConfig(;
        AxisType, figsize, gridsize, xlims, ylims, zlims, rowgap, colgap, auto
    )
    
    play_config = PlayConfig(;
        do_slide, do_record, speedup, figname, at_times, at_steps
    )

    @debug "Configurations initialized" vis_config axis_config play_config

    (;structure,traj) = bot
    ndim = get_num_of_dims(structure)
    structure.state.system.q .= traj.q[begin]
    update!(structure)
    st_ini = deepcopy(structure)

    # Create subgrids.
    grid1 = fig[1,1] = Makie.GridLayout(;tellheight=false)
    subgrids = [
        grid1[i,j] = Makie.GridLayout(;tellheight=false)
        for i in 1:gridsize[1], j = 1:gridsize[2]
    ] |> permutedims

    @debug "Subgrids created" size(subgrids) length(subgrids)

    parsed_steps = @match (play_config.at_times, play_config.at_steps) begin
        (a::Nothing,b::Nothing) => [1 for _ in eachindex(subgrids)]
        (a,b::Nothing) => [time2step(a[i],traj.t) for i in eachindex(subgrids)]
        (a::Nothing,b) => b
        (a,b) => begin
            @warn "Ignoring `at_times` because `at_steps` is provided. (忽略 `at_times`，因为已提供 `at_steps`)"
            b
        end
    end

    @debug "Steps parsed" parsed_steps

    for subgrid_idx in eachindex(parsed_steps)
        @debug "Processing subgrid" subgrid_idx total=length(parsed_steps)

        subgrid = subgrid_idx > length(subgrids) ? subgrids[end] : subgrids[subgrid_idx]
        this_step = parsed_steps[subgrid_idx]
        this_time = Observable(traj.t[this_step])
        
        goto_step!(bot,this_step;policy,)
        viz_st = VizStructure(structure)
        
        axtitle = map(this_time) do tt
            vis_config.title_format_func(subgrid_idx,tt)
        end

        # Set up the axis.
        ax = setup_axis!(subgrid, axis_config, ndim, vis_config.show_mesh)

        if vis_config.show_ground
            @debug "Drawing ground"
            draw_ground!(ax, ground, axis_config)
        end

        if vis_config.show_wire || vis_config.show_mesh || vis_config.show_cables || vis_config.show_labels || vis_config.show_loci
            @debug "Visualizing robot components"
            if vis_config.show_init
                @debug "Showing initial state"
                # Show the initial state.
                vis!(ax,st_ini;
                    show_mesh=vis_config.show_mesh,
                    show_wire=vis_config.show_wire,
                    isref=true,
                    show_cables=false,
                    show_loci=vis_config.show_loci,
                    kargs...
                )
            end
            # Show the current state.
            vis!(ax,viz_st;
                show_mesh=vis_config.show_mesh,
                show_wire=vis_config.show_wire,
                show_labels=vis_config.show_labels,
                show_loci=vis_config.show_loci,
                show_cables=vis_config.show_cables,
                kargs...
            )
        end

        sup!(ax,viz_st,subgrid_idx)

        if vis_config.show_title    
            @debug "Showing title"
            Makie.Label(
                subgrid[1, 1, Makie.Top()],
                axtitle,
                padding = (0,0,0,0),
                justification = :right,
                lineheight = 1.0,
                halign = :center,
                valign = :bottom,
            )
        end

        if vis_config.show_background
            @debug "Drawing background"
            draw_background!(ax, axis_config)
        end

        if play_config.do_slide || play_config.do_record
            @debug "Setting up interaction/recording"
            setup_interaction!(subgrid, traj, bot, structure, viz_st, ax, this_time, this_step, policy, play_config, axis_config, vis_config)
            
            if play_config.do_record
                 @info "Starting recording... (开始录制...)"
                 record_trajectory!(fig, traj, bot, structure, viz_st, this_time, policy, play_config)
                 @info "Recording finished. (录制完成)"
            end
        end
        
        Makie.colgap!(grid1, axis_config.colgap)
        Makie.rowgap!(grid1, axis_config.rowgap)
    end

    if (fig isa Makie.Figure) && !play_config.do_record
        @debug "Saving figure" play_config.figname
        savefig(fig, play_config.figname)
        # Makie.DataInspector(fig)
    end
    fig
end

function savefig(fig,figname=nothing;backend=Makie.current_backend())
    if !isa(figname,Nothing)
        if isdefined(Main,:figdir)
            figpath = joinpath(Main.figdir,figname)
        else
            figpath = figname
        end
        if nameof(backend) == :CairoMakie
            @info "Saving to $figpath.pdf"
            Makie.save(figpath*".pdf",fig;backend,pt_per_unit=1)
        else
            px_per_unit = 600 / 72 # 300dpi, approximately 8.33
            @info "Saving to $figpath.png"
            Makie.save(figpath*".png",fig;backend,px_per_unit)
        end
    end
    fig
end
