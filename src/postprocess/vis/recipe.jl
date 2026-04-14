
Makie.@recipe Vis  (object,) begin
    id = 1
    fontsize = 12
    linewidth = 1
    isref = false
    show_labels = false
    label_prefix = ""
    show_loci = false
    point_color = :black
    show_mass_centers = false
    show_mass_center_labels = false
    loci_idx = nothing
    loci_colors = :black
    loci_markers_sizes = 12
    show_loci_labels = false
    show_axes = false
    show_curvature = false
    axes_arrows_scale = 1.0
    show_mesh = true
    show_wire = false
    rigid_label_color = :darkblue
    mesh_color = nothing
    ref_color = :lightgrey
    show_cables = true
    show_cable_labels = false
    cable_color = :deepskyblue
    cable_label_color = :darkgreen
    slack_linestyle = :dash
    cable_width = 0.6
	clip_planes = Makie.Plane3f[]
end

function Makie.args_preferred_axis(::Type{<: Vis}, object)
    N = get_num_of_dims(object)
    return N == 3 ? Makie.LScene : Makie.Axis
end

# Type recipes
Makie.convert_arguments(::Makie.PointBased, x::VectorOfArray) = (Matrix(x),)

function scale_axes_arrows_callback(vis)
    on(Makie.events(vis).keyboardbutton) do event
        if event.action == Makie.Keyboard.press || event.action == Makie.Keyboard.repeat
            if event.key == Makie.Keyboard.period
                vis.axes_arrows_scale[] *= 1.1
            elseif event.key == Makie.Keyboard.comma
                vis.axes_arrows_scale[] /= 1.1
            end
        end
    end
end

function parse_loci_body(num_of_loci, loci_idx_input, loci_colors_input, loci_markers_sizes_input) 
    if loci_idx_input isa Nothing
        loci_idx = collect(1:num_of_loci)
    else
        loci_idx = loci_idx_input
    end
    if loci_colors_input isa Symbol
        loci_colors = fill(loci_colors_input,num_of_loci)
    else
        loci_colors = loci_colors_input
    end
    if loci_markers_sizes_input isa Number
        loci_markers_sizes = fill(loci_markers_sizes_input,num_of_loci)
    else
        loci_markers_sizes = loci_markers_sizes_input
    end
    @eponymtuple(loci_idx, loci_colors, loci_markers_sizes)
end

function parse_nodes(bodyid, num_of_loci, loci_idx_input, loci_colors_input, loci_markers_sizes_input) 
    if loci_idx_input isa Nothing
        loci_idx = collect(1:num_of_loci)
    else
        loci_idx = loci_idx_input[bodyid]
    end
    if loci_colors_input isa Symbol
        loci_colors = fill(loci_colors_input,num_of_loci)
    elseif loci_colors_input isa Makie.Colorant
        loci_colors = fill(loci_colors_input,num_of_loci)
    else
        loci_colors = loci_colors_input[bodyid]
    end
    if loci_markers_sizes_input isa Number
        loci_markers_sizes = fill(loci_markers_sizes_input,num_of_loci)
    else
        loci_markers_sizes = loci_markers_sizes_input[bodyid]
    end
    @eponymtuple(loci_idx, loci_colors, loci_markers_sizes)
end

# Apparatus Empty
function Makie.data_limits(vis::Vis{Tuple{S}};) where {S<:Apparatus}
    N = get_num_of_dims(vis[:object][])
    if N == 3
        ret = Makie.Rect3d((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    else
        ret = Makie.Rect2d((0.0, 0.0), (1.0, 1.0))
    end
    ret
end

# # Apparatus
function Makie.plot!(vis::Vis{Tuple{S}};) where {S<:Apparatus}
    Makie.scatter!(vis, Int[])
    vis
end

# LocusState
function Makie.plot!(vis::Vis{Tuple{S}};) where S <: LocusState{N} where N

    map!(vis.attributes, [:object, ], [:locus_position, :locus_curvature]) do args...
        locus_position = [args[1].frame.position |> GB.Point]
        locus_curvature = [args[1].curvature |> GB.Vec]
        (locus_position, locus_curvature)
    end
   
    
    #note avoid MethodError: no method matching draw_atomic
    Makie.scatter!(vis, vis.locus_position; 
        color = vis.loci_colors, 
        markersize = vis.loci_markers_sizes,
        visible = vis.show_loci, 
    )

    Makie.text!(vis,
        vis.label_prefix[]*"p$(vis.id[])" ,
        position = vis.locus_position,
        color = vis.rigid_label_color,
        align = (:left, :top),
        offset = (-5, -10),
        visible = vis.show_labels
    )

    map!(vis.attributes, [:object, :axes_arrows_scale, :show_loci, :show_axes], 
        [:axes_positions, :axes_arrows, :axes_lengthscale, :axes_visible]) do args...
        locus_state = args[1]
        axes_arrows_scale = args[2]
        show_loci = args[3]
        show_axes = args[4]
        axes_positions = [
            locus_state.frame.position |> GB.Point
            for _ = 1:N
        ]

        axes_arrows = [
            locus_state.frame.axes.X[:,i] |> GB.Vec
            for i = 1:N
        ]

        # axes_texts = N == 2 ? ["t","n"] : ["t","b","n"]

        # arrows_ends_ob = lift(locus_state_ob) do locus_state_ob
        #     [
        #         (locus_state_ob.frame.position + 1.6 .*locus_state_ob.frame.axes.X[:,i]) |> GB.Point
        #         for i = 1:N
        #     ]
        # end

        axes_lengthscale = axes_arrows_scale
        axes_visible = show_loci && show_axes
        (axes_positions, axes_arrows, axes_lengthscale, axes_visible)
    end

    #the first two components scale the radius (in x/y direction) 
    #the last scales the length of the cone
    if N == 3
        Makie.arrows3d!(vis, vis.axes_positions, vis.axes_arrows; 
            color = [:red,:green,:blue], 
            lengthscale = vis.axes_lengthscale, # master‐length of each arrow
            visible = vis.axes_visible,
            # normalize    = true,    # normalize directions → unit vectors
            align        = :tail,   # base of arrow sits at O
            markerscale = 1,
            # taillength   = 0,     # no “back‐facing” tail
            # shaftlength  = 0.7,       # 75% of s is shaft
            tiplength    = 0.25,       # 25% is tip
            shaftradius  = 0.02,       # thickness of shaft
            tipradius    = 0.06,       # ﬂared cone tip
        )
         # Makie.text!(vis,
        #     texts,
        #     position = arrows_ends_ob,
        #     color = colors,
        #     align = (:left, :center),
        #     # offset = (0,2viz.fontsize[]*(rand()-0.5))
        # )
        Makie.arrows3d!(vis, 
            vis.locus_position, 
            vis.locus_curvature;
            lengthscale = vis.axes_lengthscale,
            visible = vis.show_curvature
        ) 
    else
        # Skip arrows for 2D structures.
    end
    
    vis

end

# AbstractRigidBody
function Makie.plot!(vis::Vis{Tuple{S}};) where S <:AbstractRigidBody
    body_ob = vis[:object]
    # body decorations
    #-- mass center
    map!(vis.attributes, [:object,], [:mass_center, :mass_center_position]) do args...
        mass_center = args[1].state.mass_locus_state
        mass_center_position = mass_center.frame.position |> GB.Point
        (mass_center, mass_center_position)
    end

    #note avoid MethodError: no method matching draw_atomic
    Makie.scatter!(vis, vis.mass_center_position; markersize = 1e-6, color=color = (:red, 0.0))

    (;show_axes, show_loci_labels, show_mass_centers, show_mass_center_labels, axes_arrows_scale) = vis

    vis!(vis,vis.mass_center; 
        label_prefix = "r$(body_ob[].prop.id)", id="g", 
        show_axes,
        show_labels = show_mass_center_labels,
        show_loci = show_mass_centers,
        axes_arrows_scale,
        loci_colors = :black
    )

    #-- nodes
    loci_idx, loci_colors, loci_markers_sizes = parse_loci_body(
        length(body_ob[].prop.loci),
        vis.loci_idx[], vis.loci_colors[], vis.loci_markers_sizes[]
    )
    
    # loci_graph = Makie.ComputeGraph()
    # Makie.add_input!(loci_graph, :body, vis.attributes[:object])

    for (j,i) in enumerate(loci_idx)
        locus_state_sym = Symbol("locus_state_$i")
        map!(vis.attributes, [:object, ], 
                locus_state_sym
            ) do args...
            # @show "doing $locus_state_sym computations"
            args[1].state.loci_states[i]
        end
        vis!(vis, vis.attributes[locus_state_sym]; 
            label_prefix = "r$(body_ob[].prop.id)", 
            id = i, 
            show_axes,
            show_curvature = vis.show_curvature,
            show_labels = show_loci_labels,
            axes_arrows_scale,
            show_loci = vis.show_loci,
            loci_colors = loci_colors[j],   
            loci_markers_sizes = loci_markers_sizes[j]
        )
    end

    #-- body mesh
    map!(vis.attributes, [:object, :mesh_color,:show_mesh,:show_wire,:linewidth], 
        [:body_mesh, :body_mesh_visible, :body_mesh_strokewidth]) do args...
        body_mesh = build_mesh(args[1],color=args[2])
        body_mesh_visible = args[3] || args[4]
        body_mesh_strokewidth = args[4] ? args[5] : 0
        (body_mesh, body_mesh_visible, body_mesh_strokewidth)
    end

    Makie.poly!(vis, vis.body_mesh; shading = Makie.Automatic(), 
        visible = vis.body_mesh_visible,
        strokewidth = vis.body_mesh_strokewidth
    )
    vis
end


# VizStructure
function Makie.plot!(vis::Vis{Tuple{S}}; ) where S <:VizStructure
    mesh_color = vis.mesh_color
    cable_color = vis.cable_color
    cable_label_color = vis.cable_label_color
    if vis.show_labels[]
        show_cable_labels = 
        show_mass_center_labels = 
        show_loci_labels = true
    else
        show_cable_labels = vis.show_cable_labels
        show_mass_center_labels = vis.show_mass_center_labels
        show_loci_labels = vis.show_loci_labels
    end
    if vis.show_loci[]
        show_mass_centers = true
    else
        show_mass_centers = vis.show_mass_centers
    end

    object_ob = vis[:object]

    (;connectivity) = object_ob[]
    (;num_of_bodies,num_of_apparatuses) = connectivity

    for appar_ob in object_ob[].apparatuses
        vis!(vis,appar_ob;
            cable_color,
            cable_width = vis.cable_width,
            cable_label_color = vis.cable_label_color,
            slack_linestyle = vis.slack_linestyle,
            show_cable_labels,
        )
    end
    
    for body_ob in object_ob[].bodies
        loci_idx, loci_colors, loci_markers_sizes = parse_nodes(
            body_ob[].prop.id, length(body_ob[].prop.loci),
            vis.loci_idx[], vis.loci_colors[], vis.loci_markers_sizes[]
        )
        # @show loci_idx
        vis!(vis,body_ob;
            show_mass_centers,
            show_mass_center_labels,
            show_loci = vis.show_loci,
            show_loci_labels,
            mesh_color,
            show_axes = vis.show_axes,
            show_curvature = vis.show_curvature,
            axes_arrows_scale = vis.axes_arrows_scale,
            point_color = vis.point_color,
            show_mesh = vis.show_mesh,
            show_wire = vis.show_wire,
            loci_idx, 
            loci_colors,
            loci_markers_sizes,
        )
    end
    # Makie.DataInspector(vis)
    scale_axes_arrows_callback(vis)
    vis
end

# Structure
function Makie.plot!(vis::Vis{Tuple{S}},; kwargs...) where S <:Structure

    on(Makie.events(vis).keyboardbutton) do event
        if event.action == Makie.Keyboard.press || event.action == Makie.Keyboard.repeat
            if event.key == Makie.Keyboard.period
                vis.axes_arrows_scale[] *= 1.1
            elseif event.key == Makie.Keyboard.comma
                vis.axes_arrows_scale[] /= 1.1
            end
        end
    end

    viz_st = VizStructure(vis.object[])
    vis!(vis,viz_st; 
        show_axes = vis.show_axes,
        axes_arrows_scale = vis.axes_arrows_scale,
        show_mass_centers = vis.show_mass_centers,
        show_mass_center_labels = vis.show_mass_center_labels,
        show_loci_labels = vis.show_loci_labels,
        show_mesh = vis.show_mesh,
        show_wire = vis.show_wire,
        loci_idx = vis.loci_idx,
        loci_colors = vis.loci_colors,
        loci_markers_sizes = vis.loci_markers_sizes,
        mesh_color = vis.mesh_color,
        cable_color = vis.cable_color,
        cable_width = vis.cable_width,
        cable_label_color = vis.cable_label_color,
        slack_linestyle = vis.slack_linestyle,
        show_cable_labels = vis.show_cable_labels,
        show_labels = vis.show_labels,
        show_loci = vis.show_loci,
        point_color = vis.point_color,
    )
    vis
end

function Makie.plot!(vis::Vis{Tuple{S}};) where S <:Robot
    @error "Vis recipe for Robot is not implemented yet. Use Structure instead."
end
