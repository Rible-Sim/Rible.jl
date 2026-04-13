
# CableJoint
function Makie.plot!(vis::Vis{Tuple{S}};) where {S<:Apparatus{<:CableJoint}}

    map!(vis.attributes, [:object, :slack_linestyle], 
        [:cable_segs, :cable_linestyle, :cable_mid_point, :cable_id]) do args...
        cab = args[1]
        slack_linestyle = args[2]
        cable_id = "c$(cab.id)"
        point_start = cab.force.state.start
        point_stop = cab.force.state.stop 
        cable_segs = [(GB.Point(point_start), GB.Point(point_stop))]
        cable_mid_point = (point_start .+ point_stop)./2 |> GB.Point
        
        cable_linestyle = if cab.force.state.length > cab.force.state.restlen
            Makie.Linestyle([0.0,1e7])
        else
            slack_linestyle == :solid ? Makie.Linestyle([0.0, 1e7]) : Makie.Linestyle([0.0, 5.0, 10.0])
        end
        
        (cable_segs, cable_linestyle, cable_mid_point, cable_id)
    end

    # slackonly=false,
    # noslackonly=true

    Makie.linesegments!(
        vis, vis.cable_segs, 
        color = vis.cable_color, 
        linewidth = vis.cable_width, 
        linestyle = vis.cable_linestyle
    )
  
    # show cable labels
    Makie.text!(vis,
        vis.cable_id,
        position = vis.cable_mid_point,
        color = vis.cable_label_color,
        align = (:left, :top),
        offset = (-5, -10),
        visible = vis.show_cable_labels
    )
    vis
end


# ClusterJoint
function Makie.plot!(vis::Vis{Tuple{S}};) where {S<:Apparatus{<:ClusterJoint}}
    cluster_cable_ob = vis[:object]
    n_segs = length(cluster_cable_ob[].force.segments)
    map!(vis.attributes, [:object,], [:cluster_cable_segs,]) do args...
        cluster_cable = args[1]
        cluster_cable_segs = [
            begin
                iseg = cluster_cable.force.segments[i].force
                point_start = iseg.state.start |> GB.Point
                point_stop = iseg.state.stop |> GB.Point
                (point_start, point_stop)
            end
            for i = 1:n_segs
        ]
        (cluster_cable_segs,)
    end
    Makie.linesegments!(
        vis, vis.cluster_cable_segs,
        color=:red,
        linewidth=0.5,
        linestyle=:solid
    )
end
