
function plot_tower2d!(ax,
                    st::RB.Structure;
                    cablecolor = :dodgerblue,
                    barcolor = :darkred,
                    markercolor = :darkblue,
                    refcolor = :darkgrey,
                    isref = false,
                    markit = false)
    if isref
        cablecolor = refcolor
        barcolor = refcolor
        markercolor = refcolor
    end
    (;num_of_bodies) = st
    (;tensioned) = st.connectivity
    (;cables) = st.apparatuses
    ncables = length(cables)
    ndim = RB.get_num_of_dims(st)
    T = RB.get_numbertype(st)

    rbs = RB.get_bodies(st)
    linesegs_bars = Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
    linesegs_tris =  Vector{Tuple{Point{ndim,T},Point{ndim,T}}}()
    ploys_tris = Vector{Vector{Point{ndim,T}}}()
    for rb in rbs
        if body.state.cache.funcs.nmcs isa RB.NCF.NC2D2P
            push!(
                linesegs_bars,
                (Point(body.state.loci_states[1]),Point(body.state.loci_states[2])),
            )
        elseif body.state.cache.funcs.nmcs isa RB.NCF.NCMP
        else
            append!(
                linesegs_tris,
                [
                    (Point(body.state.loci_states[1]), Point(body.state.loci_states[2])),
                    (Point(body.state.loci_states[2]), Point(body.state.loci_states[3])),
                    (Point(body.state.loci_states[3]), Point(body.state.loci_states[1])),
                ]
            )
            push!(
                ploys_tris,
                [Point(body.state.loci_states[i]) for i = 1:3]
            )
        end
    end
    linesegs_rbs = vcat(linesegs_bars,linesegs_tris)
    linesegments!(ax, linesegs_bars, color = barcolor, linewidth = barwidth)
    poly!(ax,ploys_tris, color = barcolor, transparency = false)
    if markit
        scatter!(ax,[rbs[3].state.rg], color = markercolor, markersize = (2 |> pt2px))
    end

    linesegs_noslack_cables = get_linesegs_cables(st;noslackonly=true)
    linesegs_slack_cables = get_linesegs_cables(st;slackonly=true)
    linesegments!(ax, linesegs_noslack_cables, color = cablecolor, linewidth = cablewidth)
    linesegments!(ax, linesegs_slack_cables, color = cablecolor, linewidth = cablewidth, linestyle = :dash)

end

function plot_one_bar_one_tri!(ax,
        tgob,
        sgi;
        cablecolor = :dodgerblue,
        barcolor = :darkred,
        markercolor = :darkblue,
        refcolor = :darkgrey,
        isref = false,
        markit = false)
    if isref
        cablecolor = refcolor
        barcolor = refcolor
        markercolor = refcolor
    end
    st = tgob[]
    (;num_of_bodies) = st
    num_of_dim = RB.get_num_of_dims(st)
    T = RB.get_numbertype(st)

    linesegs_bars = @lift begin
        num_of_dim = RB.get_num_of_dims($tgob)
        T = RB.get_numbertype($tgob)
        ret = Vector{Pair{Point{num_of_dim,T},Point{num_of_dim,T}}}()
        foreach($tgob.bodies) do body
            if body.state.cache.funcs.nmcs isa RB.NCF.NC2D2P
                push!(ret,
                    Point(body.state.loci_states[1])=>
                    Point(body.state.loci_states[2])
                )
            end
        end
        ret
    end

    ploys_tris = @lift begin
        ndim = RB.get_num_of_dims($tgob)
        T = RB.get_numbertype($tgob)
        ret = Vector{Vector{Point{ndim,T}}}()
        foreach($tgob.bodies) do body
            if body.state.cache.funcs.nmcs isa RB.NCF.NC2D2P
                nothing
            else
                push!(ret,
                    [Point(body.state.loci_states[i]) for i = 1:3]
                )
            end
        end
        ret
    end

    # linesegs_rbs = vcat(linesegs_bars,linesegs_tris)
    linesegments!(ax, linesegs_bars, color = barcolor, linewidth = barwidth)
    poly!(ax,ploys_tris, color = barcolor, transparency = false)
    # if markit
    # scatter!(ax,[rbs[3].state.rg], color = markercolor, markersize = (2 |> pt2px))
    # end

    # linesegs_noslack_cables = get_linesegs_cables(st;noslackonly=true)
    # linesegs_slack_cables = get_linesegs_cables(st;slackonly=true)
    # linesegments!(ax, linesegs_noslack_cables, color = cablecolor, linewidth = cablewidth)
    # linesegments!(ax, linesegs_slack_cables, color = cablecolor, linewidth = cablewidth, linestyle = :dash)
end