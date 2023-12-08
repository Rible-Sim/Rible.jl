function plot_decompose_tower3d(bot0,bot1;
                            cablecolor = :dodgerblue,
                            barcolor = :darkred,)
    fig_width = 1columnwidth
    fig_height = 0.8columnwidth
    fig = Figure(;resolution=(fig_width,fig_height),figure_padding=(-50,0,0,0))
    rows = (1,1,2,2)
    cols = (1,2,1,2)
    axs = [
        Axis3(fig[row,col],aspect=:data,
                # xlabel=L"x~(\mathrm{m})",
                # ylabel=L"y~(\mathrm{m})",
                # zlabel=L"z~(\mathrm{m})"
            )
        for (row,col) in zip(rows,cols)
    ]
    linesegs_noslack_cables, linesegs_slack_cables = get_linesegs_cables(bot0.st)
    linesegs_cables = vcat(linesegs_noslack_cables, linesegs_slack_cables)
    linesegments!(axs[1], linesegs_cables[10:15], color = cablecolor, linewidth = cablewidth)
    linesegments!(axs[2], linesegs_cables[16:18], color = cablecolor, linewidth = cablewidth)
    linesegments!(axs[3], linesegs_cables[ 1: 9], color = cablecolor, linewidth = cablewidth)

    linesegs_bars = get_linesegs_bars(bot0.st)
    linesegments!(axs[1], linesegs_bars[ 6+1: 6+12], color = barcolor, linewidth = barwidth)
    linesegments!(axs[2], linesegs_bars[12+1:12+12], color = barcolor, linewidth = barwidth)
    linesegments!(axs[3], linesegs_bars[   1:    9], color = barcolor, linewidth = barwidth)

    linesegments!(axs[4], reduce(vcat,get_linesegs_cables(bot1.st))[ 1: 9], color = cablecolor, linewidth = cablewidth)
    linesegments!(axs[4],              get_linesegs_bars(bot1.st)[   1: 9], color = barcolor, linewidth = barwidth)

    # zlims!(ax,0,0.5)
    # ax.xlabel = "x (m)"
    # ax.ylabel = "y (m)"
    # ax.zlabel = "z (m)"
    # ax.yticks = ax.xticks = [-0.1,0,0.1]
    # ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
    # ax.ylabeloffset = ax.xlabeloffset = 40
    # ax.alignmode =  Mixed(;left = 25, right = -25)

    for (iax,ax) in enumerate(axs)
        if iax in [1,2,3,4]
            hidedecorations!(ax)
            hidespines!(ax)
        end
        if iax in [6]
            ax.zticklabelsvisible = false
            ax.zticksvisible = false
            ax.zlabelvisible = false
        end
        xlims!(ax,-0.125,0.125)
        ylims!(ax,-0.125,0.125)
        if iax in [3,4]
            zlims!(ax,0,0.3)
        end
        ax.azimuth = 7.595530633326982
        ax.elevation = 0.12269908169872408
    end
    leftpads = [150,150,100,100]
    toppads = [20,20,-0,-0]
    subtitles = [
        "Class-1 Module",
        "Class-2 Module",
        "2-Stage Triplex (Folded)",
        "2-Stage Triplex (Unfolded)"
    ]
    for (ilabel,label) in enumerate(alphabet[1:4])
        Label(fig.layout[rows[ilabel], cols[ilabel], Bottom()],
            ("($label) "),
            font = "CMU Serif Bold",
            textsize = fontsize,
            padding = (leftpads[ilabel]-1.5fontsize, 0, 0, toppads[ilabel]),
            halign = :left,
            valign = :top
        )
    end
    for (ilabel,label) in enumerate(alphabet[1:4])
        Label(fig.layout[rows[ilabel], cols[ilabel], Bottom()],
            ("$(subtitles[ilabel])"),
            textsize = fontsize,
            padding = (leftpads[ilabel], 0, 0, toppads[ilabel]),
            halign = :left,
            valign = :top
        )
    end
    rowsize!(fig.layout, 2, Relative(0.64))
    colgr̄p!(fig.layout,0)
    rowgr̄p!(fig.layout,0)

    fig
end

function plot_compose_tower3d(bot0,bot1;
                            cablecolor = :dodgerblue,
                            barcolor = :darkred,)
    fig_width = 0.7testwidth
    fig_height = 0.6testwidth
    fig = Figure(;resolution=(fig_width,fig_height),figure_padding=(fontsize,0,0,-3fontsize))
    axs = [
        Axis3(fig[1,i],aspect=:data,
                # xlabel=L"x~(\mathrm{m})",
                # ylabel=L"y~(\mathrm{m})",
                # zlabel=L"z~(\mathrm{m})"
            )
        for i = 1:2
    ]
    plot_tower3d!(axs[1],bot0.st)
    plot_tower3d!(axs[2],bot1.st)
    # zlims!(ax,0,0.5)
    # ax.xlabel = "x (m)"
    # ax.ylabel = "y (m)"
    # ax.zlabel = "z (m)"
    # ax.yticks = ax.xticks = [-0.1,0,0.1]
    # ax.yticklabelsize = ax.xticklabelsize = ax.zticklabelsize = 27
    # ax.ylabeloffset = ax.xlabeloffset = 40
    # ax.alignmode =  Mixed(;left = 25, right = -25)

    for (iax,ax) in enumerate(axs)
        # if iax in [2]
        # 	ax.zticklabelsvisible = false
        # 	ax.zticksvisible = false
        # 	ax.zlabelvisible = false
        # end
        xlims!(ax,-0.125,0.125)
        ylims!(ax,-0.125,0.125)
        ax.azimuth = 7.595530633326982
        ax.elevation = 0.12269908169872408
    end
    leftpad = 50
    leftpadplus = -25
    toppad = 0
    subtitles = [
        "3D Tower (Initial Configuration)",
        "3D Tower (Target Configuration)",
    ]

    # ndim = RB.get_num_of_dims(bot0.st)
    # T = RB.get_numbertype(bot0.st)
    # bot0_rcs_by_cables = Vector{MVector{ndim,T}}()
    # foreach(bot0.st.connectivity.tensioned) do scnt
    # 	push!(bot0_rcs_by_cables,
    # 		(
    # 			scnt.hen.bodysig.state.loci_states[scnt.hen.pid].+
    # 			scnt.egg.bodysig.state.loci_states[scnt.egg.pid]
    # 		)./2
    # 	)
    # end
    # text!(axs[1],
    # 	["l$(i)" for (i,rc) in enumerate(bot0_rcs_by_cables)];
    # 	position = bot0_rcs_by_cables,
    # 	textsize = fontsize,
    # 	color = cablecolor,
    # 	align = (:left, :center),
    # 	# offset = (-5, -10)
    # )
    for (ilabel,label) in enumerate(alphabet[1:2])
        Label(fig.layout[1, ilabel, Bottom()],
            ("($label) "),
            font = "CMU Serif Bold",
            textsize = fontsize,
            padding = (leftpad-1.8fontsize, 0, 0, toppad),
            halign = :left,
            valign = :bottom
        )
    end
    for (ilabel,label) in enumerate(alphabet[1:2])
        Label(fig.layout[1, ilabel, Bottom()],
            ("$(subtitles[ilabel])"),
            textsize = fontsize,
            padding = (leftpad, 0, 0, toppad),
            halign = :left,
            valign = :bottom
        )
    end
    colgr̄p!(fig.layout,2fontsize)
    fig
end










