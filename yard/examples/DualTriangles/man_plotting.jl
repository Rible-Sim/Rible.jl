function get_bars(body)
    # @show body.state.loci_states[2]
    [
        Point(body.state.loci_states[1]) => Point(body.state.loci_states[3]),
        Point(body.state.loci_states[1]) => Point(body.state.loci_states[4]),
        Point(body.state.loci_states[3]) => Point(body.state.loci_states[2]),
        Point(body.state.loci_states[4]) => Point(body.state.loci_states[2])
    ]

end

function get_lower_bars(bodyid,rb)
    return [
        Point(body.state.loci_states[1]) => Point(body.state.loci_states[4]),
        Point(body.state.loci_states[4]) => Point(body.state.loci_states[2])
    ]
end

function get_upper_bars(bodyid,rb)
    return [
        Point(body.state.loci_states[1]) => Point(body.state.loci_states[3]),
        Point(body.state.loci_states[3]) => Point(body.state.loci_states[2])
    ]
end

function get_lower_bars(st)
    bars = [Observable(get_lower_bars(bodyid,rb)
            ) for (bodyid,rb) in enumerate(st.rigidbodies)]
end

function get_upper_bars(st)
    bars = [Observable(get_upper_bars(bodyid,rb)
            ) for (bodyid,rb) in enumerate(st.rigidbodies)]
end

function get_bars_and_strings(st)
    bars = [Observable(get_bars(body)
            ) for rb in st.rigidbodies]
    strings = Observable(RB.get_strings(st))
    bars, strings
end

function plotstructure!(ax,tgstruct)
    bars, strings = get_bars_and_strings(tgstruct)
    function plot!(ax,bars,strings)
        for (lineid,line) in enumerate(bars)
            linesegments!(ax, line, color = :black, linewidth = 4)
        end
        linesegments!(ax, strings, color = :deepskyblue, linewidth = 2)
    end
    plot!(ax,bars,strings)
    xlims!(ax,-0.1,0.7)
    ylims!(ax,-0.6,0.2)
    bars,strings
end

function update_scene!(st,bars,strings,q)
    cnt = st.connectivity
    RB.distribute_q_to_rbs!(st,q)
    RB.update_strings_apply_forces!(st)
    for (id,rb) in enumerate(st.rigidbodies)
        bars[id][] = get_bars(body)
    end
    strings[] = RB.get_strings(st)
    # angles = update_angles(st)
    # @show angles
end

function sliderplot(fig,ax,st,bars,strings,sol;kwargs...)
    ls_step = labelslider!(fig,"step",1:length(sol.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(st,bars,strings,sol.qs[this_step])
        # @show tstep
    end
end
function sequenceplot(fig,grid,ax,st,bars,strings,update_scene!,seq,ulabels,glabels)
    nseq = length(seq)
    ls_step = labelslider!(fig,"step",1:nseq)
    grid[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(st,bars,strings,seq.q[this_step])
        RB.set_restlen!(st,seq.u[this_step])
        analyse_slackness(st,seq.q[this_step])
        ls_step.valuelabel.text[] = string((this_step-1)/(nseq-1))
        # for (id,dlabel) in enumerate(dlabels)
        #     dlabel.text[] = "d$id=$(seq.d[this_step][id])"
        # end
        for (iu,ulabel) in enumerate(ulabels)
            ulabel.text[] = "u$iu=$(seq.u[this_step][iu])"
        end
        for (ig,glabel) in enumerate(glabels)
            glabel.text[] = "g$ig=$(seq.g[this_step][ig])"
        end

    end
end
function recordplot(scene,bars,strings,update_scene!,sol;filename="pid.mp4")
    record(scene, filename, 1:60:length(sol.ts); framerate = 30) do this_step
        update_scene!(bars,strings,sol.qs[this_step])
    end
end

function plotstructure(bot::RB.Robot)
    @unpack st, traj = bot
    fig = Figure(resolution=(1000,1000))
    ax = fig[1, 1] = Axis(fig)
    bars,strings = plotstructure!(ax,st)
    sliderplot(fig,ax,st,bars,strings,traj)
    fig
end

function plotstructure(bot::RB.Robot,seq)
    @unpack st = bot
    fig = Figure(resolution=(1280,720))
    grid_tg = fig[1, 1] = GridLayout()
    ax = grid_tg[1, 1] = Axis(fig)
    grid_para = fig[:, 2] = GridLayout(tellheight = false) # Show parameters
    ax_d = grid_para[1,1]
    ax_u = grid_para[1,2]
    ax_g = grid_para[1,3]
    # dlabels = [Label(ax_d[i, 1], "d$i=$(seq[1].d[i])", textsize = 30) for i in 1:length(seq[1].d)]
    ulabels = [Label(ax_u[i, 1], "u$i=$(seq[1].u[i])", textsize = 30) for i in 1:length(seq[1].u)]
    glabels = [Label(ax_g[i, 1], "g$i=$(seq[1].g[i])", textsize = 30) for i in 1:length(seq[1].g)]

    bars,strings = plotstructure!(ax,st)
    sequenceplot(fig,grid_tg,ax,st,bars,strings,update_scene!,seq,ulabels,glabels)
    colsize!(fig.layout, 1, Relative(1/2))
    fig
end

function bars_and_strings_segs(tgstruct;ref=false)
    bars,strings = bars_and_strings(tgstruct)

    bars_lines = Vector{Vector{Tuple{Float64,Float64}}}()
    for b in bars
        for l in b[]
            lines = [l.first.data,l.second.data]
            push!(bars_lines,lines)
        end
    end
    strings_lines = Vector{Vector{Tuple{Float64,Float64}}}()
    for s in strings[]
        lines = [s.first.data,s.second.data]
        push!(strings_lines,lines)
    end
    if ref
        bar_color_name = "orange"
        string_color_name = "orange"
    else
        bar_color_name = "black"
        string_color_name = "deepskyblue"
    end
    bar_color = plt.matplotlib.colors.to_rgba(bar_color_name)
    string_color = plt.matplotlib.colors.to_rgba(string_color_name)
    bars_segs = plt.matplotlib.collections.LineCollection(bars_lines,linewidths=2,colors=bar_color)
    strings_segs = plt.matplotlib.collections.LineCollection(strings_lines,linewidths=1,colors=string_color)
    bars_segs, strings_segs
end

function pyplotstructure(man_linear, sol_linear, tstops)

    steps = [findfirst((x)->x>i,sol_linear.ts) for i = tstops]

    fig,axs_raw = plt.subplots(2,3,figsize=(9,6))
    axs = permutedims(axs_raw)
    for (i,step) in enumerate(steps)
        RB.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
        bars_segs,strings_segs = bars_and_strings_segs(man_linear)
        ax = axs[i]
        if i <= 3
            ax.axes.xaxis.set_visible(false)
        end
        if !(i ∈[1,4])
            ax.axes.yaxis.set_visible(false)
        end
        ax.add_collection(bars_segs)
        ax.add_collection(strings_segs)
        ax.set_ylim(-0.7,0.2)
        ax.set_xlim(-0.1,0.8)
        ax.set_aspect("equal")
        ax.set_title("t=$(tstops[i])")
    end
    fig
end
