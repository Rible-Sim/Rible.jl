function sliderplot(scene,bars,strings,update_scene!,sol;kwargs...)
    step_slider,tstep = textslider(1:length(sol.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(bars,strings,sol.qs[this_step])
        # @show tstep
    end
    bigscene = hbox(step_slider,scene)
end

function recordplot(scene,bars,strings,update_scene!,sol;filename="pid.mp4")
    record(scene, filename, 1:60:length(sol.ts); framerate = 30) do this_step
        update_scene!(bars,strings,sol.qs[this_step])
    end
end

function bars_and_strings(tgstruct)
    bars = [Node([
        Point(rb.state.p[1]) => Point(rb.state.p[2]);
        Point(rb.state.p[2]) => Point(rb.state.p[3]);
        Point(rb.state.p[3]) => Point(rb.state.p[1]);
        ]) for rb in tgstruct.rigidbodies]

    strings = Node(TR.get_strings(tgstruct))
    bars, strings
end

function plotstructure(tgstruct)
    bars, strings = bars_and_strings(tgstruct)
    function plot!(scene::Scene,bars,strings)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene, line, color = :black, linewidth = 4)
        end
        linesegments!(scene, strings, color = :deepskyblue, linewidth = 2)
    end

    scene = Scene(resolution=(1000,1000))
    plot!(scene,bars,strings)
    xlims!(-0.4,1.6)
    ylims!(-1.4,0.6)
    scene,bars,strings
end
# scene,bars,strings = plotstructure(man_linear)
 # plotstructure(man_linear,sol_linear,sliderplot)
function plotstructure(tgstruct,state,func;filename="pid.mp4")

    scene,bars,strings = plotstructure(tgstruct)

    function update_scene!(bars,strings,q)
        cnt = tgstruct.connectivity
        # analyse_slackness(tgstruct,q)
        for id in tgstruct.mvbodyindex
            rb = tgstruct.rigidbodies[id]
            pindex = cnt.body2q[id]
            cache = rb.state.cache
            p = rb.state.p
            for (i,ap) in enumerate(p)
                ap .= cache.Cp[i]*q[pindex]
            end
            bars[id][] = [
                    Point(rb.state.p[1]) => Point(rb.state.p[2]);
                    Point(rb.state.p[2]) => Point(rb.state.p[3]);
                    Point(rb.state.p[3]) => Point(rb.state.p[1]);
                ]
        end
        strings[] = TR.get_strings(tgstruct)
        # angles = update_angles(tgstruct)
        # @show angles
    end
    func(scene,bars,strings,update_scene!,state;filename)
end

# plotstructure(manipulator,sol,sliderplot)
# plotstructure(manipulator,sol,recordplot)

# plotstructure(manipulator)

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
        bar_color_name = "lightgrey"
        string_color_name = "lightgrey"
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
        TR.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
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
