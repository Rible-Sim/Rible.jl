function get_cables(tg)
    (;connected) = tg.connectivity
    ret = Vector{Pair{Point2{Float64},Point2{Float64}}}(undef,tg.ncables)
    map!(ret,connected) do scnt
        Point(scnt.end1.rbsig.state.rps[scnt.end1.pid]) =>
        Point(scnt.end2.rbsig.state.rps[scnt.end2.pid])
    end
    ret
end

function rb_bars(rb)
    (;lncs) = rb.state.cache.funcs
    if lncs isa TR.NCF.LocalNCF2D2P
        return [
            Point(rb.state.rps[1]) => Point(rb.state.rps[2]);
        ]
    else
        return [
            Point(rb.state.rps[1]) => Point(rb.state.rps[2]);
            Point(rb.state.rps[2]) => Point(rb.state.rps[3]);
            Point(rb.state.rps[3]) => Point(rb.state.rps[1]);
        ]
    end
end

function update_scene!(tg,bars,cables,q)
    cnt = tg.connectivity
    TR.update_rigids!(tg,q)
    foreach(tg.rigidbodies) do rb
        bars[rbid][] = rb_bars(rb)
    end
    cables[] = TR.get_cables(tg)
    # angles = update_angles(tg)
    # @show angles
end

function sliderplot(fig,ax,tg,bars,cables,traj;kwargs...)
    ls_step = labelslider!(fig,"step",1:length(traj))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(tg,bars,cables,traj.q[this_step])
        # @show tstep
    end
end

function recordplot(scene,bars,cables,update_scene!,sol;filename="pid.mp4")
    record(scene, filename, 1:60:length(sol.ts); framerate = 30) do this_step
        update_scene!(bars,cables,sol.qs[this_step])
    end
end

function bars_and_cables(tg)
    bars = Vector{Observable}(undef,tg.nrigids)
    map!(bars,tg.rigidbodies) do rb
        Observable(rb_bars(rb))
    end
    # @show bars
    cables = Observable(get_cables(tg))
    bars, cables
end

function plotstructure!(ax,tg)
    bars, cables = bars_and_cables(tg)
    function plot!(ax,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(ax, line, color = :black, linewidth = 4)
        end
        linesegments!(ax, cables, color = :deepskyblue, linewidth = 2)
    end
    plot!(ax,bars,cables)
    bars,cables
end

function plotstructure(bot::TR.TensegrityRobot)
    @unpack tg, traj = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = Axis(fig)
    bars,cables = plotstructure!(ax,tg)
    sliderplot(fig,ax,tg,bars,cables,traj)
    ax.aspect = DataAspect()
    fig
end


function bars_and_cables_segs(tgstruct;ref=false)
    bars,cables = bars_and_cables(tgstruct)

    bars_lines = Vector{Vector{Tuple{Float64,Float64}}}()
    for b in bars
        for l in b[]
            lines = [l.first.data,l.second.data]
            push!(bars_lines,lines)
        end
    end
    cables_lines = Vector{Vector{Tuple{Float64,Float64}}}()
    for s in cables[]
        lines = [s.first.data,s.second.data]
        push!(cables_lines,lines)
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
    cables_segs = plt.matplotlib.collections.LineCollection(cables_lines,linewidths=1,colors=string_color)
    bars_segs, cables_segs
end

function pyplotstructure(man_linear, sol_linear, tstops)

    steps = [findfirst((x)->x>i,sol_linear.ts) for i = tstops]

    fig,axs_raw = plt.subplots(2,3,figsize=(9,6))
    axs = permutedims(axs_raw)
    for (i,step) in enumerate(steps)
        TR.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
        bars_segs,cables_segs = bars_and_cables_segs(man_linear)
        ax = axs[i]
        if i <= 3
            ax.axes.xaxis.set_visible(false)
        end
        if !(i ∈[1,4])
            ax.axes.yaxis.set_visible(false)
        end
        ax.add_collection(bars_segs)
        ax.add_collection(cables_segs)
        ax.set_ylim(-0.7,0.2)
        ax.set_xlim(-0.1,0.8)
        ax.set_aspect("equal")
        ax.set_title("t=$(tstops[i])")
    end
    fig
end
