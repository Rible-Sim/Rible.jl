function get_cables(st)
    (;connected) = st.connectivity
    ret = Vector{Pair{Point2{Float64},Point2{Float64}}}(undef,st.ncables)
    map!(ret,connected) do scnt
        Point(scnt.hen.rbsig.state.loci_states[scnt.hen.pid]) =>
        Point(scnt.egg.rbsig.state.loci_states[scnt.egg.pid])
    end
    ret
end

function rb_bars(body)
    (;lncs) = body.state.cache.funcs
    if lncs isa RB.NCF.LocalNCF2D2P
        return [
            Point(body.state.loci_states[1]) => Point(body.state.loci_states[2]);
        ]
    else
        return [
            Point(body.state.loci_states[1]) => Point(body.state.loci_states[2]);
            Point(body.state.loci_states[2]) => Point(body.state.loci_states[3]);
            Point(body.state.loci_states[3]) => Point(body.state.loci_states[1]);
        ]
    end
end

function update_scene!(st,bars,cables,q)
    cnt = st.connectivity
    RB.update_rigids!(st,q)
    foreach(st.rigidbodies) do body
        bars[bodyid][] = rb_bars(body)
    end
    cables[] = RB.get_cables(st)
    # angles = update_angles(st)
    # @show angles
end

function sliderplot(fig,ax,st,bars,cables,traj;kwargs...)
    ls_step = labelslider!(fig,"step",1:length(traj))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(st,bars,cables,traj.q[this_step])
        # @show tstep
    end
end

function recordplot(scene,bars,cables,update_scene!,sol;filename="pid.mp4")
    record(scene, filename, 1:60:length(sol.ts); framerate = 30) do this_step
        update_scene!(bars,cables,sol.qs[this_step])
    end
end

function bars_and_cables(st)
    bars = Vector{Observable}(undef,st.nrigids)
    map!(bars,st.rigidbodies) do body
        Observable(rb_bars(body))
    end
    # @show bars
    cables = Observable(get_cables(st))
    bars, cables
end

function plotstructure!(ax,st)
    bars, cables = bars_and_cables(st)
    function plot!(ax,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(ax, line, color = :black, linewidth = 4)
        end
        linesegments!(ax, cables, color = :deepskyblue, linewidth = 2)
    end
    plot!(ax,bars,cables)
    bars,cables
end

function plotstructure(bot::RB.Robot)
    @unpack st, traj = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = Axis(fig)
    bars,cables = plotstructure!(ax,st)
    sliderplot(fig,ax,st,bars,cables,traj)
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
        RB.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
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
