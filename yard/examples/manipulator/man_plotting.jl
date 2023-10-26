function get_bars(rb)
    [
    Point(rb.state.rps[1]) => Point(rb.state.rps[2]);
    Point(rb.state.rps[2]) => Point(rb.state.rps[3]);
    Point(rb.state.rps[3]) => Point(rb.state.rps[1]);
    ]
end

function get_lower_bars(rbid,rb)
    if isodd(rbid)
        return [
            Point(rb.state.rps[2]) => Point(rb.state.rps[3]),
            Point(rb.state.rps[3]) => Point(rb.state.rps[1])
        ]
    else
        return [
            Point(rb.state.rps[1]) => Point(rb.state.rps[2])
        ]
    end
end

function get_upper_bars(rbid,rb)
    if iseven(rbid)
        return [
            Point(rb.state.rps[2]) => Point(rb.state.rps[3]),
            Point(rb.state.rps[3]) => Point(rb.state.rps[1])
        ]
    else
        return [
            Point(rb.state.rps[1]) => Point(rb.state.rps[2])
        ]
    end
end

function get_lower_bars(st)
    bars = [Observable(get_lower_bars(rbid,rb)
            ) for (rbid,rb) in enumerate(st.rigidbodies)]
end

function get_upper_bars(st)
    bars = [Observable(get_upper_bars(rbid,rb)
            ) for (rbid,rb) in enumerate(st.rigidbodies)]
end

function get_cables(st)
    (;tensioned) = st.connectivity
    ndim = RB.get_ndim(st)
    T = RB.get_numbertype(st)
    ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
    mapreduce(
        (scnt)->
        Point(scnt.hen.rbsig.state.rps[scnt.hen.pid]) =>
        Point(scnt.egg.rbsig.state.rps[scnt.egg.pid]),
        vcat,
        tensioned.connected
        ;init=ret
    )
end

function get_bars_and_cables(st)
    rbs = RB.get_bodies(st)
    bars = Observable(
                reduce(
                    vcat,[
                        get_bars(rb)
                        for rb in rbs
                    ]
                )
            )
    cables = Observable(get_cables(st))
    bars, cables
end

function plotstructure!(ax,st)
    bars, cables = get_bars_and_cables(st)
    function plot!(ax,bars,cables)
        linesegments!(ax, bars, color = :black, linewidth = 4)
        linesegments!(ax, cables, color = :deepskyblue, linewidth = 2)
    end
    plot!(ax,bars,cables)
    xlims!(ax,-0.1,0.7)
    ylims!(ax,-0.6,0.2)
    bars,cables
end

function update_scene!(st,bars,cables,q)
    cnt = st.connectivity
    RB.update_rigids!(st,q)
    rbs = RB.get_bodies(st)
    bars[] = reduce(
        vcat,[
            get_bars(rb)
            for rb in rbs
        ]
    )
    cables[] = get_cables(st)
    # angles = update_angles(st)
    # @show angles
end

function sliderplot(fig,ax,st,bars,cables,traj;kwargs...)
    ls_step = labelslider!(fig,"step",1:length(traj.t))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(st,bars,cables,traj.q[this_step])
        # @show tstep
    end
end

function sequenceplot(fig,grid,ax,st,bars,cables,update_scene!,seq,ulabels,glabels)
    nseq = length(seq)
    ls_step = labelslider!(fig,"step",1:nseq)
    grid[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        # analyse_slackness(st,sol.qs[this_step])
        update_scene!(st,bars,cables,seq.q[this_step])
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

function recordplot(scene,bars,cables,update_scene!,sol;filename="pid.mp4")
    record(scene, filename, 1:60:length(sol.ts); framerate = 30) do this_step
        update_scene!(bars,cables,sol.qs[this_step])
    end
end

function plotstructure(bot::RB.Robot)
    (;st, traj) = bot
    fig = Figure(resolution=(1280,720))
    ax = Axis(fig[1, 1])
    ax.aspect = DataAspect()
    bars,cables = plotstructure!(ax,st)
    sliderplot(fig,ax,st,bars,cables,traj)
    fig
end

function plotstructure(bot::RB.Robot,seq)
    (;st) = bot
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

    bars,cables = plotstructure!(ax,st)
    sequenceplot(fig,grid_tg,ax,st,bars,cables,update_scene!,seq,ulabels,glabels)
    colsize!(fig.layout, 1, Relative(1/2))
    fig
end
