using Makie
AbstractPlotting.__init__()

function plotstructure(tg)

    function rb_bars(rb)
        [
            Point(rb.state.r) => Point(p)
            for p in rb.state.p[1:3]
        ]
    end
    bars = [Node(rb_bars(rb)
            ) for rb in tg.rigidbodies]
    rbs = tg.rigidbodies
    cables = TR.get_strings(tg)
    #return bars,cables
    function plot!(scene::Scene,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 4)
        end
        linesegments!(scene,cables, color = :deepskyblue, linewidth = 2)
    end

    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars,cables)

    xlims!(-0.1,0.8)
    ylims!(-0.2,0.3)
    # cameracontrols(scene).rotationspeed[] = 0.01
    scene
end

plotstructure(spine)
function plotstructure(tg,sol)

    function rb_bars(rb)
        [
            Point(rb.state.r) => Point(p)
            for p in rb.state.p[1:3]
        ]
    end
    bars = [Node(rb_bars(rb)
            ) for rb in tg.rigidbodies]
    rbs = tg.rigidbodies
    cables = TR.get_strings(tg)
    #return bars,cables
    function plot!(scene::Scene,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 4)
        end
        #linesegments!(scene,cables, color = :deepskyblue, linewidth = 2)
    end

    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars,cables)

    xlims!(-0.1,0.8)
    ylims!(-0.2,0.3)

    function update_scene!(scene,q)
        cnt = tg.connectivity
        TR.distribute_q_to_rbs!(tg,q,zero(q))
        for (id,rb) in enumerate(tg.rigidbodies)
            bars[id][] = rb_bars(rb)
        end
        #cables[] = TR.get_strings(tg)
        # angles = update_angles(tg)
        # @show angles
    end
    function sliderplot(scene,update_scene!,sol)
        step_slider,tstep = textslider(1:length(sol.ts),"step",start=1)
        on(tstep) do this_step
            update_scene!(scene,sol.qs[this_step])
        end
        bigscene = hbox(step_slider,scene)
    end
    sliderplot(scene,update_scene!,sol)
end

plotstructure(spine,sol)
