using Makie
AbstractPlotting.__init__()
function plotstructure(tg,sol)

    function rb_bars(rb)
        vcat([
            Point(rb.state.p[1]) => Point(p)
            for p in rb.state.p[2:4]],
            Point(rb.state.p[2]) => Point(rb.state.p[3]),
            Point(rb.state.p[3]) => Point(rb.state.p[4]),
            Point(rb.state.p[4]) => Point(rb.state.p[2]))
    end
    function tg_cables(rb1,rb2)
        vcat(
            [Point(rb1.state.p[j]) => Point(rb2.state.p[j])
             for j = 1:4],
            [Point(rb1.state.p[j]) => Point(rb2.state.p[1])
            for j = 2:4]
        )
    end

    bars = [Node(rb_bars(rb)
            ) for rb in tg.rigidbodies]
    rbs = tg.rigidbodies
    cables = [Node(tg_cables(rbs[i],rbs[i+1])
            ) for i in 1:tg.nbody-1]
    #return bars,cables
    function plot!(scene::Scene,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 10)
        end
        for (lineid,line) in enumerate(cables)
            linesegments!(scene,line, color = :deepskyblue, linewidth = 4)
        end
    end

    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars,cables)
    cameracontrols(scene).rotationspeed[] = 0.01
    # xlims!(-0.1,1.0)
    # ylims!(-0.5,0.6)

    function update_scene!(scene,q)
        cnt = tg.connectivity
        TR.distribute_q_to_rbs!(tg,q,zero(q))
        for (id,rb) in enumerate(tg.rigidbodies)
            bars[id][] = rb_bars(rb)
        end
        for id in 1:tg.nbody-1
            cables[id][] = tg_cables(rbs[id],rbs[id+1])
        end
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
plotstructure(linkn,sol)

b,c = plotstructure(link3)

energys = [TR.energy(linkn,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
plot(energys)
penergys = [TR.potential_energy(linkn,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
plot(penergys)
kenerygs = [TR.kinetic_energy_coords(linkn,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
plot(kenerygs)
plot(kenerygs.+penergys)
