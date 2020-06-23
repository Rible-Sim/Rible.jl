using Makie
AbstractPlotting.__init__()
function plotstructure(tgstruct,state)
    bars = [Node([
        Point(rb.state.p[1]) => Point(rb.state.p[2]);
        ]) for rb in tgstruct.rigidbodies]

    cnt = tgstruct.connectivity
    rbs = tgstruct.rigidbodies
    function tg_cables(is,s)
        a,b = cnt.string2ap[is]
        state1 = rbs[a.rbid].state
        p1 = state1.p[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.p[b.apid]
        Point(p1) => Point(p2)
    end
    cables = Node([tg_cables(is,s)
             for (is,s) in enumerate(tgstruct.strings)])
    # return  bars,cables
    function plot!(scene::Scene,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 8)
        end
        linesegments!(scene,cables, color = :deepskyblue, linewith = 4)
    end

    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars,cables)
    # return scene
    # xlims!(scene, -0.1,0.1)
    # ylims!(scene, -0.2,0.0)


    function update_scene!(scene,q,tgstruct,bars,cables)
        cnt = tgstruct.connectivity
        TR.distribute_q_to_rbs!(tgstruct,q,zero(q))
        for (id,rb) in enumerate(tgstruct.rigidbodies)
            bars[id][] = [
                    Point(rb.state.p[1]) => Point(rb.state.p[2]);
                ]
        end
        cables[] = [tg_cables(is,s)
                 for (is,s) in enumerate(tgstruct.strings)]
    end


    step_slider,tstep = textslider(1:length(state.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(scene,state.qs[this_step],tgstruct,bars,cables)
    end
    bigscene = hbox(step_slider,scene)
end
plotstructure(tail,sol)

energys = [TR.kinetic_energy_coords(tail,sol.qs[it],sol.q̇s[it]) for it = 1:length(sol.ts)]
plot(energys)
