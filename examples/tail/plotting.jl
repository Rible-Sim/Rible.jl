using Makie
AbstractPlotting.__init__()
function plotstructure(tgstruct,state)
    bars = [Node([
        Point(rb.state.p[1]) => Point(rb.state.p[2]);
        ]) for rb in tgstruct.rigidbodies]

    function get_strings()
        string2ap = tgstruct.connectivity.string2ap
        rbs = tgstruct.rigidbodies
        [Point(rbs[s[1].rbid].state.p[s[1].apid]) =>
         Point(rbs[s[2].rbid].state.p[s[2].apid])
         for s in string2ap]
    end
    strings = Node(get_strings())
    function plot!(scene::Scene,bars,strings)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene, line, color = :black, linewidth = 4)
        end
        linesegments!(scene, strings, color = :deepskyblue, linewith = 2)
    end
    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars,strings)

    xlims!(scene, -0.1,0.1)
    ylims!(scene, -0.2,0.0)

    function update_scene!(scene,q,tgstruct,bars,strings)
        cnt = tgstruct.connectivity
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
                ]
        end
        strings[] = get_strings()
    end


    step_slider,tstep = textslider(1:length(state.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(scene,state.qs[this_step],tgstruct,bars,strings)
    end
    bigscene = hbox(step_slider,scene)

    # record(scene, "2DVert.mp4", 1:10:length(state.ts); framerate = 20) do istep
    #     update_scene!(scene,state.qs[istep],tgstruct,bars,strings)
    # end

end
plotstructure(tail,sol)

energys = [TR.energy(sol.qs[it],sol.q̇s[it],tail) for it = 1:length(sol.ts)]
plot(energys)
