
function sliderplot(scene,bars,strings,update_scene!,sol)
    step_slider,tstep = textslider(1:length(sol.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(bars,strings,sol.qs[this_step])
    end
    bigscene = hbox(step_slider,scene)
end

function recordplot(scene,bars,strings,update_scene!,sol)
    record(scene, "swing.mp4", 1:10:length(sol.ts); framerate = 20) do istep
        update_scene!(bars,strings,sol.qs[istep])
    end
end

function bars_and_strings(tgstruct)
    bars = [Node([
        Point(rb.state.p[1]) => Point(rb.state.p[2]);
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
        linesegments!(scene, strings, color = :deepskyblue, linewith = 2)
    end
    scene = Scene(resolution=(632,951))
    plot!(scene,bars,strings)

    xlims!(scene, -0.2,0.2)
    ylims!(scene, -0.6,0.0)
    bars,strings,scene
end

function plotstructure(tgstruct,sol,func)

    bars,strings,scene = plotstructure(tgstruct)

    function update_scene!(bars,strings,q)
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
        strings[] = TR.get_strings(tgstruct)
    end
    func(scene,bars,strings,update_scene!,sol)
end
