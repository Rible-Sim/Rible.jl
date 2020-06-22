using Makie
AbstractPlotting.__init__()
function plotstructure(tgstruct,state)
    bars = [Node([
        Point(rb.state.p[1]) => Point(rb.state.p[2]);
        ]) for rb in tgstruct.rigidbodies]
    function plot!(scene::Scene,bars)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 4)
        end
    end

    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars)

    xlims!(scene, -0.1,0.1)
    ylims!(scene, -0.2,0.0)


    function update_scene!(scene,q,tgstruct,bars)
        cnt = tgstruct.connectivity
        for (id,rb) in enumerate(tgstruct.rigidbodies)
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
    end


    step_slider,tstep = textslider(1:length(state.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(scene,state.qs[this_step],tgstruct,bars)
    end
    bigscene = hbox(step_slider,scene)
end
plotstructure(tail,sol)

energys = [TR.energy(sol.qs[it],sol.q̇s[it],tail) for it = 1:length(sol.ts)]
plot(energys)
