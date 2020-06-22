using Makie
AbstractPlotting.__init__()

function sliderplot(scene,update_scene!,state)
    step_slider,tstep = textslider(1:length(state.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(scene,state.qs[this_step])
    end
    bigscene = hbox(step_slider,scene)
end

function recordplot(scene,update_scene!,state)
    record(scene, "pid.mp4", 1:20:length(state.ts); framerate = 20) do istep
        update_scene!(scene,state.qs[istep])
    end
end

function plotstructure(tgstruct,state,func)
    bars = [Node([
        Point(rb.state.p[1]) => Point(rb.state.p[2]);
        Point(rb.state.p[2]) => Point(rb.state.p[3]);
        Point(rb.state.p[3]) => Point(rb.state.p[1]);
        ]) for rb in tgstruct.rigidbodies]
    function plot!(scene::Scene,bars)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 4)
        end
    end

    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars)
    xlims!(-0.1,1.0)
    ylims!(-0.5,0.6)

    function update_scene!(scene,q)
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
                    Point(rb.state.p[2]) => Point(rb.state.p[3]);
                    Point(rb.state.p[3]) => Point(rb.state.p[1]);
                ]
        end
        # angles = update_angles(tgstruct)
        # @show angles
    end
    func(scene,update_scene!,state)
end

plotstructure(manipulator,sol,sliderplot)
plotstructure(manipulator,sol,recordplot)

energys = [R2.energy(state.qs[it],state.q̇s[it],manipulator) for it = 1:length(state.ts)]
plot(energys)
ylims!(0.0,0.002)
manipulator.rigidbodies[1]
