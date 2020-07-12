using Makie
AbstractPlotting.__init__()
function plotstructure(tgstruct,sol,confuncs)
    rb = tgstruct.rigidbodies[1]
    bars = Node(vcat(
            [Point(rb.state.p[i]) => Point(rb.state.p[i+4]) for i = 1:4],
            [Point(rb.state.p[i]) => Point(rb.state.p[i+1]) for i = vcat(1:3,5:7)],
            Point(rb.state.p[4]) => Point(rb.state.p[1]),
            Point(rb.state.p[8]) => Point(rb.state.p[5])
            ))

    function plot!(scene::Scene,bars)
        linesegments!(scene, bars, color = :black, linewidth = 4)
    end
    scene = Scene(resolution=(1200,1200))
    plot!(scene,bars)
    cameracontrols(scene).rotationspeed[] = 0.01
    zlims!(scene, 0.0,12)
    xlims!(scene, -6,6)
    ylims!(scene, -6,6)

    function update_scene!(scene,q,tgstruct,bars)
        TR.distribute_q_to_rbs!(tgstruct,q,zero(q))
        cnt = tgstruct.connectivity
        rb = tgstruct.rigidbodies[1]
        bars[] = vcat(
                [Point(rb.state.p[i]) => Point(rb.state.p[i+4]) for i = 1:4],
                [Point(rb.state.p[i]) => Point(rb.state.p[i+1]) for i = vcat(1:3,5:7)],
                Point(rb.state.p[4]) => Point(rb.state.p[1]),
                Point(rb.state.p[8]) => Point(rb.state.p[5])
                )

        afcs,get_gaps,update_contacts = confuncs(tgstruct)
        update_contacts(q,get_gaps(q))
        activecontacts = [afc.impulse.active for afc in afcs]
        @show activecontacts
    end

    step_slider,tstep = textslider(1:length(sol.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(scene,sol.qs[this_step],tgstruct,bars)
    end
    bigscene = hbox(step_slider,scene)

    # record(scene, "2DVert.mp4", 1:10:length(sol.ts); framerate = 20) do istep
    #     update_scene!(scene,sol.qs[istep],tgstruct,bars)
    # end
end
plotstructure(tgrb1,sol,confuncs)
