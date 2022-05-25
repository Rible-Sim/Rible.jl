rb_bars(rb) = [Point(rb.state.rps[1]) => Point(rb.state.rps[2]);]

function bars_and_cables(tgstruct)
    bars = [Node(rb_bars(rb)) for rb in tgstruct.rigidbodies]
    cables = Node(TR.get_cables(tgstruct))
    bars, cables
end

function plotstructure!(ax,tgstruct)
    bars, cables = bars_and_cables(tgstruct)
    function plot!(ax,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(ax, line, color = :black, linewidth = 4)
        end
        linesegments!(ax, cables, color = :deepskyblue, linewidth = 2)
    end
    plot!(ax,bars,cables)
    xlims!(ax,[-0.2,0.2])
    ylims!(ax,[-0.6,0.05])
    bars,cables
end

function update_scene!(tg,bars,cables,q)
    cnt = tg.connectivity
    TR.distribute_q_to_rbs!(tg,q)
    for (id,rb) in enumerate(tg.rigidbodies)
        bars[id][] = rb_bars(rb)
    end
    cables[] = TR.get_cables(tg)
    # angles = update_angles(tg)
    # @show angles
end

function sliderplot(fig,tr,ax,cables)
    @unpack tg, traj = tr
    ls_step = labelslider!(fig,"step",1:length(traj.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        # analyse_slackness(tg,sol.qs[this_step])
        update_scene!(tg,ax,cables,traj.qs[this_step])
        # show_camera(ax.scene)
    end
end
function recordplot(fig,tr,ax,cables;record_name)
    @unpack tg, traj = tr
    record(fig, record_name, 1:length(traj.ts); framerate = 30) do this_step
        update_scene!(tg,ax,cables,traj.qs[this_step])
        # update_cam!(ax.scene, eyepos,lookat)
    end
end


function plotstructure(tr::TR.TensegrityRobot;record_name="")
    @unpack tg, traj = tr
    fig = Figure(resolution=(1000,1000))
    ax = fig[1, 1] = Axis(fig)
    bars,cables = plotstructure!(ax,tg)
    ax.aspect = DataAspect()
    if record_name==""
        sliderplot(fig,tr,ax,cables)
    else
        recordplot(fig,tr,ax,cables;record_name)
    end
    fig
end
