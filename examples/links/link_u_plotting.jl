# using Makie
# AbstractPlotting.__init__()

function rb_bars(rb)
    vcat([
        Point(rb.state.p[1]) => Point(p)
        for p in rb.state.p[2:4]],
        Point(rb.state.p[2]) => Point(rb.state.p[3]),
        Point(rb.state.p[3]) => Point(rb.state.p[4]),
        Point(rb.state.p[4]) => Point(rb.state.p[2]))
end

function bars_and_strings(tg)
    bars = [Node(rb_bars(rb)
            ) for rb in tg.rigidbodies]
    cables = Node(TR.get_strings(tg))
    bars,cables
end

function plotstructure!(scene,tg)
    bars,cables = bars_and_strings(tg)
    #return bars,cables
    function plot!(scene,bars,cables)
        for (lineid,line) in enumerate(bars)
            linesegments!(scene,line, color = :black, linewidth = 2)
        end
        linesegments!(scene,cables, color = :deepskyblue, linewidth = 2)
    end

    p = plot!(scene,bars,cables)
    # cameracontrols(scene.scene).rotationspeed[] = 0.01
    # xlims!(-0.1,1.0)
    # ylims!(-0.5,0.6)
    bars,cables,p
end

function plotstructure(tg)
    scene = Scene(resolution=(1600,1000))
    bars,cables,p = plotstructure!(scene,tg)
    cameracontrols(scene).rotationspeed[] = 0.01
    scene
end

function recordplot(scene,update_scene!,state)
    record(scene, "links_damped.mp4", 1:20:length(state.ts); framerate = 30) do this_step
        update_scene!(tg,bars,cables,sol.qs[this_step])
        update_cam!(scene, eyepos)
    end
end

function sliderplot(scene,tg,bars,cables,update_scene!,sol)
    step_slider,tstep = textslider(1:length(sol.ts),"step",start=1)
    on(tstep) do this_step
        analyse_slackness(tg,sol.qs[this_step])
        update_scene!(bars,cables,sol.qs[this_step])
    end
    bigscene = hbox(step_slider,scene)
end

function update_scene!(tg,bars,cables,q)
    cnt = tg.connectivity
    TR.distribute_q_to_rbs!(tg,q)
    for (id,rb) in enumerate(tg.rigidbodies)
        bars[id][] = rb_bars(rb)
    end
    cables[] = TR.get_strings(tg)
    # angles = update_angles(tg)
    # @show angles
end

function plotstructure(tg,sol,func)
    scene = Scene(resolution=(1600,1000))
    bars,cables,p = plotstructure!(scene,tg)
    cameracontrols(scene).rotationspeed[] = 0.01
    eyepos = Vec3f0(0.028828546, 0.39417213, 0.20687535)
    lookat = Vec3f0(0.035752997, 0.0071086287, 0.09679556)
    update_cam!(scene, eyepos, lookat)
    function update_scene!(bars,cables,q)
        cnt = tg.connectivity
        TR.distribute_q_to_rbs!(tg,q)
        for (id,rb) in enumerate(tg.rigidbodies)
            bars[id][] = rb_bars(rb)
        end
        cables[] = TR.get_strings(tg)
        # angles = update_angles(tg)
        # @show angles
    end
    func(scene,tg,bars,cables,update_scene!,sol)
end

# plotstructure(linkn)
#
# plotstructure(linkn,sol)
function bars_and_strings_segs_3D(tgstruct)
    bars,strings = bars_and_strings(tgstruct)

    bars_lines = Vector{Vector{Tuple{Float64,Float64,Float64}}}()
    for b in bars
        for l in b[]
            lines = [l.first.data,l.second.data]
            push!(bars_lines,lines)
        end
    end
    strings_lines = Vector{Vector{Tuple{Float64,Float64,Float64}}}()
    for s in strings[]
        lines = [s.first.data,s.second.data]
        push!(strings_lines,lines)
    end
    color_black = plt.matplotlib.colors.to_rgba("black")
    color_blue = plt.matplotlib.colors.to_rgba("deepskyblue")
    bars_segs = plt.art3D.Line3DCollection(bars_lines,linewidths=2,colors=color_black)
    strings_segs = plt.art3D.Line3DCollection(strings_lines,linewidths=1,colors=color_blue)
    bars_segs, strings_segs
end
