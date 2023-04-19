# using Makie
# AbstractPlotting.__init__()

function rb_bars(rb)
    vcat([
        Point(rb.state.rps[1]) => Point(p)
        for p in rb.state.rps[2:4]],
        Point(rb.state.rps[2]) => Point(rb.state.rps[3]),
        Point(rb.state.rps[3]) => Point(rb.state.rps[4]),
        Point(rb.state.rps[4]) => Point(rb.state.rps[2]))
    # vcat([
    #     Point(rb.state.rps[17]) => Point(p)
    #     for p in rb.state.rps[7:9]],
    #     Point(rb.state.rps[7]) => Point(rb.state.rps[8]),
    #     Point(rb.state.rps[8]) => Point(rb.state.rps[9]),
    #     Point(rb.state.rps[9]) => Point(rb.state.rps[7]))
end

function bars_and_strings(tg)
    bars = [Observable(rb_bars(rb)
            ) for rb in tg.rigidbodies]
    strings = Observable(TR.get_strings(tg))
    bars,strings
end

function plotstructure!(ax,tg)
    bars,strings = bars_and_strings(tg)
    #return bars,strings
    function plot!(ax,bars,strings)
        for (lineid,line) in enumerate(bars)
            linesegments!(ax,line, color = :black, linewidth = 6)
        end
        linesegments!(ax,strings, color = RGB(0x4A86E8), linewidth = 3)
    end

    plot!(ax,bars,strings)
    # cameracontrols(ax.scene).rotationspeed[] = 0.01
    # xlims!(ax,-0.15,0.15)
    # ylims!(ax,-0.15,0.15)
    # zlims!(ax,-0.15,0.10)
    bars,strings
end

function plotstructure(tg)
    fig = Figure(resolution=(1000,1000))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,strings = plotstructure!(ax,tg)
    cameracontrols(ax.scene).rotationspeed[] = 0.01
    fig
end

function recordplot(fig,ax,tg,bars,strings,update_scene!,sol,eyepos,lookat)
    framerate = 30
    speed = 4
    h = 1e-2
    record(fig, "links.mp4", 1:round(Int,1/(framerate/speed)/h):length(sol.qs); framerate) do this_step
        analyse_slackness(tg,sol.qs[this_step])
        update_scene!(tg,bars,strings,sol.qs[this_step])
        update_cam!(ax.scene, eyepos, lookat)
    end
end

function sliderplot(fig,ax,tg,bars,strings,update_scene!,sol)
    ls_step = labelslider!(fig,"step",1:length(sol.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        analyse_slackness(tg,sol.qs[this_step])
        update_scene!(tg,bars,strings,sol.qs[this_step])
        camera = cameracontrols(ax.scene)
        lookat = camera.lookat[]
        eyepos = camera.eyeposition[]
        @show eyepos, lookat
    end
end

function sequenceplot(fig,grid,ax,tg,bars,strings,update_scene!,seq,dlabels,ulabels,glabels)
    nseq = length(seq)
    ls_step = labelslider!(fig,"step",1:nseq)
    grid[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        # analyse_slackness(tg,sol.qs[this_step])
        update_scene!(tg,bars,strings,seq.q[this_step])
        ls_step.valuelabel.text[] = string((this_step-1)/(nseq-1))
        for (id,dlabel) in enumerate(dlabels)
            dlabel.text[] = "d$id=$(seq.d[this_step][id])"
        end
        for (iu,ulabel) in enumerate(ulabels)
            ulabel.text[] = "u$iu=$(seq.u[this_step][iu])"
        end
        for (ig,glabel) in enumerate(glabels)
            glabel.text[] = "g$ig=$(seq.g[this_step][ig])"
        end
        # camera = cameracontrols(ax.scene)
        # lookat = camera.lookat[]
        # eyepos = camera.eyeposition[]
        # @show eyepos, lookat

    end
end

function update_scene!(tg,bars,strings,q)
    cnt = tg.connectivity
    TR.distribute_q_to_rbs!(tg,q)
    for (id,rb) in enumerate(tg.rigidbodies)
        bars[id][] = rb_bars(rb)
    end
    strings[] = TR.get_strings(tg)
    # angles = update_angles(tg)
    # @show angles
end

function plotstructure(bot::TR.TensegrityRobot;do_record=false)
    @unpack tg,traj = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,strings = plotstructure!(ax,tg)
    # cameracontrols(ax.scene).mouse_rotationspeed[] = 0.01
    eyepos = Vec3f(0.21623512, 0.5631353, 0.21988256)
    lookat = Vec3f(0.057202004, 0.018000102, 0.12657055)
    update_cam!(ax.scene, eyepos, lookat)
    if do_record
        recordplot(fig,ax,tg,bars,strings,update_scene!,traj,eyepos,lookat)
    else
        sliderplot(fig,ax,tg,bars,strings,update_scene!,traj)
    end
    fig
end

function plotstructure(bot::TR.TensegrityRobot,seq)
    @unpack tg = bot
    fig = Figure(resolution=(1280,720))
    grid_tg = fig[1, 1] = GridLayout()
    ax = grid_tg[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    grid_para = fig[:, 2] = GridLayout(tellheight = false) # Show parameters
    ax_d = grid_para[1,1]
    ax_u = grid_para[1,2]
    ax_g = grid_para[1,3]
    dlabels = [Label(ax_d[i, 1], "d$i=$(seq[1].d[i])", textsize = 30) for i in 1:length(seq[1].d)]
    ulabels = [Label(ax_u[i, 1], "u$i=$(seq[1].u[i])", textsize = 30) for i in 1:length(seq[1].u)]
    glabels = [Label(ax_g[i, 1], "g$i=$(seq[1].g[i])", textsize = 30) for i in 1:length(seq[1].g)]

    bars,strings = plotstructure!(ax,tg)
    # cameracontrols(ax.scene).mouse_rotationspeed[] = 0.01
    eyepos = Vec3f(0.21623512, 0.5631353, 0.21988256)
    lookat = Vec3f(0.057202004, 0.018000102, 0.12657055)
    update_cam!(ax.scene, eyepos, lookat)
    sequenceplot(fig,grid_tg,ax,tg,bars,strings,update_scene!,seq,dlabels,ulabels,glabels)
    colsize!(fig.layout, 1, Relative(1/2))
    fig
end
# plotstructure(linkn)
#
# plotstructure(linkn,sol)
function bars_and_strings_segs_3D(tgstruct;ref=false)
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
    if ref
        bar_color_name = string_color_name = "orange"
        bar_ls = string_ls = "--"
        bar_lw = string_lw = 1.0
    else
        bar_color_name = "black"
        string_color_name = "deepskyblue"
        bar_ls = string_ls = "-"
        bar_lw = 1.5
        string_lw = 0.8
    end
    bar_color = plt.matplotlib.colors.to_rgba(bar_color_name)
    string_color = plt.matplotlib.colors.to_rgba(string_color_name)
    bars_segs = plt.art3D.Line3DCollection(bars_lines,linestyles=bar_ls,linewidths=bar_lw,colors=bar_color)
    strings_segs = plt.art3D.Line3DCollection(strings_lines,linestyles=string_ls,linewidths=string_lw,colors=string_color)
    bars_segs, strings_segs
end
