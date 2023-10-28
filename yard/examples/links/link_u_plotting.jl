# using Makie
# AbstractPlotting.__init__()

function rb_bars(body)
    vcat([
        Point(body.state.loci_states[1]) => Point(p)
        for p in body.state.loci_states[2:4]],
        Point(body.state.loci_states[2]) => Point(body.state.loci_states[3]),
        Point(body.state.loci_states[3]) => Point(body.state.loci_states[4]),
        Point(body.state.loci_states[4]) => Point(body.state.loci_states[2]))
    # vcat([
    #     Point(body.state.loci_states[17]) => Point(p)
    #     for p in body.state.loci_states[7:9]],
    #     Point(body.state.loci_states[7]) => Point(body.state.loci_states[8]),
    #     Point(body.state.loci_states[8]) => Point(body.state.loci_states[9]),
    #     Point(body.state.loci_states[9]) => Point(body.state.loci_states[7]))
end

function bars_and_strings(st)
    bars = [Observable(rb_bars(body)
            ) for rb in st.rigidbodies]
    strings = Observable(RB.get_strings(st))
    bars,strings
end

function plotstructure!(ax,st)
    bars,strings = bars_and_strings(st)
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

function plotstructure(st)
    fig = Figure(resolution=(1000,1000))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,strings = plotstructure!(ax,st)
    cameracontrols(ax.scene).rotationspeed[] = 0.01
    fig
end

function recordplot(fig,ax,st,bars,strings,update_scene!,sol,eyepos,lookat)
    framerate = 30
    speed = 4
    h = 1e-2
    record(fig, "links.mp4", 1:round(Int,1/(framerate/speed)/h):length(sol.qs); framerate) do this_step
        analyse_slackness(st,sol.qs[this_step])
        update_scene!(st,bars,strings,sol.qs[this_step])
        update_cam!(ax.scene, eyepos, lookat)
    end
end

function sliderplot(fig,ax,st,bars,strings,update_scene!,sol)
    ls_step = labelslider!(fig,"step",1:length(sol.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        analyse_slackness(st,sol.qs[this_step])
        update_scene!(st,bars,strings,sol.qs[this_step])
        camera = cameracontrols(ax.scene)
        lookat = camera.lookat[]
        eyepos = camera.eyeposition[]
        @show eyepos, lookat
    end
end

function sequenceplot(fig,grid,ax,st,bars,strings,update_scene!,seq,dlabels,ulabels,glabels)
    nseq = length(seq)
    ls_step = labelslider!(fig,"step",1:nseq)
    grid[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        # analyse_slackness(st,sol.qs[this_step])
        update_scene!(st,bars,strings,seq.q[this_step])
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

function update_scene!(st,bars,strings,q)
    cnt = st.connectivity
    RB.distribute_q_to_rbs!(st,q)
    for (id,rb) in enumerate(st.rigidbodies)
        bars[id][] = rb_bars(body)
    end
    strings[] = RB.get_strings(st)
    # angles = update_angles(st)
    # @show angles
end

function plotstructure(bot::RB.Robot;do_record=false)
    @unpack st,traj = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,strings = plotstructure!(ax,st)
    # cameracontrols(ax.scene).mouse_rotationspeed[] = 0.01
    eyepos = Vec3f(0.21623512, 0.5631353, 0.21988256)
    lookat = Vec3f(0.057202004, 0.018000102, 0.12657055)
    update_cam!(ax.scene, eyepos, lookat)
    if do_record
        recordplot(fig,ax,st,bars,strings,update_scene!,traj,eyepos,lookat)
    else
        sliderplot(fig,ax,st,bars,strings,update_scene!,traj)
    end
    fig
end

function plotstructure(bot::RB.Robot,seq)
    @unpack st = bot
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

    bars,strings = plotstructure!(ax,st)
    # cameracontrols(ax.scene).mouse_rotationspeed[] = 0.01
    eyepos = Vec3f(0.21623512, 0.5631353, 0.21988256)
    lookat = Vec3f(0.057202004, 0.018000102, 0.12657055)
    update_cam!(ax.scene, eyepos, lookat)
    sequenceplot(fig,grid_tg,ax,st,bars,strings,update_scene!,seq,dlabels,ulabels,glabels)
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
