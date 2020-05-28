using Makie
bars = [Node([
    Point(rb.state.p[1]) => Point(rb.state.p[2]);
    ]) for rb in rbs]
function AbstractPlotting.plot!(scene::Scene,bars)
    for (lineid,line) in enumerate(bars)
        linesegments!(scene,line, color = :black, linewidth = 4)
    end
end

scene = Scene(resolution=(1200,1200))
plot!(scene,bars)

xlims!(scene, -0.1,0.1)
ylims!(scene, -0.2,0.0)

function update_scene!(scene,q,st2d,bars)
    cnt = st2d.connectivity
    for (id,rb) in enumerate(rbs)
        pindex = cnt[id]
        auxs = rb.state.auxs
        p = rb.state.p
        for (i,ap) in enumerate(p)
            ap .= auxs.Cp[i]*q[pindex]
        end
        bars[id][] = [
            Point(rb.state.p[1]) => Point(rb.state.p[2]);
            ]
    end
end
update_scene!(scene,state.qs[5],tailstruct,bars)

step_slider,tstep = textslider(1:length(state.ts),"step",start=1)
on(tstep) do this_step
    update_scene!(scene,state.qs[this_step],tailstruct,bars)
end
bigscene = hbox(step_slider,scene)

record(scene, "four_vert.mp4", 1:10:length(state.ts); framerate = 10) do istep
    update_scene!(scene,state.qs[istep],rbs,bars)
end

PackageCompiler.create_sysimage(
    :Makie;
    sysimage_path="MakieSys.so",
    precompile_execution_file="precompile_file.jl")
)
