function rb_bars(rb)
    [Point(rb.state.r) => Point(p)
        for p in rb.state.p]
end

function rb_sphere(rb)
    Sphere(Point3f0(rb.state.r), 0.006f0)
end

function bars_and_strings(tg)
    bars = [Node(rb_bars(rb)) for rb in tg.rigidbodies]
    cables = Node(TR.get_strings(tg))
    spheres = [Node(rb_sphere(rb)) for rb in tg.rigidbodies]
    bars, cables, spheres
end

function plotstructure!(scene,tg)
    bars, cables, spheres = bars_and_strings(tg)
    #return bars,cables
    function plot!(scene,bars,cables,spheres)
        for bar in bars
            linesegments!(scene, bar, color = :grey, linewidth = 10)
        end
        for sphere in spheres
            mesh!(sphere, color = :grey)
        end
        linesegments!(scene,cables, color = :deepskyblue, linewidth = 4)

    end

    p = plot!(scene,bars,cables,spheres)
    # cameracontrols(scene.scene).rotationspeed[] = 0.01
    # xlims!(-0.1,1.0)
    # ylims!(-0.5,0.6)
    bars,cables,spheres,p
end

function plotstructure(tg)
    scene = Scene(resolution=(1600,800))
    bars,cables,spheres,p = plotstructure!(scene,tg)
    cameracontrols(scene).rotationspeed[] = 0.01
    scene
end

function sliderplot(scene,bars,cables,spheres,update_scene!,sol)
    step_slider,tstep = textslider(1:length(sol.ts),"step",start=1)
    on(tstep) do this_step
        update_scene!(bars,cables,spheres,sol.qs[this_step])
    end
    bigscene = hbox(step_slider,scene)
end

function update_eye!(scene)
    lookat = Vec3f0(0.00580111, -0.05311403, 0.0038683226)
    eyepos = Vec3f0(0.051770717, -0.052829295, 0.1323204)
    update_cam!(scene, eyepos, lookat)
end

function recordplot(scene,bars,cables,spheres,update_scene!,sol)
    record(scene, "spine3d.mp4", 1:length(sol.ts); framerate = 30) do this_step
        update_scene!(bars,cables,spheres,sol.qs[this_step])
        update_eye!(scene)
    end
end

function plotstructure(tg,sol,func)
    scene = Scene(resolution=(1600,800))
    bars,cables,spheres,p = plotstructure!(scene,tg)
    cameracontrols(scene).rotationspeed[] = 0.01

    function update_scene!(bars,cables,spheres,q)
        cnt = tg.connectivity
        TR.distribute_q_to_rbs!(tg,q)
        for (id,rb) in enumerate(tg.rigidbodies)
            bars[id][] = rb_bars(rb)
            spheres[id][] = rb_sphere(rb)
        end
        cables[] = TR.get_strings(tg)
        # angles = update_angles(tg)
        # @show angles
    end

    func(scene,bars,cables,spheres,update_scene!,sol)
end
