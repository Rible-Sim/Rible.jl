function rb_bars(rb)
    Cg = rb.state.cache.Cg
    q = rb.state.coords.q
    rg = Cg*q
    [Point(rg) => Point(p)
        for p in rb.state.rps]
end

function rb_sphere(rb)
    Cg = rb.state.cache.Cg
    q = rb.state.coords.q
    rg = Cg*q
    Sphere(Point3f(rg), 0.006f0)
end

function bars_and_cables(tg)
    bars = [Observable(rb_bars(rb)) for rb in tg.rigidbodies]
    cables = Observable(TR.get_cables(tg))
    spheres = [Observable(rb_sphere(rb)) for rb in tg.rigidbodies]
    bars, cables, spheres
end

function plotstructure!(ax,tg)
    bars, cables, spheres = bars_and_cables(tg)
    #return bars,cables
    function plot!(ax,bars,cables,spheres)
        for bar in bars
            linesegments!(ax, bar, color = :grey, linewidth = 10)
        end
        for sphere in spheres
            mesh!(ax, sphere, color = :grey)
        end
        linesegments!(ax,cables, color = :deepskyblue, linewidth = 4)

    end

    plot!(ax,bars,cables,spheres)
    # cameracontrols(scene.scene).rotationspeed[] = 0.01
    # xlims!(-0.1,1.0)
    # ylims!(-0.5,0.6)
    bars,cables,spheres
end

function plotstructure(tg::TR.TensegrityStructure)
    scene = Scene(resolution=(1280,720))
    bars,cables,spheres,p = plotstructure!(scene,tg)
    # cameracontrols(scene).rotationspeed[] = 0.01
    scene
end

function sliderplot(fig,ax,tg,bars,cables,spheres,traj)
    ls_step = labelslider!(fig,"step",1:length(traj.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(tg,bars,cables,spheres,traj.qs[this_step])
    end
    fig
end

function update_eye!(scene)
    lookat = Vec3f(0.00580111, -0.05311403, 0.0038683226)
    eyepos = Vec3f(0.051770717, -0.052829295, 0.1323204)
    update_cam!(scene, eyepos, lookat)
end

function recordplot(fig,ax,tg,bars,cables,spheres,traj)
    framerate = 30
    speed = 4
    h = 1e-2
    record(fig, "spine3d.mp4", 1:round(Int,1/(framerate/speed)/h):length(traj.ts); framerate = 30) do this_step
        update_scene!(tg,bars,cables,spheres,traj.qs[this_step])
        update_eye!(ax.scene)
    end
end

function update_scene!(tg,bars,cables,spheres,q)
    cnt = tg.connectivity
    TR.distribute_q_to_rbs!(tg,q)
    for (id,rb) in enumerate(tg.rigidbodies)
        bars[id][] = rb_bars(rb)
        spheres[id][] = rb_sphere(rb)
    end
    cables[] = TR.get_cables(tg)
    # angles = update_angles(tg)
    # @show angles
end

function plotstructure(bot::TR.TensegrityRobot;do_record=false)
    @unpack tg,traj = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,cables,spheres = plotstructure!(ax,tg)
    # cameracontrols(ax.scene).rotationspeed[] = 0.01
    # update_eye!(ax.scene)
    if do_record
        recordplot(fig,ax,tg,bars,cables,spheres,traj)
    else
        sliderplot(fig,ax,tg,bars,cables,spheres,traj)
    end
    fig
end
