function rb_bars(body)
    Cg = body.state.cache.Cg
    q = body.state.coords.q
    rg = Cg*q
    [Point(rg) => Point(p)
        for p in body.state.loci_states]
end

function rb_sphere(body)
    Cg = body.state.cache.Cg
    q = body.state.coords.q
    rg = Cg*q
    Sphere(Point3f(rg), 0.006f0)
end

function bars_and_cables(st)
    bars = [Observable(rb_bars(body)) for rb in st.rigidbodies]
    cables = Observable(RB.get_cables(st))
    spheres = [Observable(rb_sphere(body)) for rb in st.rigidbodies]
    bars, cables, spheres
end

function plotstructure!(ax,st)
    bars, cables, spheres = bars_and_cables(st)
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

function plotstructure(st::RB.Structure)
    scene = Scene(resolution=(1280,720))
    bars,cables,spheres,p = plotstructure!(scene,st)
    # cameracontrols(scene).rotationspeed[] = 0.01
    scene
end

function sliderplot(fig,ax,st,bars,cables,spheres,traj)
    ls_step = labelslider!(fig,"step",1:length(traj.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        update_scene!(st,bars,cables,spheres,traj.qs[this_step])
    end
    fig
end

function update_eye!(scene)
    lookat = Vec3f(0.00580111, -0.05311403, 0.0038683226)
    eyepos = Vec3f(0.051770717, -0.052829295, 0.1323204)
    update_cam!(scene, eyepos, lookat)
end

function recordplot(fig,ax,st,bars,cables,spheres,traj)
    framerate = 30
    speed = 4
    h = 1e-2
    record(fig, "spine3d.mp4", 1:round(Int,1/(framerate/speed)/h):length(traj.ts); framerate = 30) do this_step
        update_scene!(st,bars,cables,spheres,traj.qs[this_step])
        update_eye!(ax.scene)
    end
end

function update_scene!(st,bars,cables,spheres,q)
    cnt = st.connectivity
    RB.distribute_q_to_rbs!(st,q)
    for (id,rb) in enumerate(st.rigidbodies)
        bars[id][] = rb_bars(body)
        spheres[id][] = rb_sphere(body)
    end
    cables[] = RB.get_cables(st)
    # angles = update_angles(st)
    # @show angles
end

function plotstructure(bot::RB.Robot;do_record=false)
    @unpack st,traj = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    bars,cables,spheres = plotstructure!(ax,st)
    # cameracontrols(ax.scene).rotationspeed[] = 0.01
    # update_eye!(ax.scene)
    if do_record
        recordplot(fig,ax,st,bars,cables,spheres,traj)
    else
        sliderplot(fig,ax,st,bars,cables,spheres,traj)
    end
    fig
end
