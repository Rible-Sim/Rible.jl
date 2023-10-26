
function get_cube_bars(st)
    rb = RB.get_bodies(st)[1]
    ps = Point.(rb.state.rps)
    fs = [
        GeometryBasics.QuadFace(1,3,4,2),
        GeometryBasics.QuadFace(5,6,8,7),
        GeometryBasics.QuadFace(3,7,8,4),
        GeometryBasics.QuadFace(1,2,6,5),
        GeometryBasics.QuadFace(2,4,8,6),
        GeometryBasics.QuadFace(3,1,5,7)
    ]
    bars = vcat(
        [ps[i]=>ps[i+4] for i = 1:4],
        [ps[i]=>ps[i+2] for i in [1,2,5,6]],
        [ps[i]=>ps[i+1] for i in [1,3,5,7]],
    )
    GeometryBasics.mesh(GeometryBasics.Mesh(ps,fs),normaltype=Vec3f),bars
end

function update_scene!(ax,st,cube_bars,q)
    # RB.reset_forces!(st)
    RB.update_rigids!(st,q)
    newcube, newbars = get_cube_bars(st)
    cube, bars = cube_bars
    cube[] = newcube
    bars[] = newbars

    # afcs,get_gaps,update_contacts = confuncs(st)
    # update_contacts(q,get_gaps(q))
    # activecontacts = [afc.impulse.active for afc in afcs]
    # @show activecontacts
end

function sliderplot(fig,st,bot,ax,cube_bars,contact_dynamics)
    # 𝒈,get_indices,,get_contacts = contact_dynamics(st)
    (;traj) = bot
    sg = SliderGrid(fig[2,1],
            (label="step", range=1:length(traj.t))
        )
    on(sg.sliders[1].value) do this_step
        # analyse_slackness(st,sol.qs[this_step])
        q = traj.q[this_step]
        update_scene!(ax,st,cube_bars,q)
        # @show get_contacts(q)
        camera = cameracontrols(ax.scene)
        lookat = camera.lookat[]
        eyepos = camera.eyeposition[]
        # @show lookat, eyepos
    end
end

# record(scene, "2DVert.mp4", 1:10:length(sol.ts); framerate = 20) do istep
#     update_scene!(scene,sol.qs[istep],st,cube)
# end

function vis!(ax,bot,cube_bars)
    cube, bars = Observable.(cube_bars)
    xmin, xmax = -10, 10
    ymin, ymax = -10, 10
    zmin, zmax =  -1, 19
    points = [
        Point3(xmin,ymin,zmin),Point3(xmin,ymin,zmin),
        Point3(xmin,ymin,zmin),Point3(xmax,ymax,zmin),
        Point3(xmin,ymin,zmin),Point3(xmin,ymax,zmin),
        Point3(xmin,ymin,zmin),Point3(xmax,ymax,zmin),
        Point3(xmin,ymin,zmax),Point3(xmin,ymin,zmax),
        Point3(xmin,ymin,zmax),Point3(xmax,ymax,zmax),
        Point3(xmin,ymin,zmax),Point3(xmin,ymax,zmax),
        Point3(xmin,ymin,zmax),Point3(xmax,ymax,zmax),
    ]
    ground = RB.Plane([0,0,1.0],[0.0,0.0,0.0])
    ground_mesh = GeometryBasics.Mesh(Rect(points),
                    MarchingCubes(eps=1e-14)) do v
        RB.signed_distance(v,ground)
    end
    mesh!(ax, ground_mesh, color = :white, ambient = Vec3f(0.9, 0.9, 0.9),)
    mesh!(ax, cube, color = :lightblue,
                            linewidth = 4,
                            lightposition = Vec3f(0, 0, 15))
    linesegments!(ax, bars, linewidth = 2)
    # zlims!(ax.scene, 0.0,10)
    # xlims!(ax.scene, -12,12)
    # ylims!(ax.scene, -12,12)
    # ax.aspect = DataAspect()
    cube, bars
end

function vis(bot,contact_dynamics;do_record=false)
    (;st, traj) = bot
    fig = Figure(resolution=(1280,720))
    ax = fig[1,1] = LScene(fig, scenekw = (show_axis = false, camera = cam3d!, raw = false))
    cube_bars = vis!(ax,st,get_cube_bars(st))
    # cameracontrols(ax.scene).rotationspeed[] = 0.01
    if do_record
        lookat = [5.520309, -10.742862, 0.6962326]
        eyepos = [17.465916, -22.283669, 3.495206]
        framerate = 30
        h = 1e-3
        record(fig, "cube.mp4", 1:round(Int,1/(30*1)/h):length(qs); framerate) do this_step
            q = traj.qs[this_step]
            update_scene!(ax,st,cube_bars,q)
            update_cam!(ax.scene, eyepos, lookat)
        end
    else
        sliderplot(fig,st,bot,ax,cube_bars,contact_dynamics)
    end
    fig
end
