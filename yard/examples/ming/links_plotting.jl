function rb_bars(body)
    ps = body.state.p
    cs = [[1,2,3],[4,7,10],[5,8,11],[6,9,12],[10,11,12]]
    vcat(collect(Base.Iterators.flatten([[Point(ps[c[1]]) => Point(ps[c[2]]),
            Point(ps[c[2]]) => Point(ps[c[3]]),
            Point(ps[c[3]]) => Point(ps[c[1]])] for c in cs])),
            [Point(ps[i]) => Point(ps[i+9]) for i in 1:3])
end

function read_rb_mesh()
    pyramid = load("zlzt_1_1.obj")
    pyramid_coordinates = coordinates(pyramid)
    pyramid_pure_coordinates = metafree(pyramid_coordinates)
    pyramid_normals = normals(pyramid)
    pyramid_pure_coordinates_new = map((x)->RotX(π)*x*0.001,pyramid_pure_coordinates)
    pyramid_normals_new = Ref(RotX(π)).*pyramid_normals
    pyramid_faces = faces(pyramid)
    new_pyramid = GeometryBasics.Mesh(meta(pyramid_pure_coordinates_new;normals=pyramid_normals_new),pyramid_faces)
end

function bars_and_strings(st)
    bars = [Node(rb_bars(body)
            ) for rb in st.rigidbodies]
    strings = Node(RB.get_strings(st))
    bars,strings
end

function pyramids_and_strings(st,pyramid)
    pyramids = [pyramid for rb in st.rigidbodies]
    strings = Node(RB.get_strings(st))
    pyramids,strings
end

function get_spheres(st)
    r = 0.002
    [Makie.Sphere(Point3f0(rp), r) for rb in st.rigidbodies for rp in body.state.loci_states]
end

function pyramids_strings_and_spheres(st,pyramid)
    pyramids,strings = pyramids_and_strings(st,pyramid)
    spheres = Node.(get_spheres(st))
    pyramids,strings,spheres
end

function update_rb_mesh!(meshobj,rb)
    o = [0,0,0]
    n1 = [1,0,0]
    n2 = [0,1,0]
    n3 = [0,0,1]
    @unpack c,C = body.state.cache.funcs
    q = body.state.coords.q
    r = C(c(o))*q
    R1 = C(c(n1))*q-r
    R2 = C(c(n2))*q-r
    R3 = C(c(n3))*q-r
    R = hcat(R1,R2,R3)
    Makie.translate!(meshobj,r...)
    uq = UnitQuaternion(RotMatrix(R))
    Makie.rotate!(meshobj,Quaternion(uq.x,uq.y,uq.z,uq.w))
end

function plotstructure!(ax,st,pyramid)
    # pyramids,strings = pyramids_and_strings(st,pyramid)
    pyramids,strings,spheres = pyramids_strings_and_spheres(st,pyramid)
    #return bars,strings

    function plot!(ax,pyramids,strings)
        for (lineid,pyramid) in enumerate(pyramids)
            meshobj = mesh!(ax,pyramid, color = :grey)
            update_rb_mesh!(meshobj,st.rigidbodies[lineid])
        end
        for sphere in spheres
            mesh!(ax,sphere,color=:orange)
        end
        linesegments!(ax,strings, color = :deepskyblue, linewidth = 2)
    end

    plot!(ax,pyramids,strings)
    # cameracontrols(ax.scene).rotationspeed[] = 0.01
    # xlims!(ax.scene,-0.1,0.1)
    # ylims!(ax.scene,-0.1,0.1)
    # zlims!(ax.scene,-0.05,0.25)
    pyramids,strings,spheres
end

function plotstructure(st::RB.Structure,pyramid::GeometryBasics.Mesh)
    fig = Figure(resolution=(1000,1000))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    pyramids,strings,spheres = plotstructure!(ax,st,pyramid)
    cameracontrols(ax.scene).rotationspeed[] = 0.01
    fig,ax
end

function recordplot(fig,tr,ax,strings,eyepos,lookat;record_name)
    @unpack st, traj = tr
    record(fig, record_name, 1:length(traj.ts); framerate = 30) do this_step
        update_scene!(st,ax,strings,traj.qs[this_step])
        update_cam!(ax.scene, eyepos,lookat)
    end
end
function show_camera(scene)
    camera = cameracontrols(scene)
    lookat = camera.lookat[]
    eyepos = camera.eyeposition[]
    @show eyepos, lookat
end

function sliderplot(fig,tr,ax,strings,spheres)
    @unpack st, traj = tr
    ls_step = labelslider!(fig,"step",1:length(traj.ts))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        # analyse_slackness(st,sol.qs[this_step])
        update_scene!(st,ax,strings,spheres,traj.qs[this_step])
        show_camera(ax.scene)
    end
end

function update_scene!(st,ax,strings,spheres,q)
    cnt = st.connectivity
    RB.reset_forces!(st)
    RB.distribute_q_to_rbs!(st,q)
    RB.update_strings_apply_forces!(st)
    for (id,rb) in enumerate(st.rigidbodies)
        update_rb_mesh!(ax.scene[1+id],rb)
    end
    strings[] = RB.get_strings(st)
    new_spheres = get_spheres(st)
    foreach((sphere,new_sphere)->sphere[]=new_sphere,spheres,new_spheres)
    # angles = update_angles(st)
    # @show angles
end

function plotstructure(tr::RB.Robot,pyramid;record_name="")
    @unpack st = tr
    fig = Figure(resolution=(1280,720))
    ax = fig[1, 1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    _,strings,spheres = plotstructure!(ax,st,pyramid)
    # cameracontrols(ax.scene).rotationspeed[] = 0.01
    eyepos = Vec3f(0.20397526, 0.5129121, 0.17975415)
    lookat = Vec3f(0.034443498, 0.0022780984, 0.0027566003)
    update_cam!(ax.scene, eyepos, lookat)
    if record_name==""
        sliderplot(fig,tr,ax,strings,spheres)
    else
        recordplot(fig,tr,ax,strings,eyepos,lookat;record_name)
    end
    fig
end

function plotstructure(trs::Vector{<:RB.Robot},nx,ny,pyramid::GeometryBasics.Mesh)
    fig = Figure(resolution=(1280,720))
    eyepos = Vec3f(0.20397526, 0.5129121, 0.17975415)
    lookat = Vec3f(0.034443498, 0.0022780984, 0.0027566003)
    ntr = length(trs)
    @assert ntr <= nx*ny
    i = 0
    for ix = 1:nx
        for iy = 1:ny
            ax = fig[ix, iy] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
            i += 1
            st = trs[i].st
            plotstructure!(ax,st,pyramid)
            cameracontrols(ax.scene).rotationspeed[] = 0.01
            update_cam!(ax.scene, eyepos, lookat)
        end
    end
    fig
end
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
        bar_color_name = "orange"
        string_color_name = "orange"
    else
        bar_color_name = "black"
        string_color_name = "deepskyblue"
    end
    bar_color = plt.matplotlib.colors.to_rgba(bar_color_name)
    string_color = plt.matplotlib.colors.to_rgba(string_color_name)
    bars_segs = plt.art3D.Line3DCollection(bars_lines,linewidths=2,colors=bar_color)
    strings_segs = plt.art3D.Line3DCollection(strings_lines,linewidths=1,colors=string_color)
    bars_segs, strings_segs
end

function get_trajectory(bodyid,pid,st,sol,step_range=:)
    rp = VectorOfArray(Vector{Vector{Float64}}())
    for q in sol.qs
        RB.distribute_q_to_rbs!(st,q)
        push!(rp,st.rigidbodies[bodyid].state.p[pid])
    end
    rp
end
