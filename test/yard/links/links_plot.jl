function torus(R=0.5,r=0.02,ns=36,ne=8)
    coords = Matrix{Float64}(undef,3,ne*ns)
    θ = 2π/ne
    α = 2π/ns
    for j = 0:ns-1
        for i = 0:ne-1
            secxyz = [r*cos(i*θ) + R, 0.0, r*sin(i*θ)]
            coords[:,j*ne+i+1] = RotZ(j*α)*(secxyz)
        end
    end
    function connect2sec(a,b)
        around = Matrix{Int64}(undef,3,2length(a))
        apair = hcat([a[1],b[1],a[2]],[a[2],b[1],b[2]])
        around[:,1:2] = apair
        for i = 1:length(a)-2
            nextpair = apair .+ i
            around[:,2i+1:2i+2] = nextpair
        end
        around[:,end-1:end] = hcat([a[end],b[end],a[1]],[a[1],b[end],b[1]])
        around
    end
    faces = Matrix{Int64}(undef,3,2ne*ns)
    for i = (0:ns-2)
        faces[:,i*2ne+1:(i+1)*2ne] = connect2sec((1:ne).+i*ne,(ne+1:2ne).+i*ne)
    end
    faces[:,end-2ne+1:end] = connect2sec((1:ne).+ne*(ns-1),1:ne)
    [Point3(coords[:,i]) for i = 1:ne*ns], [Face{3}(faces[:,i]) for i = 1:2ne*ns]
end

function linkmesh(r=0.5,h=1.0,d=r/10)
    toruspoints,torusfaces = torus(r,d/2)
    offset_toruspoints = [Point(p[1],p[2],p[3]+h) for p in toruspoints]
    torusmesh = GLNormalMesh(offset_toruspoints,torusfaces)
    αs = [0.0,2π/3,-2π/3]
    ri = Point3f0(0., 0., 0.)
    rk = [Point3f0(r*cos(α), r*sin(α), h) for α in αs]
    meshbar = [GLNormalMesh(Cylinder{3, Float32}(rk[i],ri,Float32(d/2)),20) for i = 1:3]
    linkmesh = merge(meshbar...,torusmesh)
end
linkmesh(0.5,1.0)
set_theme!(size = (1280,720),center = false)
scene = Scene()
function update_rb_position!(scene,bodyid,rb)
    qtn = Rotations.Quat(body.state.R)
    rotate!(scene[bodyid],AbstractPlotting.Quaternion(qtn.x,qtn.y,qtn.z,qtn.w))
    translate!(scene[bodyid],body.state.loci_states)
end
function draw_rbs!(scene,rbs)
    for (bodyid,rb) in enumerate(rbs)
        mesh!(scene,linkmesh(0.5,1.0),color=:grey, show_axis=false)
        # mesh!(scene, Sphere(Point3(tglink.NC.r̄1), 0.04), color=:cyan, show_axis=false)
        # mesh!(scene, Sphere(Point3(tglink.NC.r̄2), 0.04), color=:cyan, show_axis=false)
        # mesh!(scene, Sphere(Point3(tglink.NC.r̄3), 0.04), color=:cyan, show_axis=false)
        # mesh!(scene, Sphere(Point3(tglink.NC.r̄4), 0.04), color=:cyan, show_axis=false)

        update_rb_position!(scene,bodyid,rb)
    end
end
draw_rbs!(scene,rbs)

function update_rbs_position!(scene,rbs)
    for (bodyid,rb) in enumerate(rbs)
        update_rb_position!(scene,bodyid,rb)
    end
end
update_rbs_position!(scene,rbs)
function string2linesegment(st)
    rbida,pida = connectivity[st.ida]
    rbidb,pidb = connectivity[st.idb]
    rba = rbs[rbida]
    rbb = rbs[rbidb]
    ra = rba.state.p[pida]
    body = rbb.state.p[pidb]
    ls = Point3(ra...) => Point3(body...)
end


function draw_strings!(scene,sts)
    lss = [string2linesegment(st) for st in sts]
    linesegments!(scene, lss, linewidth = 4.0, color = :magenta, show_axis = false)
end

draw_strings!(scene,sts)

<<<<<<< Updated upstream
function update_strings!(scene,sts)
    lss = [string2linesegment(st) for st in sts]
    scene.plots[end].input_args[1][] = lss
end
update_strings!(scene,sts)
=======
function update_apparatuses!(scene,sts)
    lss = [string2linesegment(st) for st in sts]
    scene.plots[end].input_args[1][] = lss
end
update_apparatuses!(scene,sts)
>>>>>>> Stashed changes

function update_scene!(scene,tgsys,q,q̇)
    TRS.reset_forces!.(tgsys.rigidbodies)
    TRS.distribute_q_to_rbs!(tgsys.rigidbodies.movables,q,q̇)
    TRS.compute_string_forces!(tgsys)
    update_rbs_position!(scene,tgsys.rigidbodies)
<<<<<<< Updated upstream
    update_strings!(scene,tgsys.strings)
=======
    update_apparatuses!(scene,tgsys.cables)
>>>>>>> Stashed changes
end
update_scene!(scene,tgsys,state.qs[1],state.q̇s[1])

step_slider,tstep = textslider(1:length(state.ts),"step",start=1)
on(tstep) do this_step
    update_scene!(scene,tgsys,state.qs[this_step],state.q̇s[this_step])
end
bigscene = vbox(step_slider,scene)
rotate_cam!(scene,0.0,0.0,0.0)
record(vbox(step_slider,scene), "output.mp4", 1:10:length(state.ts)) do i
     tstep[] = i
end
