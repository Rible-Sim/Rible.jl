rbs = [cube]
set_theme!(resolution = (1280,720),center = false)
scene = Scene()
cube.state.r .= 0.0
function update_rb_position!(scene,rbid,rb)
    qtn = Rotations.Quat(rb.state.R)
    rotate!(scene[rbid],AbstractPlotting.Quaternion(qtn.x,qtn.y,qtn.z,qtn.w))
    translate!(scene[rbid],rb.state.r)
end
function draw_rbs!(scene,rbs)
    for (rbid,rb) in enumerate(rbs)
        mesh!(scene,HyperRectangle(Vec3f0(-0.5), Vec3f0(1)),color=:grey, show_axis=false)
        update_rb_position!(scene,rbid,rb)
    end
end
draw_rbs!(scene,rbs)

function update_scene!(scene,rbs,q,q̇)
    TRS.distribute_q_to_rbs!(rbs,q,q̇)
    update_rb_position!(scene,1,cube)
end
update_scene!(scene,rbs,qs[1],q̇s[1])

step_slider,step = textslider(1:length(ts),"step",start=1)
on(step) do this_step
    update_scene!(scene,rbs,qs[this_step],q̇s[this_step])
    #ke,pe,en = energy(cube,qs[this_step],q̇s[this_step])
    #@show ke,pe,en
end
box = hbox(step_slider,scene)

record(scene, "output.mp4", 1:length(ts), framerate = 60) do i
     step[] = i
end
