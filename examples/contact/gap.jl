
function gap2(rb)
    @unpack state = rb
    @unpack R,r = state
    a = 34
    h = 69
    i = @view R[:,1]
    j = @view R[:,2]
    k = @view R[:,3]

    iz = i[3]
    jz = j[3]
    dom = √(iz^2+jz^2)
    cosϕ = iz/(dom+eps(iz))
    sinϕ = jz/(dom+eps(jz))
    ϕ = atan(sinϕ,cosϕ)
    #ϕ = acos(cosϕ)
    θ = ϕ + π
    @show cosϕ
    rc2 = r + h*k + a*i*cos(θ) + a*j*sin(θ)
    ret = rc2[3]
    rc2,ret
end

rbs = [rb1]

gap2(rb1)
scene = Scene()
function update_rb_position!(scene,rbid,rb)
    for (id,point) = enumerate(rb.state.p)
        point .= rb.state.r + rb.state.R*rb.prop.anchorpoints[id].p
    end
    qtn = Rotations.Quat(rb.state.R)
    rotate!(scene[rbid],AbstractPlotting.Quaternion(qtn.x,qtn.y,qtn.z,qtn.w))
    translate!(scene[rbid],rb.state.r)
    translate!(scene[rbid+1],rb.state.p[4])
    rc2,g2 = gap2(rb)
    translate!(scene[rbid+2],rc2)
end
function draw_rbs!(scene,rbs)
    for (rbid,rb) in enumerate(rbs)
        mesh!(scene,linkmesh(34.0,69.0),color=:grey)
        mesh!(scene,Sphere(Point3(rb.state.p[4]...), 4.0), color=:cyan, show_axis=false)
        mesh!(scene,Sphere(Point3(rb.state.p[4]...), 4.0), color=:cyan, show_axis=false)

        update_rb_position!(scene,rbid,rb)
    end
end
draw_rbs!(scene,rbs)
scene
x_slider,x = textslider(-50e-1:1e-4:50e-1,"x",start=0.0)
y_slider,y = textslider(-50e-1:1e-4:50e-1,"y",start=0.0)
z_slider,z = textslider(-69e-1:1e-4:69e-1,"z",start=0.0)
a_slider,a = textslider(0:0.1:2π,"a",start=0.0)
b_slider,b = textslider(0:0.1:2π,"b",start=0.0)
c_slider,c = textslider(0:0.1:2π,"c",start=0.0)

rbid = 2
on(x) do this_x
    rb1.state.r[1] = this_x
    update_rb_position!(scene,rbid,rb1)
end
on(y) do this_y
    rb1.state.r[2] = this_y
    update_rb_position!(scene,rbid,rb1)
end
on(z) do this_z
    rb1.state.r[3] = this_z
    update_rb_position!(scene,rbid,rb1)
end
on(a) do this_a
    R = RotXYZ(rb1.state.R)
    R_a = R.theta1
    R_b = R.theta2
    R_c = R.theta3
    rb1.state.R .= RotXYZ(this_a,R_b,R_c)
    update_rb_position!(scene,rbid,rb1)
end
on(b) do this_b
    R = RotXYZ(rb1.state.R)
    R_a = R.theta1
    R_b = R.theta2
    R_c = R.theta3
    @show R_a,R_b,R_c
    rb1.state.R .= RotXYZ(R_a,this_b,R_c)
    update_rb_position!(scene,rbid,rb1)
end
on(c) do this_c
    R = RotXYZ(rb1.state.R)
    R_a = R.theta1
    R_b = R.theta2
    R_c = R.theta3
    rb1.state.R .= RotXYZ(R_a,R_b,this_c)
    update_rb_position!(scene,rbid,rb1)
end
vbox(scene,hbox(z_slider,y_slider,x_slider,c_slider,b_slider,a_slider))
