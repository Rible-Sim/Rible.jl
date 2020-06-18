using Makie
using GeometryTypes

function coords2state(rb)
    @unpack mass,inertia = rb.prop
    @unpack r,ṙ,R,ω,coords = rb.state
    @unpack x,q,p,L = coords
    r .= x
    ṙ .= p/mass
    R .= SMatrix(Quat(q[1],q[2],q[3],q[4]))
    Jt = R*inertia*transpose(R)
    ω .= inv(Jt)*L
end

function sol2state!(tgsys,sol)
    rbs = tgsys.rigidbodies
    sts = tgsys.strings
    cnt = tgsys.connectivity
    nc = 13
    for i in eachindex(rbs)
        @unpack x,q,p,L = rbs[i].state.coords
        coords_sol = @view sol[(i-1)*nc+1:(i-1)*nc+13]
        x .= coords_sol[1:3]
        q .= coords_sol[4:7]
        p .= coords_sol[8:10]
        L .= coords_sol[11:13]
        coords2state(rbs[i])
        for ip in eachindex(rbs[i].prop.anchorpoints)
            rbs[i].state.p[ip] .= point_position(rbs[i].prop,ip,
                                    rbs[i].state.r,rbs[i].state.R)
        end
    end
    for i in eachindex(sts)
        st = sts[i]
        rbid1,pid1 = cnt[st.pid1]
        rbid2,pid2 = cnt[st.pid2]
        rb1 = rbs[rbid1]
        rb2 = rbs[rbid2]
        st.𝐬 .= rb2.state.p[pid2] - rb1.state.p[pid1]
    end
end
#sol2state!(tgsys,sol[:,1])
sol = solution
sol2state!(tgsys,sol[1])
color_list = [RGBAf0(1,0,0,1),RGBAf0(0,1,0,1),RGBAf0(0,0,1,1)]
scene = Scene(show_axis=false,limits=FRect3D((-100., -100.,-100.),(200., 200., 200.)))
function initialplot(s,rbs)
    for i in eachindex(rbs)
        box = HyperRectangle(Vec3f0(-0.25,-0.25,-1.0), Vec3f0(0.5,0.5,2.0))
        #ball = HyperSphere(Point(0.0,0.0,0.0),rbs[i].object.shape.radius)
        #rbmesh = mesh!(scene,box,color = RGBAf0(1,0,0,1),show_axis = false, center=false)[end]
        rbmesh = mesh!(scene,box,color = color_list[i],show_axis = false, center=false)[end]
        x = rbs[i].state.r
        translate!(rbmesh,x...)
        q = rbs[i].state.coords.q
        rotate!(rbmesh,Quaternion(q[2],q[3],q[4],q[1]))
    end
    #ground = HyperRectangle(Vec(-10.,-10.,-0.1), Vec(20.,20.,0.1))
    #gmesh = mesh!(scene,ground,color = RGBAf0(0,1,0,1),show_axis = false, center=false)[end]
end

initialplot(scene,tgsys.rigidbodies)
scene
using AbstractPlotting
tslider,tt = textslider(1:length(sol), "t",start = 1)
on(tt) do t
    func(t)
end

function func(t)
    sol2state!(tgsys,sol[t])
    rbs = tgsys.rigidbodies
    for i in eachindex(rbs)
        rbmesh = scene.plots[i]
        x = rbs[i].state.r
        translate!(rbmesh,x...)
        q = rbs[i].state.coords.q
        rotate!(rbmesh,Quaternion(q[2],q[3],q[4],q[1]))
    end
end

record(scene, "video.mp4", 1:200) do i
    func(i) # or some other manipulation of the Scene
end
hbox(scene,tslider)



boxmesh = GLNormalMesh((box,RGBAf0(1,0,0,1)))

scene = Scene(show_axis=false)
box = HyperRectangle(Vec3f0(-0.25,-0.25,-1.0), Vec3f0(0.5,0.5,2.0))

boxscene1 = mesh!(scene,box,color = RGBAf0(1,0,0,1),show_axis = false, center=false)[end]

boxscene2 = mesh!(scene,box,color = RGBAf0(0,1,0,1),show_axis = false, center=false)[end]

xslider,x = textslider(0.0:0.1:1.0, "x",start = 0.0)
yslider,y = textslider(0.0:0.1:1.0, "y",start = 0.0)
zslider,z = textslider(0.0:0.1:1.0, "z",start = 0.0)
on(x) do xx
    translate!(boxscene1,xx,0.0,0.0)
end
on(y) do yy
    translate!(boxscene2,0.0,yy,0.0)
end
on(z) do zz
    translate!(boxscene1,0.0,0.0,zz)
end
rxslider,rx = textslider(0.0:0.1:π, "rx",start = 0.0)
ryslider,ry = textslider(0.0:0.1:π, "ry",start = 0.0)
rzslider,rz = textslider(0.0:0.1:π, "rz",start = 0.0)
on(rx) do rxx
    q = Quat(RotX(rxx))
    rotate!(boxscene1,Quaternion(q.x,q.y,q.z,q.w))
end
on(ry) do ryy
    q = Quat(RotY(ryy))
    rotate!(boxscene2,Quaternion(q.x,q.y,q.z,q.w))
end
on(rz) do rzz
    q = Quat(RotZ(rzz))
    rotate!(boxscene1,Quaternion(q.x,q.y,q.z,q.w))
end
hbox(scene,vbox(xslider,yslider,zslider,rxslider,ryslider,rzslider),
    parent = Scene(resolution = (1000, 500)))

scene
