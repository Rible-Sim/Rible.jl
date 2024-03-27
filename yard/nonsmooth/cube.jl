
#-- cube
cubebot = make_cube()
plot_traj!(cubebot)

function cube_contact_dynfuncs(bot)
contact_dynfuncs(bot;
    flatplane = RB.Plane([0,0,1.0],[0,0,0.0])
)
end

ros = [[0,0,z] for z = 10.0:-0.1:0.0]
qs = [
begin
    cube = make_cube(origin_position,RotXY(π/4,π/4))
    q = RB.get_coords(cube.st)
end
for origin_position in ros
]
# p = RB.generalized_alpha(1.0)
Ro = [
0.368869  -0.824063   0.429949;
0.769524   0.011314  -0.638518;
0.521314   0.566385   0.63831
]
origin_velocity = [0.0,0.0,0.0]
ωo = [3.0,4.0,5.0]
# cube = make_cube()
# cube = make_cube(ros[1])
cube = make_cube(ros[2],Ro,origin_velocity,0.1.*ωo)
tspan = (0.0,1.91)
tspan = (0.0,5.0)
dt = 1e-3

prob = RB.DynamicsProblem(cube,cube_contact_dynfuncs)
RB.solve!(prob,
RB.ZhongCCP();tspan,dt,ftol=1e-10,maxiters=50,exception=false
)

plot_traj!(cube;
rigidcolor=:white,
xlims = [-10,10],
ylims = [-10,10],
zlims = [-1e-3,10],
showinfo=false,
showwire=true
)

ME = RB.mechanical_energy!(cube)
lines(cube.traj.t,ME.E)

ts,qs,vs,as,ṽ̇s = RB.NSSFC.nssfc(nq,nλ,q0,v0,p,h,contact_dynfuncs(cube),tspan;tol=1e-6,imax=10)

ts,cs,qs,vs,ps,λs,friction_coefficients = RB.NSSFC.nhsolve(nq,nλ,nμ,q0,v0,contact_dynfuncs(cube),tspan;
        dt=h,ftol=1e-6,imax=50,exception=true)

ts,cs,qs,vs,ps,λs,friction_coefficients = RB.NSSFC.ipsolve(nq,nλ,nμ,q0,v0,contact_dynfuncs(cube),tspan;
        dt=h,ftol=1e-10,imax=50,exception=false)

ts,qs,vs,as,ṽ̇s = RB.NSSFC.nssfc(n,b,q26026,v26026,p,h,contact_dynfuncs(cube),tspan;tol=1e-10,imax=100)

cuberef = deepcopy(cube)
prob = RB.DynamicsProblem(cuberef,dynfuncs,tspan)
RB.solve!(prob,RB.Zhong06();dt=h,ftol=1e-14)

vis(cuberef,contact_dynamics)

cube = make_cube(ros[1],Ro,origin_velocity,ωo)
# cube = make_cube(ros[1])
# cube = make_cube()

RB.prepare_traj!(cube.traj;tspan,dt=h,restart=true)

cube.traj.t[2:end] = ts[2:end]
cube.traj.q[2:end] = qs[2:end]
cube.traj.q̇[2:end] = vs[2:end]

vis(cube,contact_dynamics)

vis(cube,contact_dynamics;do_record=true)

plot(cube.traj.ts,VectorOfArray(cube.traj.qs)[7,:],color="yellow")
plot!(cuberef.traj.ts,VectorOfArray(cuberef.traj.qs)[7,:],color="blue")

ME = RB.mechanical_energy!(cube;gravity=true,)

ME.E
scatter(ME.E)
kes,epes,gpes,restitution_coefficients,es_err = analyse_energy(cube;gravity=true,elasticity=false)
es_off = OffsetArray(restitution_coefficients,0:length(restitution_coefficients)-1)
kes_off = OffsetArray(kes,0:length(restitution_coefficients)-1)
gpes_off = OffsetArray(gpes,0:length(restitution_coefficients)-1)
steprange = 0:length(restitution_coefficients)-1
steprange = 558:560
steprange = 5574:5600
steprange = 55740:56000
es1 = es_off[steprange[1]]
plot(steprange,(es_off[steprange].-es1)./es1)
plot(steprange,kes_off[steprange])
plot(steprange,gpes_off[steprange])
smoothsteps = findall((x)->x==0,cs)
plotsmoothsteps = [x for x in smoothsteps if x∈steprange]
plot!(plotsmoothsteps,es_off[plotsmoothsteps],color=:blue)

-0.0016049537616469219-(-0.0016048102389008441)

plot(epes[6400:6500])
ylims!(78,80)

plot(ts,kes[begin+1:end])
plot(ts,gpes[begin+1:end])

ros = [q[1:3] for q in qs]
ṙos = [q̇[1:3] for q̇ in vs]
r̈os = [ṽ̇[1:3] for ṽ̇ in ṽ̇s]
aos = [a[1:3] for a in as]

plot(ts,norm.(ros))

plot(ts,norm.(ṙos))

plot!(ts,norm.(ros))

plot(ts,norm.(r̈os))


plot(ts,norm.(aos))

plot!(ts,gpes[begin+1:end])

norm(ṙos[1393])/norm(ṙos[1392])

y_split = [[0.09510319720143898, -0.008672570657933804, 0.09470694080223288], [0.23134216351905035, 0.1781735286559687, 0.10655619736781423], [0.25492807263830397, 0.18779000020080971, -0.14810414657700796]]
[yi[1] - norm(yi[2:3]) for yi in y_split]
[norm(yi[2:3]) for yi in y_split]
Λ_split = [[6.077211920457062e-5, 5.541879908804313e-6, -6.051890646476349e-5], [1.6851978930926635e-17, -1.4431551427742453e-17, -7.662097843905011e-18], [1.0182657945568905e-17, -7.963039923870281e-18, 5.867035717865434e-18]]
[Λi[1] - norm(Λi[2:3]) for Λi in Λ_split]
[Λi[1] for Λi in Λ_split]
[norm(Λi[2:3]) for Λi in Λ_split]
