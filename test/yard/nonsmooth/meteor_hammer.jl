
#-- meteor hammer

cablemesh = sample(ancs,e,1000)

mesh(cablemesh,transparency=false)

function make_hammer(id,r̄ijkl,origin_position,R,ri,rj=nothing,rk=nothing,rl=nothing;
        μ,e,
        ω = zero(origin_position),
        contactable = true,
        visible = true,
        pres_idx = Int[],
        cstr_idx = collect(1:6),
    )
    # free_idx = collect(1:6)
    m = 4.58794901
    Īg = SMatrix{3,3}(
        Matrix(Diagonal([
            0.00522676,
            0.00512073,
            0.00512073,
            ])
        )
    )
    mass_locus = SVector(-0.20,0,0)
    # loci = deepcopy(r̄ijkl)
    loci = vcat(
        [
            SVector(0.22*cos(i*π/2),0.0,0.22*sin(i*π/2)) + mass_locus
            for i = -2:1
        ],
        [
            RotY( π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + mass_locus
            for i = 1:6
        ],
        [
            RotY(-π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + mass_locus
            for i = [1,2,4,5]
        ]
    )
    axes = [
        SVector(1.0,0,0) 
        for _ in eachindex(loci)
    ]
    friction_coefficients = fill(μ,length(loci))
    restitution_coefficients = fill(e,length(loci))
    @show m,diag(Īg),mass_locus,length(loci),friction_coefficients,restitution_coefficients
    prop = RB.RigidBodyProperty(
                id,contactable,m,Īg,
                mass_locus,loci,axes,
                friction_coefficients,restitution_coefficients;
                visible=visible
                )
    origin_velocity = zero(origin_position)
    if rj isa Nothing
        nmcs = RB.NCF.NC1P3V(ri,origin_position,R,origin_velocity,ω)
    elseif rk isa Nothing
        nmcs = RB.NCF.NC2P2V(ri,rj,origin_position,R,origin_velocity,ω)
    elseif rl isa Nothing
        nmcs = RB.NCF.NC3P1V(ri,rj,rk,origin_position,R,origin_velocity,ω)
    else
        nmcs = RB.NCF.NC4P(ri,rj,rk,rl,origin_position,R,origin_velocity,ω)
    end
    # cf = RB.NCF.CoordinateFunctions(nmcs,q0,pres_idx,free_idx)
    # @show typeof(nmcs)
       # trimesh = Meshes.discretize(box,) |> simple2mesh
    meteormesh = load(joinpath(RB.assetpath("流星锤.STL"))) |> make_patch(;
        scale = 1/400,
        trans = [-0.20,0,0],
    )
    state = RB.RigidBodyState(prop,nmcs,origin_position,R,origin_velocity,ω,pres_idx,cstr_idx)
    body = RB.RigidBody(prop,state,meteormesh)
end

function cable_ancf(pres_idx, 𝐞, L = 1.0) 
    radius = 2.0e-3
    # ancs = RB.ANCF.ANC3DRURU(8.96e3;E=110e6,L,radius)
    # mat_cable = filter(
    #     row->row.name == "Nylon 66",
    #     material_properties
    # )[1]
    # ρ = ustrip(Unitful.kg/Unitful.m^3,mat_cable.density)
    # E = ustrip(Unitful.Pa,mat_cable.modulus_elas)
    ρ = 1.03e3
    E = 0.2e9
    @myshow ρ, E
    ancs = RB.ANCF.ANC3DRURU(ρ;E,L,radius)
    mass = RB.ANCF.build_mass(ancs)
    @show mass
    T = typeof(L)
    mass_locus = SVector(L/2,T(0),T(0))
    r̄p1 = SVector(T(0),T(0),T(0))
    # r̄p2 = SVector(L,T(0),T(0))
    loci = [
        r̄p1,
    ]
    prop = RB.FlexibleBodyProperty(
        1,
        :cable,
        mass,
        mass_locus,
        # length(loci),
        loci
    )
    # cache = RB.BodyCache(prop,ancs,𝐞)
    state = RB.FlexibleBodyState(prop,ancs,𝐞;pres_idx)
    fb = RB.FlexibleBody(prop,state)
end

function make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.5*√3, 0.0, 0.5-1e-2],
        rjx = [ 0.0,-1.0, 0.0],
        doDR=false,
        μ = 0.5,
        e = 0.9,
        L = 1.0,
        nx = 2,
        R = RotY(deg2rad(-60)),
        ω = zero(ri)
    )
    if doDR
        fb_pres_idx = [7,8,9]
        rb_pres_idx = [1,2,3]
    else
        fb_pres_idx = Int[]
        rb_pres_idx = Int[]
    end
    𝐞 = [ri;rix;rj;rjx]
    fb1 = cable_ancf(fb_pres_idx, 𝐞, L)
    # @show fb1.prop.mass
    subfbs,subsm = RB.subdivide(fb1,nx)
    r̄ijkl = SVector{3,Float64}.(
        0.1 .*[
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
        ]
    )
    rigidbody = make_hammer(
        nx+1,
        r̄ijkl,
        rj,
        R,
        rj;
        μ,e,
        ω,
        pres_idx=rb_pres_idx,
        visible=ifelse(!isempty(rb_pres_idx),true,false)
    )
    fbs = TypeSortedCollection(vcat(subfbs,rigidbody))
    # fbs = TypeSortedCollection([fb1,])
    numbered = RB.number(fbs)
    # indexed = RB.index(fbs,)
    sm = zeros(Int,size(subsm,1)+3,size(subsm,2)+1)
    sm[1:size(subsm,1),1:size(subsm,2)] = subsm
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)  ] = 7:9
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)+1] = 1:3
    # display(sm)
    indexed = RB.index(fbs,sm)
    ss = Int[]
    apparatuses = (cables = ss,)
    connected = RB.connect(fbs,)
    tensioned = @eponymtuple(connected,)
    cst1 = RB.FixedIndicesConstraint(1,[1,2,3],ri)
    jointed = RB.join((cst1,),indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(fbs,apparatuses,cnt,)
    bot = RB.Robot(st)
end

function flexcable_contact_dynfuncs(bot,ground_plane)
    contact_dynfuncs(bot;
        flatplane = ground_plane
    )
end

# don't bother
R = RotXY(deg2rad(0),deg2rad(-0))
inclined_plane = RB.Plane(R*[0,0,1.0],[-1.039230,-0.0,0.781])
p1 =  [-1.039230,0,0.780]
p2 = [-1.173952,-0.11,0.8646278]
p3 = [-0.901447,-0.112500,0.762216]
n = (p2-p1) × (p3 - p1 )

inclined_plane = RB.Plane(n,p3.+[0,0,0.001])


flexcable_DR = make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
RB.GDR!(flexcable_DR;β=2e-4,maxiters=2e5,verbose=false)
flexcable = make_flexcable(;
    ri  = [ 0.0, 0.0, 1.5],
    rix = [ 0.0, 0.0,-1.0],
    rj  = [-0.6*√3, 0.0, 1.2],
    rjx = [ 0.0,-1.0, 0.0],
    e = 0.0,
    μ=0.1,L=1.2,nx=5,R=RotY(-π/2),
    ω=-200n,
)

# flexcable_DR = make_flexcable(;R = RotY(-π/2),L=1.5,nx=5,doDR=true)
# RB.GDR!(flexcable_DR;β=1e-3,maxiters=1e5)
# flexcable = make_flexcable(;R = RotY(-π/2),e=0.0,μ=0.01,L=1.5,nx=5)

flexcable_DR.traj.q[end][end-8:end] .= flexcable.traj.q[end][end-8:end]
flexcable_DR.traj.q̇[end][end-8:end] .= flexcable.traj.q̇[end][end-8:end]

RB.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])


with_theme(theme_pub;
        Poly= (
            transparency = true,
        ) 
    ) do 
    plot_traj!(
        flexcable;
        ground=inclined_plane,
        xlims=(-1.2,1.0),
        ylims=(-0.5,0.5),
        zlims=(-0.4,1.8),
    )
end

# sliding and avg err 
tspan = (0.0,1.5)
h = 2e-4

RB.solve!(
    RB.DynamicsProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
    RB.ZhongCCP();
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

plotsave_contactpoints(flexcable,)

rp2 = RB.get_trajectory!(flexcable,6,10)
lines(rp2[2,:])

me = RB.mechanical_energy!(flexcable)
lines((me.E.-me.E[begin])./me.E[begin])

#dt
# dts = [1e-1,3e-2,1e-2,3e-3,1e-3,1e-4]
dts = [1e-2,5e-3,2e-3,0.99e-3,5e-4,2e-4,2e-5]
flexcables_dt = [
    begin
        flexcable_dt = deepcopy(flexcable)
        RB.solve!(
            RB.DynamicsProblem(flexcable_dt,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
                RB.ZhongCCP();
                tspan=(0.0,0.12),dt,ftol=1e-14,maxiters=50,exception=false)
    end
    for dt in dts
]

_, err_avg = get_err_avg(flexcables_dt;bid=6,pid=10,di=1)

with_theme(theme_pub;
        size = (0.9tw,0.45tw),
        figure_padding = (fontsize,fontsize,0,0),
        Axis3=(
            azimuth = 7.595530633326987,
            elevation = 0.14269908169872403
        )
    ) do
    bot = flexcable
    (;t) = bot.traj
    bots = flexcables_dt
    rp1 = RB.get_trajectory!(bot,6,1)
    rp12 = RB.get_trajectory!(bot,6,12)
    rp10 = RB.get_trajectory!(bot,6,10)
    fig = Figure()
    gd2 = fig[1,2] = GridLayout()
    gd3 = fig[2,2] = GridLayout()
    gd1 = fig[:,1] = GridLayout()
    steps = 1:1000:length(bot.traj)
    laststep = last(steps)
    cg = cgrad(:winter, length(steps), categorical = true, rev = true)
    nstep = length(steps)
    alphas = fill(0.2,nstep)
    alphas[1:3] = [1, 0.2, 0.2]
    alphas[end] = 1
    plot_traj!(
        bot;
        AxisType=Axis3,
        fig = gd1,
        # dorecord=true,
        showtitle=false,
        showpoints=false,
        showlabels=false,
        showarrows=false,
        showcables=false,
        showwire=false,
        showmesh=false,
        xlims=(-1.2,1.0),
        ylims=(-0.3,0.8),
        zlims=(-0.0,1.8),
        doslide=false,
        showinfo=false,
        ground=inclined_plane,
        sup! = (ax,_,_) -> begin
            for (i,step) in enumerate(steps)
                RB.goto_step!(bot,step)
                tgvis = deepcopy(bot.structure)
                (;r,g,b) = cg[i]
                viz!(ax,tgvis;
                    meshcolor=Makie.RGBAf(r,g,b,alphas[i])
                )
            end
            lines!(ax,rp1[begin:laststep],color=:red)
            lines!(ax,rp12[begin:laststep],color=:blue)
            lines!(ax,rp10[begin:laststep],color=:purple)
            # lines!(ax,rp2[:,506:end],color=:blue)
            handlemesh = load(joinpath(RB.assetpath(),"把柄.STL")) |> make_patch(;
                scale = 1/400,
                trans = [ 0.0, 0.0, 1.55],
                rot = RotY(-π/2)
            )
            hidey(ax)
            mesh!(ax,handlemesh)
        end
    )
    ax2 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = "Energy (J)",
        # aspect = DataAspect()
    )
    lines!(ax2,t,me.E,label="E")
    lines!(ax2,t,me.T,label="T")
    lines!(ax2,t,me.V,label="V")
    # axislegend(ax2)
    Legend(gd2[1,2],ax2)
    ax3 = Axis(gd3[1,1])
    plot_convergence_order!(ax3,dts[begin:end-1],err_avg;show_orders=true)
    Legend(gd3[1,2],ax3)
    colsize!(fig.layout,1,0.45tw)
    rowgap!(fig.layout,0)
    Label(
        gd1[1,1,TopLeft()],"($(alphabet[1]))",font = :bold,
        padding = (0,0,0,0),
        justification = :right,
        lineheight = 1.0,
        halign = :center,
        valign = :bottom,        
    )
    Label(
        gd2[1,1,TopLeft()],"($(alphabet[2]))",font = :bold
    )
    Label(
        gd3[1,1,TopLeft()],"($(alphabet[3]))",font = :bold
    )
    # savefig(fig,"meteor_sliding")
    fig
end

# swing
flexcable_DR = make_flexcable(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.6*√3, 0.0, 1.2],
        rjx = [ 0.0,-1.0, 0.0],
        L=1.2,nx=5,doDR=true
)
RB.GDR!(flexcable_DR;β=2e-4,maxiters=2e5,verbose=true)
flexcable = make_flexcable(;
    ri  = [ 0.0, 0.0, 1.5],
    rix = [ 0.0, 0.0,-1.0],
    rj  = [-0.6*√3, 0.0, 1.2],
    rjx = [ 0.0,-1.0, 0.0],
    e = 0.5,
    μ = 0.9,
    L=1.2,nx=5,R=RotY(-π/2)
)
# flexcable_DR.traj.q[end][end-8:end] .= flexcable.traj.q[end][end-8:end]

RB.set_new_initial!(flexcable,flexcable_DR.traj.q[end],flexcable_DR.traj.q̇[end])

plane_normal = RotZY(deg2rad(-7.5),deg2rad(-90))*[0,0,1.0]
inclined_plane = RB.Plane(plane_normal,[0.05,0,0])

plane_normal'*rp+inclined_plane.d

GM.activate!();with_theme(theme_pub;
        Poly= (
            transparency = true,
        ) 
    ) do 
    plot_traj!(flexcable;
    ground=inclined_plane
    )
end

h = 2e-4
tspan = (0.0,1.5)
RB.solve!(
    RB.DynamicsProblem(flexcable,(x)->flexcable_contact_dynfuncs(x,inclined_plane)),
    RB.ZhongCCP();
    tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=false
)

me = RB.mechanical_energy!(flexcable)
lines(flexcable.traj.t,(me.E.-me.E[begin])./me.E[begin])
rp2 = RB.get_trajectory!(flexcable,6,2)
vp2 = RB.get_velocity!(flexcable,6,2)
pnvp2 = [plane_normal'*v for v in vp2 ]
lines(pnvp2)

plotsave_contact_persistent(flexcable,tol=0)

with_theme(theme_pub;
        size = (0.95tw,0.5tw),
        figure_padding = (0,fontsize,0,0),
        Axis3=(
            azimuth = 8.125530633326981,
            elevation = 0.18269908169872404
        ),
        Poly=(
            transparency=true,
        )
    ) do
    bot = flexcable
    (;t) = bot.traj
    fig = Figure()
    gd2 = fig[1,2] = GridLayout()
    gd3 = fig[2,2] = GridLayout()
    gd1 = fig[:,1] = GridLayout()
    steps = 1:500:4001
    alphas = [1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, 0.2, 1]
    cg = cgrad(:winter, length(steps), categorical = true, rev = true)
    plot_traj!(
        bot;
        AxisType=Axis3,
        fig = gd1,
        # dorecord=true,
        showtitle=false,
        showpoints=false,
        showlabels=false,
        showarrows=false,
        showcables=false,
        showmesh=false,
        showwire=false,
        xlims=(-1.2,0.1),
        ylims=(-0.5,0.5),
        zlims=(-0.1,1.8),
        doslide=false,
        showinfo=false,
        ground=inclined_plane,
        sup! = (ax,_,_) -> begin
            for (i,step) in enumerate(steps)
                RB.goto_step!(bot,step)
                tgvis = deepcopy(bot.structure)
                (;r,g,b) = cg[i]
                viz!(ax,tgvis;meshcolor=Makie.RGBA(r,g,b,alphas[i]))
            end
            lines!(ax,rp2[:,begin:4001],color=:blue)
            handlemesh = load(joinpath(RB.assetpath(),"把柄.STL")) |> make_patch(;
                scale = 1/400,
                trans = [ 0.0, 0.0, 1.55],
                rot = RotY(-π/2)
            )
            mesh!(ax,handlemesh)
        end
    )
    ax2 = Axis(gd2[1,1],
        xlabel = tlabel,
        ylabel = L"\acute{v}_n~(\mathrm{m/s})",
        # aspect = DataAspect()
    )
    impact_time = 2995*h
    step_before_impact = time2step(impact_time,t)
    @myshow step_before_impact
    vminus, vplus = pnvp2[step_before_impact:step_before_impact+1]
    @myshow vminus, vplus, vplus/(-vminus)
    vlines!(ax2,[impact_time];linestyle=:dash)
    lines!(ax2,t,pnvp2,color=:red)
    ylims!(ax2,-7,7)
    xlims!(ax2,extrema(t)...)
    poly!(ax2,Point2f.([[0.45,-7],[0.45,7],[0.75,7],[0.75,-7]]),
        color=Makie.RGBAf(0,0,0,0),
        strokecolor=:black,
        strokewidth=1,
        )
    ax22 = Axis(gd2[1,2],
        xlabel = tlabel,
        ylabel = L"\acute{v}_n~(\mathrm{m/s})",
        # aspect = DataAspect()
    )
    vlines!(ax22,[impact_time];linestyle=:dash)
    lines!(ax22,t,pnvp2,color=:red)
    ylims!(ax22,-7,7)
    xlims!(ax22,0.45,0.75)
    # axislegend(ax2)
    ax3 = Axis(gd3[1,1],
        xlabel = tlabel,
        ylabel = "Energy (J)"
    )
    vlines!(ax3,[impact_time];linestyle=:dash)
    lines!(ax3,t,me.E,label="E")
    lines!(ax3,t,me.T,label="T")
    lines!(ax3,t,me.V,label="V")
    Legend(gd3[1,2],ax3,orientation=:horizontal,tellheight=false)
    xlims!(ax3,extrema(t)...)
    # axislegend(ax3,position=:lc)
    colsize!(fig.layout,1,0.35tw)
    colsize!(gd3,1,0.3tw)
    rowgap!(fig.layout,0)
    Label(
        gd1[1,1,TopLeft()],"($(alphabet[1]))",font = :bold,
    )
    Label(
        gd2[1,1,TopLeft()],"($(alphabet[2]))",font = :bold
    )
    Label(
        gd2[1,2,TopLeft()],"($(alphabet[3]))",font = :bold
    )
    Label(
        gd3[1,1,TopLeft()],"($(alphabet[4]))",font = :bold
    )
    savefig(fig,"meteor_swing")
    fig
end

#-- meteor end