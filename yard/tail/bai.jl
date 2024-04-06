using Revise #jl
import Rible as RB

include(joinpath(pathof(RB),"../../yard/stability_stiffness.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = ""

include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl
#-- preamble

mesh_head = load(RB.assetpath("tail/装配体头部-1.2.3.STL")) |> RB.make_patch(;trans = [10, 0, 0],rot=RotZYX(π/2,-π/2,π/2))
mesh_rib3 = load(RB.assetpath("tail/肋片-3.STL")) |> RB.make_patch(;trans = [43+10, 0, 0], rot=RotZYX(π/2,-π/2,π/2))
mesh_headrib = GB.merge([mesh_head,mesh_rib3])
mesh_rib4 = load(RB.assetpath("tail/肋片-4.STL")) |> RB.make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_rib5 = load(RB.assetpath("tail/肋片-5 - 副本.STL")) |> RB.make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_rib7 = load(RB.assetpath("tail/肋片-7 - 副本.STL")) |> RB.make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_fin = load(RB.assetpath("tail/尾鳍-v1-软.STL")) |> RB.make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_bar = load(RB.assetpath("tail/摆杆-v1.2.2.STL")) |> RB.make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh(mesh_fin,color=:slategrey)
meshes = [
    mesh_headrib,
    mesh_rib4,
    mesh_rib5,
    mesh_rib7,
    mesh_fin,
    mesh_bar
]

质心 = [
    [-1.05,0], #装配体头部 + 肋片-3
    [-1.42,0], #肋片-4
    [-0.71,0], #肋片-5
    [-1.0,0], #肋片-7 
    [141.191,0], #尾鳍-v1 - 硬
    [50,0] #摆杆-v1.2.2
]

质量 = [
    55.31, #装配体头部 + 肋片-3
    33.34, #肋片-4
    26.79, #肋片-5
    14.86, #肋片-7 
    0.095e-3, #尾鳍-v1 - 硬
    14.86 #摆杆-v1.2.2
]

惯量 = [
    83897.99, #装配体头部 + 肋片-3
    38489.59, #肋片-4
    23523.93, #肋片-5
    6176.25, #肋片-7 
    308.454, #尾鳍-v1 - 硬
    308.454 #摆杆-v1.2.2
]

# 定义一个函数来新建机器人
function make_bai(meshes,质心,质量,惯量)

    nb = 6
    ver_lengths_raw = [
        #  0   43;
        # 18   18;
        0    43+18+18-10;
        18   16;
        19.1  9;
        18    8;
        18    8;
         0   16.29;
         0  150
    ]
    ver_lengths = sum(ver_lengths_raw,dims=2)
    hor_lengths = [
        32.61,
        26.21,
        20.64,
        14.71,
        13.24,
        # 13.24,
         2.5
    ]
    nhor = length(hor_lengths)
    P = [
        zeros(2,)
        for i = 1:nhor+1
    ]
    for i = 2:nhor+1
        P[i] = P[i-1] + [0.0,-ver_lengths[i-1]]
    end

    function rigidbody(i,P)
        if i == 1
            contactable = false
            visible = true
            ci = collect(1:6)
            cstr_idx = Int[]
        elseif i in [2,nb]
            contactable = true
            visible = true
            ci = collect(1:2)
            cstr_idx = collect(1:3)
        else
            contactable = true
            visible = true
            ci = Int[]
            cstr_idx = collect(1:3)
        end
        α = -π/2
        if i == 1
            ri = SVector{2}(0.0,  0.0)
            rj = SVector{2}(P[i+1] )
        elseif i < nb
            ri = SVector{2}(P[i])
            rj = SVector{2}(P[i+1])
        else            
            ri = SVector{2}(0.0,  0.0)
            rj = SVector{2}(0.0,-150.0)
        end
        b = hor_lengths[i]
        offset_fix = begin
            if i == 1
                ret = 43 + 10
            elseif i == nb
                ret = 150
            else
                ret = 0
            end
            ret 
        end
        offset = ver_lengths_raw[i,1]
        ro = ri + SVector{2}(0.0,-offset)
        mass_locus = SVector{2}(质心[i])
        m = 质量[i]
        inertia = 惯量[i]/2
        Ī = SMatrix{2,2}([
            inertia 0
            0 inertia
        ])
        r̄p1 = SVector{2}(offset_fix,-b,) 
        r̄p2 = SVector{2}(offset_fix, b,)
        r̄p3 = SVector{2}(-offset, 0.0)
        r̄p4 = SVector{2}(norm(rj-ri)-offset, 0.0)
        loci = [r̄p1,r̄p2,r̄p3,r̄p4]
        nr̄p = length(loci)
        ṙo = zeros(2); ω = 0.0
        prop = RB.RigidBodyProperty(
            i,
            contactable,
            m,
            Ī,
            SVector{2}(mass_locus),
            loci;
            visible = visible,
        )
        @show i, ri, rj
        nmcs = RB.NCF.NC2P1V(ri, rj, ro, RB.rotation_matrix(α))
        state = RB.RigidBodyState(prop, ro, α, ṙo, ω, )
        coords = RB.NonminimalCoordinates(nmcs, ci, cstr_idx)
        mesh_rigid = meshes[i]
        body = RB.RigidBody(prop, state, coords, mesh_rigid)
    end
    rbs = [
        rigidbody(i, P) for i = 1:nb
    ]
    rigdibodies = TypeSortedCollection(rbs)

    ncables = 2(nb-1)
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for i = 1:ncables
        original_restlens[i] = 22.0
        ks[i] = 1000.0
    end
    spring_dampers = [RB.DistanceSpringDamper2D( original_restlens[i], ks[i], 0.0;slack=true) for i = 1:ncables]  

    matrix_cnt_raw = Vector{Matrix{Int}}()
    for i = 1:nb-1
        s = zeros(2, nb)
        s[1, i  ] =  1
        s[1, i+1] = -1
        s[2, i  ] =  2
        s[2, i+1] = -2
        push!(matrix_cnt_raw, s)
    end
    connecting_matrix = reduce(vcat, matrix_cnt_raw)
    # display(connecting_matrix)
    cables = RB.connect(rigdibodies, spring_dampers; connecting_matrix)
    apparatuses = TypeSortedCollection(
        cables
    )

    numbered = RB.number(rigdibodies,apparatuses)
    matrix_sharing_raw = Vector{Matrix{Int}}()
    for i = 1:nb-2
        s = zeros(2, nb)
        s[1:2, i] = 3:4
        s[1:2, i+1] = 1:2
        push!(matrix_sharing_raw, s)
    end
    s = zeros(2, nb)
    s[1:2,  1] = 1:2
    s[1:2, nb] = 1:2
    push!(matrix_sharing_raw, s)
    sharing_matrix = reduce(vcat, matrix_sharing_raw)
    # display(sharing_matrix)
    indexed = RB.index(rigdibodies, apparatuses; sharing_matrix)

    cnt = RB.Connectivity(numbered, indexed,)

    st = RB.Structure(rigdibodies, apparatuses, cnt)

    gauges = Int[]
    actuators = Int[]
    #
    ## acs = [
    ##     RB.RegisterActuator(1,
    ##         [1:2(nb-1)],
    ##         original_restlens,
    ##     )
    ## ]
    
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )

    bot = RB.Robot(st, hub)
end

tail = make_bai(meshes,质心,质量,惯量)
bot = tail
botvis = deepcopy(bot)
plot_traj!(
    botvis;
    fontsize = 14 |> pt2px,
    show_axis = false,
    showmesh = true,
    showpoints = false,
    showlabels = false,
    xlims = (-50,50),
    ylims = (-300,40),
    showground = false
)

# 静力平衡
RB.check_static_equilibrium_output_multipliers(bot.structure;)

# 刚度、稳定性分析
ω²,δq̌ = RB.undamped_eigen(bot.structure)

## δq̌ = [Ň*orthovm[:,i] for i in axes(orthovm,2)]
scaling=0.2
nk = 5
RB.reset!(botvis)
for i = 1:nk
    push!(botvis.traj,deepcopy(botvis.traj[end]))
    botvis.traj.t[end] = i
    δq̌i = δq̌[i]
    ratio = norm(δq̌i)/norm(botvis.traj.q̌[begin])
    botvis.traj.q̌[end] .= botvis.traj.q̌[begin] .+ scaling.*δq̌i/ratio
end
with_theme(theme_pub;
        fontsize = 16,
        Axis3 = (
            azimuth = 4.7078743833269865,
            elevation = 1.307620956698724
        )
    ) do
    fig = Figure()
    gd2 = fig[1,1:nk] = GridLayout()
    plot_traj!(
        botvis;
        fig = gd2,
        AxisType=Axis3,
        showinfo = true,
        gridsize=(1, nk),
        atsteps=2:nk+1,
        doslide=false,
        showmesh = false,
        showwire = false,
        showlabels=true,
        showpoints=true,
        # showcables = false,
        showground = false,
        xlims = (-150,150),
        ## ylims = (-400,300),
        ylims = (-400,0),
        zlims = (-1,0),
        slack_linestyle = :solid,
        ## showinit = true,
        refcolor = Makie.RGBA{Float32}(0.5,0.5,0.5,0.5),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi+1])) ", font=:bold),
                    "Mechanism Mode $sgi"
                )
        end,
        sup! = (ax,tgob,sgi)->begin
            bodies = RB.get_bodies(tgob[])
            r1p3 = bodies[1].state.loci_states[3].frame.position
            r1p4 = bodies[1].state.loci_states[4].frame.position
            r2p4 = bodies[2].state.loci_states[4].frame.position
            r3p4 = bodies[3].state.loci_states[4].frame.position
            r4p4 = bodies[4].state.loci_states[4].frame.position
            r5p4 = bodies[5].state.loci_states[4].frame.position
            r5g = bodies[5].state.mass_locus_state.frame.position
            @show r1p3, r1p4
            lines!(ax,
                [
                    r1p3,r1p4,r2p4,r3p4,r4p4,r5p4,r5g
                ]
            )
            ## hidedecorations!(ax)
            ## xlims!(ax,-50,50)
            ## RB.hidez(ax)
            ## ylims!(ax,-300,40)
        end,
    )
    fig
end
