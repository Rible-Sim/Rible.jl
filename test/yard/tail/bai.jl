using LinearAlgebra
using StaticArrays, SparseArrays, TypeSortedCollections
using EponymTuples
using Unitful, Match, Printf
using GeometryBasics, Meshing
using Rotations, CoordinateTransformations
using Unitful
using LaTeXStrings
using Makie
import GLMakie as GM
GM.activate!()
import GeometryBasics as GB
using FileIO, MeshIO
using Revise
import Rible as RB

cd("examples/tail")
include("../vis.jl")
includet("../vis.jl")
#-- preamble

mesh_head = load("装配体头部-1.2.3.STL") |> make_patch(;trans = [10, 0, 0],rot=RotZYX(π/2,-π/2,π/2))
mesh_rib3 = load("肋片-3.STL") |> make_patch(;trans = [43+10, 0, 0], rot=RotZYX(π/2,-π/2,π/2))
mesh_headrib = GB.merge([mesh_head,mesh_rib3])
mesh_rib4 = load("肋片-4.STL") |> make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_rib5 = load("肋片-5 - 副本.STL") |> make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_rib7 = load("肋片-7 - 副本.STL") |> make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_fin = load("尾鳍-v1-软.STL") |> make_patch(;rot=RotZYX(π/2,-π/2,π/2))
mesh_bar = load("摆杆-v1.2.2.STL") |> make_patch(;rot=RotZYX(π/2,-π/2,π/2))
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
        nmcs = RB.NCF.NC2P1V(ri, rj, ro, α)
        state = RB.RigidBodyState(prop, nmcs, ro, α, ṙo, ω, ci, cstr_idx)
        mesh_rigid = meshes[i]
        body = RB.RigidBody(prop, state, mesh_rigid)
        rb
    end
    rbs = [
        rigidbody(i, P) for i = 1:nb
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
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
    matrix_sharing = reduce(vcat, matrix_sharing_raw)
    # display(matrix_sharing)
    indexed = RB.index(rigdibodies, matrix_sharing)
    # indexed = RB.index(rigdibodies, )

    ncables = 2(nb-1)
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for i = 1:ncables
        original_restlens[i] = 22.0
        ks[i] = 1000.0
    end
    cables =
        [RB.DistanceSpringDamper2D( original_restlens[i], ks[i], 0.0;slack=true) for i = 1:ncables]  #
    apparatuses = (cables = cables,)
    acs = [
        RB.ManualActuator(1,
            [1:2(nb-1)],
            original_restlens,
        )
    ]
    hub = (actuators = acs,)

    matrix_cnt_raw = Vector{Matrix{Int}}()
    for i = 1:nb-1
        s = zeros(2, nb)
        s[1, i  ] =  1
        s[1, i+1] = -1
        s[2, i  ] =  2
        s[2, i+1] = -2
        push!(matrix_cnt_raw, s)
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)
    # display(matrix_cnt)
    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    cnt = RB.Connectivity(numbered, indexed, tensioned)

    st = RB.Structure(rigdibodies, apparatuses, cnt)
    bot = RB.Robot(st, hub)
end

tail = make_bai(meshes,质心,质量,惯量)

plot_traj!(
    tail;
    fontsize = 14 |> pt2px,
    showmesh = false,
    xlims = (-50,50),
    ylims = (-300,40),
    showground = false
)

# 静力平衡
RB.check_static_equilibrium_output_multipliers(tail.st;)

# 刚度、稳定性分析
RB.undamped_eigen(tail.st)

# 动力学分析
