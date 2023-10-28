# # 二维类脊椎张拉整体动力学仿真
## using Literate #hide
## Literate.markdown("examples/tail/dynamics.jl", "src/";name="tail") #hide

# 加载所需程序包
using LinearAlgebra
using StaticArrays, SparseArrays, TypeSortedCollections
using EponymTuples
using Unitful, Match, Printf
using GeometryBasics, Meshing
using Rotations, CoordinateTransformations
using Makie
import GLMakie as GM
GM.activate!()
using Revise
import Rible as RB

include("examples/vis.jl")

# 定义一个函数来新建机器人
function make_new_tail(n)

    nver = n
    nhor = n + 1
    nb = nver + nhor
    ver_lengths = fill(0.2,nver)
    hor_lengths = fill(0.1,nhor)
    O = zeros(2, 2, nhor)
    for i = 2:nhor
        O[:, 1, i] = O[:, 1, i-1] + [0.0,-ver_lengths[i-1]]
    end
    for i = 1:nhor
        O[:, 2, i] = O[:, 1, i] + [hor_lengths[i], 0.0]
    end

    function rigidbody(i,O)
        lev,pos = divrem(i,2)
        if i == 1
            movable = false
            constrained = true
            ci = collect(1:6)
            Φi = Int[]
        else
            movable = true
            if i == 2
                constrained = true
                ci = collect(1:2)
                Φi = [1]
            else
                constrained = false
                ci = Int[]
                if isodd(pos)
                    Φi = collect(1:3)
                else
                    Φi = [1]
                end
            end
        end

        if isodd(pos)
            α = 0.0
            ri = SVector{2}(O[:,1, lev+1])
            rj = SVector{2}(O[:,2, lev+1])
            b = norm(rj-ri)
            mass_locus = SVector{2}(0.0,0.0)
            r̄p1 = SVector{2}(-b,0.0)
            r̄p2 = SVector{2}( b,0.0)
            m = inertia = 0.1
            Ī = SMatrix{2,2}([
                inertia 0
                0 inertia
            ])
        else
            α = -π / 2
            ri = SVector{2}(O[:,1, lev])
            rj = SVector{2}(O[:,1, lev+1])
            b = norm(rj-ri)
            mass_locus = SVector{2}(b/2,0.0)
            r̄p1 = SVector{2}(0.0,0.0)
            r̄p2 = SVector{2}(  b,0.0)
            m = inertia = 0.1
            Ī = SMatrix{2,2}([
                inertia 0
                0 0
            ])
        end
        loci = [r̄p1,r̄p2]
        nr̄p = length(loci)
        ṙo = zeros(2); ω = 0.0
        prop = RB.RigidBodyProperty(
            i,
            movable,
            m,
            Ī,
            SVector{2}(mass_locus),
            loci;
            constrained = constrained,
        )
        ro = ri
        if isodd(pos)
            lncs, _ = RB.NCF.NC2P1V(ri, rj, ro, α, ṙo, ω)
        else
            lncs, _ = RB.NCF.NC2D2P(ri, rj, ro, α, ṙo, ω)
        end
        state = RB.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)
        body = RB.RigidBody(prop, state)
    end
    rbs = [
        rigidbody(i, O) for i = 1:nb
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = RB.number(rigdibodies)
    matrix_sharing_raw = Vector{Matrix{Int}}()
    for i = 1:2n
        s = zeros(2, nb)
        if isodd(i)
            s[1:2, i] = 1:2
            s[1:2, i+1] = 1:2
        else
            s[1:2, i] = 3:4
            s[1:2, i+1] = 1:2
        end
        push!(matrix_sharing_raw, s)
    end
    matrix_sharing = reduce(vcat, matrix_sharing_raw)
    display(matrix_sharing)
    indexedcoords = RB.index(rigdibodies, matrix_sharing)

    ncables = 4n
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for i = 1:ncables
        lev,j = divrem(i-1,4)
        original_restlens[i] = ifelse(j ∈ [1, 0], ver_lengths[lev+1], ver_lengths[lev+1])
        ks[i] = ifelse(j ∈ [1, 0], 100, 100.0)
    end
    cables =
        [RB.Cable2D(i, 0.5original_restlens[i], ks[i], 0.0) for i = 1:ncables]  #
    tensiles = (cables = cables,)
    acs = [
        RB.ManualActuator(1,
            [4(i-1)+j for i = 1:n for j = 1:4],
            5original_restlens[[4(i-1)+j for i = 1:n for j = 1:4]],
        )
    ]
    hub = (actuators = acs,)

    matrix_cnt_raw = Vector{Matrix{Int}}()
    for i = 1:n
        s = zeros(4, nb)
        s[1, 2i-1] = 1
        s[1, 2i+1] = -1
        s[2, 2i] = 1
        s[2, 2i+1] = -1
        s[3, 2i] = 1
        s[3, 2i+1] = -2
        s[4, 2i-1] = 2
        s[4, 2i+1] = -2
        push!(matrix_cnt_raw, s)
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)
    display(matrix_cnt)
    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    cnt = RB.Connectivity(numberedpoints, indexedcoords, tensioned)

    st = RB.Structure(rigdibodies, tensiles, cnt)
    bot = RB.Robot(st, hub)
end

# 节数
n = 4

# 张拉整体脊椎
tail = make_new_tail(n)

# 检查连接点
rbs = RB.get_bodies(tail.st)
plot_rigid(rbs[2])


# 设置初始条件
tail.traj.q̇[begin][tail.st.connectivity.indexed.mem2sysfull[end][1:4]] .= [0.1,0.0,0.1,0.0]

# 定义广义力函数
function dynfuncs(bot)
    (;st) = bot
    function F!(F,q,q̇,t)
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.update_tensiles!(st)
        ## RB.apply_gravity!(st)
        RB.generate_forces!(st)
        RB.get_force!(F,st)
    end
    Jac_F! = nothing
    @eponymtuple(F!)
end

# 动力学仿真问题
prob = RB.SimProblem(tail,dynfuncs)

# 动力学仿真求解
RB.solve!(prob,RB.Zhong06();dt=0.01,tspan=(0.0,10.0),ftol=1e-13,verbose=true)

# 可视化
plot_traj!(tail;showmesh=false,showinfo=false,showground=false)

# 系统机械能
me = RB.mechanical_energy!(tail)
scatter(me.E)
