```@meta
EditURL = "<unknown>/examples/tail/dynamics2.jl"
```

# 二维类脊椎张拉整体动力学仿真

````@example tail
# using Literate #hide
# Literate.markdown("examples/tail/dynamics2.jl", "src/";name="tail") #hide
````

加载所需程序包

````@example tail
using LinearAlgebra
using SparseArrays
using StaticArrays
using Makie
import GLMakie
GLMakie.activate!()
# using BenchmarkTools
# using LaTeXStrings
# using NLsolve
using EponymTuples
using TypeSortedCollections
using Revise
import TensegrityRobots as TR
# includet("tail_define.jl")
include("examples/tail/make_new_tail.jl")
# includet("plotting.jl")
# includet("../analysis.jl")
include("examples/vis.jl")
````

节数

````@example tail
n = 4
````

张拉整体脊椎

````@example tail
tail = make_new_tail(n)
````

检查连接点

````@example tail
rbs = TR.get_rigidbodies(tail.tg)
plot_local_points(rbs[2])
````

设置初始条件

````@example tail
tail.traj.q̇[begin][tail.tg.connectivity.indexed.mem2sysfull[end][1:4]] .= [0.1,0.0,0.1,0.0]
````

定义广义力函数

````@example tail
function dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        # TR.apply_gravity!(tg)
        TR.generate_forces!(tg)
        TR.get_force!(F,tg)
        # F .= 0
    end
    Jac_F! = nothing
    @eponymtuple(F!)
end
````

动力学仿真问题

````@example tail
# dynfuncs(tail)
prob = TR.SimProblem(tail,dynfuncs)
# @code_warntype TR.SimProblem(twobaronetri,dynfuncs)
````

动力学仿真求解

````@example tail
TR.solve!(prob,TR.Zhong06();dt=0.01,tspan=(0.0,10.0),ftol=1e-13,verbose=true)
# sol = TR.solve(prob,TR.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
##@code_warntype TR.solve(prob,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
````

可视化

````@example tail
plot_traj!(tail)
````

系统机械能

````@example tail
me = TR.mechanical_energy!(tail)
scatter(me.E)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

