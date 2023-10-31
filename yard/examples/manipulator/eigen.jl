using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
using Makie
AbstractPlotting.__init__()

using Dierckx

using Revise

using TensegritySolvers; const TS = TensegritySolvers
using Robot
const TR = Robot

cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")
function eigenanalysis(n,α=π/6,k=1.e4)
    num_of_dof = 6
    manipulator = man_ndof(num_of_dof,k=k,c=0.0)
    as = Vector{Vector{Float64}}()

    Y = build_Y(manipulator)

    q0,_,_ = RB.get_initial(manipulator)

    λ0,_,a= RB.inverse(manipulator,deepcopy(manipulator),Y)
    push!(as,a)

    ωs = Vector{Vector{Float64}}()
    Zs = Vector{Matrix{Float64}}()

    ω0,Z0 = RB.undamped_eigen(manipulator,q0,λ0)

    push!(ωs,ω0)
    push!(Zs,Z0)

    for i = 1:n
        iθ = -i*α/n
        refman = man_ndof(num_of_dof,θ=iθ,k=k,c=0.0) # reference
        refqi,_,_ = RB.get_initial(refman)
        refλi,Δu,a= RB.inverse(manipulator,refman,Y)
        push!(as,a)
        actmani = deepcopy(manipulator)

        ωi,Zi = RB.undamped_eigen(actmani,refqi,refλi)
        push!(ωs,ωi)
        push!(Zs,Zi)
    end
    ωs,Zs,as
end

n = 10
ωs,Zs,as = eigenanalysis(n,π/12)
x = 0:(30/n):30
y = [a[1] for a in as]
spl = Spline1D(z,y,k=1)
plt.plot(0:0.1:30,spl(0:0.1:30))
at = 30
f(x) = at*(x/at)^5
z = f.(x)
spl = Spline1D(z,y)
plt.plot(0:0.1:at,f.(0:0.1:at))


function simulate_interpolation(spl,num_of_dof=6;k=4.e3,c=0.0)
    manipulator = man_ndof(num_of_dof,k=k,c=c)

    q0,q̇0,λ0 = RB.get_initial(manipulator)

    function interpolation_actuate(tgstruct,q0,target_t)

        M = RB.build_massmatrix(tgstruct)
        Φ = RB.build_Φ(tgstruct,q0)
        A = RB.build_A(tgstruct)

        function F!(F,q,q̇,t)
            if t < target_t
                # RB.actuate!(tgstruct,fill(spl(t),6))
                RB.actuate!(tgstruct,t/target_t*fill(spl(target_t),6))
            end
            RB.reset_forces!(tgstruct)
            RB.distribute_q_to_rbs!(tgstruct,q,q̇)
            RB.update_cables_apply_forces!(tgstruct)
            RB.assemble_forces!(F,tgstruct)
        end

        M,Φ,A,F!,nothing

        #A,Φ,∂T∂q̇!,F!,M!,nothing
    end



    dt = 0.01 # Same dt used for PID AND Dynamics solver
    prob = TS.DyProblem(interpolation_actuate(manipulator,q0,30.0),q0,q̇0,λ0,(0.0,100.0))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-10,verbose=true)
    manipulator, sol
end

# bars_and_cables_segs(man_linear)
man_interpolation, sol_interpolation = simulate_interpolation(spl)
plotstructure(man_interpolation,sol_interpolation,sliderplot)

for i = 1:n+1
    plt.plot(ωs[i],label="$i,k=300")
end
plt.legend()


ωs6,Zs6 = eigenanalysis(n,π/4,6.e2)

plt.plot(ωs6)
plt.legend()


x = [0., 1., 2., 3., 4.]
y = [-1., 0., 7., 26., 63.]  # x.^3 - 1.
spl = Spline1D(x, y; k = 3)
plt.plot(0:0.1:4,spl(0:0.1:4))

function freevibra(tgstruct,q0)

    M = RB.build_massmatrix(tgstruct)
    Φ = RB.build_Φ(tgstruct,q0)
    A = RB.build_A(tgstruct)

    function F!(F,q,q̇,t)

        RB.reset_forces!(tgstruct)
        RB.distribute_q_to_rbs!(tgstruct,q,q̇)
        RB.update_cables_apply_forces!(tgstruct)
        # F .= Q̃*RB.fvector(tgstruct)
        # RB.apply_gravity!(tgstruct,factor=0.01)
        # F .= G
        RB.assemble_forces!(F,tgstruct)
        # @show isapprox(F,Q̃*RB.fvector(tgstruct))
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = freevibra(manipulator,q0)

prob = TS.DyProblem(freevibra(manipulator,q0),q0,q̇0,λ,(0.0,50.0))


function perturbation!(q0,q̇0)
    q̇0[end] = 0.0001
end
perturbation!(q0,q̇0)

sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-12,verbose=true)
δq = [q-q0 for q in sol.qs]

plt.plot(transpose(hcat(δq...)))

q = RB.undamped_modal_solve!(manipulator,q0,q̇0,λ,50.0,dt)
plt.plot(transpose(q))
