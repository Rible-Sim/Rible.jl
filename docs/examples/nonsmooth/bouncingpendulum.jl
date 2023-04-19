using LinearAlgebra
using Parameters
using GLMakie
using RecursiveArrayTools
using Revise
import TensegrityRobots as TR
function pendulum(θ0 = -π/12)
    m = 1
    l = 1
    J = 0.1
    ag = 10.0
    e = 0.8
    q0 = [l*cos(θ0), l*sin(θ0), θ0]
    v0 = zeros(3)
    𝒞 = Set([1,2,3])
    𝒰c = Set([1,2])
    𝐞 = vcat(zeros(length(𝒰c)),[e])

    function 𝒈(q)
        x,y,θ = q
        [
        x-l*cos(θ),
        y-l*sin(θ),
        x - √2/2
        ]
    end

    function 𝒈𝒒(q)
        x,y,θ = q
        [
        1 0  l*sin(θ);
        0 1 -l*cos(θ);
        1 0  0
        ]
    end

    function ∂𝒈𝒒T𝛌∂𝒒(q,λ)
        x,y,θ = q
        λ1,λ2,λ3 = λ
        [
        0 0  0;
        0 0  0;
        0 0  λ1*l*cos(θ) + λ2*l*sin(θ)
        ]
    end

    function ∂𝒈𝒒𝐯∂𝒒(q,v)
        x,y,θ = q
        v1,v2,v3 = v
        [
        0 0  v3*l*cos(θ);
        0 0  v3*l*sin(θ);
        0 0  0
        ]
    end

    function 𝐟(q,v,t)
        [
        0, -m*ag, 0
        ]
    end

    function 𝐌(q)
        Diagonal([m,m,J])
    end

    function ∂𝐌𝐚∂𝐪(q,a)
        zero(𝐌(q))
    end

    function ∂𝐟∂𝐪(q,v,t)
        zero(𝐌(q))
    end

    function ∂𝐟∂𝐯(q,v,t)
        zero(𝐌(q))
    end

    Jacobians = (∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯)
    dyfuncs = 𝐌,𝒈,𝒈𝒒,𝐟,Jacobians
    q0,v0,𝒞,𝒰c,𝐞,dyfuncs
end

q0,v0,𝒞,𝒰c,𝐞,dyfuncs = pendulum()
𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians


# θs = collect(-π/12:-0.01:-π/3)
# qs = [
#     begin
#         q0,_ = pendulum(θ)
#         q0
#     end for θ in θs
# ]

function simulate(;h = 1e-3,tspan=(0.0,4.0))
    ρ∞ = 0.8
    p = TR.generalized_α(ρ∞)
    n = 3
    c = 3
    ū = 2
    t = 0.0
    q0,v0,𝒞,𝒰c,𝐞,dyfuncs = pendulum()
    qs = TR.NSGA.initialize_St(n,c,ū,q0,v0,t,p,h,𝒞,𝒰c,𝐞,dyfuncs,tspan)
end

function plot_pendulum(qs;do_record=false,h=1e-3)
    x,y,θ = qs[1]
    xs = Node([0.0,x])
    ys = Node([0.0,y])
    fig = Figure()
    ax = fig[1,1] = Axis(fig)
    ls_step = labelslider!(fig,"step",1:length(qs))
    fig[2,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        qi = qs[this_step]
        xi,yi,θi = qi
        xs[] = [0.0,xi]
        ys[] = [0.0,yi]
        # active_sets = TR.initialize_St(n,c,ū,qi,v0,t,p,h,𝒞,𝒰c,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
        # @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
        # @show 𝒜,𝒜c,ℬ,ℬc
        # 𝐫ˢ,𝐫ᵖ,𝐫ᵛ =
    end
    vlines!(ax,[0])
    hlines!(ax,[0])
    vlines!(ax,[√2/2])
    scatterlines!(ax,xs,ys)
    ylims!(ax,-1.5,0.5)
    xlims!(ax,-0,1.5)
    ax.aspect = DataAspect()
    framerate = 30
    if do_record
        record(fig, "bouncingpendulum.mp4", 1:Int(ceil(1/framerate/h)):length(qs); framerate) do this_step
            qi = qs[this_step]
            xi,yi,θi = qi
            xs[] = [0.0,xi]
            ys[] = [0.0,yi]
        end
    end
    fig
end

qs = simulate(;tspan=(0.0,4.0))
plot_pendulum(qs)
plot_pendulum(qs;do_record=true)
g3 = VectorOfArray(𝒈.(qs))[3,:]
lines(collect(0:1e-3:4),g3)
xlims!(0,4)
ylims!(0,0.3)
qsref = simulate(;h=1e-4,tspan=(0.0,4.0))
g3ref = VectorOfArray(𝒈.(qsref))[3,:]
lines(collect(0:1e-3:4),g3)
lines!(collect(0:1e-4:4),g3ref)
