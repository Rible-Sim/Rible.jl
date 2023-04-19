using LinearAlgebra
using StaticArrays
# using Parameters
using GLMakie
using Cthulhu
using RecursiveArrayTools
using BenchmarkTools
using GeometryBasics
using Meshing
using Rotations
using Revise
import TensegrityRobots as TR
cd("examples/nonsmooth")
# includet("plotting.jl")
includet("../analysis.jl")

function pointmass_contact_dynfuncs_one(m=1.0,μ=0.3,θ=π/12,e=0.9)
    o = zero(m)
    T = typeof(m)
    nq = 3
    M = Matrix{T}(m*I,nq,nq)
    zeroM = zero(M)
    𝐌(q) = M
    𝚽(q) = Vector{eltype(q)}()
    𝚽𝐪(q) = Matrix{eltype(q)}(undef,0,nq)
    𝚿(q,q̇) = Vector{eltype(q)}()
    𝚿𝐪(q,q̇) = Matrix{eltype(q)}(undef,0,nq)
    𝚿𝐪̇(q) = Matrix{eltype(q)}(undef,0,nq)
    gravity = 9.81
    G = [o,o,-m*gravity]
    function 𝐅!(F,q,q̇,t)
        F .= G
    end
    ∂𝐌𝐚∂𝐪(q,a) = zeros(eltype(q),nq,nq)
    ∂𝚽𝐪T𝛌∂𝒒(q,λ) = zeros(eltype(q),nq,nq)
    ∂𝚿𝐪̇T𝛍∂𝒒(q,μ) = zeros(eltype(q),nq,nq)
    function Jac_F!(∂F∂q,∂F∂v,q,q̇,t)
        ∂F∂q .= 0.0
        ∂F∂v .= 0.0
    end
    ∂𝐅∂𝐪(q,q̇,t) = zeroM
    ∂𝐅∂𝐯(q,q̇,t) = zeroM
    a = tan(θ)
    inclined_plane = (nv=[-a,0,1]./norm([-a,0,1]),origin=zeros(3))
    function 𝐠(q)
        nv = inclined_plane.nv
        d = transpose(nv)*q/norm(nv)
        gaps = [d]
    end

    function get_indices(q)
        gaps = 𝐠(q)
        active_indices = findall((x)->x≤0, gaps)
        length(active_indices),active_indices,gaps
    end

    function get_FCs(active_indices,q)
        n = inclined_plane.nv
        function make_FC(pid)
            cf = TR.ContactFrame(n)
            C = I
            gd = TR.GeneralizedDirections{3}(cf,C)
            TR.FrictionCone(μ,e,cf,gd)
        end
        [make_FC(i) for i in active_indices]
    end

    function get_D(active_indices,q)
        FCs = get_FCs(active_indices,q)
        Ds = [[transpose(FC.gd.Dn);transpose(FC.gd.Du);transpose(FC.gd.Dw)] for FC in FCs]
        D = reduce(vcat,Ds)
        μs = [one(FC.μ) for FC in FCs]
        es = [FC.e for FC in FCs]
        H = Diagonal(reduce(vcat,[[inv(FC.μ),one(FC.μ),one(FC.μ)] for FC in FCs]))
        D, μs, es, H
    end

    Jacobians = (Jac_F!,𝚿𝐪,∂𝚽𝐪T𝛌∂𝒒,∂𝚿𝐪̇T𝛍∂𝒒)
    contact_funcs = (𝐠,get_indices,get_FCs,get_D)
    M,𝚽,𝚽𝐪,𝚿,𝚿𝐪̇,𝐅!,Jacobians,contact_funcs
end

function pointmass_contact_dynfuncs(m=1.0,μ=0.3,θ=π/12,e=0.9)
    o = zero(m)
    T = typeof(m)
    nq = 3
    M = Matrix{T}(m*I,nq,nq)
    zeroM = zero(M)
    𝐌(q) = M
    𝚽(q) = Vector{eltype(q)}()
    𝚽𝐪(q) = Matrix{eltype(q)}(undef,0,nq)
    𝚿(q,q̇) = Vector{eltype(q)}()
    𝚿𝐪(q,q̇) = Matrix{eltype(q)}(undef,0,nq)
    𝚿𝐪̇(q) = Matrix{eltype(q)}(undef,0,nq)
    gravity = 9.81
    G = [o,o,-m*gravity]
    function 𝐅!(F,q,q̇,t)
        F .= G
    end
    ∂𝐌𝐚∂𝐪(q,a) = zeros(eltype(q),nq,nq)
    ∂𝚽𝐪T𝛌∂𝒒(q,λ) = zeros(eltype(q),nq,nq)
    ∂𝚿𝐪̇T𝛍∂𝒒(q,μ) = zeros(eltype(q),nq,nq)
    function Jac_F!(∂F∂q,∂F∂v,q,q̇,t)
        ∂F∂q .= 0.0
        ∂F∂v .= 0.0
    end
    ∂𝐅∂𝐪(q,q̇,t) = zeroM
    ∂𝐅∂𝐯(q,q̇,t) = zeroM
    a = tan(θ)
    inclined_plane = (nv=[-a,0,1]./norm([-a,0,1]),origin=zeros(3))
    function 𝐠(q)
        nv = inclined_plane.nv
        d = transpose(nv)*q/norm(nv)
        gaps = [d]
    end

    function get_indices(q)
        gaps = 𝐠(q)
        active_indices = findall((x)->x≤0, gaps)
        length(active_indices),active_indices,gaps
    end

    function get_FCs(active_indices,q)
        n = inclined_plane.nv
        function make_FC(pid)
            cf = TR.ContactFrame(n)
            C = I
            gd = TR.GeneralizedDirections{3}(cf,C)
            TR.FrictionCone(μ,e,cf,gd)
        end
        [make_FC(i) for i in active_indices]
    end

    function get_D(active_indices,q)
        FCs = get_FCs(active_indices,q)
        Ds = [[transpose(FC.gd.Dn);transpose(FC.gd.Du);transpose(FC.gd.Dw)] for FC in FCs]
        D = reduce(vcat,Ds)
        μs = [FC.μ for FC in FCs]
        es = [FC.e for FC in FCs]
        D, μs, es
    end

    Jacobians = (Jac_F!,𝚿𝐪,∂𝚽𝐪T𝛌∂𝒒,∂𝚿𝐪̇T𝛍∂𝒒)
    contact_funcs = (𝐠,get_indices,get_FCs,get_D)
    M,𝚽,𝚽𝐪,𝚿,𝚿𝐪̇,𝐅!,Jacobians,contact_funcs
end

μ = 0.29
θ = 0.0
m = 5.0
M,𝚽,𝚽𝐪,𝚿,𝚿𝐪̇,𝐅!,Jacobians,contact_funcs = pointmass_contact_dynfuncs(m,μ,θ)
Jac_F!,𝚿𝐪,∂𝚽𝐪T𝛌∂𝒒,∂𝚿𝐪̇T𝛍∂𝒒 = Jacobians
𝐠,get_indices,get_FCs,get_D = contact_funcs
q0 = [0.0,0.0,0.5]
get_indices(q0)
@code_warntype get_indices(q0)
# θ = 0.0
# a = tan(θ)
# inclined_plane = (nv=[-a,0,1]./norm([-a,0,1]),origin=zeros(3))
#
# x = 21.3
# y = 0.0
# z = tan(θ)*x
# q = [x,y,z] + -0.5inclined_plane.nv
# 𝐌(q)
# 𝚽(q)
# 𝚽𝐪(q)
# 𝐅(q,q,0)
# λ = ones(1)
# ∂𝚽𝐪T𝛌∂𝒒,∂𝚽𝐪𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐅∂𝐪,∂𝐅∂𝐯 = Jacobians
# ∂𝚽𝐪T𝛌∂𝒒(q,λ)
# ∂𝚽𝐪𝐯∂𝒒(q,q)
# ∂𝐌𝐚∂𝐪(q,q)
# ∂𝐅∂𝐪(q,q,0)
# ∂𝐅∂𝐯(q,q,0)
# 𝐠,get_indices,get_FCs,get_D = contact_funcs
# 𝐠(q)
# u,active_indices,gap = get_indices(q)
# get_FCs(active_indices,q)

μ = 0.5
θ = π/12
m = 5.0
a = tan(θ)
tspan = (0.0,0.322)
tspan = (0.0,0.145)
tspan = (0.0,4.0)
h = 1e-2
p = TR.generalized_α(0.8,h)
nq = 3
nλ = 0
nμ = 0

# ts,qs,vs = TR.NSSFC.nssfc(nq,nλ,q0,q̇0,p,h,pointmass_contact_dynfuncs(;m,μ,θ),tspan;tol=1e-11,imax=10)
ts,cs,qs,vs,ps,λs,μs = TR.NSSFC.nhsolve(nq,nλ,nμ,q0,q̇0,pointmass_contact_dynfuncs(m,μ,θ),tspan;dt=h,ftol=1e-11,exception=false)

ts,cs,qs,vs,ps,λs,μs = TR.NSSFC.nhsolve(nq,nλ,nμ,q0,q̇0,pointmass_contact_dynfuncs_one(m,μ,θ),tspan;dt=h,ftol=1e-11,exception=false)

ts,cs,qs,vs,ps,λs,μs = TR.NSSFC.ipsolve(nq,nλ,nμ,q0,q̇0,
    pointmass_contact_dynfuncs_one(m,μ,θ),tspan;dt=h,ftol=1e-12,imax=50,exception=false)

TR.NSSFC.NTScale([21.86239409424172, 0.0, 18.033828064662508], [29.060492941155868, 0.0, 29.060492941155868])

@descend_code_warntype TR.NSSFC.nhsolve(nq,nλ,nμ,q0,q̇0,pointmass_contact_dynfuncs(m,μ,θ),tspan;dt=h,ftol=1e-11)

plot(ts,norm.(qs))
plot(ts,[1/2*m*v⋅v + m*q[3]*9.81 for (q,v) in zip(qs,vs)])
plot(ts,abs.(norm.(vs)))

plane_n = [-a,0,1]
plane_r = zeros(3)
plane = TR.Plane(plane_n,plane_r)
TR.signed_distance(q̇0,plane)

function vis_pointmass(qs,plane;do_record=false)
    fig = Figure(resolution=(1920,1080))
    ax = fig[1,1] = LScene(fig, scenekw = (show_axis=false, camera = cam3d!, raw = false))
    points = [
        Point3f(-0.25,0.0,-0.25),
        Point3f(0.0,0.25,-1.0),
        Point3f(0.0,-0.25,-1.0),
        Point3f(0.5,0.0,0.5)
    ]
    plane_mesh = GeometryBasics.Mesh(Rect(points), NaiveSurfaceNets()) do v
        TR.signed_distance(v,plane)
    end
    mesh!(ax,plane_mesh,color="white",ambient = Vec3f(0.9, 0.9, 0.9),diffuse = Vec3f(0.4, 0.4, 0.4))
    p = Observable([Point3(qs[1])])
    scatter!(p,color="blue")

    if do_record
        lookat = [0.11257553, -7.635099f-7, 0.079256475]
        eyepos = [-0.012789218, -0.68866956, 0.13630177]
        framerate = 30
        h = 1e-3
        record(fig, "pointmass.mp4", 1:round(Int,1/(30*4)/h):length(qs); framerate) do this_step
            p[] = [Point3(qs[this_step])]
            update_cam!(ax.scene, eyepos, lookat)
        end
    else
        sg = SliderGrid(fig[2,1],
                (label="step",range=collect(eachindex(qs)))
            )
        on(sg.sliders[1].value) do this_step
            p[] = [Point3(qs[this_step])]
            camera = cameracontrols(ax.scene)
            lookat = camera.lookat[]
            eyepos = camera.eyeposition[]
            # @show lookat, eyepos
        end
    end
    fig
end

vis_pointmass(qs,plane)

vis_pointmass(qs,plane;do_record=true)

function random_socv(;p = 3, l = 1.0, μ=1)
    x1 = l*rand(p-1)
    x = vcat(μ*norm(x1)+l*rand()+l*1,x1)
end

x_split = [random_socv() for i = 1:10]
y_split = [random_socv() for i = 1:10]

TR.NSSFC.:⊙(x_split,y_split)

TR.NSSFC.:⊘(x_split,y_split)

TR.NSSFC.:⊙(x_split,TR.NSSFC.:⊘(x_split,y_split)) .- y_split

TR.NSSFC.:⊘(x_split,TR.NSSFC.:⊙(x_split,y_split)) .- y_split

TR.NSSFC.:⊙(y_split,TR.NSSFC.:⊘(y_split,x_split)) .- x_split

W = TR.NSSFC.NTScale(x_split,y_split)


x = reduce(vcat,x_split)
y = reduce(vcat,y_split)
z = W*x
z - inv(W)*y

W_blocks = TR.NSSFC.NTScale.(x_split,y_split)
ΘW_blocks = TR.NSSFC.NTScale_Anderson.(x_split,y_split)

@btime TR.NSSFC.NTScale(x_split[1],y_split[1])
@btime TR.NSSFC.NTScale_Anderson(x_split[1],y_split[1])
W_blocks .- ΘW_blocks


z_split = W_blocks.*x_split
z_split .- inv.(W_blocks).*y_split

TR.NSSFC.:⊙(z_split,z_split) .- TR.NSSFC.:⊙(x_split,y_split)

x_split = [[3.8760483208829286, 0.0, 3.876048320487105]]
y_split = [[0.8591722005249399, 0.0, -0.859172200010235]]
Δx_split = [[-6.352136087299624e-9, -0.0, -2.8367571879860825e-9]]
Δy_split = [[-6.318927479303759e-9, 0.0, -2.4428789769374963e-9]]

x_split = [random_socv() for i = 1:100]
y_split = [random_socv() for i = 1:100]
Δy_split = [rand(3) for i = 1:100]
Δx_split = [rand(3) for i = 1:100]
# W_blocks = TR.NSSFC.NTScale.(x_split,y_split)
# z_split = W_blocks.*x_split
J = Diagonal(vcat(one(x_split[1][begin]),-one.(x_split[1][begin+1:end])))
αmax_x = TR.NSSFC.find_cone_step_length(x_split,Δx_split,J)
αmax_y = TR.NSSFC.find_cone_step_length(y_split,Δy_split,J)
α_max = min(αmax_x,αmax_y)
αmax = TR.NSSFC.find_cone_step_length(z_split,W_blocks,Δy_split,Δx_split,J)
α = 0.99αmax
x2_split = x_split.+α.*Δx_split
y2_split = y_split.+α.*Δy_split
minimum([transpose(xi)*J*xi for xi in x2_split])
minimum([transpose(yi)*J*yi for yi in y2_split])

x = [3.8760483201787226, -1.1339090909894236, 3.706480912365776]
Δx = [-1.1620053987098517e-15, 2.0031919628115605e-9, 6.128284079931474e-10]
J = Diagonal(vcat(one(x[begin]),-one.(x[begin+1:end])))
TR.NSSFC.find_cone_step_length(x,Δx,J)
