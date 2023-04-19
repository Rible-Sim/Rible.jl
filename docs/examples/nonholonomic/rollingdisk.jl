using LinearAlgebra
using StaticArrays
using Parameters
using GLMakie
using RecursiveArrayTools
using TypeSortedCollections
using ForwardDiff
using Revise
import TensegrityRobots as TR

function make_rollingdisk()
    m = 1.0
    R = 1.0
    r̄g = zeros(2)
    Īg = SMatrix{2,2}(
		[
			1/4*m*R^2       0.0;
			      0.0 1/4*m*R^2;
		]
	)
    ro = [0.0,R]
    α = 0.0
    ṙo = [0.0,0.0]
    ω = 4.0
    aps = [r̄g]
    movable = true
	constrained = false
    prop = TR.RigidBodyProperty(1,true,m,Īg,r̄g,aps)
    ri = copy(ro)
    lncs,q,q̇ = TR.NCF.NC1P2V(ri,ro,α,ṙo,ω)
    state = TR.RigidBodyState(prop,lncs,ro,α,ṙo,ω)
    rb1 = TR.RigidBody(prop,state)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = TR.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = TR.index(rbs,matrix_sharing)
	#
	ss = Int[]
	tensiles = (cables = ss,)
	hub = nothing
	#
	connections = TR.connect(rbs,zeros(Int,0,0))

	jointedmembers = TR.unjoin()

	cnt = TR.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers)
	tg = TR.TensegrityStructure(rbs,tensiles,cnt)
    bot = TR.TensegrityRobot(tg,nothing)
end

function dynfuncs(bot)
    @unpack tg = bot
    nq = 6
    nμ = 2
    nu = 2
    M = Matrix(TR.build_M(tg))
    Φ = TR.make_Φ(tg)
    A = TR.make_A(tg)
    R = 1.0
    function Cc(q)
        x,y,u1,u2,v1,v2 = q
        c1 = -R*v1/(v1*u2-u1*v2)
        c2 = -R*u1/(u1*v2-v1*u2)
        C = [1 0 c1  0 c2  0;
             0 1  0 c1  0 c2]
    end

    Ψ(q,q̇) = Cc(q)*q̇
    Ψq̇(q) = Cc(q)

    function F!(F,q,q̇,t)
        F .= 0
    end

    function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
        ∂F∂q .= 0
        ∂F∂q̇ .= 0
    end

    E = Diagonal([0.4,0.4])

    wall1_x = -2.0
    wall2_x =  2.0
    function g(q)
        x,y,u1,u2,v1,v2 = q
        [x-wall1_x, wall2_x-x]
    end

    function gq(q)
        T = eltype(q)
        o = zero(T)
        i = one(T)
        [ i o o o o o;
         -i o o o o o;]
    end

    function ∂gqᵀΛ∂q(q,Λ)
        # function gqᵀΛ(q)
        #     transpose(gq(q))*Λ
        # end
        ret = zeros(eltype(q),nq,nq)
        # ForwardDiff.jacobian!(ret,gqᵀΛ,q)
        # ret
    end

    function ∂gqq̇∂q(q,q̇)
        function gqq̇(q)
            gq(q)*q̇
        end
        ret = zeros(eltype(q),nu,nq)
        # ForwardDiff.jacobian!(ret,gqq̇,q)
        # ret
    end

    function ∂Aᵀλ∂q(q,λ)
        TR.∂Aᵀλ∂q̌(tg,λ)
    end

    function Ψq(q,q̇)
        function inner_Ψ(q)
            Ψ(q,q̇)
        end
        ret = zeros(eltype(q),nμ,nq)
        ForwardDiff.jacobian!(ret,inner_Ψ,q)
        ret
    end

    function ∂Bᵀμ∂q(q,μ)
        function Bᵀμ(q)
            transpose(Ψq̇(q))*μ
        end
        ret = zeros(eltype(q),nq,nq)
        ForwardDiff.jacobian!(ret,Bᵀμ,q)
        ret
    end

    jacobians = Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q
    contact_funcs = E,g,gq,∂gqᵀΛ∂q,∂gqq̇∂q
    M,Φ,A,Ψ,Ψq̇,F!,jacobians,contact_funcs
end

disk = make_rollingdisk()

M,Φ,A,Ψ,Ψq̇,F!,jacobians,contact_funcs = dynfuncs(disk)


h = 1e-3
tspan = (0.0,10.0)
q0 = TR.get_q(disk.tg)
q̇0 = TR.get_q̇(disk.tg)
prob = TR.SimProblem(disk,dynfuncs)
ts,qs,q̇s,ps,λs,μs = TR.nhsolve(prob,6,3,2,2,q0,q̇0;tspan,dt=h,ftol=1e-10,exception=false)
ts,qs,q̇s,ps,λs,μs = TR.snhsolve(prob,6,3,2,2,q0,q̇0;tspan,dt=h,ftol=1e-12,exception=true)
ts,qs,q̇s,ps,λs,μs = TR.ipsolve(prob,6,3,2,q0,q̇0;dt=h,ftol=1e-10,exception=false)

TR.prepare_traj!(disk.traj;tspan,dt=h,restart=true)

disk.traj.t[2:end] = ts[2:end]
disk.traj.q[2:end] = qs[2:end]
disk.traj.q̇[2:end] = q̇s[2:end]
disk.traj.q̇
lines([q̇[1] for q̇ in disk.traj.q̇])
function vis(bot;do_record=false)
    (;tg, traj) = bot
    # (;q, q̇) = traj
    fig = Figure(resolution=(1920,1080))
    ax = fig[1,1] = Axis(fig)
    centerpoint = Point2f(traj.q[1][1:2])
    direction = Point2f(traj.q[1][3:4])
    points = Observable([centerpoint])
    directions = Observable([direction])
    radius = 1.0
    circle = Observable(GLMakie.Circle(centerpoint,radius))
    mesh!(ax,circle,color=:white)
    arrows!(ax,points,directions)
    vlines!(ax,[-3,3],linewidth=2)
    ax.aspect = DataAspect()
    xlims!(ax,-4,4)
    ylims!(ax,-0.5,5)
    if do_record
        framerate = 30
        h = 1e-3
        record(fig, "rollingdisk.mp4", 1:round(Int,1/(30/2)/h):length(traj.q); framerate) do this_step
            qi = traj.q[this_step]
            centerpoint = Point2f(qi[1:2])
            direction = Point2f(qi[3:4])
            points[] = [centerpoint]
            directions[] = [direction]
            circle[] = GLMakie.Circle(centerpoint,radius)
        end
    else
		sg = SliderGrid(fig[2,1],
					(label = "step",range = 1:length(traj.q), startvalue = 1),
				)
        on(sg.sliders[1].value) do this_step
            qi = traj.q[this_step]
            centerpoint = Point2f(qi[1:2])
            direction = Point2f(qi[3:4])
            points[] = [centerpoint]
            directions[] = [direction]
            circle[] = GLMakie.Circle(centerpoint,radius)
        end
    end
    fig
end

vis(disk)
vis(disk;do_record=true)
