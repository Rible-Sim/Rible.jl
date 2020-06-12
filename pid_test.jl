using Test
@testset "trigonometric identities" begin
           θ = 2/3*π
           @test sin(-θ) ≈ -sin(θ)
           @test cos(-θ) ≈ cos(θ)
           @test sin(2θ) ≈ 2*sin(θ)*cos(θ)
           @test cos(2θ) ≈ cos(θ)^2 - sin(θ)^2
end;

pid = PID(1, 1., 1., setpoint=0.)
@test pid(0) == 0

pid = PID(1., 0., 0., setpoint=10.)
@test pid(0) == 10
@test pid(5) == 5
@test pid(-5) == 15

pid = PID(1., 0., 0., setpoint=-10.)
@test pid(0) == -10
@test pid(5) == -15
@test pid(-5) == -5
@test pid(-15) == 5

using DifferentialEquations
using Plots
pyplot()
x0 = 1.0
ẋ0 = 0.0
# m*ẍ + k*x + c*ẋ = p
# p = Kp*(setpoint-x) + Kd*(-ẋ)
# m*ẍ + k*x + Kp*(x-setpoint) + c*ẋ + Kd*ẋ = 0
# m*ẍ + k*x + Kp*x + c*ẋ + Kd*ẋ = Kp*setpoint

function fmaker()
    pid = PID(1.0,1.0,1.0,setpoint=2.0)
    function innerf(ẋ,x,p,t)
        k = 1.0
        c = 1.0
        m = 1.0
        λ = 1.0
        p = update!(pid,x,t)
        1/m*(p - k*x - c*ẋ - λ*x^3)
    end
end

function pid_maker()
    function pid_affect!(intor)
        #@show typeof(intor.u[1]),intor.t
        output = update!(pid,intor.u[2];dt=intor.t-intor.tprev)
        #intor.p = output
    end
end
pid_callback = DiscreteCallback((u,t,integrator)->true,pid_maker())
prob = SecondOrderODEProblem(fmaker(),ẋ0,x0,(0.0,100.0),0.0,callback=pid_callback)
sol = solve(prob,DPRKN6(),dt=0.01,adaptive=false)

plot(sol,label=["ẋ" "x"])
ylims!(-1,3)
