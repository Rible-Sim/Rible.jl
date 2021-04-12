
uniform_1 = Uniform(-1.0,1.0)

function sample_ball_Muller(uniform_1,d)
    u = rand(uniform_1,d)
    r = rand()^(1/d)
    r./norm(u).*u
end
sample_ball_Muller(uniform_1,2)
fig,ax,plotscatter = scatter([Point2(sample_ball_Muller(uniform_1,2)) for i = 1:4000])
ax.aspect = DataAspect()
fig
exp0_5 = Exponential(0.5)
function sample_ball_exp(uniform_1,exp0_5,d)
    x = rand(uniform_1,d)
    denom = sqrt(rand(exp0_5) + sum(x.^2))
    x./denom
end

sample_ball_exp(uniform_1,exp0_5,2)
fig,ax,plotscatter = scatter([Point2(sample_ball_exp(uniform_1,exp0_5,2)) for i = 1:4000])
ax.aspect = DataAspect()
fig
function sample_ball_drop(uniform_1,d)
    u = rand(uniform_1,d+2)
    u ./= norm(u)
    u[1:d]
end
sample_ball_drop(uniform_1,2)
fig,ax,plotscatter = scatter([Point2(sample_ball_drop(uniform_1,2)) for i = 1:4000])
ax.aspect = DataAspect()
fig
n = 100
@code_warntype sample_ball_Muller(uniform_1,n)
@code_warntype sample_ball_exp(uniform_1,exp0_5,n)
@code_warntype sample_ball_drop(uniform_1,n)

@btime sample_ball_Muller(uniform_1,n)
@btime sample_ball_exp(uniform_1,exp0_5,n)
@btime sample_ball_drop(uniform_1,n)
