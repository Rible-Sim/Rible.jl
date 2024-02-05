using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/adjoint.jl"))

using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\IMSD 2025\LaTex_Abstract" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/IMSD 2024/LaTex_Abstract" #src
end #src
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl
tw = 455.8843 #pt |> pt2px
scalefactor = 4

#--  slider crank 
include(joinpath(pathof(RB),"../../examples/robots/cartpole.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/cartpole.jl"))#jl

coordsType=RB.QCF.QC
coordsType=RB.NCF.NC

cp_terminal = cart_pole(;coordsType)
f(t) = [1.0]

x = rand(ns)
jac = Zygote.jacobian(x -> Lux.apply(actor, x, ps, st)[1], rand(ns))[1]

zero(ComponentArray(ps))
axes(ComponentArray(ps))

ComponentArray(ps)[:]
Lux.apply(actor, x, ps, st)[1]
pa = ComponentArray(ps) 
pa += ComponentArray(ps)[:]
ps, st = Lux.setup(rng,actor)
q = RB.get_coords(cp_terminal.structure)
q̇ = RB.get_velocs(cp_terminal.structure)

control = x -> Lux.apply(actor,x,ps,st)[1]
control_jac = (s) -> Zygote.jacobian(control, s)[1]
control(vcat(q,q̇))
control_jac(vcat(q,q̇))

π_θ = (p) -> Lux.apply(actor,x,ps,st)[1]

jac = (p) -> Zygote.jacobian(
    control, 
    s
)[1]

cp_sim = cart_pole(;θ0=π/4,y0 = -1.0,coordsType,)
rng = Random.default_rng()
Random.seed!(rng, 0)
nq = RB.get_num_of_free_coords(cp_sim.structure)
ns = 2nq
na = 1
actor = Chain(
    Dense(ns, 4, elu;),
    Dense(4, na, ;),
)
ps,st = Lux.setup(rng,actor)
pa = deepcopy(ComponentArray(ps))
cp_policy = RB.ActorPolicy(
    @eponymtuple(
        actor,
        ps = deepcopy(pa),
        st,
    )
)

function forward_dyn(θ=ones(10),params=nothing;
        tspan = (0.0,0.1),
        dt = 1e-2
    )
    ## times = tspan[1]+dt/2:dt:tspan[2]
    ## f(t) = [extrapolate(scale(interpolate(u, BSpline(Constant())), times), Line())(t)]
    bot = cart_pole(;θ0=π/4,y0 = -1.0,coordsType,)
    RB.solve!(
        RB.DynamicsProblem(bot,actor),
        RB.DynamicsSolver(
            RB.Zhong06()
        );
        tspan,
        dt
    )
    bot
end

function get_J(bot,policy)
    function inner_J(θ,params=nothing)
        policy.nt.ps[:] .= θ
        tspan = (0.0,1.0)
        dt = 1e-2
        @show θ
        RB.solve!(
            RB.DynamicsProblem(bot,policy),
            RB.DynamicsSolver(
                RB.Zhong06()
            );
            tspan,
            dt,
            ftol=1e-7
        )
        RB.cost!(bot,bot.traj.q[end],bot.traj.q̇[end],bot.traj.t[end])
    end
end

function get_dJ!(bot,policy)
    function inner_dJ!(G,θ,params=nothing)
        policy.nt.ps[:] .= θ
        @show θ
        tspan = (0.0,1.0)
        dt = 1e-2
        dsprob = RB.DynamicsSensitivityProblem(bot,policy)
        dsolver = RB.DynamicsSolver(
            RB.Zhong06()
        )
        adsolver = RB.DiscreteAdjointDynamicsSolver(
            dsolver
        )
        dssolver = RB.AdjointDynamicsSensitivitySolver(
            dsolver,
            adsolver
        )
        ## adprob  = RB.AdjointDynamicsProblem(bot,nothing)  )
        _,solvercache = RB.solve!(
            dsprob,
            dssolver;
            tspan,
            dt,
            ftol=1e-7
        )
        G .= sum(solvercache.cache.∂J∂θᵀ)
    end
end

get_J(cp_sim,cp_policy)(pa)

G_fd = FiniteDiff.finite_difference_gradient(get_J(cp_sim,cp_policy), pa)

G = zeros(axes(pa))
get_dJ!(cp_sim,cp_policy)(G,pa)
G

rel_err = (G .- G_fd .|> abs) ./ abs.(G) 
gidx = findall(isfinite,rel_err[:]) 
rel_err[gidx]
(G .- G_fd)[gidx]
G[gidx]
 ((x)->findall(isfinite,x)) |> maximum

optprob = OptimizationFunction(get_J(cp_sim,cp_policy), Optimization.AutoFiniteDiff())
prob = Optimization.OptimizationProblem(optprob, deepcopy(pa), )
sol = solve(prob, Optim.LBFGS())
sol.u

optprob = OptimizationFunction(get_J(cp_sim,cp_policy), grad = get_dJ!(cp_sim,cp_policy))
prob = Optimization.OptimizationProblem(optprob, deepcopy(pa),)
sol = solve(prob, Optim.BFGS(alphaguess = Optim.InitialStatic(;alpha=0.01), linesearch = Optim.HagerZhang()))
sol.u
sol
get_J(cp_sim,cp_policy)(ComponentArray(θ))
get_dJ!(cp_sim,cp_policy)(G,ComponentArray(θ))


tspan = (0.0,1.0)
dt = 1e-2
cp_policy.nt.ps[:] .= sol.u
RB.solve!(
    RB.DynamicsProblem(cp_sim,cp_policy),
    RB.DynamicsSolver(
        RB.Zhong06()
    );
    tspan,
    dt
)
RB.cost!(cp_sim,cp_sim.traj.q[end],cp_sim.traj.q̇[end],cp_sim.traj.t[end])

plot_traj!(cp_sim;showmesh=false,showground=false)
# gauges terminal
m1 = RB.measure(cp_terminal.structure,cp_terminal.hub.gauges.data[1][1])
## m2 = RB.measure(cp_terminal,cp_terminal.hub.gauges.data[2][1])

# gauges
m1 = RB.measure(cp_sim.structure,cp_sim.hub.gauges.data[1][1])
## m2 = RB.measure(cp_sim,cp_sim.hub.gauges.data[2][1])

# cost in errors
RB.cost!(cp_sim,cp_sim.traj.q[end],cp_sim.traj.q̇[end],cp_sim.traj.t[end])

RB.error_jacobian(cp_sim)
## RB.measure_jacobian(cp_sim,cp_sim.hub.gauges.data[2][1])

RB.cost_jacobian!(cp_sim,)

# cost in actuators
RB.get_num_of_actions(cp_sim.hub.actuators.data[1][1])
cp_sim.hub.actuators.data[1][1]
RB.generalized_force(cp_sim,cp_sim.hub.actuators.data[1][1])
RB.generalized_force(cp_sim,cp_sim.hub.actuators.data[2][1])
RB.actions_jacobian(cp_sim,cp_sim.hub.actuators.data[1][1])
RB.actions_jacobian(cp_sim,cp_sim.hub.actuators.data[2][1])

RB.actions_jacobian(cp_terminal,cp_terminal.hub.actuators.data[2][1])

jac = Zygote.jacobian(p -> Lux.apply(actor, x, p, state_in)[1], ComponentArray(ps) )[1]

plot_traj!(bot;showmesh=false,showground=false)


# Construct the layer
model = Chain(BatchNorm(128), Dense(128, 256, tanh), BatchNorm(256),
              Chain(Dense(256, 1, tanh), Dense(1, 10)))

# Get the device determined by Lux
device = gpu_device()

# Parameter and State Variables
ps, st = Lux.setup(rng, model) .|> device

# Dummy Input
x = rand(rng, Float32, 128, 2) |> device

# Run the model
y, st = Lux.apply(model, x, ps, st)
y
# Gradients
jac = Zygote.jacobian(p -> Lux.apply(model, x, p, st)[1], caps)[1]

using ComponentArrays

caps = ps |> ComponentArray

Lux.apply(model, x, caps, st)

import ReinforcementLearningCore as RLCore
import ReinforcementLearningBase as RLBase
using ReinforcementLearningCore
using ReinforcementLearningBase
using StableRNGs
using Distributions
using Flux: glorot_uniform
using Flux
seed = 123
rng = StableRNG(seed)
ns = 1
na = 1
actor = Chain(
    Dense(ns, 256, relu; init = glorot_uniform(rng)),
    Dense(256, na; init = glorot_uniform(rng)),
)
critic = Chain(
    Dense(ns, 256, relu; init = glorot_uniform(rng)),
    Dense(256, 1; init = glorot_uniform(rng)),
)
optimizer = ADAM(1e-3)
approximator = RLCore.ActorCritic(;
    actor,
    critic,
    optimizer,
) |> gpu
policy = PPOPolicy(
    approximator,
    γ = 0.99f0,
    λ = 0.95f0,
    clip_range = 0.1f0,
    max_grad_norm = 0.5f0,
    n_epochs = 4,
    n_microbatches = 4,
    actor_loss_weight = 1.0f0,
    critic_loss_weight = 0.5f0,
    entropy_loss_weight = 0.001f0,
    update_freq = UPDATE_FREQ,
)

using Lux, Random, Optimisers, Zygote
# using LuxCUDA, LuxAMDGPU # Optional packages for GPU support

# Seeding
rng = Random.default_rng()
Random.seed!(rng, 0)

# Construct the layer
model = Chain(BatchNorm(128), Dense(128, 256, tanh), BatchNorm(256),
              Chain(Dense(256, 1, tanh), Dense(1, 10)))

# Get the device determined by Lux
device = gpu_device()

# Parameter and State Variables
ps, st = Lux.setup(rng, model) .|> device

# Dummy Input
x = rand(rng, Float32, 128, 2) |> device

# Run the model
y, st = Lux.apply(model, x, ps, st)

# Gradients
gs = gradient(p -> sum(Lux.apply(model, x, p, st)[1]), ps)[1]

# Optimization
st_opt = Optimisers.setup(Optimisers.Adam(0.0001), ps)
st_opt, ps = Optimisers.update(st_opt, ps, gs)

rng = Random.default_rng()
Random.seed!(rng, 0)
x = [1, 2, 3]

"""
    ActorCritic(;actor, critic, optimizer=Adam())
The `actor` part must return logits (*Do not use softmax in the last layer!*), and the `critic` part must return a state value.
"""
Base.@kwdef struct ActorCritic{A,C,O}
    actor::A
    critic::C
end

@functor ActorCritic

#####
# GaussianNetwork
#####

export GaussianNetwork

"""
    GaussianNetwork(;pre=identity, μ, σ, min_σ=0f0, max_σ=Inf32)

Returns `μ` and `σ` when called.  Create a distribution to sample from using
`Normal.(μ, σ)`. `min_σ` and `max_σ` are used to clip the output from
`σ`. `pre` is a shared body before the two heads of the NN. σ should be > 0. 
You may enforce this using a `softplus` output activation. The `squash` function is
applied elementwise to the action. If squash is `tanh`, a correction is applied to the
logpdf. Other squashing functions are not supported except for identity.
"""

logpdfcorrection(z, ::F) where F <: typeof(tanh) = -sum(log.(1 .- tanh.(z).^2), dims = 1)
logpdfcorrection(s, f) = 0
inversesquash(::F) where F <: typeof(tanh) = atanh
inversesquash(::F) where F <: typeof(identity) = identity

Base.@kwdef struct GaussianNetwork{P,U,S,F}
    pre::P = identity
    μ::U
    σ::S
    min_σ::Float32 = 0.0f0
    max_σ::Float32 = Inf32
    squash::F
end

GaussianNetwork(pre, μ, σ; squash = identity) = GaussianNetwork(pre, μ, σ, 0.0f0, Inf32, squash)

ns = 72
pre = Chain(
    Dense(ns, 64, relu;),
    Dense(64, 64, relu;),
)
μ = Chain(Dense(64, 1, tanh;), vec)
σ = Chain(Dense(64, 1, softplus;), vec)

model =  GaussianNetwork(pre, μ, σ,)

pre_ps, pre_st = Lux.setup(rng,model.pre)
μ_ps, μ_st = Lux.setup(rng,model.μ)
σ_ps, σ_st = Lux.setup(rng,model.σ)
s = rand(Float32,72)
x,_ = Lux.apply(model.pre,s,pre_ps,pre_st)
μ,_ = Lux.apply(model.μ,x,μ_ps,μ_st)
raw_σ,_ = Lux.apply(model.σ,x,σ_ps,σ_st)
σ = clamp.(raw_σ, model.min_σ, model.max_σ)

z = ignore_derivatives() do
    noise = randn(rng, Float32, size(μ))
    μ .+ σ .* noise
end
model.squash.(z)

noise = randn(rng, Float32, size(μ))
μ .+ σ .* noise

@functor GaussianNetwork

"""
This function is compatible with a multidimensional action space.

- `rng::AbstractRNG=Random.default_rng()`
- `is_sampling::Bool=false`, whether to sample from the obtained normal distribution. 
- `is_return_log_prob::Bool=false`, whether to calculate the conditional probability of getting actions in the given state.
"""
function (model::GaussianNetwork)(rng::AbstractRNG, s; is_sampling::Bool=false, is_return_log_prob::Bool=false)
    x = model.pre(s)
    μ, raw_σ = model.μ(x), model.σ(x)
    σ = clamp.(raw_σ, model.min_σ, model.max_σ)
    if is_sampling
        z = ignore_derivatives() do
            noise = randn(rng, Float32, size(μ))
            μ .+ σ .* noise
        end
        if is_return_log_prob
            logp_π = diagnormlogpdf(μ, σ, z) .+ logpdfcorrection(z, typeof(model.squash))
            return model.squash.(z), logp_π
        else
            return model.squash.(z)
        end
    else
        return μ, σ
    end
end

model(rng,[1.0])


"""
    (model::GaussianNetwork)(rng::AbstractRNG, state::AbstractArray{<:Any, 3}, action_samples::Int)

Sample `action_samples` actions from each state. Returns a 3D tensor with dimensions `(action_size x action_samples x batchsize)`.
`state` must be 3D tensor with dimensions `(state_size x 1 x batchsize)`. Always returns the logpdf of each action along.
"""
function (model::GaussianNetwork)(rng::AbstractRNG, s::AbstractArray{<:Any, 3}, action_samples::Int)
    x = model.pre(s)
    μ, raw_σ = model.μ(x), model.σ(x)
    σ = clamp.(raw_σ, model.min_σ, model.max_σ)
    z = ignore_derivatives() do
        noise = randn(rng, Float32, (size(μ, 1), action_samples, size(μ, 3))...)
        μ .+ σ .* noise
    end
    logp_π = diagnormlogpdf(μ, σ, z) .+ logpdfcorrection(z, typeof(model.squash))
    return model.squash.(z), logp_π
end

function (model::GaussianNetwork)(state; is_sampling::Bool=false, is_return_log_prob::Bool=false)
    model(Random.default_rng(), state; is_sampling=is_sampling, is_return_log_prob=is_return_log_prob)
end

function (model::GaussianNetwork)(state, action_samples::Int)
    model(Random.default_rng(), state, action_samples)
end

function (model::GaussianNetwork)(state, action)
    x = model.pre(state)
    μ, raw_σ = model.μ(x), model.σ(x)
    σ = clamp.(raw_σ, model.min_σ, model.max_σ)
    logp_π = diagnormlogpdf(μ, σ, inversesquash(model.squash).(action)) .+ logpdfcorrection(inversesquash(model.squash).(action), typeof(model.squash))
    return logp_π
end


using Lux
using Random

x_dim = 72
num_actions = 1

policy_network = Chain(
    Dense(x_dim, 128, relu),
    Dense(128, num_actions)
)
# Define the optimizer
## optimizer = ADAM(1e-3)

rng = Random.default_rng()
Random.seed!(rng, 0)
ps, st = Lux.setup(rng,policy_network)

# Function to sample an action from the policy
function sample_action(x)
    probs = softmax(Lux.apply(policy_network,x,ps,st)[1])
    action = rand(Categorical(probs))
    return action
end

sample_action(rand(72))

# Loss function for policy gradient
function policy_gradient_loss(states, actions, rewards)
    log_probs = logsoftmax(policy_network(states))
    selected_log_probs = sum(log_probs .* Flux.onehotbatch(actions, 1:num_actions), dims=1)
    loss = -sum(rewards .* selected_log_probs)
    return loss
end
# Function to run one episode and collect training data
function run_episode(env)
    states = []
    actions = []
    rewards = []
    state = reset!(env)
    done = false
    while !done
        action = sample_action(state)
        next_state, reward, done, _ = step!(env, action)
        push!(states, state)
        push!(actions, action)
        push!(rewards, reward)
        state = next_state
    end
    return states, actions, rewards
end
# Function to calculate the discounted rewards
function discount_rewards(rewards, gamma=0.99)
    discounted_rewards = zeros(length(rewards))
    running_add = 0
    for t in reverse(1:length(rewards))
        running_add = running_add * gamma + rewards[t]
        discounted_rewards[t] = running_add
    end
    return discounted_rewards
end
# Training loop
num_episodes = 1000
for i in 1:num_episodes
    states, actions, rewards = run_episode(env)
    discounted_rewards = discount_rewards(rewards)
    # Normalize rewards for better stability
    discounted_rewards = (discounted_rewards .- mean(discounted_rewards)) ./ std(discounted_rewards)
    # Calculate loss and perform a gradient update
    loss = policy_gradient_loss(hcat(states...), actions, discounted_rewards)
    grads = gradient(() -> loss, Flux.params(policy_network))
    Flux.Optimise.update!(optimizer, Flux.params(policy_network), grads)
    
    # Logging
    println("Episode: $i, Total Reward: $(sum(rewards))")
end