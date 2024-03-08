prob = RB.DynamicsProblem(
    bot_small, policy,
    rigid_contacts,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    ),
    RB.EulerEytelwein(),
)
solver = RB.DynamicsSolver(
    RB.Zhong06(),
    RB.InnerLayerContactSolver(
        RB.InteriorPointMethod()
    ),
    RB.MonolithicApparatusSolver(
        RB.SmoothedFischerBurmeister()
    ),
)
controller = ((prescribe!)=nothing, (actuate!)=nothing)
simulator = RB.Simulator(prob, solver, controller; tspan, dt)
solver_cache = RB.generate_cache(
    simulator, solver; dt
)

(; env, bot) = prob
(; traj, structure) = bot

qₖ₋₁ = traj.q[1]
q̇ₖ₋₁ = traj.q̇[1]
qˣ = qₖ₋₁ .+ dt./2 .*q̇ₖ₋₁
T = eltype(qₖ₋₁)

contact_cache = RB.activate_frictional_contacts!(structure,env,solver_cache,qˣ)
(; na) = contact_cache.cache
nΛ = 3na
Λₖ = zeros(T,nΛ)
RB.get_frictional_directions_and_positions!(structure, contact_cache, qₖ₋₁, q̇ₖ₋₁, Λₖ)

GeometryBasics.Mesh.Sphere((0.,0.,0.), 1.)
s1 = Meshes.Sphere((0.,0.,0.),1.)
mesh = Meshes.discretize(s1, Meshes.RegularDiscretization(30, 30))

m1 = mesh |> RB.simple2mesh