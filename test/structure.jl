
q = rand(28)
λ = rand(25)
∂Aᵀλ∂q(q,λ) = RB.cstr_forces_jacobian(sc.structure,q,λ)
A(q) = RB.cstr_jacobian(sc.structure,q)
A(q)
my = ∂Aᵀλ∂q(q,λ)
using FiniteDiff

ref = FiniteDiff.finite_difference_jacobian(
    (q) -> A(q)'*λ,
    q,
    Val(:central),
)
my ≈ ref