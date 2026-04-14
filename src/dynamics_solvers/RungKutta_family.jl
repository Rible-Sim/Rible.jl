
struct RKIntegrator{tableauType} <: AbstractIntegrator
    tableau::tableauType
end

struct RungeKuttaConstants{T, TableauType}
    nq::Int
    nλ::Int
    ns::Int
    nx::Int
    ny::Int
    h::T
    tableau::TableauType
end

struct RungeKuttaJacobianWorkspace{SparseMT,MT, VT, VVT}
    M::SparseMT
    M⁻¹::SparseMT
    F::VT
    ∂F∂q::MT
    ∂F∂q̇::MT
    ∂F∂u::MT
    ∂F∂s::MT
    ∂S∂q::MT
    ∂S∂s::MT
    A::MT
    ∂Aᵀλ∂q::MT
    x::VT
    qpₖ₊₁::VT
    qpₖ::VT
    qp_rk_stages::VVT
    Res::VT
    Jac::MT
    JK::MT
end

struct RungeKutta_Constant_Mass_Cache{solT, T, TableauType, SparseMT, MT, VT, VVT, optionsType}
    solver::solT
    consts::RungeKuttaConstants{T, TableauType}
    jacobian_workspace::RungeKuttaJacobianWorkspace{SparseMT, MT, VT, VVT}
    options::optionsType
end

struct RungeKuttaSolverState{CoordinatesStateType, T}
    state_k::CoordinatesStateType
    state_kp1::CoordinatesStateType
    inst_state::CoordinatesStateType
    h::T
end

function compute_rungekutta_residual! end

function compute_rungekutta_jacobian! end