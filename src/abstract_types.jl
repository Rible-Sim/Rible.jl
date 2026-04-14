
# Coordinate system abstract types
abstract type AbstractCoordinates{N,T} end
abstract type AbstractMinimalCoordinates{N,T} <: AbstractCoordinates{N,T} end
abstract type AbstractNonminimalCoordinates{N,T} <: AbstractCoordinates{N,T} end


"""
Abstract interface for body properties.
"""
abstract type AbstractBodyProperty{N,T} end

"""
Abstract interface for body states.
$(TYPEDEF)
"""
abstract type AbstractBodyState{N,T} end

"""
Abstract interface for body caches.
"""
abstract type AbstractBodyCache end

"""
Abstract interface for bodies.
"""
abstract type AbstractBody{N,T} end



# Body abstract types
abstract type AbstractRigidBody{N,T} <: AbstractBody{N,T} end
abstract type AbstractRigidBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} <: AbstractBodyState{N,T} end

abstract type AbstractFlexibleBody{N,T} <: AbstractBody{N,T} end
abstract type AbstractFlexibleBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractFlexibleBodyState{N,T} <: AbstractBodyState{N,T} end

abstract type AbstractConnectivity end

abstract type AbstractCoordinatesState end

abstract type AbstractStructure{bType,stype,cType} end

abstract type AbstractStructureState end

# Common structs and types

"""
Inertia cache for storing mass matrices and related quantities.
$(TYPEDEF)
"""
mutable struct InertiaCache{MType,JType,GType}
    M::MType
    M⁻¹::MType
    ∂Mq̇∂q::JType
    ∂M⁻¹p∂q::JType
    "Coriolis force"
    Ṁq̇::GType
    "centrifugal force"
    ∂T∂qᵀ::GType
    dirty::Bool
end

function InertiaCache(bodies, cnt)
    (; num_of_full_coords, bodyid2sys_full_coords) = cnt
    T = get_numbertype(bodies)
    M = spzeros(T, num_of_full_coords, num_of_full_coords)
    
    foreach(bodies) do body
        memfull = bodyid2sys_full_coords[body.prop.id]
        M[memfull, memfull] .+= body.cache.inertia_cache.M
    end
    
    M⁻¹ = sparse(inv(Matrix(M)))

    ∂Mq̇∂q = spzeros(T, num_of_full_coords, num_of_full_coords)
    ∂M⁻¹p∂q = spzeros(T, num_of_full_coords, num_of_full_coords)
    Ṁq̇ = spzeros(T, num_of_full_coords)
    ∂T∂qᵀ = spzeros(T, num_of_full_coords)
    return InertiaCache(M, M⁻¹, ∂Mq̇∂q, ∂M⁻¹p∂q, Ṁq̇, ∂T∂qᵀ, true)
end


function Base.sort(tsc::TypeSortedCollection)
    sort!(reduce(vcat,tsc.data))
end



# --- Interfaces ---

"""
    AbstractStateExtractor

Extracts specific quantities of interest (x) from the robot state.
"""
abstract type AbstractStateExtractor end

"""
    AbstractFeaturizer

Transforms the extracted state x into a feature vector φ.
"""
abstract type AbstractFeaturizer end

"""
    AbstractTimeBasis

Abstract type for time-dependent basis functions.
"""
abstract type AbstractTimeBasis end

"""
    AbstractParamFun

Maps features φ to control actions u.
"""
abstract type AbstractParamFun end

abstract type AbstractPolicy end

abstract type AbstractContinuousPolicy <: AbstractPolicy end

abstract type AbstractDiscretePolicy <: AbstractPolicy end

abstract type OpenloopContinuousPolicy <: AbstractContinuousPolicy end
abstract type FeedbackContinuousPolicy <: AbstractContinuousPolicy end

abstract type AbstractTimePolicy <: OpenloopContinuousPolicy end

abstract type OpenloopDiscretePolicy <: AbstractDiscretePolicy end
abstract type FeedbackDiscretePolicy <: AbstractDiscretePolicy end

struct NoPolicy <: AbstractPolicy end

abstract type AbstractContactModel end
struct Contactless <: AbstractContactModel end
struct RestitutionFrictionCombined{restitutionType,frictionType} <: AbstractContactModel 
    restitution::restitutionType
    friction::frictionType
end

abstract type AbstractFrictionModel <: AbstractContactModel end
struct Frictionless <: AbstractFrictionModel end
struct CoulombFriction <: AbstractFrictionModel end
struct PolyhedralCoulombFriction <: AbstractFrictionModel end
struct MaximumDissipation <: AbstractFrictionModel end

abstract type AbstractRestitutionModel <: AbstractContactModel end
struct Inelastic <: AbstractRestitutionModel end
struct NewtonRestitution <: AbstractRestitutionModel end
struct PoissonRestitution <: AbstractRestitutionModel end
struct StrangeRestitution <: AbstractRestitutionModel end

abstract type AbstractApparatusModel end
struct Naive <: AbstractApparatusModel end

abstract type AbstractSolver end
abstract type AbstractDynamicsSolver <: AbstractSolver end
abstract type AbstractForwardDynamicsSolver <: AbstractDynamicsSolver end
abstract type AbstractAdjointDynamicsSolver <: AbstractDynamicsSolver end
abstract type AbstractDynamicsSensitivitySolver <: AbstractDynamicsSolver end

abstract type AbstractComplementaritySolver <: AbstractSolver end
struct InteriorPointMethod <: AbstractComplementaritySolver end
const IPM = InteriorPointMethod
struct AcceleratedProjectedGradientDescent <: AbstractComplementaritySolver end
const APGD = AcceleratedProjectedGradientDescent
struct AlternatingDirectionMethodofMultipliers <: AbstractComplementaritySolver end
const ADMM = AlternatingDirectionMethodofMultipliers
struct ProjectedGaussSeidel <: AbstractComplementaritySolver end
const PGS = ProjectedGaussSeidel
struct ProjectedGaussJacobi <: AbstractComplementaritySolver end
const PGJ = ProjectedGaussJacobi
struct SmoothedFischerBurmeister  <: AbstractComplementaritySolver end
struct SemismoothNewton  <: AbstractComplementaritySolver end
struct CableHole <: AbstractComplementaritySolver end

abstract type AbstractContactSolver <: AbstractSolver end
struct NoContactSolver <: AbstractContactSolver end

struct MonolithicContactSolver{complementarity_solverType} <: AbstractContactSolver
    complementarity_solver::complementarity_solverType
end
struct InnerLayerContactSolver{complementarity_solverType} <: AbstractContactSolver
    complementarity_solver::complementarity_solverType
end

abstract type AbstractBodySolver <: AbstractSolver end
struct NoBodySolver <: AbstractBodySolver end

abstract type AbstractApparatusSolver <: AbstractSolver end
struct NoApparatusSolver <: AbstractApparatusSolver end

struct MonolithicApparatusSolver{complementarity_solverType} <: AbstractApparatusSolver
    complementarity_solver::complementarity_solverType
end

struct InnerLayerApparatusSolver{complementarity_solverType} <: AbstractApparatusSolver
    complementarity_solver::complementarity_solverType
end

abstract type AbstractIntegrator end



abstract type AbstractMechanicsProblem end

abstract type AbstractDynamicsProblem <: AbstractMechanicsProblem end

# Constraint methods

"""
    cstr_function(coords, q, d)

Evaluate constraint functions.
"""
function cstr_function end

# Constraint methods

"""
    cstr_function(ret, coords, q, d)

Evaluate constraint functions.
"""
function cstr_function! end

"""
    cstr_jacobian(coords, q)

Evaluate constraint jacobian.
"""
function cstr_jacobian end

"""
    cstr_jacobian!(ret, coords, q)

Evaluate constraint jacobian.
"""
function cstr_jacobian! end

"""
    cstr_hessians(coords)

Get constraint hessians.
"""
function cstr_hessians end

"""
    cstr_forces_jacobian(coords, λ)

Get constraint forces jacobian.
"""
function cstr_forces_jacobian end

"""
    add_cstr_forces_jacobian!(ret, coords, λ)

Get constraint forces jacobian.
"""
function add_cstr_forces_jacobian! end

"""
    cstr_forces_jacobian!(ret, coords, λ)

Get constraint forces jacobian.
"""
function cstr_forces_jacobian! end

"""
    cstr_velocity_jacobian(coords, q̇)

Get constraint velocity jacobian.
"""
function cstr_velocity_jacobian end

"""
    cstr_velocity_jacobian!(ret, coords, q̇)

Get constraint velocity jacobian.
"""
function cstr_velocity_jacobian! end

"""
    make_cstr_function(coords, ...)

Make constraint function.
"""
function make_cstr_function end

"""
    make_cstr_jacobian(coords, ...)

Make constraint jacobian.
"""
function make_cstr_jacobian end


# Mass matrix methods

"""
Check if coordinate system has constant mass matrix.
$(TYPEDSIGNATURES)
"""
function has_constant_mass_matrix end


# Force methods

"""
    apply_field!(F, body, field, q)

Apply field forces to the body.
"""
function apply_field! end

"""
    field_jacobian!(∂F̌∂q̌, body, field, q)

Calculate field force jacobian.
"""
function field_jacobian! end

"""
    clear_forces!(body)

Clear all forces and torques on the body.
"""
function clear_forces! end

"""
    reset!(body)

Reset body state including contact states.
"""
function reset! end



"""
    get_num_of_intrinsic_cstr(coords)

Return the number of intrinsic constraints.
"""
function get_num_of_intrinsic_cstr end

"""
    get_num_of_dims(coords)

Return the spatial dimension of the coordinate system.
"""
function get_num_of_dims end

"""
    get_num_of_coords(coords)

Return the total number of coordinates.
"""
function get_num_of_coords end


const get_num_of_full_coords = get_num_of_coords

"""
    get_num_of_dof(coords)

Return the number of degrees of freedom.
"""
function get_num_of_dof end

"""
    get_num_of_cstr(coords)

Return the number of constraints.
"""
function get_num_of_cstr end
