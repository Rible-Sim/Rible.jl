# Abstract types for environment fields
abstract type AbstractField end
abstract type AbstractGravity <: AbstractField end
struct NoField <: AbstractField end

# Abstract types for geometries
abstract type AbstractGeometry end
struct NoGeometry <: AbstractGeometry end

abstract type AbstractEnvironment end

struct StaticEnvironment{geometryType,fieldType} <: AbstractEnvironment
    geometry::geometryType
    field::fieldType
end

EmptyEnv() = StaticEnvironment(
    NoGeometry(),
    NoField(),
)

abstract type AbstractContactEnvironment <: AbstractEnvironment end
abstract type RigidBodyContactEnvironment <: AbstractContactEnvironment end

has_gravity(::NoField) = false
has_gravity(env::AbstractEnvironment) = has_gravity(env.field)

iscontact(env::StaticEnvironment) = true
iscontact(env::StaticEnvironment{NoGeometry}) = false
iscontact(env::AbstractContactEnvironment) = true


struct Gravity{T} <: AbstractGravity
    acceleration::T
    function Gravity(g=9.81)
        new{typeof(g)}(g)
    end
end

GravityEnv(g=9.81) = StaticEnvironment(
    NoGeometry(),
    Gravity(g),
)



has_gravity(::Gravity) = true
