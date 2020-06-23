module TensegrityRobot

using SparseArrays
using LinearAlgebra
using Parameters
using StaticArrays
using NLsolve

include("naturalcoordinates.jl")
using .NaturalCoordinates
include("string.jl")

include("rigidbody.jl")
#include("rigidbody3d.jl")

include("tensegrity.jl")

include("inverse.jl")

include("control.jl")

struct TGRobot2D{ST,CT}
    tgstruct::ST
    hub::CT
end
# Write your package code here.

end
