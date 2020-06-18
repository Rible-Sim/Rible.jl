module Robot2D
using SparseArrays
using LinearAlgebra
using Parameters
using StaticArrays
using NLsolve
include("NC.jl")
include("string.jl")
include("rigidbody.jl")

include("tensegrity.jl")

include("inverse.jl")

include("control.jl")

struct TGRobot2D{ST,CT}
    st2d::ST
    hub::CT
end


end # module
