
using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using Test
using SafeTestsets

const GROUP = get(ENV, "GROUP", "ALL")

group_enabled(groups::AbstractString...) = GROUP == "ALL" || GROUP in groups

@time begin
    if group_enabled("ADJOINT")
        @time @safetestset "Tensegrity Adjoint Tests" include("adjoint.jl")
    end

    if group_enabled("STIFFNESS")
        @time @safetestset "Tensegrity Stiffness Tests" include("stiffness.jl")
    end
    
    # Add more test groups as needed
    if group_enabled("DYNAMICS", "DYNAMICS_A")
        @time @safetestset "Tensegrity Dynamics Tests" include("superball.jl")
    end

    if group_enabled("DYNAMICS", "DYNAMICS_B")
        @time @safetestset "Cluster Tensegrity Dynamics Tests" include("cluster.jl")
    end

    if group_enabled("ES")
         @time @safetestset "ES Examples" include("es_examples.jl")
    end
end
