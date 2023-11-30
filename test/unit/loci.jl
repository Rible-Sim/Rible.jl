using Revise
import Rible as RB
using Test, Random
using StaticArrays
using Rotations
Random.seed!(100)

@testset for N = 2:3
    M = 2N-3
    T = Float64
    position = @SVector rand(N)
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = RB.Axes(normal)
    velocity = @SVector rand(N)
    ω = @SVector rand(M)
    @time "CartesianFrame $(N)D" frame1 = RB.CartesianFrame(position,axes,velocity,ω)
    frame2 = RB.CartesianFrame(position,velocity,axes,ω)
    @test  frame1 == frame2
    force  = @MVector rand(N)
    torque = @MVector rand(M)
    @time "ContactState $(N)D" cs = RB.ContactState(normal)
    @time "Locus $(N)D" locus = RB.Locus(position,normal)
    @time "LocusState $(N)D" locus_state1 = RB.LocusState(position,velocity,force,torque) 
    @time "LocusState $(N)D from frame $(N)D" locus_state2 = RB.LocusState(locus,frame1,force,torque)
end


