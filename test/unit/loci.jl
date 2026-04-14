import Rible as RB
using Test, Random
using StaticArrays
using Rotations
Random.seed!(100)

@testset for N = 2:3
    M = 2N-3
    T = Float64
    position = @MVector rand(N)
    normal = SVector{N}(ifelse(i==1,one(T),zero(T)) for i = 1:N)
    axes = RB.Axes(normal)
    velocity = @MVector rand(N)
    ω = @MVector rand(M)
    @time "CartesianFrame $(N)D" frame1 = RB.CartesianFrame(position,axes,velocity,ω)
    frame2 = RB.CartesianFrame(position,velocity,axes,ω)
    @test  frame1.position == frame2.position
    @test  frame1.axes == frame2.axes
    @test  frame1.velocity == frame2.velocity
    @test  frame1.angular_velocity == frame2.angular_velocity
    @test  frame1.local_angular_velocity == frame2.local_angular_velocity
    force  = @MVector rand(N)
    torque = @MVector rand(M)
    @time "ContactState $(N)D" cs = RB.ContactState(normal)
    @time "Locus $(N)D" locus = RB.Locus(SVector(position),normal)
    @time "LocusState $(N)D" locus_state1 = RB.LocusState(position,velocity,force,torque) 
    @time "LocusState $(N)D from frame $(N)D" locus_state2 = RB.LocusState(locus,frame1,force,torque)
end

@testset "Anchor locus-index constructors" begin
    p1 = SVector(1.0, 0.0, 0.0)
    p2 = SVector(0.0, 1.0, 0.0)
    n1 = SVector(0.0, 0.0, 1.0)
    n2 = SVector(0.0, 1.0, 0.0)
    l1 = RB.Locus(p1, n1)
    l2 = RB.Locus(p2, n2)
    body = (prop = (loci = [l1, l2],),)

    a = RB.Anchor(body, 1, 2)
    @test a.position == l1.position
    @test a.trl_axes == l2.axes
    @test a.rot_axes == l2.axes
end

@testset "NC3D1P1V intrinsic constraint count" begin
    ri = SVector(0.0, 0.0, 0.0)
    u = SVector(1.0, 0.0, 0.0)
    nmcs = RB.NCF.NC3D1P1V(ri, u)
    @test RB.get_num_of_cstr(nmcs) == 1

    q = vcat(ri, u)
    ret = zeros(RB.get_num_of_cstr(nmcs))
    RB.cstr_function!(ret, nmcs, q)
    @test ret ≈ [0.0]

    nmcs_default = RB.NCF.NC3D1P1V(ri, u; cstr_idx=:)
    @test RB.get_cstr_idx(nmcs_default) == [1]

    nmcs_empty = RB.NCF.NC3D1P1V(ri, u; cstr_idx=Int[])
    @test RB.get_num_of_cstr(nmcs_empty) == 0
    @test RB.get_cstr_idx(nmcs_empty) == Int[]

    @test_throws ArgumentError RB.NCF.NC3D1P1V(ri, u; cstr_idx=[2])
end
