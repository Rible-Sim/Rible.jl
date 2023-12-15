function make_hammer(id,r̄ijkl,origin_position,R,ri,rj=nothing,rk=nothing,rl=nothing;
        μ,e,
        ω = zero(origin_position),
        contactable = true,
        visible = true,
        pres_idx = Int[],
        cstr_idx = collect(1:6),
    )
    # free_idx = collect(1:6)
    m = 4.58794901
    Īg = SMatrix{3,3}(
        Matrix(Diagonal([
            0.00522676,
            0.00512073,
            0.00512073,
            ])
        )
    )
    mass_locus = SVector(-0.20,0,0)
    # loci = deepcopy(r̄ijkl)
    loci = vcat(
        [
            SVector(0.22*cos(i*π/2),0.0,0.22*sin(i*π/2)) + mass_locus
            for i = -2:1
        ],
        [
            RotY( π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + mass_locus
            for i = 1:6
        ],
        [
            RotY(-π/4)*SVector(0.0,0.22*cos(i*π/3),0.22*sin(i*π/3)) + mass_locus
            for i = [1,2,4,5]
        ]
    )
    axes = [
        SVector(1.0,0,0) 
        for _ in eachindex(loci)
    ]
    friction_coefficients = fill(μ,length(loci))
    restitution_coefficients = fill(e,length(loci))
    @show m,diag(Īg),mass_locus,length(loci),friction_coefficients,restitution_coefficients
    prop = RB.RigidBodyProperty(
                id,contactable,m,Īg,
                mass_locus,loci,axes,
                friction_coefficients,restitution_coefficients;
                visible=visible
                )
    origin_velocity = zero(origin_position)
    if rj isa Nothing
        nmcs = RB.NCF.NC1P3V(ri,origin_position,R)
    elseif rk isa Nothing
        nmcs = RB.NCF.NC2P2V(ri,rj,origin_position,R)
    elseif rl isa Nothing
        nmcs = RB.NCF.NC3P1V(ri,rj,rk,origin_position,R)
    else
        nmcs = RB.NCF.NC4P(ri,rj,rk,rl,origin_position,R)
    end
    # cf = RB.NCF.CoordinateFunctions(nmcs,q0,pres_idx,free_idx)
    # @show typeof(nmcs)
    # trimesh = Meshes.discretize(box,) |> simple2mesh
    meteormesh = load(joinpath(RB.assetpath("流星锤.STL"))) |> RB.make_patch(;
        scale = 1/400,
        trans = [-0.20,0,0],
    )
    coords = RB.NonminimalCoordinates(
        nmcs, pres_idx,cstr_idx
    )
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)
    body = RB.RigidBody(prop,state,coords,meteormesh)
end