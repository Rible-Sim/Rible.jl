function make_meteor_hammer(;
        ri  = [ 0.0, 0.0, 1.5],
        rix = [ 0.0, 0.0,-1.0],
        rj  = [-0.5*√3, 0.0, 0.5-1e-2],
        rjx = [ 0.0,-1.0, 0.0],
        doDR=false,
        μ = 0.5,
        e = 0.9,
        L = 1.0,
        nx = 2,
        R = RotY(deg2rad(-60)),
        ω = zero(ri)
    )
    if doDR
        fb_pres_idx = [7,8,9] # the end that attached to the hammer
        rb_pres_idx = [1,2,3]
    else
        fb_pres_idx = Int[]
        rb_pres_idx = Int[]
    end
    𝐞 = [ri;rix;rj;rjx]
    fb1 = cable_ancf(fb_pres_idx, 𝐞, L)
    # @show fb1.prop.mass
    subfbs,subsm = RB.subdivide(fb1,nx)
    r̄ijkl = SVector{3,Float64}.(
        0.1 .*[
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
        ]
    )
    rigidbody = make_hammer(
        nx+1,
        r̄ijkl,
        rj,
        R,
        rj;
        μ,e,
        ω,
        pres_idx=rb_pres_idx,
        visible=ifelse(!isempty(rb_pres_idx),true,false)
    )
    bodies = TypeSortedCollection(vcat(subfbs,rigidbody))
  
    cst1 = RB.FixedIndicesConstraint(1,subfbs[1],[1,2,3],ri)
    apparatuses = TypeSortedCollection(
        [cst1]
    )

    sm = zeros(Int,size(subsm,1)+3,size(subsm,2)+1)
    sm[1:size(subsm,1),1:size(subsm,2)] = subsm
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)  ] = 7:9
    sm[size(subsm,1)+1:size(subsm,1)+3,size(subsm,2)+1] = 1:3
    # display(sm)
    indexed = RB.index(bodies,apparatuses,sm)
    numbered = RB.number(bodies,apparatuses)

    cnt = RB.Connectivity(numbered,indexed,)
    st = RB.Structure(bodies,apparatuses,cnt,)
    bot = RB.Robot(st)
end
