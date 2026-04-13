

function build_mesh(body::AbstractRigidBody;update=true,color=nothing)
    (;mesh) = body
    if mesh isa GB.Mesh
        frame = to_3D(body.state.origin_frame)
        # if update
            origin_position = frame.position
            R = frame.axes.X
            # R = one(RotMatrix{3, Float64})
        # else
        #     origin_position = zero.(frame.position)
        #     R = one(frame.axes.X)
        # end
        trans = Translation(origin_position)
        rot = LinearMap(R)
        ct = trans ∘ rot
        updated_pos = GB.Point3f.(ct.(mesh.position))
        fac = GB.faces(mesh)
        nls = GB.normals(updated_pos,fac)
        
        colors = parse_mesh_color(mesh,color)
        
        return GB.Mesh(fac; position=updated_pos, normal=nls, color=colors)
    else
        # Return a simple default mesh for non-GB.Mesh objects
        # Create a simple cube mesh
        positions = [
            GB.Point3f(0.0, 0.0, 0.0),
            GB.Point3f(1.0, 0.0, 0.0),
            GB.Point3f(1.0, 1.0, 0.0),
        ]
        
        faces = [
            GB.TriangleFace(1, 2, 3), 
        ]
        
        normals = GB.normals(positions, faces)
        colors = parse_mesh_color(GB.Mesh(faces; position=positions), color)
        
        return GB.Mesh(faces; position=positions, normal=normals, color=colors)
    end
end

function parse_mesh_color(mesh::GB.Mesh,color=nothing)
    num_pos = mesh.position |> length
    if color isa CT.Colorant
        colors = fill(color,num_pos)
    elseif !(color isa Nothing)
        parsedcolor = parse(CT.RGB{Float32},color)
        colors = fill(parsedcolor,num_pos)
    elseif hasproperty(mesh,:color)
        colors = mesh.color
    else
        parsedcolor = parse(CT.RGB{Float32},:slategrey) 
        colors = fill(parsedcolor,num_pos)
    end
    colors
end

function make_patch(;trans=[0.0,0,0],rot=RotX(0.0),scale=1,color=nothing)
    function patch(mesh::GB.Mesh{Dim,T}) where {Dim,T}
        ct = Translation(trans) ∘ LinearMap(rot)
        updated_pos = ct.(mesh.position.*scale)
        fac = GB.faces(mesh)
        nls = GB.normals(updated_pos,fac;normaltype=GB.Vec{Dim,Float32})
        colors = parse_mesh_color(mesh,color)
        GB.Mesh(fac; position=updated_pos,normal=nls,color=colors)
    end
end


function getIsosurface(plane, rect::GB.HyperRectangle{3}; nSteps = 100, level=0.0)
     # Use this to set the density 
    (;origin, widths) = rect
    xr,yr,zr = ntuple(i->range(origin[i],origin[i]+widths[i],nSteps),3)
    A = [ 
        # cos(x)*sin(y)+cos(y)*sin(z)+cos(z)*sin(x)
        signed_distance(SVector(x,y,z),plane)
        for x in xr, y in yr, z in zr
    ]
    # @show extrema(xr) extrema(yr) extrema(zr)
    mc = MarchingCubes.MC(A,Int)
    MarchingCubes.march(mc,level)
    F = [GB.TriangleFace{Int64}(f) for f in mc.triangles]
    V = [GB.Point{3,Float64}(origin .+ widths .*p./nSteps) for p in mc.vertices]
    # F,
    GB.Mesh(V,F)
end

function get_groundmesh(f::Function,rect)
    GB.Mesh(f, rect,  MarchingCubes()) |> make_patch(;color = :snow)
end

function get_groundmesh(plane::Plane,rect)
    getIsosurface(plane, rect)  |> make_patch(;color = :snow)
end


function get_groundmesh(env_geo::Vector{<:ContactGeometry},rect)
    map(env_geo) do surface
        getIsosurface(surface, rect)
    end |> GB.merge |> make_patch(;color = :snow)
end

function get_groundmesh(static_env::StaticEnvironment,rect)
    get_groundmesh(static_env.geometry,rect)
end
function get_groundmesh(::Nothing,rect)
    plane_n = [0,0,1.0]
    plane_r = zeros(3)
    plane = Plane(plane_n,plane_r)
    get_groundmesh(plane,rect)
end


function endpoints2mesh(
        p1,
        p2;
        radius=norm(p2-p1)/40,
        n1=10,n2=2,
        color=:slategrey
    )
    
    start = Meshes.Point(p1...)
    finish = Meshes.Point(p2...)
    cyl_bar = Meshes.Cylinder(
        start,
        finish,
        radius
    )
    cylsurf_bar = Meshes.boundary(cyl_bar)
    # Meshes.sample(cylsurf_bar,Meshes.RegularSampling(10,3))
    cyl_bar_simple = Meshes.discretize(cylsurf_bar,Meshes.RegularDiscretization(n1,n2))
    simple2mesh(cyl_bar_simple,color)
end

function spbasis(n)
    a = abs.(n)
    if (a[1]≥0 && a[2]≥0) || (a[1]≤0 && a[2]≤0) 
        v = SVector(n[1]+1,n[2]-1,n[3])
    else
        v = SVector(n[1]-1,n[2]-1,n[3])
    end
    t = n×v |> normalize
    b = n×t |> normalize
    t,b
end

function build_mesh(fb::FlexibleBody{<:FlexibleBodyProperty},nsegs=100;color=:slategrey)
    (;state,coords,cache) = fb
    (;e) = cache
    ancs = coords
    (;L,radius) = ancs
    T = typeof(L)
    V = T <: AbstractFloat ? T : Float64
    _r = ANCF.make_r(ancs,e)
    _rₓ = ANCF.make_rₓ(ancs,e)

    sz = (10,nsegs)
    φmin, φmax = V(0), V(2π)
    xmin, xmax = V(0), V(L)
    δφ = (φmax - φmin) / sz[1]
    φrange = range(φmin, stop=φmax-δφ, length=sz[1])
    xrange = range(xmin, stop=xmax,    length=sz[2])

    function point(φ, x)
        o = _r(x) |> Meshes.Point
        n = _rₓ(x) |> normalize
        u,v = spbasis(n)
        R = [u v n]
        R*Meshes.Vec(radius*cos(φ), radius*sin(φ), 0.0) + o
    end

    points = IterTools.ivec(point(φ, x) for φ in φrange, x in xrange) |> collect
     # connect regular samples with quadrangles
    nx, ny = sz
    topo   = Meshes.GridTopology((nx-1, ny-1))
    middle = collect(Meshes.elements(topo))
    for j in 1:ny-1
        u = (j  )*nx
        v = (j-1)*nx + 1
        w = (j  )*nx + 1
        z = (j+1)*nx
        quad = Meshes.connect((u, v, w, z))
        push!(middle, quad)
    end

    connec = middle
    Meshes.SimpleMesh(points, connec) |> (x)->simple2mesh(x,color)
end


function simple2mesh(sp,color=:slategrey)
    dim   = Meshes.embeddim(sp)
    nvert = Meshes.nvertices(sp)
    nelem = Meshes.nelements(sp)
    verts = Meshes.vertices(sp)
    topo  = Meshes.topology(sp)
    elems = Meshes.elements(topo)

    # coords of vertices
    coords = Meshes.coords.(verts) .|> Meshes.CoordRefSystems.values .|> (x) -> ustrip.(x)
    # fan triangulation (assume convexity)
    tris4elem = map(elems) do elem
      I = Meshes.indices(elem)
      [[I[1], I[i], I[i+1]] for i in 2:length(I)-1]
    end

    # flatten vector of triangles
    tris = [tri for tris in tris4elem for tri in tris]
    points  = GB.Point.(coords)
    faces  = GB.TriangleFace{UInt64}.(tris)
    # T = Float64 # Bug
    T = Float32 # OK
    nls = GB.normals(points,faces;normaltype=GB.Vec{3,T})
    parsedcolor = parse(CT.RGB{Float32},color)
    colors = fill(parsedcolor,length(points))
    
    GB.Mesh(faces; position = points,normal=nls,color=colors)
end
