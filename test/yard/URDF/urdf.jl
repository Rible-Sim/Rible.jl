using LinearAlgebra
using StaticArrays
using EzXML
using ColorTypes
using Makie
import GLMakie as GM
GM.activate!()
import GeometryBasics as GB
using CoordinateTransformations
using Rotations
using EponymTuples
using RowEchelon
using TypeSortedCollections
using FileIO, MeshIO
using Meshing
using Unitful, Match
using LaTeXStrings
using Cthulhu
using JET
using Printf
using AbbreviatedStackTraces
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true
ENV["JULIA_STACKTRACE_MINIMAL"] = true
using Revise
import Rible as RB
include("../vis.jl"); includet("../vis.jl")
#-- preamble

function parse_vector(T,v)
    parse.(T,split(v))
end

function parse_origin(T,origin)
    default = @SVector zeros(T,3)
    # roll-pitch-yaw
    rpy = default
    xyz = default
    if origin !== nothing
        if haskey(origin,"rpy")
            rpy = parse_vector(T,origin["rpy"]) |> SVector{3}
        end
        if haskey(origin,"xyz")
            xyz = parse_vector(T,origin["xyz"]) |> SVector{3}
        end
    end
    RB.URDF.Origin(xyz,rpy)
end

function parse_material(T,material)    
    default = @SVector zeros(T,3)
    name = nothing
    color = nothing
    texture = nothing
    if material !== nothing
        name = material["name"]
        color_element = findfirst("color",material)
        if color_element !== nothing
            if haskey(color_element,"rgba")
                color_vector = parse_vector(T,color_element["rgba"]) |> SVector{4}
                color = Makie.RGBAf(color_vector...)
            end
        end

        texture_element = findfirst("texture",material)
        if texture_element !== nothing
            texture = texture_element["filename"]
        end
    end
    RB.URDF.Material(name,color,texture)
end

function parse_links(T,urdf)    
    xroot = root(urdf)
    @assert nodename(xroot) == "robot"    
    robot_materials = [   
        begin
            material = parse_material(T,material_element)
            material.name => material
        end
        for material_element in findall("/robot/material",urdf,)
    ] |> Dict
    links_elements = findall("/robot/link",urdf,)
    linksdict = [
        begin
            linkname = link["name"]
            inertial = findfirst("inertial",link)                
            inertial_origin = parse_origin(Float64,nothing)
            mass_value = 1.0
            im = SMatrix{3,3}(Matrix(0.01*I,3,3))
            if inertial !== nothing
                inertial_origin = parse_origin(Float64,
                    findfirst("origin",inertial)
                )
                mass = findfirst("mass",inertial)
                mass_value_parsed = parse(Float64,mass["value"])
                if mass_value_parsed != 0.0
                    mass_value = mass_value_parsed
                    inertia = findfirst("inertia",inertial)
                    ixx = parse(Float64,inertia["ixx"])
                    ixy = parse(Float64,inertia["ixy"])
                    ixz = parse(Float64,inertia["ixz"])
                    iyy = parse(Float64,inertia["iyy"])
                    iyz = parse(Float64,inertia["iyz"])
                    izz = parse(Float64,inertia["izz"])
                    im = SMatrix{3,3}([
                        ixx ixy ixz;
                        ixy iyy iyz;
                        ixz iyz izz
                    ])
                end
            end
            inertial = RB.URDF.Inertial(
                inertial_origin,
                mass_value,
                im
            )
            # origin::Origin{T}
            # mass::T
            # inertia::SMatrix{3,3,T}
            # @show inerital_origin
            # @show mass_value
            # @show im
            visual = findfirst("visual",link)
            if visual !== nothing
                visual_name = findfirst("name",visual)
                visual_origin = parse_origin(Float64,
                    findfirst("origin",visual)
                )
                visual_geometry = findfirst("geometry",visual)
                g = firstelement(visual_geometry)
                if nodename(g) == "mesh"
                    filename = g["filename"]
                    filename_striped = replace(filename,r"^package:\/\/"=>s"")
                    filepath,ext = filename_striped |> splitext
                    fileobj = filepath*".obj"
                    if !isfile(fileobj)
                        error("$fileobj no found")
                    end
                    if haskey(g,"scale")
                        scale = parse_vector(Float64,g["scale"]) |> SVector{3}
                    else
                        scale = SVector{3}(1.0,1.0,1.0)
                    end
                    visual_mesh_geometry = RB.URDF.MeshGeometry(
                        fileobj,
                        scale
                    )
                else
                    @warn "Found $(nodename(g)) geometry"
                    visual_mesh_geometry = nothing
                end
                if visual_name isa Nothing
                    visual_name = "NA"
                end
                material_element = findfirst("material",visual)
                if material_element isa Nothing
                    visual_material = nothing
                else
                    material_parsed = parse_material(T,material_element)
                    if material_parsed.name === nothing
                        visual_material = nothing
                    elseif material_parsed.color === material_parsed.texture === nothing
                        visual_material = robot_materials[material_parsed.name]
                    else
                        visual_material = material_parsed
                    end
                end
                visual = RB.URDF.Visual(                
                    visual_name,
                    visual_origin,
                    visual_mesh_geometry,
                    visual_material
                )
            end
            # collisions = findall("collision",link)
            # foreach(collisions) do collision
            #     if collision !== nothing
            #         collision_origin = parse_origin(Float64,
            #             findfirst("origin",collision)
            #         )
            #         collision_geometry = findfirst("geometry",collision)
            #         for g in eachelement(collision_geometry)
            #             if nodename(g) == "box"
            #                 size = parse_vector(Float64,g["size"])
            #                 # @show size
            #             elseif nodename(g) == "cylinder"
            #                 length = parse(Float64,g["length"])
            #                 radius = parse(Float64,g["radius"])
            #                 # @show length, radius
            #             elseif nodename(g) == "sphere"
            #                 radius = parse(Float64,g["radius"])
            #                 # @show radius
            #             end

            #         end

            #         # @show collision_origin
            #     end
            # end
            linkname => RB.URDF.Link{Float64}(
                linkname,
                inertial,
                visual,
                nothing
            )
        end
        for (ilink,link) in enumerate(links_elements)
    ] |> Dict
end

function parse_joints(T,urdf)
    xml_joints = findall("/robot/joint",urdf,)
    joints = [
        begin
            jointname = joint["name"]
            jointtype = joint["type"]
            # eles = elements(joint)
            # @show nodecontent(joint)
            parent = findfirst("parent",joint,)["link"]
            child = findfirst("child",joint,)["link"]
            origin = parse_origin(T,findfirst("origin",joint,))
            if jointtype in ["revolute","prismatic"]
                axis = findfirst("axis",joint,)
                axis_parsed = parse_vector(T,axis["xyz"],) |> SVector{3}
                axis_parent = inv(RotXYZ(origin.rpy...))*axis_parsed
                limit = findfirst("limit",joint,)
                limit_parsed = (
                    # command_effort = parse(T,limit["command_effort"]),
                    # current = parse(T,limit["current"]),
                    effort = parse(T,limit["effort"]),
                    # gear_velocity = parse(T,limit["gear_velocity"]),
                    lower = parse(T,limit["lower"]),
                    upper = parse(T,limit["upper"]),
                    velocity = parse(T,limit["velocity"]),
                )
                dynamics = findfirst("dynamics",joint,)
                damping = 0.0
                friction = 0.0
                if dynamics !== nothing                    
                    if haskey(dynamics,"damping")
                        damping = parse(T,dynamics["damping"])
                    end                  
                    if haskey(dynamics,"friction")
                        friction = parse(T,dynamics["friction"])
                    end
                end
                dynamics_parsed = @eponymtuple(damping,friction)
                mobi = RB.URDF.Mobility(
                    axis_parent,
                    nothing,
                    dynamics_parsed,
                    limit_parsed
                )
            else
                mobi = nothing
            end
            # @show jointname,jointtype
            RB.URDF.Joint(
                jointname,
                jointtype,
                origin,
                parent,
                child,
                mobi
            )
        end

        for joint in xml_joints
    ]
end

function link_to_rigidbody(i,link::RB.URDF.Link, 
        ro = SVector(0.0,0.0,0.0),
        R = RotX(0.0),
        cT = RB.QBF.QC;
        constrained=false,
        ci = Int[],
        Φi = collect(1:6)
    )
    (;inertial,visual) = link
    if inertial isa Nothing
        mass = 0.0
        inertia = @SMatrix zeros(3,3)
        mass_locus = @SVector zeros(3)
        R̄ = RotX(0.0)
    else
        (;mass,inertia) = inertial
        mass_locus = inertial.origin.xyz
        R̄ = RotXYZ(inertial.origin.rpy...)
    end
    ṙo = zero(ro)
    Ω = zero(ro)
    ω = R*Ω
    # μ = 0.5,
    # e = 0.9,
    movable = true
    
    loci = [zero(ro)]
    prop = RB.RigidBodyProperty(
        i,
        movable,
        mass,inertia,
        mass_locus,loci;
        constrained,
        type=Symbol(link.name)
    )
    if cT == RB.QBF.QC
        nmcs = RB.QBF.QC(mass,inertia)
    else
        ri = ro
        nmcs, _ = RB.NCF.NC1P3V(ri,ro,R,ṙo,ω)
    end
    state = RB.RigidBodyState(prop,nmcs,ro,R,ṙo,ω,ci,Φi)
    
    rigidmesh = nothing
    if visual !== nothing
        (;geometry,material) = visual
        color_parsed = parse(Makie.RGBAf,:slategrey) 
        if material !== nothing
            if material.color !== nothing
                color_parsed = material.color
            end
        end
        if geometry !== nothing
            (;filename,scale) = geometry
            mesh_loaded = load(filename)
            scaled_position = map(mesh_loaded.position) do p
                p .* scale |> 
                (Translation(visual.origin.xyz...) ∘ LinearMap(RotXYZ(visual.origin.rpy...)))
            end
            colors = fill(color_parsed,length(scaled_position))
            uv = GB.texturecoordinates(mesh_loaded)
            if uv isa Nothing
                rigidmesh = GB.Mesh(
                    GB.meta(scaled_position;
                    normals=GB.normals(mesh_loaded),
                    color=colors
                    ), 
                    GB.faces(mesh_loaded)
                )
            else
                rigidmesh = GB.Mesh(
                    GB.meta(scaled_position;
                    normals=GB.normals(mesh_loaded),
                    uv,
                    color=colors
                    ), 
                    GB.faces(mesh_loaded)
                )
            end
        end
    end
    # base_uv_texture = load("anymal_b_simple_description/meshes/base_uv_texture.jpg")
    # base_uv_texture_sampled = GM.Sampler(base_uv_texture; x_repeat=:repeat,y_repeat=:repeat)
    # mesh(rigidmesh,color=base_uv_texture_sampled)
    
    @show link.name, ro, R
    RB.RigidBody(prop,state,rigidmesh)
end

function find_base_link(joints)
    parents = [joint.parent for joint in joints] |> unique
    childs = [joint.child for joint in joints] |> unique
    base_candidates = [
        parent 
        for parent in parents
        if !(parent in childs)
    ]
    # @show base_candidates
    @assert length(base_candidates) == 1
    base_candidates[1]
end

function parse_joints_links!(
        joints_input::AbstractVector{<:RB.URDF.Joint{T}},
        linksdict
    ) where T
    joints = deepcopy(joints_input)
    linknames_all = keys(linksdict) |> collect
    o = zero(T)
    id_trans = Translation(o,o,o) ∘ LinearMap(one(RotMatrix{3, T}))
    base_link = find_base_link(joints)
    i = 1
    rb_base = link_to_rigidbody(
        i,
        linksdict[base_link],
        id_trans.translation,
        id_trans.linear,
        1;        
        constrained=true,
        ci = collect(1:12),
        Φi = Int[]
    )
    jointid = 0
    transformsdict = Dict([base_link => id_trans])
    bodiesdict = Dict{String,Any}([base_link => rb_base])
    jointcsts = []
    while !isempty(joints)
        # @show keys(transformsdict)
        for (ijoint,joint) in enumerate(joints)
            (;name,type,parent,child,origin) = joint
            if !(haskey(transformsdict,parent))
                "Skip"
            else # has parent  
                # has no this child
                if !(haskey(transformsdict,child))
                    (;rpy,xyz) = joint.origin
                    # @show rpy,xyz
                    parentlink = linksdict[parent]
                    childlink = linksdict[child]
                    if childlink.inertial isa Nothing
                        @warn "Removing a zero mass leaf"
                        delete!(linksdict,child)
                        deleteat!(joints,ijoint)
                    else
                        if (parentlink.inertial === nothing) && (type == "fixed")
                            #@warn "$(parentlink.name) has zero mass"
                            @warn "Replacing $parent with $child"

                            parentlink.inertial = childlink.inertial
                            parentlink.visual = childlink.visual
                            foreach(joints) do otherjoint
                                if otherjoint.parent == child
                                    otherjoint.parent = parent
                                end
                            end
                            delete!(linksdict,child)
                            #@warn "$name to $(childlink.name) is fixed"
                            # replace link
                        else
                            i += 1
                            jointtrans = Translation(xyz...) ∘ LinearMap(RotXYZ(rpy...))
                            tf_child = transformsdict[parent] ∘ jointtrans
                            push!(transformsdict, child => tf_child)                            
                            rb_child = link_to_rigidbody(
                                i,
                                linksdict[child],
                                tf_child.translation,
                                tf_child.linear,
                                1
                            )
                            push!(bodiesdict, child => rb_child)
                            jointid += 1
                            if joint.type == "fixed"
                                @warn "fixed joint: $(joint.name)"
                                push!(bodiesdict[parent].prop.loci,joint.origin.xyz)
                                (;C,c) = bodiesdict[parent].state.cache.funcs
                                push!(bodiesdict[parent].state.cache.Cps,C(c(joint.origin.xyz)))
                                fixedaxis = [1.0,0,0]
                                push!(bodiesdict[parent].prop.axes,fixedaxis)
                                push!(bodiesdict[parent].state.as,bodiesdict[parent].state.R*fixedaxis)
                                fixedjoint = RB.FixedJoint(
                                    jointid,
                                    RB.End2End(
                                        ijoint,
                                        RB.ID(bodiesdict[parent],
                                        length(bodiesdict[parent].prop.loci),
                                        length(bodiesdict[parent].prop.axes)),
                                        RB.ID(bodiesdict[child],1),
                                    )
                                )
                                push!(jointcsts,fixedjoint)
                            elseif joint.type == "revolute"
                                push!(bodiesdict[parent].prop.loci,joint.origin.xyz)
                                (;C,c) = bodiesdict[parent].state.cache.funcs
                                push!(bodiesdict[parent].state.cache.Cps,C(c(joint.origin.xyz)))
                                push!(bodiesdict[parent].prop.axes,joint.mobility.axis)
                                push!(bodiesdict[parent].state.as,bodiesdict[parent].state.R*joint.mobility.axis)
                                revjoint = RB.RevoluteJoint(
                                    jointid,
                                    RB.End2End(
                                        ijoint,
                                        RB.ID(bodiesdict[parent],
                                        length(bodiesdict[parent].prop.loci),
                                        length(bodiesdict[parent].prop.axes)),
                                        RB.ID(bodiesdict[child],1),
                                    )
                                )
                                push!(jointcsts,revjoint)
                            elseif joint.type == "prismatic"
                                push!(bodiesdict[parent].prop.loci,joint.origin.xyz)
                                (;C,c) = bodiesdict[parent].state.cache.funcs
                                push!(bodiesdict[parent].state.cache.Cps,C(c(joint.origin.xyz)))
                                push!(bodiesdict[parent].prop.axes,joint.mobility.axis)
                                push!(bodiesdict[parent].state.as,bodiesdict[parent].state.R*joint.mobility.axis)
                                prmjoint = RB.PrismaticJoint(
                                    jointid,
                                    RB.End2End(
                                        ijoint,
                                        RB.ID(bodiesdict[parent],
                                        length(bodiesdict[parent].prop.loci),
                                        length(bodiesdict[parent].prop.axes)),
                                        RB.ID(bodiesdict[child],1),
                                    )
                                )
                                push!(jointcsts,prmjoint)
                            end
                        end
                        deleteat!(joints,ijoint)
                    end
                    break
                end
            end
        end
    end
    transformsdict,bodiesdict,jointcsts
end

function urdf_to_rigidrobot(urdf,T=Float64)
    linksdict = parse_links(T,urdf)
    joints = parse_joints(T,urdf)
    transformsdict,bodiesdict,jointcsts = parse_joints_links!(joints,linksdict)
    # @show collect(values(bodiesdict))
    rbs = TypeSortedCollection(collect(values(bodiesdict)))
	numbered = RB.number(rbs)
	indexed = RB.index(rbs,)
    ss = Int[]
	tensiles = (cables = ss,)
	connected = RB.connect(rbs,zeros(Int,0,0))
	tensioned = @eponymtuple(connected,)
    jointed = RB.join(TypeSortedCollection(jointcsts),indexed)
	cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
	# contacts = [RB.Contact(i,μ,e) for i = [5]]
	st = RB.Structure(rbs,tensiles,cnt,)
    bot = RB.Robot(st)
end
#-- functions

# urdf = readxml(raw"spot_ros\spot_description\urdf\spot.urdf.xacro")

cd(raw"examples/URDF")

urdf = readxml(raw"anymal_b_simple_description/urdf/anymal.urdf")

cd(raw"examples/URDF/unitree_ros/robots")
urdf = readxml(raw"a1_description/urdf/a1.urdf")

urdf = readxml(raw"example-robot-data/robots/double_pendulum_description/urdf/double_pendulum.urdf")

urdf = readxml(raw"example-robot-data/robots/anymal_b_simple_description/robots/anymal.urdf")

urdf = readxml(raw"example-robot-data/robots/panda_description/urdf/panda.urdf")

newbot = urdf_to_rigidrobot(urdf)

RB.ASSETS_DIR

starsimage = load(RB.assetpath("stars.jpg")) |> rotr90


q = RB.get_q(newbot.st)
Φ = RB.make_Φ(newbot)
findall((x)->abs(x)>1e-10,Φ(q))
Φ(q)
A = RB.make_A(newbot)

_, pidx = rref_with_pivots(transpose(A(q)))
ridx = collect(1:size(A(q),1))
deleteat!(ridx,pidx)
ridx
mass_matrices = RB.build_MassMatrices(newbot)
mass_matrices.M̌ |> size

newbot.st.nconstraints
λ = rand(newbot.st.nconstraints)
RB.∂Aᵀλ∂q̌(newbot.st,λ)
@descend_code_warntype RB.∂Aᵀλ∂q̌(newbot.st,λ)
@report_opt RB.∂Aᵀλ∂q̌(newbot.st,λ)
Makie.inline!(false)
plot_traj!(
    newbot;
    # AxisType=Axis3,
    showmesh=true,
    showpoints=true,
    showlabels=false,
    showground=false,    
    sup! = (ax,tgob,sgi)->begin
        # hidedecorations!(ax)
        nothing
    end
)

function dynfuncs(bot)
    (;st) = bot
    function F!(F,q,q̇,t)
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        # RB.update_tensiles!(st)
        RB.apply_gravity!(st)
        F .= RB.generate_forces!(st)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        # RB.update_tensiles!(st)
        # RB.build_tangent_stiffness_matrix!(∂F∂q̌,st)
        # RB.build_tangent_stiffness_matriẋ!(∂F∂q̌̇,st)
    end

	@eponymtuple(F!,Jac_F!)
end

dt = 1e-3
tspan = (0.0,1.1)
prob = RB.SimProblem(newbot,dynfuncs)

RB.solve!(prob,RB.Zhong06();tspan,dt,ftol=1e-10,maxiters=50,verbose=true,exception=true)
RB.solve!(prob,RB.Alpha(0.7);tspan,dt,ftol=1e-10,verbose=true,exception=true)

newbot.traj
with_theme(theme_pub;
    figure_padding = (0,0,0,0),) do 
    plot_traj!(
        newbot;
        figsize=:FHD,
        showpoints=false,
        showlabels=false,
        showground=false,
        doslide=true
    )
end
ME = RB.mechanical_energy!(newbot;gravity=true)
ME.V |> lines
ME.T |> lines
ME.E |> lines
ME.E
jie = load("末节.STL")
mesh(jie;color=:white)
scatter!([Point3f(0.0,0.0,0.0)])