struct Plane{T}
    n::SArray{Tuple{3},T,1,3}
    d::T
    r::SArray{Tuple{3},T,1,3}
end

function Plane(n::AbstractVector{T},r::AbstractVector{T}) where T
    n /= norm(n)
    a = n[1]
    b = n[2]
    c = n[3]
    d = -(a*r[1]+b*r[2]+c*r[3])
    Plane(SVector{3}(n),d,SVector{3}(r))
end

function Plane(n::AbstractVector{T},d::T) where T
    n /= norm(n)
    a = n[1]
    b = n[2]
    c = n[3]
    o = zero(d)
    x = o
    y = o
    z = d/(-c)
    Plane(SVector{3}(n),d,SVector{3}(x,y,z))
end

function Plane(a::T,b::T,c::T,d::T) where T
    n = [a,b,c]
    n /= norm(n)
    Plane(n,d)
end

function signed_distance(x::AbstractVector{T},p::Plane) where T
    @unpack n, d = p
    transpose(n)*x + d
end

function distance(x::AbstractVector{T},p::Plane) where T
    abs(signed_distance(x,p))
end

function ison(x::AbstractVector{T},p::Plane;tol=eps(T)) where T
    distance(x,p) < tol
end

function show_camera(scene)
    camera = cameracontrols(scene)
    lookat = camera.lookat[]
    eyepos = camera.eyeposition[]
    @show eyepos, lookat
end

struct AdamsResults
    names::Vector{String}
    steps::Vector{Vector{Float64}}
end

function AdamsResults(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    StepMap = findfirst("//adams:Analysis/adams:StepMap", adams_xml.root, ns)
    variable_names = [ele["name"] for ele in eachelement(StepMap)]
    Data =  findall("//adams:Analysis/adams:Data", adams_xml.root, ns)
    dynamics_Data = Data[findfirst((x)-> occursin("dynamic",x["name"]), Data)]
    steps = Vector{Vector{Float64}}()
    for step_node in eachelement(dynamics_Data)
        step_data = parse.(Float64,split(step_node.content,"\n")[2:end-1])
        @assert length(variable_names) == length(step_data)
        push!(steps,step_data)
    end
    AdamsResults(variable_names, steps)
end

function (ar::AdamsResults)(name::String,i=1)
    index = findall((x)->x==name,ar.names)
    [step[index[i]] for step in ar.steps]
end
