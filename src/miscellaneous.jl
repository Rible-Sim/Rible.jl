
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

function (ar::AdamsResults)(name::String)
    index = findfirst((x)->x==name,ar.names)
    [step[index] for step in ar.steps]
end
