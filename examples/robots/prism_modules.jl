function prism_modules(;
    κ = nothing,
    b = 1.010,
    r = 0.5463,
    m = 6,
    α = 2π/m,
    h = 0.35,
    n = 1,
    θ = 2π/3,
    p = 3,
)
d, l = divrem(n,2)
bps = [
    begin
        if j |> iseven
            ps = [
                [r*cos(i*α),r*sin(i*α),j*h]
                for i = 1:m
            ]
        else
            ps = [
                [r*cos(i*α+θ),r*sin(i*α+θ),j*h]
                for i = 1:m
            ]
        end
        ps .|> SVector{3} |> CircularArray
    end
    for j = 0:n
]

ys = [0,-√3*r/2,√3*r/2]
xs = [0,-3r/2,-3r/2]


# bar
A_bar = ((5e-3)^2-(4e-3)^2)*π
ρ_bar = 1800.0
b_bar = 1.01
m_bar = A_bar*b_bar*ρ_bar
bars = [
    begin
        ks = m*n*(k-1)
        js = m*(j-1)
        if j |> iseven
            [
                make_3d_bar(
                    ks+js+i,
                    bps[j  ][i] + SVector(xs[k],ys[k],0),
                    bps[j+1][i] + SVector(xs[k],ys[k],0); 
                    #ci = ifelse(j==1,[1,2,3],Int[]),
                    m = m_bar,
                    radius_ratio = 1/60,
                ) for i = 1:m
            ]
        else
            [
                make_3d_bar(
                    ks+js+i,
                    bps[j  ][i-1] + SVector(xs[k],ys[k],0),
                    bps[j+1][i-1] + SVector(xs[k],ys[k],0);
                    m = m_bar,
                    radius_ratio = 1/60,
                ) for i = 1:m
            ]
        end
    end
    for j = 1:n, k = 1:p
]

nbars = length.(bars) |> sum
nb = nbars


rbs = vcat(
    reduce(vcat,bars),
) |> TypeSortedCollection
numbered = RB.number(rbs)
cm = CircularArray(collect(1:m))
# sharing_elas = ElasticArray{Int}(undef, nb, 0)
# for j = 1:n
#     is = 2m*(j-1)
#     for i = 1:m
#         for k = 1:p
#             row = zeros(Int,nb)
#             row[is+cm[i]  ]   = k+3
#             row[is+m+cm[i+1]] = k
#             append!(sharing_elas,row)		
#         end
#     end
# end	
# for j = 2:n
#     is = 2m*(j-2)+m
#     for i = 1:m
#         for k = 1:p
#             row = zeros(Int,nb)
#             row[is+cm[i]  ]   = k+3
#             row[is+m+cm[i+2]] = k
#             append!(sharing_elas,row)		
#         end
#     end
# end
# sharing = Matrix(sharing_elas')
# display(sharing)
# indexed = RB.index(rbs,sharing)
indexed = RB.index(rbs,)
    
@myshow nb

mat_cable = filter(
    row->row.name == "Kevlar",
    material_properties
)[1]
radius = 0.32e-3Unitful.m
cable_hor_length = 0.5463Unitful.m
cable_dia_length = 0.6488017339680894Unitful.m
κ_hor = ustrip(Unitful.N/Unitful.m,(mat_cable.modulus_elas)*π*(radius)^2/cable_hor_length)
κ_dia = ustrip(Unitful.N/Unitful.m,(mat_cable.modulus_elas)*π*(radius)^2/cable_dia_length)
@myshow κ_hor, κ_dia
cable_hor_restlength = ustrip(Unitful.m,cable_hor_length)-17/κ_hor
cable_dia_restlength = ustrip(Unitful.m,cable_dia_length)-17*1.1876290206261921/κ_dia
# @show uconvert(Unitful.N/Unitful.m,κ),ustrip(Unitful.N/Unitful.m,κ)
κs = typeof(κ_hor)[]
friction_coefficients = typeof(cable_hor_restlength)[]
connecting_elas = ElasticArray{Int}(undef, nb, 0)
for k = 1:p
    # # lower
    ks = m*n*(k-1)
    idx = collect(1:m)
    # if k == 2 
    #     deleteat!(idx,[1,])
    # elseif k == 3 
    #     deleteat!(idx,[5,6])
    # end
    for i = idx
        row = zeros(Int,nb)
        row[ks+cm[i  ]] =  1
        row[ks+cm[i+1]] = -1
        append!(connecting_elas,row)
        push!(κs,κ_hor)
        push!(friction_coefficients,cable_hor_restlength)
    end

    for j = 0:n-1
        is = m*j

        # cross
        if j |> iseven
            for i = idx
                row = zeros(Int,nb)
                row[ks+is+cm[i  ]] =  1
                row[ks+is+cm[i-1]] = -2
                append!(connecting_elas,row)
                push!(κs,κ_dia)
                push!(friction_coefficients,cable_dia_restlength)
            end
        else
            for i = idx
                row = zeros(Int,nb)
                row[ks+is+cm[i+1]] = -2
                row[ks+is+cm[i  ]] =  1
                append!(connecting_elas,row)
                push!(κs,κ_dia)
                push!(friction_coefficients,cable_dia_restlength)
            end
        end

        # upper 
        for i = idx
            row = zeros(Int,nb)
            row[ks+is+cm[i-2]] =  2
            row[ks+is+cm[i-1]] = -2
            append!(connecting_elas,row)
            push!(κs,κ_hor)
            push!(friction_coefficients,cable_hor_restlength)
        end
    end
end
ncables_prism = size(connecting_elas,2)
connecting = Matrix(connecting_elas')
# display(connecting)
connected = RB.connect(rbs,connecting)

ncables = size(connecting,1)
# @assert ncables == ncables_prism + ncables_outer

cables_prism = [
    RB.DistanceSpringDamper3D(
    i,
    friction_coefficients[i], #restlength  
    κs[i],
    0.0; #damping
    slack=true) 
    for i = 1:ncables_prism
]
cables = cables_prism
acs = [
    RB.RegisterActuator(
        1,
        collect(1:ncables_prism),
        zeros(ncables_prism)
    ),
]
apparatuses = (cables = cables,)
hub = (actuators = acs,)


# csts_bar2bar = [
#     RB.PinJoint(
#         i+m*(j-1),
#         RB.Hen2Egg(
#             i+m*(j-1),
#             RB.Signifier(bars[j][cm[i+1]],1,1),
#             RB.Signifier(bars[j][cm[i]]  ,2,1)
#         )
#     )
#     for j = 1:d for i = 1:m
# ]
if p == 3
    csts_bar2bar = [
        # one to two
        RB.PinJoint(1,RB.Hen2Egg(
                RB.Signifier(bars[    1][cm[4]],1),
                RB.Signifier(bars[n+1][cm[8]],1)
        )),
        RB.PinJoint(2,RB.Hen2Egg(
                RB.Signifier(bars[    1][cm[5]],1),
                RB.Signifier(bars[n+1][cm[7]],1)
        )),
        RB.PinJoint(3,RB.Hen2Egg(
                RB.Signifier(bars[    1][cm[3]],2),
                RB.Signifier(bars[n+1][cm[11]],2)
        )),
        RB.PinJoint(4,RB.Hen2Egg(
                RB.Signifier(bars[    1][cm[2]],2),
                RB.Signifier(bars[n+1][cm[12]],2)
        )),
        # one to three
        RB.PinJoint(5,RB.Hen2Egg(
                RB.Signifier(bars[     1][cm[ 3]],1),
                RB.Signifier(bars[2n+1][cm[13]],1)
        )),
        RB.PinJoint(6,RB.Hen2Egg(
                RB.Signifier(bars[     1][cm[ 4]],1),
                RB.Signifier(bars[2n+1][cm[18]],1)
        )),
        RB.PinJoint(7,RB.Hen2Egg(
                RB.Signifier(bars[     1][cm[ 1]],2),
                RB.Signifier(bars[2n+1][cm[17]],2)
        )),
        RB.PinJoint(8,RB.Hen2Egg(
                RB.Signifier(bars[     1][cm[ 2]],2),
                RB.Signifier(bars[2n+1][cm[16]],2)
        )),
        # two to three
        RB.PinJoint(9,RB.Hen2Egg(
                RB.Signifier(bars[ n+1][cm[ 9]],1),
                RB.Signifier(bars[2n+1][cm[17]],1)
        )),
        RB.PinJoint(10,RB.Hen2Egg(10,
                RB.Signifier(bars[ n+1][cm[ 7]],2),
                RB.Signifier(bars[2n+1][cm[15]],2)
        )),
    ]
    csts = csts_bar2bar

    jointedmembers = RB.join(csts,indexed)
else
    jointedmembers = RB.unjoin()
end

cnt = RB.Connectivity(numbered,indexed,@eponymtuple(connected,),jointedmembers)

st = RB.Structure(rbs,apparatuses,cnt)
bot = RB.Robot(st,hub)
end