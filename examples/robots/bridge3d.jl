function bridge3d(;
    k = nothing,
    n = 2,		
    d = 10.0,
    hw = 9.0,
    h = 4.5,
    b = 6.0,
    hww = 3.0,
    l =  (hw + hww )^2+(h+b)^2 |> sqrt
)

hw = sqrt(l^2 - (h+b)^2 ) - hww

up = [
    [ hw*i,6+3*i+12*j, h+b]
    for i in [-1,1], j = 0:n-1
] .|> SVector{3}
reverse!(@view up[2,:])
lo = [
    [-hww*i,6+3*i+12*j,0]
    for i in [-1,1], j = 0:n-1
] .|> SVector{3}
reverse!(@view lo[2,:])
mid = begin
    ret = [
        [3.0*i,6*j,b]
        for i in [-1,1], j = 0:2n
    ]
    reverse!(@view ret[2,:])
    # @show ret
    ret
end |> transpose |> vec .|> SVector{3}
# display(mid)
lbars = [
    make_3d_bar(
    i,
    lo[1,i],
    up[1,i]; ci = [1,2,3]) 
    for i = 1:n
]
rbars = [
    make_3d_bar(
    n+i,
    lo[2,i],
    up[2,i]; ci = [1,2,3]) 
    for i = 1:n
]
deckcenter = SVector(0.0,n/2*12,b)
nodes_deck = vcat(mid,[SVector(0.0,0.0,b)]).-Ref(deckcenter)
hei = 0.25
box = Meshes.Box(
    Meshes.Point(nodes_deck[begin]), 
    Meshes.Point(nodes_deck[2n+2]+SVector(0,0,hei))
)
deck = new_deck(
    2n+1,nodes_deck, # loci
    deckcenter,#ro
    RotX(0.0),#R
    deckcenter,#ri
    box;
)
# ;k=500.0,k1=1000.0,c=0.0,
# 			d = 0.1*√2/2, r2 = 0.11,
# 			ijkl=1)

rbs = TypeSortedCollection(vcat(lbars,rbars,[deck]))
numberedpoints = RB.number(rbs)
indexedcoords = RB.index(rbs)
# # #
cnt_matrix_elas = ElasticArray{Int}(undef, 2n+1, 0)
# left
for i = 1:n
    if i == 1
        js, je = 2i-1, 2i+2
    elseif i == n
        js, je = 2i-3, 2i+1
    else
        js, je = 2i-3, 2i+2
    end
    for j = js:je
        row = zeros(Int,2n+1)
        row[i] = 2
        row[2n+1] = -j
        append!(cnt_matrix_elas,row)
    end
    for j = js:je
        row = zeros(Int,2n+1)
        row[i] = 1
        row[2n+1] = -j
        append!(cnt_matrix_elas,row)
    end
    for j = js:je
        row = zeros(Int,2n+1)
        row[i] = 1
        row[2n+1] =-( 4n+3-j)
        append!(cnt_matrix_elas,row)
    end
    for j = 1:n
        row = zeros(Int,2n+1)
        row[i] = 1
        row[n+j] = -2
        append!(cnt_matrix_elas,row)
    end
    # if i != 1
    # 	row = zeros(Int,2n+1)
    # 	row[[i,i+n-1]] .= [2,-1]
    # 	append!(cnt_matrix_elas,row)
    # end
    # if i != n
    # 	row = zeros(Int,2n+1)
    # 	row[[i,i+1]] .= [2,-2]
    # 	append!(cnt_matrix_elas,row)
    # end
end
# Right
for i = n+1:2n
    if i == n+1
        js, je = 2i, 2i+3
    elseif i == 2n
        js, je = 2i-2, 2i+2
    else
        js, je = 2i-1, 2i+3
    end
    for j = js:je
        row = zeros(Int,2n+1)
        row[i] = 2
        row[2n+1] = -j
        append!(cnt_matrix_elas,row)
    end
    
    for j = js:je
        row = zeros(Int,2n+1)
        row[i] = 1
        row[2n+1] = -j
        append!(cnt_matrix_elas,row)
    end
    for j = js:je
        row = zeros(Int,2n+1)
        row[i] = 1
        row[2n+1] = -( 4n+3 - j)
        append!(cnt_matrix_elas,row)
    end
    for j = 1:n
        row = zeros(Int,2n+1)
        row[i] = 1
        row[j] = -2
        append!(cnt_matrix_elas,row)
    end
# 	if i != 1
# 		row = zeros(Int,2n+1)
# 		row[[i,i+n-1]] .= [2,-1]
# 		append!(cnt_matrix_elas,row)
# 	end
    # if i != 2n
    # 	row = zeros(Int,2n+1)
    # 	row[[i,i+1]] .= [2,-2]
    # 	append!(cnt_matrix_elas,row)
    # end
end
# cross
for i = 1:n
    if  i == 1
        js, je = 2n+1-i, 2n+1-i
    else
        js, je = 2n+1-i, 2n+1-i+1
    end
    for j = js:je		
        row = zeros(Int,2n+1)
        row[i] = 2
        row[j] = -2
        append!(cnt_matrix_elas,row)
    end
end
# display(cnt_matrix_elas)
cnt_matrix = Matrix(transpose(cnt_matrix_elas))
ncables = size(cnt_matrix,1)
hncables = div(ncables,2)
mat_cable = filter(
    row->row.name == "Nylon 66",
    material_properties
)[1]
diameter = 1e-2
κ = (mat_cable.modulus_elas |> ustrip)*1e9*π*(diameter/2)^2/10
@show κ
if k isa Nothing
    cables = [RB.DistanceSpringDamper3D(0.0,κ,0.0;slack=false) for i = 1:ncables]
else
    cables = [RB.DistanceSpringDamper3D(0.0,k[i],0.0;slack=false) for i = 1:ncables]
end
acs = [
    RB.ManualActuator(i,[i,i+hncables],zeros(2))
    for i = 1:hncables
]
apparatuses = (cables = cables,)
hub = (actuators = acs,)
# 	# triplex 1
# 	1 0 0 -1  0  0  0 0 0 0;
# 	0 1 0  0 -1  0  0 0 0 0;
# 	0 0 1  0  0 -1  0 0 0 0;
# 	# triplex intersect
# 	0 0 0  1 -1  0 0 0 0 0;
# 	0 0 0  0  1 -1 0 0 0 0;
# 	0 0 0 -1  0  1 0 0 0 0;
# 	# triplex 2
# 	0 0 0  1  0  0 -2 0 0 0;
# 	0 0 0  0  1  0 -3 0 0 0;
# 	0 0 0  0  0  1 -4 0 0 0;
# 	# Inner
# 	0 0 0 0 0 0 1 -2 0 0;
# 	0 0 0 0 0 0 1 -3 0 0;
# 	0 0 0 0 0 0 1 -4 0 0;
# 	# Outer
# 	0 0 0 0 0 0 2 -2 0 0;
# 	0 0 0 0 0 0 3 -3 0 0;
# 	0 0 0 0 0 0 4 -4 0 0;
# 	#
# 	0 0 0 0 0 0 0 2 -2 0;
# 	0 0 0 0 0 0 0 3 -4 0;
# 	0 0 0 0 0 0 0 4 -3 0;
# 	# Outer
# 	0 0 0 0 0 0 0 0 2 -2;
# 	0 0 0 0 0 0 0 0 3 -3;
# 	0 0 0 0 0 0 0 0 4 -4;
# 	# Inner
# 	0 0 0 0 0 0 0 0 2 -1;
# 	0 0 0 0 0 0 0 0 3 -1;
# 	0 0 0 0 0 0 0 0 4 -1;
# 	]
connected = RB.connect(rbs,cnt_matrix)
# ss = Int[]
# apparatuses = (cables = ss,)
# connected = RB.connect(rbs,zeros(Int,0,0))
# #
#
# cst1 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb1_to_3[1],2),RB.ID(rb4,1)))
# cst2 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb1_to_3[2],2),RB.ID(rb4,2)))
# cst3 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb1_to_3[3],2),RB.ID(rb4,3)))
# jointedmembers = RB.join((cst1,cst2,cst3),indexedcoords)
#
cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

st = RB.Structure(rbs,apparatuses,cnt)
bot = RB.Robot(st,hub)
end