function sample_ball_Muller(d,σ=1.0)
    i = one(σ)
    u = rand(Uniform(-i,i),d)
    r = rand(typeof(σ))^(i/d)
    σ.*r./norm(u).*u
end

function sample_chart(chart)
    @unpack nX, σ, Uc, xc, neighbor_vectors = chart
    local y_rand
    while true
        y_rand = sample_ball_Muller(nX,σ)
        inside = true
        for y_nb in neighbor_vectors
            hp = transpose(y_rand.-y_nb./2)*y_nb
            if hp > 0
                inside = false
                break
            end
        end
        if inside
            break
        end
    end
    x_rand = xc + Uc*y_rand
end


function sample_atlas(atlas)
    nchart = length(atlas)
    ic = rand(DiscreteUniform(1,nchart))
    chart = atlas[ic]
    x_rand = sample_chart(chart)
end

function find_nearest_state(states,x)
    tree = KDTree(states,reorder=false)
    i, dis = nn(tree,x)
    i
end


function sample_action(tr,mag=0.1)
    @unpack hub = tr
    @unpack actuators = hub
    nact = length(actuators)
    u_dis = Uniform(-mag,mag)
    u = rand(u_dis,nact)
    function random_action(t)
        u
    end
end
# draw(PNG("wheel10.png", 16cm, 16cm), gplot(mg,nodelabel=1:nv(mg),edgelabel=1:ne(mg)))

# y0 = chart0.ϕc(chart0.xc)
# y0 += 0.1rand(length(y0))
# # @code_warntype chart0.ϕc(chart0.xc)
# ψc!,Jac_ψc! = chart0.make_ψc(y0)
# nlsolve(ψc!,Jac_ψc!,chart0.xc;method=:newton,ftol=1e-14,iterations=50)
# fig,ax,plotscatter = scatter([Point2(sample_chart(chart0)) for i = 1:4000])
# ax.aspect = DataAspect()
# fig
