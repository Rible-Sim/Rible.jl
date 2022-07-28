

# coζ = zeros(Float64, 8, 8)
# dζdq = zeros(Float64, 8, 78)
# b = zeros(Float64, 8, 6)
# T = zeros(Float64, 8, 8)
# ∂l∂q = zeros(Float64, 6, 78)
# q = zeros(Float64, 78)
# dζdq2 = zeros(Float64, 8, 78)
# XLSX.openxlsx("t1.xlsx", mode="r") do xf
#     sh = xf["Temp"]
#     s1 = xf["前向差分"]
#     sh2 = xf["t1"]
#     for j in 1:78
#         for i in 8+1:8+8
#             dζdq[i-8, j] = sh[i, j]
#         end
#         for i in 8+8+1
#             q[j] = sh[i, j]
#         end
#     end
#     for i in 1:8
#         for j in 1:8
#             coζ[i, j] = sh[i, j]
#         end
#     end
#     for i in 133+1:133+8
#         for j in 1:78
#             dζdq2[i-133, j] = s1[i, j]
#         end
#     end
#     for i in 1:8
#         for j in 1:6
#             b[i, j] = sh2[i, j]
#         end
#         for j in 1:8
#             T[i, j] = sh2[i+8, j]
#         end
#     end
#     for i in 1:6
#         for j in 1:78
#             ∂l∂q[i, j] = sh2[i+16, j]
#         end
#     end
# end
# display(dζdq2)
# data = 1/2 * coζ * dζdq
# display(data)
# display(1/2 * coζ * T*b*∂l∂q)

# t1 = deepcopy(∂l∂q)

# # t1[2,:] = ∂l∂q[4,:]
# display(1/2 * coζ * T*b*t1)
# display(coζ*T*b)
# display(dζdq2 * 2)
# a1 = (-0.0072904+0.10001*0.0492836)/-0.113402
# a2 = (-0.213266+0.10001*0.998785)/-0.113402
# a3 = (0.00587082-0.0999971*0.041649)/-0.113387
# a4 = (0.426718-0.0999971*1.99957)/-0.113387
# t1[2, 1] = a1; t1[2, 2] = a2
# t1[3, 5] = a3; t1[3, 6] = a4
# display(t1)
# display(1/2 * coζ * T*b*t1)
# display(t1)
# TR.Record_build_∂ζ∂q(tg, q, "t1.xlsx", "t1")
# id = 201
# ∂ζ∂q = TR.build_∂ζ∂q(bot.tg, bot.traj[id].q̌)
# # sᵏ = bot.traj[id].s
# # ζ = TR.build_ζ(bot.tg)
# # n = length(ζ)
# # κ₁ = 10; κ₂ = 10
# # coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sᵏ[i])^2)^(-1/2) for i in 1:n])
# # coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
# # data = 1/2 * coζ * ∂ζ∂q
#
# # t1 = zeros(Float64, 8, 78)
# # t1[1,1] = -0.00188; t1[1,2]
#
# # XLSX.openxlsx("t1.xlsx", mode="rw") do xf
# #     sh = xf["Temp"]
# #     length, width = size(data)
# #     for i in 1:length
# #         sh[i, 1:width] = data[i, :]
# #     end
# # end