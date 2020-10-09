### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 48217560-f3f9-11ea-27c7-ad6b4997aa59
using PlutoUI

# ╔═╡ 0d64b9e0-f3fa-11ea-2472-2bc283193ca5
md"""
### A Tensegrity Manipulator

$k = 1e4 g/s^2 = 10kg\cdot m/s^2/m = 10N/m$
"""

# ╔═╡ 7dd09690-f404-11ea-13d2-3b3c5a241fde
md"""$(LocalResource("man_freefall_energy.png",:width => 500))"""

# ╔═╡ f731ddd0-f3fc-11ea-19f6-3f86fa5f7a15
md"""
Without cables.

Free vibration in the gravitational field. Essentially a double pendulum. 
"""

# ╔═╡ 3a0083e0-f3f9-11ea-0293-bf817a1f9fcf
md"""$(LocalResource("man_freevibra_energy.png",:width => 500))"""

# ╔═╡ a3684e70-f404-11ea-3e6d-cb13153c89d3
md"""
With Cables.

Free vibraion under the gravity and tension forces. Damping set to zero. 
"""

# ╔═╡ 1eebe730-f402-11ea-3e74-b5fed3ccf530
md"""
#### Linear Actuation
"""

# ╔═╡ 2fba6c80-f402-11ea-24ee-a38aef67194a
md"""
$(LocalResource("linearactuating_c=0.0.png", :width => 800))
"""

# ╔═╡ 06aa8100-f406-11ea-3fa4-090a9d9baf68
md"""
Zero damping.
"""

# ╔═╡ fa86e620-f405-11ea-1917-c10e4c4b6ca6
md"""
$(LocalResource("linearactuating_c=100000.0.png", :width => 800))
"""

# ╔═╡ 1a6020b0-f406-11ea-12f7-55ea2d351393
md"""Damping coefficient $c=1e5$"""

# ╔═╡ 27f2e580-f412-11ea-1b21-9d6fb21f5b9a
md"""
#### Inverse Kinematics
"""

# ╔═╡ e51b147e-f411-11ea-3b1d-dfb814fb8d8b
md"""
$(LocalResource("man_inverse.png", :width => 500))
"""

# ╔═╡ 21c93920-f4a8-11ea-04fe-7d124d71263d
md"""
#### Modal Analysis
"""

# ╔═╡ 3e95f702-f4a8-11ea-2c95-99c1f4048590
md"""
$(LocalResource("man_mode.png",:width => 800))
"""

# ╔═╡ Cell order:
# ╠═48217560-f3f9-11ea-27c7-ad6b4997aa59
# ╟─0d64b9e0-f3fa-11ea-2472-2bc283193ca5
# ╠═7dd09690-f404-11ea-13d2-3b3c5a241fde
# ╟─f731ddd0-f3fc-11ea-19f6-3f86fa5f7a15
# ╠═3a0083e0-f3f9-11ea-0293-bf817a1f9fcf
# ╟─a3684e70-f404-11ea-3e6d-cb13153c89d3
# ╟─1eebe730-f402-11ea-3e74-b5fed3ccf530
# ╟─2fba6c80-f402-11ea-24ee-a38aef67194a
# ╟─06aa8100-f406-11ea-3fa4-090a9d9baf68
# ╟─fa86e620-f405-11ea-1917-c10e4c4b6ca6
# ╟─1a6020b0-f406-11ea-12f7-55ea2d351393
# ╟─27f2e580-f412-11ea-1b21-9d6fb21f5b9a
# ╟─e51b147e-f411-11ea-3b1d-dfb814fb8d8b
# ╟─21c93920-f4a8-11ea-04fe-7d124d71263d
# ╟─3e95f702-f4a8-11ea-2c95-99c1f4048590
