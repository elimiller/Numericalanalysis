### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 918fde30-9882-11ef-33e6-1358e287ed87
md"""
# MATH566 - Lesson 20

## Adaptive quadrature: recursive quadrature
"""

# ╔═╡ 2a109ea7-d3b6-46ad-9fa0-4d91f75a12a9
md"""
### Quadrature rules
"""

# ╔═╡ c38e9521-de2d-4170-870a-34ca0daff86e
function trapezoid(a,b,f)
	return 0.5*(b-a)*(f(a)+f(b))
end

# ╔═╡ 65a2824c-0004-4959-8b37-de59e9e6e2e1
function Simpson(a,b,f)
	c = a + 0.5*(b-a)
	return (b-a)*(f(a)+4*f(c)+f(b))/6
end

# ╔═╡ 9dd230fd-f086-49f1-8629-2c0e0daf5966
md"""
## Recursive quadrature
"""

# ╔═╡ 9cd46b80-95ce-4363-8f9d-6392156e076d
md"""
Approximate

$\int_a^b f(t) dt  = \sum_i w_if_i$

"""

# ╔═╡ c876be8b-465c-4331-a163-0ef3edf87963
function RecQ(a,b,f,quadR,relerr,lvl)
	MaxLevel = 6
	Qab = quadR(a,b,f)
	c = a + 0.5*(b-a)
	Qac = quadR(a,c,f); Qcb = quadR(c,b,f)
	Qacb = Qac + Qab
	if (abs(Qab-Qacb) <= relerr*abs(Qacb) || (lvl >= MaxLevel))
		return Qacb
	else
		return RecQ(a,c,f,quadR,relerr,lvl+1) + RecQ(c,b,f,quadR,relerr,lvl+1)
	end
end

# ╔═╡ 60e98df4-325e-4ad6-8fad-bf673f6540f2
md"""
Test on first few members of the monomial basis.
"""

# ╔═╡ 5322ae12-8c7d-401f-a40e-ce1bd0de0c5d
f0(t) = 1; f1(t) = t; f2(t) = t^2

# ╔═╡ f2106196-1b88-40f5-b58a-18e73e27426a
RecQ(0,1,f0,trapezoid,0.001,0)

# ╔═╡ 27a6d621-6dea-4650-bdd3-be9e90d03314
md"""
Extra Credit opportunity. Fix the quadrature above.  Plot the function and lines corresponding to the evaluaiton intervals
"""

# ╔═╡ Cell order:
# ╠═918fde30-9882-11ef-33e6-1358e287ed87
# ╟─2a109ea7-d3b6-46ad-9fa0-4d91f75a12a9
# ╠═c38e9521-de2d-4170-870a-34ca0daff86e
# ╠═65a2824c-0004-4959-8b37-de59e9e6e2e1
# ╟─9dd230fd-f086-49f1-8629-2c0e0daf5966
# ╠═9cd46b80-95ce-4363-8f9d-6392156e076d
# ╠═c876be8b-465c-4331-a163-0ef3edf87963
# ╟─60e98df4-325e-4ad6-8fad-bf673f6540f2
# ╠═5322ae12-8c7d-401f-a40e-ce1bd0de0c5d
# ╠═f2106196-1b88-40f5-b58a-18e73e27426a
# ╟─27a6d621-6dea-4650-bdd3-be9e90d03314
