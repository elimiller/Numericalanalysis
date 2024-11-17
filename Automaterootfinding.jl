### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 3849b625-12b8-4b0c-b893-10bea5bec05f
using LinearAlgebra

# ╔═╡ 0cc786f0-a50f-11ef-1124-fbdbbb7c69b3
function polynomialfoil(p1,p2)
    # p1 and p2 are arrays of coefficients, where the i-th element is the coefficient for t^i
    result_degree = length(p1) + length(p2) - 2
    result = zeros(Float64, result_degree + 1)
    
    # Multiply each term of p1 by each term of p2 and accumulate the results
    for i in 1:length(p1)
        for j in 1:length(p2)
            result[i + j - 1] += p1[i] * p2[j]
        end
    end
    
    return result
end

# ╔═╡ 624f4b45-c518-47fa-9f4e-fcc2155dac1d
function doublefactorial(n)
	list = [i for i = 1:n]
	product = 1
	for j = 0:2:n-1
		product *= list[n-j]
	end
	return product
end

# ╔═╡ ccb7b862-df6c-43fb-a42e-e8308509019e
function innerproduct(f,g)
	h = polynomialfoil(f,g)
	sum = 0 ; a = sqrt(pi)
	 n = length(h)
	for i = 0:2:n-1 #Degree is n-1, so for odd it's zero. Even are i-1
		sum += h[i+1]*doublefactorial(i-1)/2^(i/2)
	end
	return (a/sqrt(4*pi))*sum
end
	

# ╔═╡ 8ecd0fac-c9db-4996-bb38-2d12689e5da4
function Basis(n)
	z = []
	for i = 1:n+1 #ChatGPT's code
        row = [j == i ? 1 : 0 for j in 1:n+1] 
        push!(z, row)
    end
	for i = 1:(n+1)
		v = z[i]
		sum = zeros(n+1)
		if  i != 1
			for j = 1:(i-1)
				sum += innerproduct(v,z[j]).*z[j] 
			end
		end
		x = v .- sum 
		h = innerproduct(x,x)
		J = x/sqrt(h)
		z[i] = J
	end
	return z
end

# ╔═╡ 26b34eb1-6e76-4b61-98ae-f94508871810
B = Basis(11)

# ╔═╡ b9b88780-34e8-44bd-962c-6c27b7a5bdc7
B[10]

# ╔═╡ 1b63c404-15de-42e0-9e11-23f7e8845fa0
function Newton1(x,x1,f,fp)
	fpx = fp(x)
	if abs(fpx) > 100 *  eps(Float64)
		return x - f(x)/fpx
	else
		throw(DomainError(x, "Derivative too close to zero"))
	end
end

# ╔═╡ eae93a04-2838-4e76-87df-c991ef4d0d87
function FindRoot1(a,b,f,fp,g,eps,delta,niter)
	for n = 1:niter
		if (abs(a-b) <= eps*abs(a)) && (abs(f(a)) <= delta)
			return a,n
		end
		x = g(a,b,f,fp)
		b = a; a = x
	end
	throw(ConvergenceException("Failed to converge in $niter iterations"))
end

# ╔═╡ 5d17ba68-6d1e-495f-ad1d-15c38b02b438
function bisect(a,b,f,fp)
  if (a>b) a,b=b,a end
  fa=f(a); fb=f(b)
  δ=b-a; c=(a+b)/2
  while ((δ>0.001) && (fa*fb<=0))
    δ=δ/2; c=a+δ; fc=f(c)
    if (fa*fc<=0)
      b,fb=c,fc
    else
      a,fa=c,fc
    end
  end
  return c
end

# ╔═╡ 67c25763-1bb5-4f02-bf61-4ed2d5dc79c5
function nthdegreepolynomial(B,n,x)
	p = B[n+1]
	sum = 0
	for i = 1:length(p)
		sum+= p[i]*x^(i-1)
	end
	return sum
end

# ╔═╡ f4a8480d-d62e-44ab-acb8-2d2c779caea5
function nthdegreepolynomialp(B,n,x)
	p = B[n+1]
	sum = 0
	for i = 1:length(p)
		sum += p[i]*(i-1)*x^(i-2)
	end
	return sum
end

# ╔═╡ f408185a-2e45-44cb-af25-501855aff113
#R2_1 = FindRoot1(-1.5,0.5,Phi2,Phi2p,Newton1,1e-10,1e-10,20)

# ╔═╡ f14daa4a-4b6b-4438-af62-06e5f918c318
B[10]

# ╔═╡ d89e9a1c-e39e-4e0b-8b6e-cc0b3c7a2f8e
function secant(x0,x1,f,fp)
	fx0 = f(x0);fx1 = f(x1)
	s = (fx1-fx0)/(x1-x0)
	return x1 - fx1/s
end

# ╔═╡ c55038f8-a8f2-406e-90d7-013f06c29ae7
function find_rootsbasis(B,n)
	function polynomial(x)
		return nthdegreepolynomial(B,n,x)
	end
	function polynomialp(x)
		return nthdegreepolynomialp(B,n,x)
	end
	roots = []
	a = -(n)/2; b = a+1
	for i in 1:n
		root, iterations = FindRoot1(a,b,polynomial,polynomialp,secant,100*eps(Float64),100*eps(Float64),40)
		push!(roots,[root,iterations])
		a = root+0.8; b = a + 0.6
	end
	return roots
end

# ╔═╡ 049d014c-dc60-4bf3-9c6e-96be466e8775
find_rootsbasis(B,4)

# ╔═╡ ee7b1f74-d4e6-48d3-a5d4-e70294d7e6b0
find_rootsbasis(B,6)

# ╔═╡ 8d09cbcf-e421-4f32-911c-e28565afff3b
find_rootsbasis(B,8)

# ╔═╡ be2d3404-270d-4d7b-ab2d-15acde813b28
find_rootsbasis(B,2)

# ╔═╡ b999cb05-7454-4bb9-bfaa-0cc9dca2fd82
find_rootsbasis(B,3)

# ╔═╡ 1a3e608e-71e4-4d9c-92d9-ecc08bd03b0c
find_rootsbasis(B,5)

# ╔═╡ 9c07b430-3199-4be5-9cef-94773193fa85
find_rootsbasis(B,8)

# ╔═╡ 52ef7b59-4e08-490d-a89f-c7fc9c0c39c8
find_rootsbasis(B,7)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╠═3849b625-12b8-4b0c-b893-10bea5bec05f
# ╠═0cc786f0-a50f-11ef-1124-fbdbbb7c69b3
# ╠═624f4b45-c518-47fa-9f4e-fcc2155dac1d
# ╠═ccb7b862-df6c-43fb-a42e-e8308509019e
# ╠═8ecd0fac-c9db-4996-bb38-2d12689e5da4
# ╠═26b34eb1-6e76-4b61-98ae-f94508871810
# ╠═b9b88780-34e8-44bd-962c-6c27b7a5bdc7
# ╠═1b63c404-15de-42e0-9e11-23f7e8845fa0
# ╠═eae93a04-2838-4e76-87df-c991ef4d0d87
# ╠═5d17ba68-6d1e-495f-ad1d-15c38b02b438
# ╠═67c25763-1bb5-4f02-bf61-4ed2d5dc79c5
# ╠═f4a8480d-d62e-44ab-acb8-2d2c779caea5
# ╠═f408185a-2e45-44cb-af25-501855aff113
# ╠═f14daa4a-4b6b-4438-af62-06e5f918c318
# ╠═d89e9a1c-e39e-4e0b-8b6e-cc0b3c7a2f8e
# ╠═c55038f8-a8f2-406e-90d7-013f06c29ae7
# ╠═049d014c-dc60-4bf3-9c6e-96be466e8775
# ╠═ee7b1f74-d4e6-48d3-a5d4-e70294d7e6b0
# ╠═8d09cbcf-e421-4f32-911c-e28565afff3b
# ╠═be2d3404-270d-4d7b-ab2d-15acde813b28
# ╠═b999cb05-7454-4bb9-bfaa-0cc9dca2fd82
# ╠═1a3e608e-71e4-4d9c-92d9-ecc08bd03b0c
# ╠═9c07b430-3199-4be5-9cef-94773193fa85
# ╠═52ef7b59-4e08-490d-a89f-c7fc9c0c39c8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
