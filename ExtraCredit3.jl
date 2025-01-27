### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3849b625-12b8-4b0c-b893-10bea5bec05f
using LinearAlgebra

# ╔═╡ 6c843526-c610-4606-9db1-5469fad46e9a
md"""
# Heat equation extra credit
"""

# ╔═╡ 1e5a55a7-0c64-41d7-b7ca-25356be4a0f8
md"""
###### $F(y) = \dfrac{1}{\sqrt{4\pi}}\int_{-\infty}^{\infty}e^{\dfrac{-(x-y)^{2}}{4}}f(x)dx$
"""

# ╔═╡ 9be5da64-88f5-45f3-be43-74c598a1ff6a
md"""
Let $t = \dfrac{(x-y)}{2}$
"""

# ╔═╡ e3c3a47e-e77e-4555-a8e3-1320ea9101a4
md"""
Then the integral becomes
###### $F(y) = \dfrac{1}{\sqrt{4\pi}}\int_{-\infty}^{\infty}e^{-t^{2}}f^{*}(t)dt$
"""

# ╔═╡ 45c803aa-3905-4f3c-9fd8-203f3574ef9f
md"""
Earlier the substituition was made $t=\dfrac{x-y}{2}$ was made. Since the initial heat distribution is $\cos(\sin(x)) + \sin(\cos(x))$, we must solve for x in terms of t, and y to find what to plug in for the funciton of $F(y)$, which is the heat disttribution after one unit time. random points y will be chosen to compare iteration values. Also, solving for x is $2t+y$. Subbing back into the funciton 2t + y, for an array of y's, we get our quadrature formula. Passing k through the function, k will be the array of y values of interest. And we can't forget, the differential substitution puts us off by a factor of 1/2, so must multiply the entire formula by 2 
"""

# ╔═╡ be211c5b-3eab-466f-abeb-c900ceaec05b
md"""
#### Ask ChatGPT: Write a script to foil polynomials in julia
"""

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

# ╔═╡ 765f0ef5-d1a8-497d-9110-46d854ada83b
md"""
The formula derived for the integral $I_k = \int_{-\infty}^{\infty}e^{-t^{2}}t^{k}dt$ is as follows
$I_k = \begin{cases}
\dfrac{(k-1)!!}{2^{(\dfrac{k}{2})}}\sqrt{\pi} & \text{if k is even}\\
0 & \text{if k is odd}
\end{cases}$

Where $n!! = n*(n-2)*(n-4)*...*(n-(n-1)?(n-n))$

The (n-1)?(n-n) means the last term depends on whether n is even or odd
"""

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
function Basis(n) #Gram_schmidt orthonormalization to generate basis
	z = []
	for i = 1:n+1 #ChatGPT's code wrote the row code
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

# ╔═╡ db7e86cb-265d-4026-9bb8-a2ba09af4e54
B = Basis(17);

# ╔═╡ 6e02001e-f179-47d9-a402-ee49df529bc9
innerproduct(B[13],B[13])

# ╔═╡ 0c9fab2f-d78d-4d7e-b75e-24c31a71a874
innerproduct(B[13],B[12])

# ╔═╡ 3c1fbadb-722a-4fbf-9ed3-32c2ac105c48
innerproduct(B[9],B[9])

# ╔═╡ 10e69088-95df-40a0-bbd1-7ea70809eef9
md"""
Build "identity" here
"""

# ╔═╡ ac1c22c8-eb6c-415f-a897-74f7bdf25360
md"""
Passes orthonormality w.r.t inner product
"""

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
function nthdegreepolynomialp(B,n,x) #Derivative only necessary for Newton's method
	p = B[n+1]
	sum = 0
	for i = 1:length(p)
		sum += p[i]*(i-1)*x^(i-2)
	end
	return sum
end

# ╔═╡ f14daa4a-4b6b-4438-af62-06e5f918c318
B[12]

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
		root, iterations = FindRoot1(a,b,polynomial,polynomialp,secant,1e-10,1e-10,40) #Decided on secant instead of newton, prioritizing finding roots to speed
		push!(roots,root)
		a = root+0.7; b = a + 0.5
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

# ╔═╡ 2abd0531-812e-484e-b7bc-4f045be6f777
find_rootsbasis(B,9)

# ╔═╡ 4805b26d-2ed0-46b5-9286-f160b30b4cb9
find_rootsbasis(B,10)

# ╔═╡ 1ea239a7-8263-4d12-a1e3-7f35183260ba
find_rootsbasis(B,11)

# ╔═╡ 0551f5be-d362-4204-bfeb-ca7ef8474013
R=find_rootsbasis(B,12)

# ╔═╡ 3ad28df4-5cbf-452c-abae-85241ff1ba4f
find_rootsbasis(B,13)

# ╔═╡ 640f864a-3c20-4a47-8fd5-92bc1ad23349
md"""
All the code from above is trouble shooting my guesses for where the roots will be. It seems difficult to get the roots accurate for the 13th degree polynomial, and doesn't seem worth it to try and improve the guess to go higher, as each guess at this point has a high risk of messing up the previous polynomial's roots
"""

# ╔═╡ a1dd59a4-7ceb-4384-8a5b-71c9f7aacb70
function findweights(B,n)
	R = find_rootsbasis(B,n)
	A = [R[j]^(i-1) for i in 1:n, j in 1:length(R)] #Claude A.i. helped me implement this matrix
	b = zeros(n)
	for i = 1:n
		if i % 2 == 0
			b[i] = 0
		else
		v = zeros(Int((i+1)/2))
		v[end] = 1
		b[i] = innerproduct(v,v)
		end
	end
	W = A\b
	return W		
end

# ╔═╡ 36a62fa8-f9a6-4409-bce1-dfffecb1fc6c
findweights(B,4)

# ╔═╡ c64288e3-60be-4b7d-8dd3-8f64a1de5e69
f(x) = cos(sin(x)) + sin(cos(x)) #Function that I'm trying to integrate

# ╔═╡ 616dc57f-717e-4465-81e5-a3fe0b4c5592
K = [i for i=0:3] #These are the y points

# ╔═╡ 10e386c2-ec87-4219-81b9-3c71b1c6edc7
function Gauss(B,n,f,k)
	W = findweights(B,n); R = find_rootsbasis(B,n); sum = 0
	for i = 1:length(W)
		sum += W[i]*f(2*R[i]+k) #Transformation of evaluation point based on substitution
	end
	return 2*sum #Can't forget factor of two from substitution of differential
end

# ╔═╡ 2de2a0db-45e1-46e0-acaf-a9cc951f56b3
intcalculatorvalues = (1/sqrt(4pi)).*[3.875197554120239,3.326488800239226,2.225158038987561,1.590644453177028]

# ╔═╡ 018224cb-6e74-4a5f-966a-85a6fb571b44
[Gauss(B,12,f,i) for i = 0:3] #n = 12 evaluated at all sample points

# ╔═╡ 9d924e36-b736-49f5-b3e5-456a63737973
[Gauss(B,11,f,i) for i = 0:3]

# ╔═╡ 0b3792e2-709c-42a9-adc8-22aca86d9825
[Gauss(B,10,f,i) for i = 0:3]

# ╔═╡ d1679994-b5b7-4894-8f0e-be54951065c8
[Gauss(B,9,f,i) for i = 0:3]

# ╔═╡ e2eb041c-f501-42de-ad32-0a8cfbe887a4
[Gauss(B,8,f,i) for i = 0:3]

# ╔═╡ 380cf422-e6e0-4d73-82ed-b7752647b40a
[Gauss(B,7,f,i) for i = 0:3]

# ╔═╡ 8edc181b-29c7-4cb6-9d44-c1396d8bb9e8
[Gauss(B,6,f,i) for i = 0:3]

# ╔═╡ 786b8517-2fdc-4491-a171-beeb5d330e71
[Gauss(B,5,f,i) for i = 0:3]

# ╔═╡ 43cb71e1-c5ae-4da0-aa79-83ddc912651c
[Gauss(B,4,f,i) for i = 0:3]

# ╔═╡ 9b03b118-3d18-4b2e-97a0-ec1a05c4a0cb
[Gauss(B,3,f,i) for i = 0:3]

# ╔═╡ 23e0a976-425f-4ee9-80ad-76199d0ada6b
[Gauss(B,2,f,i) for i = 0:3]

# ╔═╡ 3ab0b98e-81c3-4d5f-80cf-2be79cf5e988
[Gauss(B,1,f,i) for i = 0:3]

# ╔═╡ ab7a46d5-8d87-4ca3-8e89-5414e6cc1141
md"""
Am able to get three significant digits of accuracy using this Gauss method, any more digits would be difficult
"""

# ╔═╡ 0c658046-9a18-42db-aeac-e41ba3fcc284
md"""
Based on physical approximation, the F(y) should be 0 far from the origin
"""

# ╔═╡ 384fedb4-464e-44ed-97e4-9f753e7a59c9


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─6c843526-c610-4606-9db1-5469fad46e9a
# ╟─1e5a55a7-0c64-41d7-b7ca-25356be4a0f8
# ╠═3849b625-12b8-4b0c-b893-10bea5bec05f
# ╟─9be5da64-88f5-45f3-be43-74c598a1ff6a
# ╟─e3c3a47e-e77e-4555-a8e3-1320ea9101a4
# ╟─45c803aa-3905-4f3c-9fd8-203f3574ef9f
# ╟─be211c5b-3eab-466f-abeb-c900ceaec05b
# ╠═0cc786f0-a50f-11ef-1124-fbdbbb7c69b3
# ╠═624f4b45-c518-47fa-9f4e-fcc2155dac1d
# ╟─765f0ef5-d1a8-497d-9110-46d854ada83b
# ╠═ccb7b862-df6c-43fb-a42e-e8308509019e
# ╠═8ecd0fac-c9db-4996-bb38-2d12689e5da4
# ╠═db7e86cb-265d-4026-9bb8-a2ba09af4e54
# ╠═6e02001e-f179-47d9-a402-ee49df529bc9
# ╠═0c9fab2f-d78d-4d7e-b75e-24c31a71a874
# ╠═3c1fbadb-722a-4fbf-9ed3-32c2ac105c48
# ╠═10e69088-95df-40a0-bbd1-7ea70809eef9
# ╟─ac1c22c8-eb6c-415f-a897-74f7bdf25360
# ╠═b9b88780-34e8-44bd-962c-6c27b7a5bdc7
# ╠═1b63c404-15de-42e0-9e11-23f7e8845fa0
# ╠═eae93a04-2838-4e76-87df-c991ef4d0d87
# ╠═67c25763-1bb5-4f02-bf61-4ed2d5dc79c5
# ╠═f4a8480d-d62e-44ab-acb8-2d2c779caea5
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
# ╠═2abd0531-812e-484e-b7bc-4f045be6f777
# ╠═4805b26d-2ed0-46b5-9286-f160b30b4cb9
# ╠═1ea239a7-8263-4d12-a1e3-7f35183260ba
# ╠═0551f5be-d362-4204-bfeb-ca7ef8474013
# ╠═3ad28df4-5cbf-452c-abae-85241ff1ba4f
# ╟─640f864a-3c20-4a47-8fd5-92bc1ad23349
# ╠═a1dd59a4-7ceb-4384-8a5b-71c9f7aacb70
# ╠═36a62fa8-f9a6-4409-bce1-dfffecb1fc6c
# ╠═c64288e3-60be-4b7d-8dd3-8f64a1de5e69
# ╠═616dc57f-717e-4465-81e5-a3fe0b4c5592
# ╠═10e386c2-ec87-4219-81b9-3c71b1c6edc7
# ╠═2de2a0db-45e1-46e0-acaf-a9cc951f56b3
# ╠═018224cb-6e74-4a5f-966a-85a6fb571b44
# ╠═9d924e36-b736-49f5-b3e5-456a63737973
# ╠═0b3792e2-709c-42a9-adc8-22aca86d9825
# ╠═d1679994-b5b7-4894-8f0e-be54951065c8
# ╠═e2eb041c-f501-42de-ad32-0a8cfbe887a4
# ╠═380cf422-e6e0-4d73-82ed-b7752647b40a
# ╠═8edc181b-29c7-4cb6-9d44-c1396d8bb9e8
# ╠═786b8517-2fdc-4491-a171-beeb5d330e71
# ╠═43cb71e1-c5ae-4da0-aa79-83ddc912651c
# ╠═9b03b118-3d18-4b2e-97a0-ec1a05c4a0cb
# ╠═23e0a976-425f-4ee9-80ad-76199d0ada6b
# ╠═3ab0b98e-81c3-4d5f-80cf-2be79cf5e988
# ╟─ab7a46d5-8d87-4ca3-8e89-5414e6cc1141
# ╠═0c658046-9a18-42db-aeac-e41ba3fcc284
# ╠═384fedb4-464e-44ed-97e4-9f753e7a59c9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
