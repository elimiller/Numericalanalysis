### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 7568e14a-a5cd-4a7c-aad6-fe385b4b6ee0
using LinearAlgebra

# ╔═╡ 40aca690-a3a6-11ef-0c19-c9831e2db6e4
md"""
# Heat equation extra credit
"""

# ╔═╡ a046db6b-8eff-4b34-b4dc-96185eae9bec
md"""
###### $F(y) = \dfrac{1}{\sqrt{4\pi}}\int_{-\infty}^{\infty}e^{\dfrac{-(x-y)^{2}}{4}}f(x)dx$
"""

# ╔═╡ fc97b29b-9370-4e12-ac3d-61cc3ab831f7


# ╔═╡ febe189b-6aef-440e-93e3-c62c4e84a6eb
md"""
Let $t = \dfrac{(x-y)}{2}$
"""

# ╔═╡ 32603b51-a31d-4790-861a-3c027db655db
md"""
Then the integral becomes
###### $F(y) = \dfrac{1}{\sqrt{4\pi}}\int_{-\infty}^{\infty}e^{-t^{2}}f(t)dt$
"""

# ╔═╡ a0deba76-8a95-4e71-a87c-52283af15a97
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

# ╔═╡ f5da61bd-6154-4083-a11b-6a203a2079df
x = [0,2,1]; y = [0,1]

# ╔═╡ 0392b7bb-dd4c-4fcb-b57a-ad35311946f8
polynomialfoil(x,y)

# ╔═╡ 23d7b167-4ce3-4991-9fbd-9e48898074d0
function doublefactorial(n)
	list = [i for i = 1:n]
	product = 1
	for j = 0:2:n-1
		product *= list[n-j]
	end
	return product
end

# ╔═╡ 4eee3025-6a55-41f1-94b3-2d6791e61257
doublefactorial(5)

# ╔═╡ 6f54c03a-4ecd-47e5-ba7e-0a3d32b97003
doublefactorial(4)

# ╔═╡ fe0f6e86-bd20-42f2-a1cf-e98c292c3c5b
function innerproduct(f,g)
	h = polynomialfoil(f,g)
	sum = 0 ; a = sqrt(pi)
	 n = length(h)
	for i = 0:2:n-1 #Degree is n-1, so for odd it's zero. Even are i-1
		sum += h[i+1]*doublefactorial(i-1)/2^(i/2)
	end
	return (a/sqrt(4*pi))*sum
end
	

# ╔═╡ 9a86bebd-2eb7-4a96-93bc-6114a61ddf9d
innerproduct([0,0,1],[0,0,1])

# ╔═╡ 8c5fbdf4-4920-4a13-a40a-cf2981868b7e
sqrt(pi)

# ╔═╡ 0510381c-c5c7-4dc0-b9b8-984f497cc6e5
innerproduct([1],[0,1])

# ╔═╡ 6f08ba0a-d3e0-4944-b3c1-258c1fab8d8d
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

# ╔═╡ ec827bac-c0c4-4b69-8d90-655caa262481
B = Basis(5)

# ╔═╡ ee717841-ec9f-49e2-8e88-76a8e79db970
innerproduct(B[2],B[2])

# ╔═╡ b32fda0d-4548-415b-b1eb-953a4a0a9641
innerproduct(B[4],B[3])

# ╔═╡ d85ee79a-1bae-4c9c-850e-74d5867795c4
innerproduct(B[4],B[4])

# ╔═╡ 17404900-4e50-4472-933c-458734f77d49
function Newton1(x,x1,f,fp)
	fpx = fp(x)
	if abs(fpx) > 100 *  eps(Float64)
		return x - f(x)/fpx
	else
		exit(-1)
	end
end

# ╔═╡ 75e7ef11-0755-483f-abff-6bc72284ae9a
function FindRoot1(a,b,f,fp,g,eps,delta,niter)
	for n = 1:niter
		if (abs(a-b) <= eps*abs(a)) && (abs(f(a)) <= delta)
			return a,n
		end
		x = g(a,b,f,fp)
		b = a; a = x
	end
	exit(-2)
end

# ╔═╡ d9030718-fe65-4a1e-8cda-7af21b0ee458
md"""
Will start churning out the Gauss integrals, based on n points
"""

# ╔═╡ fb358092-3af3-49f9-8dcb-cbb7c505e830
phi2 = B[3]

# ╔═╡ daa64d26-5aae-4a4d-82a4-5f082251bafa
Phi2(x) = phi2[1] + phi2[3]*x^2

# ╔═╡ 9b35db3f-d0c9-429e-9155-07982db84b35
Phi2p(x) = 2*phi2[3]*x

# ╔═╡ d4d849e1-bee4-40f2-bae1-de299d44bfe7
#R2 = FindRoot1(0,0.2,phi3,phi3p,Newton1,1e-10,1e-10,20)

# ╔═╡ 6c47327e-8bf6-4a2c-bf97-1058a7da2dd3
R2_1 = FindRoot1(-1.5,0.5,Phi2,Phi2p,Newton1,1e-10,1e-10,20)

# ╔═╡ a35a03b4-9080-43fd-8591-1ceb3c2dbd4b
R2_2 = FindRoot1(0.5,1,Phi2,Phi2p,Newton1,1e-10,1e-10,20)

# ╔═╡ 3f3e926d-0125-4b80-9eb1-4d21f40610a3
r2_1 = R2_1[1]; r2_2 = R2_2[1]

# ╔═╡ c1177a92-9154-4634-95da-b5dfcdbca636
A2 = [1 1; r2_1 r2_2]

# ╔═╡ 792262d3-335e-41cb-b206-f18c7c858cfc
b2 = [innerproduct([1],[1]);0]

# ╔═╡ c27353e2-6ef7-4124-9e17-78fe739bd6df
W2 = A2\b2

# ╔═╡ a90bb5e0-4ea7-4eba-9c95-5b7e3d3a8618
md"""
Earlier the substituition was made $t=\dfrac{x-y}{2}$ was made. Since the initial heat distribution is $\cos(\sin(x)) + \sin(\cos(x))$, we must solve for x in terms of t, and y to find what to plug in for the funciton of $F(y)$, which is the heat disttribution after one unit time. random points y will be chosen to compare iteration values. This vector will be k, which has been passed through the gauss functions. Also, solving for x is $2t+y$. Subbing back into the funciton 2t + y, for an array of y's, we get our quadrature formula. Passing k through the function, k will be the array of y values of interest. 
"""

# ╔═╡ bc537547-e6e6-428a-9a3d-14bb63199ff3
function Gauss2test(f,k)
	w1 = W2[1]; w2 = W2[2]
	return w1*f(r2_1) + w2*f(r2_2)
end

# ╔═╡ 122a0f69-f229-4942-b2c7-4232677169cf
test(x) = 3x^3 +2x -1

# ╔═╡ 0cb70469-093f-4132-8845-14a956a13f9d
Gauss2test(test,1)

# ╔═╡ 60ced504-525f-47ee-89ca-ce8bac38825f
function Gauss2(f,k)
	w1 = W2[1]; w2 = W2[2]
	return w1*f(2*r2_1+k) + w2*f(2*r2_2+k)
end

# ╔═╡ 906c0ee5-d02a-46a9-a344-e144c82a60b2
phi3 = B[4]

# ╔═╡ e1db3fba-2ee7-4c7d-8a37-55ba219eb93f
Phi3(x) = phi3[2]*x + phi3[4]*x^3

# ╔═╡ 5d1196a3-f698-47e3-a7da-a9dcbe5be81e
Phi3p(x) = phi3[2] + 3*phi3[4]*x^2

# ╔═╡ a569e02c-a546-465f-94d2-4cbb9826281f
R3_1 = FindRoot1(-2,-1,Phi3,Phi3p,Newton1,1e-10,1e-10,20)

# ╔═╡ a6cc5fbb-3a12-4de6-ab7b-8151ddfdfe30
R3_2 = FindRoot1(-0.5,0.5,Phi3,Phi3p,Newton1,1e-10,1e-10,20)

# ╔═╡ ed6c6519-54fb-48f3-8233-53b64ae1233c
R3_3 = FindRoot1(1,2,Phi3,Phi3p,Newton1,1e-10,1e-10,20)

# ╔═╡ b1e78c65-f530-41ce-a08a-d113433568e9
r3_1 = R3_1[1]; r3_2 = R3_2[1]; r3_3 = R3_3[1]

# ╔═╡ 8e2217b3-ef3c-4f25-bf3c-9c16933090dc
A3 = [1 1 1;r3_1 r3_2 r3_3;r3_1^2 r3_2^2 r3_3^2]

# ╔═╡ c17819f7-b446-471d-90bd-79c2e2a778fe
b3 = [innerproduct([1],[1]);0;innerproduct([0,1],[0,1])]

# ╔═╡ d838413a-dc20-4df7-9bc4-a7666e972a94
W3 = A3\b3

# ╔═╡ bf3ee2e3-3637-4827-a6b5-629b564c1c25
function Gauss3(f,k)
	w1 = W3[1]; w2 = W3[2]; w3 = W3[3]
	return (1/sqrt(4*pi))*w1*f(2*r3_1+k) + w2*f(2*r3_2+k) + w3*f(2*r3_3+k)
end

# ╔═╡ 66ce8e1d-c816-49e0-9275-6332d3339e86
function Gauss3test(f,k)
	w1 = W3[1]; w2 = W3[2]; w3 = W3[3]
	return (1/sqrt(4*pi))*w1*f(r3_1*2+k) + w2*f(r3_2*2+k) + w3*f(r3_3*2+k)
end

# ╔═╡ dd911ca6-e41f-4192-979c-7dc7017c8568
phi4 = B[5]

# ╔═╡ 850db5d8-f5f6-48be-9496-51860e38f129
Phi4(x) = (phi4[5]*x^2 + phi4[3])*x^2 + phi4[1] 

# ╔═╡ b81081c0-c4b1-45f3-9533-90fa4a593267
Phi4p(x) = (4*phi4[5]*x^2 + 2*phi4[3])*x 

# ╔═╡ 6f054971-f602-4b8d-9f9c-87ebf80eb756
R4_1 = FindRoot1(-2,-1.5,Phi4,Phi4p,Newton1,1e-10,1e-10,20)

# ╔═╡ 5e65bb1c-9f30-46ae-9f70-9f8cc1b1d5a8
R4_2 = FindRoot1(-1,-0.5,Phi4,Phi4p,Newton1,1e-10,1e-10,20)

# ╔═╡ 24694e29-8f91-440a-9ddb-12b4653e72a2
R4_3 = FindRoot1(0.5,1.0,Phi4,Phi4p,Newton1,1e-10,1e-10,20)

# ╔═╡ 18e2137e-7f4e-4e22-b434-92a478a86326
R4_4 = FindRoot1(1.5,2,Phi4,Phi4p,Newton1,1e-10,1e-10,20)

# ╔═╡ 169d5c32-8fa3-46f0-bc70-1ea72b3b3a55
r4_1 = R4_1[1]; r4_2 = R4_2[1];r4_3 = R4_3[1];r4_4 = R4_4[1];

# ╔═╡ 47317edc-2229-4b66-a36c-f180fead7d99
A4 = [1 1 1 1 ;r4_1 r4_2 r4_3 r4_4;r4_1^2 r4_2^2 r4_3^2 r4_4^2; r4_1^3 r4_2^3 r4_3^3 r4_4^3]

# ╔═╡ 11c9a924-677c-48b0-a613-0c567acbc93e
[r4_1 ,r4_2, r4_3 ,r4_4]

# ╔═╡ fd3d454e-f0c4-4dfd-b488-bfc46e294cf4
b4 = [innerproduct([1],[1]);0;innerproduct([0,1],[0,1]);0]

# ╔═╡ ff4858ee-ae91-473c-8cd9-812720ec0e9d
W4 = A4\b4

# ╔═╡ ed06b821-0710-4a48-ad25-7ccb4343dacb
function Gauss4(f,k)
	w1 = W4[1];w2 = W4[2];w3 = W4[3];w4 = W4[4]
	return (1/sqrt(4*pi))*w1*f(2*r4_1+k)+w2*f(2*r4_2+k)+w3*f(2*r4_3+k)+w4*f(2*r4_4+k)
end

# ╔═╡ 484aa27e-7509-44e9-b269-86b9ab839973
phi5 = B[6]

# ╔═╡ 1667da50-0575-41a2-afaf-899b1b5256ef
Phi5(x) = (phi5[6]*x^2 + phi5[4])*x^3 + phi5[2]*x

# ╔═╡ ca6f3cfd-326d-4d55-99d8-ebc99741d9e3
Phi5p(x) = (5*phi5[6]x^2 + 3*phi5[4])*x^2 + phi5[2]

# ╔═╡ 20a4092f-dbef-4e79-ba02-58a2c2239786
R5_1 = FindRoot1(-2.5,-2,Phi5,Phi5p,Newton1,1e-10,1e-10,20);R5_2 = FindRoot1(-1.3,-1,Phi5,Phi5p,Newton1,1e-10,1e-10,20);R5_3 = FindRoot1(-0.1,0.1,Phi5,Phi5p,Newton1,1e-10,1e-10,20);R5_4 = FindRoot1(1,1.3,Phi5,Phi5p,Newton1,1e-10,1e-10,20);R5_5 = FindRoot1(2,2.5,Phi5,Phi5p,Newton1,1e-10,1e-10,20)

# ╔═╡ b88f3720-a6f0-4cb4-8caa-73a41078b57d
r5_1 = R5_1[1]; r5_2 = R5_2[1];r5_3 = R5_3[1];r5_4 = R5_4[1];r5_5 = R5_5[1]

# ╔═╡ dbc53e20-619d-4879-b46b-c30bab0fd9ef
A5 = [1 1 1 1 1;r5_1 r5_2 r5_3 r5_4 r5_5;r5_1^2 r5_2^2 r5_3^2 r5_4^2 r5_5^2;r5_1^3 r5_2^3 r5_3^3 r5_4^3 r5_5^3;r5_1^4 r5_2^4 r5_3^4 r5_4^4 r5_5^4]

# ╔═╡ b8568d67-149c-4b05-959d-da68c9cce24a
b5 = [innerproduct([1],[1]);0;innerproduct([0,1],[0,1]);0;innerproduct([0,0,1],[0,0,1])]

# ╔═╡ ee7b024b-2372-42aa-aac4-a5f4b6dfe1b8
W5 = A5\b5

# ╔═╡ b80f1681-9e80-47ad-9251-8c96a7a78217
function Gauss5(f,k)
	w1 = W5[1];w2 = W5[2];w3 = W5[3];w4 = W5[4];w5 = W5[5];
	return (1/sqrt(4*pi))*w1*f(2*r5_1+k)+w2*f(2*r5_2+k)+w3*f(2*r5_3+k)+w4*f(2*r5_4+k)+w5*f(2*r5_5 +k)
end

# ╔═╡ daf41caa-a4da-4af9-9dc8-f609e81a7ed9
function f(x)
	return sin(cos(x)) + cos(sin(x))
end

# ╔═╡ 62abf93c-a750-4119-b0db-ffbcc6d639e7
md"""
Interested in y = 0,1,2,3
"""

# ╔═╡ 3c005c45-1902-43a7-af6b-c828e26a7407
K = [i for i = 0:3]

# ╔═╡ 1b4a817b-7b66-42e9-8d33-ae25243220a6
Gauss3test.(f,K)

# ╔═╡ e04204bd-e9dc-4e0c-bf8b-d14e1081d2ec
md"""
Integral calculators estimated values for 
###### $F(y) = \dfrac{1}{\sqrt{4\pi}}\int_{-\infty}^{\infty}e^{\dfrac{-(x-y)^{2}}{4}}(\cos(\sin(x))+\sin(\cos(x)))dx$
are

[1.09317
0.938385
0.627705
0.448713]
"""

# ╔═╡ 5363e1fd-8ba2-445b-8f1f-353f94cb88a0
intcalculatorvalues = (1/sqrt(4pi)).*[3.875197554120239,3.326488800239226,2.225158038987561,1.590644453177028]

# ╔═╡ 2dbdef20-5c4a-4dc5-a009-c1382a97bc04
intcalculatorvalues

# ╔═╡ 5e9ee408-00fe-4137-89bf-5d3681499fb4
Gauss2.(f,K)

# ╔═╡ e75dfd47-2be5-4152-9e33-554686fd61ef
Gauss3.(f,K)

# ╔═╡ ac3f9948-73a6-4825-a32d-87ab1213de74
Gauss4.(f,K)

# ╔═╡ deea5c6c-68e1-4457-86b7-be0468037482
Gauss5.(f,K)

# ╔═╡ 9bc09662-419f-40df-8d32-040a90a4036a
Gauss5.(f,K)

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
# ╠═40aca690-a3a6-11ef-0c19-c9831e2db6e4
# ╠═a046db6b-8eff-4b34-b4dc-96185eae9bec
# ╠═fc97b29b-9370-4e12-ac3d-61cc3ab831f7
# ╠═febe189b-6aef-440e-93e3-c62c4e84a6eb
# ╠═32603b51-a31d-4790-861a-3c027db655db
# ╠═a0deba76-8a95-4e71-a87c-52283af15a97
# ╠═f5da61bd-6154-4083-a11b-6a203a2079df
# ╠═0392b7bb-dd4c-4fcb-b57a-ad35311946f8
# ╠═23d7b167-4ce3-4991-9fbd-9e48898074d0
# ╠═4eee3025-6a55-41f1-94b3-2d6791e61257
# ╠═6f54c03a-4ecd-47e5-ba7e-0a3d32b97003
# ╠═fe0f6e86-bd20-42f2-a1cf-e98c292c3c5b
# ╠═9a86bebd-2eb7-4a96-93bc-6114a61ddf9d
# ╠═8c5fbdf4-4920-4a13-a40a-cf2981868b7e
# ╠═0510381c-c5c7-4dc0-b9b8-984f497cc6e5
# ╠═6f08ba0a-d3e0-4944-b3c1-258c1fab8d8d
# ╠═ec827bac-c0c4-4b69-8d90-655caa262481
# ╠═ee717841-ec9f-49e2-8e88-76a8e79db970
# ╠═b32fda0d-4548-415b-b1eb-953a4a0a9641
# ╠═d85ee79a-1bae-4c9c-850e-74d5867795c4
# ╠═17404900-4e50-4472-933c-458734f77d49
# ╠═75e7ef11-0755-483f-abff-6bc72284ae9a
# ╟─d9030718-fe65-4a1e-8cda-7af21b0ee458
# ╠═fb358092-3af3-49f9-8dcb-cbb7c505e830
# ╠═daa64d26-5aae-4a4d-82a4-5f082251bafa
# ╠═9b35db3f-d0c9-429e-9155-07982db84b35
# ╠═d4d849e1-bee4-40f2-bae1-de299d44bfe7
# ╠═6c47327e-8bf6-4a2c-bf97-1058a7da2dd3
# ╠═a35a03b4-9080-43fd-8591-1ceb3c2dbd4b
# ╠═3f3e926d-0125-4b80-9eb1-4d21f40610a3
# ╠═c1177a92-9154-4634-95da-b5dfcdbca636
# ╠═792262d3-335e-41cb-b206-f18c7c858cfc
# ╠═7568e14a-a5cd-4a7c-aad6-fe385b4b6ee0
# ╠═c27353e2-6ef7-4124-9e17-78fe739bd6df
# ╟─a90bb5e0-4ea7-4eba-9c95-5b7e3d3a8618
# ╠═bc537547-e6e6-428a-9a3d-14bb63199ff3
# ╠═122a0f69-f229-4942-b2c7-4232677169cf
# ╠═0cb70469-093f-4132-8845-14a956a13f9d
# ╠═60ced504-525f-47ee-89ca-ce8bac38825f
# ╠═906c0ee5-d02a-46a9-a344-e144c82a60b2
# ╠═e1db3fba-2ee7-4c7d-8a37-55ba219eb93f
# ╠═5d1196a3-f698-47e3-a7da-a9dcbe5be81e
# ╠═a569e02c-a546-465f-94d2-4cbb9826281f
# ╠═a6cc5fbb-3a12-4de6-ab7b-8151ddfdfe30
# ╠═ed6c6519-54fb-48f3-8233-53b64ae1233c
# ╠═b1e78c65-f530-41ce-a08a-d113433568e9
# ╠═8e2217b3-ef3c-4f25-bf3c-9c16933090dc
# ╠═c17819f7-b446-471d-90bd-79c2e2a778fe
# ╠═d838413a-dc20-4df7-9bc4-a7666e972a94
# ╠═bf3ee2e3-3637-4827-a6b5-629b564c1c25
# ╠═66ce8e1d-c816-49e0-9275-6332d3339e86
# ╠═1b4a817b-7b66-42e9-8d33-ae25243220a6
# ╠═dd911ca6-e41f-4192-979c-7dc7017c8568
# ╠═850db5d8-f5f6-48be-9496-51860e38f129
# ╠═b81081c0-c4b1-45f3-9533-90fa4a593267
# ╠═6f054971-f602-4b8d-9f9c-87ebf80eb756
# ╠═5e65bb1c-9f30-46ae-9f70-9f8cc1b1d5a8
# ╠═24694e29-8f91-440a-9ddb-12b4653e72a2
# ╠═18e2137e-7f4e-4e22-b434-92a478a86326
# ╠═169d5c32-8fa3-46f0-bc70-1ea72b3b3a55
# ╠═47317edc-2229-4b66-a36c-f180fead7d99
# ╠═11c9a924-677c-48b0-a613-0c567acbc93e
# ╠═fd3d454e-f0c4-4dfd-b488-bfc46e294cf4
# ╠═ff4858ee-ae91-473c-8cd9-812720ec0e9d
# ╠═ed06b821-0710-4a48-ad25-7ccb4343dacb
# ╠═484aa27e-7509-44e9-b269-86b9ab839973
# ╠═1667da50-0575-41a2-afaf-899b1b5256ef
# ╠═ca6f3cfd-326d-4d55-99d8-ebc99741d9e3
# ╠═20a4092f-dbef-4e79-ba02-58a2c2239786
# ╠═b88f3720-a6f0-4cb4-8caa-73a41078b57d
# ╠═dbc53e20-619d-4879-b46b-c30bab0fd9ef
# ╠═b8568d67-149c-4b05-959d-da68c9cce24a
# ╠═ee7b024b-2372-42aa-aac4-a5f4b6dfe1b8
# ╠═b80f1681-9e80-47ad-9251-8c96a7a78217
# ╠═daf41caa-a4da-4af9-9dc8-f609e81a7ed9
# ╟─62abf93c-a750-4119-b0db-ffbcc6d639e7
# ╠═3c005c45-1902-43a7-af6b-c828e26a7407
# ╟─e04204bd-e9dc-4e0c-bf8b-d14e1081d2ec
# ╠═5363e1fd-8ba2-445b-8f1f-353f94cb88a0
# ╠═2dbdef20-5c4a-4dc5-a009-c1382a97bc04
# ╠═5e9ee408-00fe-4137-89bf-5d3681499fb4
# ╠═e75dfd47-2be5-4152-9e33-554686fd61ef
# ╠═ac3f9948-73a6-4825-a32d-87ab1213de74
# ╠═deea5c6c-68e1-4457-86b7-be0468037482
# ╠═9bc09662-419f-40df-8d32-040a90a4036a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
