### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ b135e352-4e49-4e0e-8876-e9c0e8965dd5
using Pkg; Pkg.add("Bessels")

# ╔═╡ 81635f89-cf57-4e4c-8ff4-2e83ee4168f1
using Bessels

# ╔═╡ 9e3a44ad-3640-45e1-9e15-f9adf8a87921
besselj(0,3)

# ╔═╡ ddbfbf50-843c-11ef-05e7-6d8faf65e014
function bisect(f,a,b,ε)
  if (a>b) a,b=b,a end
  fa=f(a); fb=f(b)
  δ=b-a; c=(a+b)/2
  while ((δ>ε) && (fa*fb<=0))
    δ=δ/2; c=a+δ; fc=f(c)
    if (fa*fc<=0)
      b,fb=c,fc
    else
      a,fa=c,fc
    end
  end
  return c
end

# ╔═╡ 1a5dc0fb-1f3d-4d38-8129-0f0c3ba17a9d
function findmachineps()
	eps = 1.0
	while 1 + eps/2 != 1
		eps /= 2
	end
	return eps
end
		

# ╔═╡ eed4c9dd-5539-48cc-9e63-cc2e4dd92c25
MaxchEps = findmachineps()

# ╔═╡ 8aaf481e-8bf5-43e1-862f-200d5e9fe9b7
function secant(xn,xnm,f)
	yn = f(xn); ynm = f(xnm)
	if abs(min(yn,ynm)) < MaxchEps
		exit
	else
		xnp = xn - yn*(xn-xnm)/(yn-ynm)
		end
end

# ╔═╡ a84ea700-f279-4c2d-b063-c721195ca859
f(x)=x^2-2; a=1; b=2; ε=0.01;

# ╔═╡ 824b257c-c6d3-40bd-aace-998140a77ac2
secant(1,2,f)

# ╔═╡ 6b8d547d-bc6e-473a-8e76-058350fd82cb
function FindRoot(a,b,f,g,eps,delta,niter)
	x0 = a + 0.5*(b-a)
	for n = 1:niter
		x1 = g(x0,f)
	end
end

# ╔═╡ Cell order:
# ╠═b135e352-4e49-4e0e-8876-e9c0e8965dd5
# ╠═81635f89-cf57-4e4c-8ff4-2e83ee4168f1
# ╠═9e3a44ad-3640-45e1-9e15-f9adf8a87921
# ╠═ddbfbf50-843c-11ef-05e7-6d8faf65e014
# ╠═8aaf481e-8bf5-43e1-862f-200d5e9fe9b7
# ╠═824b257c-c6d3-40bd-aace-998140a77ac2
# ╠═1a5dc0fb-1f3d-4d38-8129-0f0c3ba17a9d
# ╠═eed4c9dd-5539-48cc-9e63-cc2e4dd92c25
# ╠═a84ea700-f279-4c2d-b063-c721195ca859
# ╠═6b8d547d-bc6e-473a-8e76-058350fd82cb
