### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d776d4d6-825b-11ef-3d58-85ab2150ff5f
using LinearAlgebra, Colors, ForwardDiff, PlutoUI

# ╔═╡ e4a1845e-0d27-4f0c-abbb-11960b8a2452
using PlotlyLight

# ╔═╡ 812d53c7-1384-4c77-9185-cf3d6656b98b
md"""
# Dérivée à plusieurs variables

## De droite tangente à plan tangent

La dérivée univariée au point $a$ correspond à la pente de la **droite tangente** à la fonction en $a$.
Pour une fonction bivariée, la dérivée dans une direction $d$ correspond à la pente du **plan tangent** à la fonction en $a$ dans la direction $d$.
"""

# ╔═╡ b3c01426-9a70-448b-bf3e-901e3100082f
begin
	x0 = -4
	x1 = 4
	@bind cx Slider(range(x0, stop = x1 - 1, length = 101), default=(x1 + x0) / 2, show_value = true)
end

# ╔═╡ 674be307-64c1-4485-ad25-a59ff6840df8
begin
	y0 = -4
	y1 = 4
	@bind cy Slider(range(y0, stop = y1 - 1, length = 101), default=(y1 + y0) / 2, show_value = true)
end

# ╔═╡ 0356921d-c3e5-4351-9c14-0efd657be7dd
md"""
## Dérivée directionelle

Que vaut la dérivée dans la direction rouge `(1, 1)` si on connait les dérivée vertes en `x` : $\partial f/\partial x$ et en `y` : $\partial f / \partial y$ ?
Un plan est défini par un point et deux vecteurs donc le plan bleu est entièrement défini par le point $(x, f(x))$ et les deux vecteurs $(1, 0, \partial f/\partial x)$ et $(0, 1, \partial f/\partial y)$.
En utilisant la linéarité du plan, quelle est la valeur de $?$ pour que le vecteur $(1, 1, ?)$ soit dans le plan ?

$$(1, 1, \partial f/\partial x + \partial f/\partial y)$$

En général, pour une direction arbitraire, on a

$$(d_x, d_y, d_x \cdot \partial f/\partial x + d_y \cdot \partial f/\partial y)$$

La dérivée en direction $(d_x, d_y)$ est donc obtenue par produit scalaire entre le vecteur $(d_x, d_y)$ et le *gradient* $(\partial f/\partial x, \partial f/\partial y)$. Voici trois notations équivalentes pour écrire ce produit scalaire:

$$\langle (d_x, d_y),  (\partial f/\partial x, \partial f/\partial y) \rangle
=
\begin{bmatrix}d_x & d_y\end{bmatrix}
\begin{bmatrix}
  \partial f/\partial x \\ \partial f/\partial y
\end{bmatrix}
=
d_x \cdot \partial f/\partial x + d_y \cdot \partial f/\partial y$$

Le gradient est souvent dénoté avec le symbol $\nabla$:

$$\nabla f = (\partial f/\partial x, \partial f/\partial y)$$
"""

# ╔═╡ b27a8a10-c3fa-44c6-917d-14959f4f30dd
md"""
## Calculer le gradient à la main

Comment calculer les valeurs de $\partial f/\partial x$ et $\partial f/\partial y$ ?
Commençant par voir comment faire ça à la main. Nous verrons au cours suivant comment le calculer algorithmiquement.
L'astuce: pour calculer $\frac{\partial}{\partial x}f(x, y)$, il faut voir $y$ comme constant et donc voir $f$ comme fonction de $x$ uniquement. Par exemple:

$$\frac{\partial}{\partial x} \sin(x \pi) e^{-y^2 / 4}
= \pi \cos(x \pi) e^{-y^2 / 4}$$

Idem pour $y$:

$$\frac{\partial}{\partial y} \sin(x \pi) e^{-y^2 / 4}
= -\sin(x \pi) \frac{y}{2}  e^{-y^2 / 4}$$
"""

# ╔═╡ f09837e4-1c5a-4cfc-ba63-84a238d7f6c5
md"""
## Quand les dérivées sont-elles dans un plan ?

On dit dans ce cas que la fonction est *différentiable* en $(0, 0)$ et les dérivées directionelles sont donnée par produit scalaire avec le gradient.
Quels sont les conditions pour la différentiabilité ?

Voyons un exemple:
```math
\frac{x^{i_1} y^{i_2}}{x^{i_3} + y^{i_4}}
```
"""

# ╔═╡ c075846a-7e4c-4b99-ac91-3e83240fee0e
md"i\_1 = $(@bind i_1 Slider(0:6, default=2, show_value = true))"

# ╔═╡ 910d323f-583f-4e78-9374-2faa68b9e0d0
md"i\_2 = $(@bind i_2 Slider(0:6, default=2, show_value = true))"

# ╔═╡ 5f27741c-6b83-4b38-a452-caab047a2201
md"i\_3 = $(@bind i_3 Slider(2:2:6, default=2, show_value = true))"

# ╔═╡ 038a98c6-0ccf-461e-8dde-2a00ae834fe8
md"i\_4 = $(@bind i_4 Slider(2:2:6, default=2, show_value = true))"

# ╔═╡ 057b77ae-fff6-4314-bbc4-e34bf3d7d880
let
	x = y = range(-1, stop = 1, length = 100) # nombre pair de points pour ne pas évaluer en (0, 0)
	f(x, y) = x^i_1 * y^i_2 / (x^i_3 + y^i_4)
	plot.surface(; x, y, z = f.(x, y'))
end

# ╔═╡ a608649d-9f55-4337-a530-4f6ea2c27e16
md"Pour quels valeurs de `i_1`, `i_2`, `i_3`, `i_4` est ce que la fonction est continue en `(0, 0)` ? Vous y répondrez à la séance d'exercices! Ici, on se contentera de se fier au visuel."

# ╔═╡ 4ffc0380-b4da-44fd-9f63-d0f2eca4b28f
md"""
Il est nécessaire que la fonction soit continue en $(0, 0)$, donc par exemple $xy / (x^2 + y^2)$ n'est pas continue.
Y a-t-il des conditions sur $\partial f / \partial x$:
"""

# ╔═╡ 495f0638-138d-4551-b442-38cf956b2e8b
let
	x = y = range(-1, stop = 1, length = 300)
	function dfdx(x, y)
		d = -x^i_1 * y^i_2 / (x^i_3 + y^i_4)^2 * i_3 * x^(i_3 - 1)
		if i_1 > 0
			d += i_1 * x^(i_1 - 1) * y^i_2 / (x^i_3 + y^i_4)
		end
		return d
	end
	plot.surface(; x, y, z = dfdx.(x, y'))
end

# ╔═╡ 568618dc-fecb-43ba-a4a2-8284bb61b6c5
md"Et $\partial f / \partial y$:"

# ╔═╡ 6ac2de67-75d7-494a-93bf-072217d52317
let
	x = y = range(-1, stop = 1, length = 32)
	function dfdy(x, y)
		d = -x^i_1 * y^i_2 / (x^i_3 + y^i_4)^2 * i_4 * y^(i_4 - 1)
		if i_2 > 0
			d += i_2 * x^i_1 * y^(i_2 - 1) / (x^i_3 + y^i_4)
		end
		return d
	end
	plot.surface(; x, y, z = dfdy.(x, y'))
end

# ╔═╡ fd3a0cde-5ec0-4f3f-81b0-447d71971ab3
md"**Théorème**: Si $f$ est continue en $a$ et que les dérivées partielles de $f$ sont définies dans un voisinage de $a$ et continues en $a$ alors $f$ est différentiable."

# ╔═╡ 38c1fd0a-d148-47c2-bbe2-a23924073fe1
md"""
## La Hessienne

Pour obtenir la dérivée second, on dérive à nouveau le gradient.
On les notations:

$$\begin{align}
\frac{\partial}{\partial x} (\frac{\partial}{\partial x} f)
& = \frac{\partial^2}{\partial x^2} f\\
\frac{\partial}{\partial x} (\frac{\partial}{\partial y} f)
& = \frac{\partial^2}{\partial x \partial y} f\\
\frac{\partial}{\partial y} (\frac{\partial}{\partial x} f)
& = \frac{\partial^2}{\partial y \partial x} f\\
\frac{\partial}{\partial y} (\frac{\partial}{\partial y} f)
& = \frac{\partial^2}{\partial y^2} f\\
\end{align}$$

**Attention**: $\partial x^2$ ne veut **pas** dire qu'on dérive par $x^2$, c'est juste une notation pour dire qu'on dérive 2 fois par $x$.

Cette notation peut paraître mal choisie ! En effet, dans $\partial xy$ et $\partial yx$, les parties $xy$ et $yx$ ressemblent à des produits. Comme le produit est commutatif, on a à tendance à penser que dériver par $x$ puis $y$ est égal à dériver par $y$ puis $x$... À moins que ça soit vraiment égal ?

**Théorème** Si $f$ est différentiable alors
$$\frac{\partial^2}{\partial x\partial y} f = \frac{\partial^2}{\partial y\partial x} f$$.

On aime mettre ces 4 valeurs dans une matrice nommée *Hessienne*. Par ce théorème, la matrice est donc symmétrique!

$$\text{Hess}(f) = \begin{bmatrix}
  \frac{\partial^2}{\partial x^2} f & \frac{\partial^2}{\partial x \partial y} f\\
  \frac{\partial^2}{\partial x \partial y} f & \frac{\partial^2}{\partial y^2} f
\end{bmatrix}$$
"""

# ╔═╡ 942a8fe6-e911-41ae-a186-9cea58731eef
md"""
## Dérivée à plus de 2 variables

Tous les résultats se généralisent sans broncher à plus de 2 variables.
Considérons une fonction $f(x_1, x_2, \ldots, x_n)$ à $n$ variables.
Le gradient est
$$\nabla f = (\partial f/\partial x_1, \partial f/\partial x_2, \ldots, \partial f/\partial x_n)$$
La dérivée dans la direction $d = (d_1, d_2, \ldots, d_n)$ est:
$$\langle d, \nabla f \rangle
=
d^\top \nabla f
=
d_1 \cdot \partial f/\partial x_1 + \cdots + d_n \cdot \partial f/\partial x_n$$

La Hessienne est toujours symmétrique et devient:
$$\text{Hess}(f) =
\begin{bmatrix}
  \frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1 \partial x_2} &
  \frac{\partial^2 f}{\partial x_1 \partial x_3} & \cdots & \frac{\partial^2 f}{\partial x_1 \partial x_n}\\
  \frac{\partial^2 f}{\partial x_1 \partial x_2} & \frac{\partial^2 f}{\partial x_2^2} &
  \frac{\partial^2 f}{\partial x_2 \partial x_3} & \cdots & \frac{\partial^2 f}{\partial x_2 \partial x_n}\\
  \vdots & \ddots & \ddots & \ddots & \vdots\\
  \frac{\partial^2 f}{\partial x_1 \partial x_{n-1}} & \frac{\partial^2 f}{\partial x_2 \partial x_{n-1}} &
  \cdots & \frac{\partial^2 f}{\partial x_{n-1}^2} & \frac{\partial^2 f}{\partial x_{n-1} \partial x_n}\\
  \frac{\partial^2 f}{\partial x_1 \partial x_n} & \frac{\partial^2 f}{\partial x_2x_n} &
  \cdots & \frac{\partial^2 f}{\partial x_{n-1} \partial x_n} & \frac{\partial^2 f}{\partial x_n^2}\\
\end{bmatrix}$$
"""

# ╔═╡ 9bb14f1d-b7b5-43c1-a931-c5018f544806
md"## Utils"

# ╔═╡ 73248509-2b7b-4d79-9bc1-4689d843b73e
function tangent_plane(f, x0, x1, y0, y1, cx, cy; Δ = 0.5, n = 32)
	x = range(x0, stop = x1, length = n)
	y = range(y0, stop = y1, length = n)
	df(x, y) = ForwardDiff.gradient(xy -> f(xy...), [x, y])
	surface(x, y, f, label = "f(x)", zlims = (minimum(f.(x, y') + maximum.(df.(x, y'))), maximum(f.(x, y') + maximum.(df.(x, y')))), legend=:bottomleft)
	grad = df(cx, cy)
	fc = f(cx, cy)
	plane(x, y) = fc + grad' * [x - cx, y - cy]
	plot!([cx, cx + 1], [cy, cy], [fc, fc + grad[1]], linewidth = 10, color = Colors.JULIA_LOGO_COLORS.green, label = "")
	plot!([cx, cx], [cy, cy + 1], [fc, fc + grad[2]], linewidth = 10, color = Colors.JULIA_LOGO_COLORS.green, label = "")
	return plot!([cx, cx+1], [cy, cy + 1], [fc, fc + grad[1] + grad[2]], linewidth = 10, color = Colors.JULIA_LOGO_COLORS.red, label = "")
	surface!( # works with plotly but not GR
		range(cx - Δ, stop = cx + Δ, length = 3),
		range(cy - Δ, stop = cy + Δ, length = 3),
		plane,
		color = [Colors.JULIA_LOGO_COLORS.blue],
		label = "",
	)
end

# ╔═╡ fd74a3bd-e20c-43dd-badc-b1e861e69660
tangent_plane((x, y) -> sin(x * π) * exp(-y^2 / 4), x0, x1, y0, y1, cx, cy)

# ╔═╡ bf46ec06-54d8-41ba-a6b8-1ae1ac289556
function tangent_line(f, x0, x1, c; n = 100)
	x = range(x0, stop = x1, length = n)
	df(x) = ForwardDiff.derivative(f, x)
	pente = df(c)
	tangente(x) = pente * (x - c) + f(c)
	x_right = range(c, stop = c + 1, length = 2)
	plot.plot(; x, y = f.(x), color = Colors.JULIA_LOGO_COLORS.blue, label = "f(x)", ylims = (minimum(f.(x) + df.(x)), maximum(f.(x) + df.(x))), legend=:bottomleft)
	#plot(; x, y = df.(x), color = Colors.JULIA_LOGO_COLORS.purple, label = "f'(x)").
	#plot(; x = [c, c + 1], y = [tangente(x), tangent(c + 1)], color = Colors.JULIA_LOGO_COLORS.green, label = "").
	#plot(; x = [c, c + 1], y = [f(c), f(c)], color = Colors.JULIA_LOGO_COLORS.green, label = "").
	#plot(; x = [c + 1, c + 1], y = [f(c), tangente(c + 1)], color = Colors.JULIA_LOGO_COLORS.red, label = "").
	#plot.scatter(; x = [c], y = [pente], markerstrokewidth = 0, color = Colors.JULIA_LOGO_COLORS.red, label = "")
end

# ╔═╡ cdf52f5f-d541-48e7-99ca-b5d0ff2c2d17
tangent_line(x -> sin(x * π) * exp(-x^2 / 10), x0, x1, cx)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotlyLight = "ca7969ec-10b3-423e-8d99-40f33abb42bf"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Colors = "~0.13.1"
ForwardDiff = "~1.3.0"
PlotlyLight = "~0.13.1"
PlutoUI = "~0.7.75"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.2"
manifest_format = "2.0"
project_hash = "5eac52e911ef4250599510751f969203542ef6f7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Cobweb]]
deps = ["DefaultApplication", "OrderedCollections", "Scratch"]
git-tree-sha1 = "6665ec6b16446379fb76ad58a2a7b65687c77271"
uuid = "ec354790-cf28-43e8-bb59-b484409b7bad"
version = "0.7.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EasyConfig]]
deps = ["JSON3", "OrderedCollections", "StructTypes"]
git-tree-sha1 = "11fa8ecd53631b01a2af60e16795f8b4731eb391"
uuid = "acab07b0-f158-46d4-8913-50acef6d41fe"
version = "0.1.16"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cd33c7538e68650bd0ddbb3f5bd50a4a0fa95b50"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.0"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyLight]]
deps = ["Artifacts", "Cobweb", "Dates", "Downloads", "EasyConfig", "JSON3", "REPL"]
git-tree-sha1 = "ed95b3125e681e5209221a035cfc46058cfcd88f"
uuid = "ca7969ec-10b3-423e-8d99-40f33abb42bf"
version = "0.13.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "db8a06ef983af758d285665a0398703eb5bc1d66"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.75"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─812d53c7-1384-4c77-9185-cf3d6656b98b
# ╠═cdf52f5f-d541-48e7-99ca-b5d0ff2c2d17
# ╟─b3c01426-9a70-448b-bf3e-901e3100082f
# ╟─674be307-64c1-4485-ad25-a59ff6840df8
# ╠═fd74a3bd-e20c-43dd-badc-b1e861e69660
# ╟─0356921d-c3e5-4351-9c14-0efd657be7dd
# ╟─b27a8a10-c3fa-44c6-917d-14959f4f30dd
# ╟─f09837e4-1c5a-4cfc-ba63-84a238d7f6c5
# ╟─c075846a-7e4c-4b99-ac91-3e83240fee0e
# ╟─910d323f-583f-4e78-9374-2faa68b9e0d0
# ╠═057b77ae-fff6-4314-bbc4-e34bf3d7d880
# ╟─5f27741c-6b83-4b38-a452-caab047a2201
# ╟─038a98c6-0ccf-461e-8dde-2a00ae834fe8
# ╟─a608649d-9f55-4337-a530-4f6ea2c27e16
# ╟─4ffc0380-b4da-44fd-9f63-d0f2eca4b28f
# ╠═495f0638-138d-4551-b442-38cf956b2e8b
# ╟─568618dc-fecb-43ba-a4a2-8284bb61b6c5
# ╠═6ac2de67-75d7-494a-93bf-072217d52317
# ╟─fd3a0cde-5ec0-4f3f-81b0-447d71971ab3
# ╟─38c1fd0a-d148-47c2-bbe2-a23924073fe1
# ╟─942a8fe6-e911-41ae-a186-9cea58731eef
# ╟─9bb14f1d-b7b5-43c1-a931-c5018f544806
# ╠═73248509-2b7b-4d79-9bc1-4689d843b73e
# ╠═d776d4d6-825b-11ef-3d58-85ab2150ff5f
# ╠═e4a1845e-0d27-4f0c-abbb-11960b8a2452
# ╠═bf46ec06-54d8-41ba-a6b8-1ae1ac289556
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
