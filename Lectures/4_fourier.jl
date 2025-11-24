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

# ╔═╡ 583c3b75-71cf-48f6-b93d-2201c6c7d520
using Plots

# ╔═╡ 34c779d2-5a9d-4c98-ac28-e4be1a85da40
using PlutoUI, FFTW, PortAudio, Images, TestImages, ImageFiltering

# ╔═╡ e3a7dd72-f106-4865-8b24-d86056782160
using Polynomials, LinearAlgebra

# ╔═╡ b8949a96-8d97-4048-8022-ab8fabf326a1
md"# À la fréquence de Fourier"

# ╔═╡ 4309cf1d-2a7a-4a84-996b-ecc99d62dad1
md"## 1.1 La Transformée de Fourier"

# ╔═╡ df7fa000-70ce-462b-af19-fb9c4b022012
md"""
### À quoi ça sert ?

Un monde parallèle ou les choses deviennent plus simples (ou plus compliquées)

* **1.2.1 Traitement du signal** Séparation des contributions de chaque fréquence
* **1.2.2 Avantage computationel** L'opération de convolution devient une multiplication classique
"""

# ╔═╡ 58d25489-0d99-4c33-ad22-5ffb2751c3c6
md"### 1.1.1 Illustration : séparation des fréquences"

# ╔═╡ 5010b40e-86c1-4ed1-b7ae-25086ad45205
md"#### Audio"

# ╔═╡ 27d19b85-8f5b-4d53-93c2-033cca3b9006
@bind la_clicked CounterButton("play la")

# ╔═╡ 41333814-27a7-4ecf-99ce-4f972c0a37da
@bind ré_clicked Button("play ré")

# ╔═╡ 308757bd-bfbc-43d5-8e72-5903390c1458
@bind la_ré_clicked Button("play la + ré")

# ╔═╡ 8f77bda0-25b0-44fb-b235-afb6c2e725ce
md"#### Signal temporel"

# ╔═╡ 54e54a58-3fdf-43d8-a87a-22ff047a392e
md"Fréquence d'échantillonnage = $(@bind échantillonnage Slider((2).^(4:13), default=1024, show_value = true))"

# ╔═╡ 83f8f109-50fa-4174-8fcd-6bd9d8ef2d65
md"Un son est un signal continu mais on va l'[échantillonner](https://fr.wikipedia.org/wiki/%C3%89chantillonnage_(signal)). C'est à dire qu'on va prendre un nombre fini de valeur par seconde à des distance égale dans le temps. Dans cet exemple, on prend une valeur toute les $(échantillonnage)ᵉ de seconde."

# ╔═╡ 83c1cbf7-e24d-4e9a-a309-c7b60e730344
md"Ça correspond à une distance en seconde entre deux échantillons de:"

# ╔═╡ 44dbdbff-18a2-43ab-b8e5-399bf8fe9639
Δt = 1 / échantillonnage

# ╔═╡ 590c895f-d7d4-467a-90cd-90bfba59c2b5
md"nombre\_échantillons = $(@bind nombre_échantillons Slider((2).^(4:13), default=1024, show_value = true))"

# ╔═╡ 40490a0f-9a06-4283-9d3a-d2c96b6b14df
temps_total = nombre_échantillons * Δt

# ╔═╡ 41155d3b-ba20-49ee-af03-babd4ce10c68
md"En prenant $(nombre_échantillons), le signal dure $(temps_total) secondes. Comme ces nombres sont équidistants (distance de `Δt` secondes entre eux), on peut représenter cette suite de façon compact comme suit:"

# ╔═╡ f18256be-d172-4f10-9c6b-662015f535e4
temps = range(Δt, stop=temps_total, length=nombre_échantillons)

# ╔═╡ b790ddf3-37ac-45d6-87aa-ebd06c829d82
md"Le *la* est une note de musique de [fréquence 440 Hz](https://fr.wikipedia.org/wiki/La_440)"

# ╔═╡ 02084ab9-305a-4317-9fad-a4acf8cc5cd3
la = cispi.(2*440*temps)

# ╔═╡ 4a11717b-4086-44f3-9a50-b531d51c6944
md"La note de musique *ré* a une [fréquence de 293.7 Hz](https://fr.wikipedia.org/wiki/Musique_occidentale)"

# ╔═╡ 0ba2c1bc-aa82-4fbd-a704-034274349c4d
ré = cispi.(2*293.7*temps)

# ╔═╡ ef156c25-cb7d-4cd3-8412-d9f5f2bc7a3b
md"temps\_zoom = $(@bind temps_zoom Slider(range(Δt, stop=1000Δt, length=32), default=10 / 440, show_value = true))"

# ╔═╡ 4265ccd9-3bf2-4655-a0da-2d774baf32ac
begin
	index_max = something(findfirst(t -> t >= temps_zoom, temps), length(temps))
	plot(temps[1:index_max], real.(la[1:index_max]), label = "la")
	plot!(temps[1:index_max], real.(ré[1:index_max]), label = "ré")
	plot!(temps[1:index_max], real.(la[1:index_max] + ré[1:index_max]), label = "la + ré")
end

# ╔═╡ c26a3bd3-cf79-426a-b439-2b29501cee4b
md"#### Signal fréquentiel"

# ╔═╡ a315cbaf-a4b5-42ea-a0e0-f0b064b9995e
@bind freq_sel Select(["abs", "real", "imag"])

# ╔═╡ bb5eebce-e385-4c4f-b827-0ffcdead0799
let
	f = range(0, stop = échantillonnage - 1/temps_total, length=length(la))
	sel = Dict(
		"abs" => abs,
		"real" => real,
		"imag" => imag,
	)[freq_sel]
	plot(f, sel.(fft(la)), label = "la", linewidth=5)
	plot!(f, sel.(fft(ré)), label = "ré", linewidth=5, color = :green)
	plot!(f, sel.(fft(la + ré)), label = "la + ré", linewidth=2, color = :red)
end

# ╔═╡ 2488c83e-5e54-4d47-80d2-6ef41ca7ac7f
md"On a la valeur de la transformée de fourier, tous les $(1/temps_total) Hz"

# ╔═╡ 279e5f5e-4fc6-475f-8dab-c5740a470b0e
md"""
### 1.1.2 Calcul rapide de convolution

La convolution continue:
```math
(f \ast g)(x) = \int_{-\infty}^\infty f(s) g(x - s) \text{d}s
```
La convolution discrète:
```math
(f \ast g)_n = \sum_{k=-\infty}^\infty f_k g_{n - k}
```

**Propriétés algébriques:**
```math
\begin{align}
  \text{commutatif} && f \ast g & = g \ast f\\
  \text{associatif} && (f \ast g) \ast h & = f \ast (g \ast h)\\
  \text{distributif} && f \ast (g + h) & = f \ast g + f \ast h\\
  \text{conjugué} && \overline{f \ast g} & = \overline{f} \ast \overline{g}\\
\end{align}
```

La cross-corrélation:
```math
\begin{align}
  (f \star g)(x) & = \int_{-\infty}^\infty f(s) g(x + s) \text{d}s &
  (f \star g)_n & = \sum_{-\infty}^\infty f_k g_{n + k} \text{d}s
\end{align}
```

**Convolution theorem**: Notons la tranformée de Fourier de ``f`` par ``\mathcal{F}(f)``. Par la transformée de Fourier, la convolution devient un produit classique:
```math
f \ast g = \mathcal{F}^{-1}(\mathcal{F}(f) \cdot \mathcal{F}(g))
```
"""

# ╔═╡ 7cdf9587-af9f-4346-b619-ba8bdf6d6799
md"#### 1.1.2.1 Illustration: convolution d'images"

# ╔═╡ d4ac67de-9442-464c-be8b-e8ff468dd996
mandrill = testimage("mandrill")

# ╔═╡ 54ced290-6da4-4f24-a18a-e46ac116448c
md"##### Gaussian Kernel"

# ╔═╡ 93cfa0b9-dbb6-4956-bd94-575fe233c9f4
md"La fonction `gaussian(d)` donne une matrix de taille ``(2d + 1) \times (2d + 1)`` contenant la valeur d'une Gaussienne."

# ╔═╡ 5dd93a6b-f355-4927-8f3d-a3522f43934f
md"d = $(@bind d Slider(1:10, default=5, show_value = true))"

# ╔═╡ 3855c981-206d-46e6-8a9b-eda011769208
surface(-2d:2d, -2d:2d, Kernel.gaussian(d))

# ╔═╡ 53045a43-9941-4b3f-a0be-d6de969ee124
md"La fonction `imfilter` calcule la cross-corrélation."

# ╔═╡ 8f9cccf1-ec8d-448d-a86f-3955f6c238d8
imfilter(mandrill, Kernel.gaussian(d))

# ╔═╡ 5edc270e-a7c6-4e44-97b0-7697adca16a2
md"La convolution est obtenue avec `reflect`."

# ╔═╡ e1c25e06-3259-4a56-b26a-60903e57f92b
imfilter(mandrill, reflect(Kernel.gaussian(d)))

# ╔═╡ 802cbd16-2126-4d4f-aca2-df29024d5fdf
md"##### Performance"

# ╔═╡ 9fae43f1-61cd-43c6-b383-d9bbd4c2a1ac
function imfilter_time(n, alg)
	A = rand(n, n)
	B = rand(n, n)
	return @elapsed imfilter(A, B, alg)
end

# ╔═╡ 31aca097-4a85-4d63-b114-6d1828648e87
let
	n = 2 .^ (1:6)
	fir_time(n) = imfilter_time(n, Algorithm.FIR())
	fft_time(n) = imfilter_time(n, Algorithm.FFT())
    plot(n, fir_time ∘ Int, label = "Finite Impulse Response (FIR)")
	plot!(n, fft_time ∘ Int, label = "Fast Fourier Transform (FFT)")
end

# ╔═╡ e7b74b94-bf5f-4455-9fc7-922ce7d4c284
convert(AbstractArray, Kernel.Laplacian())

# ╔═╡ 101663c2-acca-4eb7-858c-20a70c2103e9
md"##### Autre kernels"

# ╔═╡ ae1e22a7-b4af-4f83-8f28-97534d7e0c20
imfilter(mandrill, Kernel.Laplacian())

# ╔═╡ eeafd828-82a7-449e-b083-f22eb0d07079
md"Différents kernels permettent d'analyser des aspects différents d'une image. Le kernel à utilisé peut aussi être **appris**, c'est la base des **Convolutional Neural Networks** (CNNs) !"

# ╔═╡ 66fd9f77-6c2a-418d-9650-39033396b584
md"#### 1.1.2.2 Illustration : Le produit de polynômes"

# ╔═╡ 0230ce41-bf49-4363-b344-3c475ef6474d
deg = 4

# ╔═╡ 680c44ba-02b9-41d1-b85f-86f49a85ca7a
p = Polynomial(rand(deg))

# ╔═╡ 14a60621-f76e-48f6-8dbe-cf4c871b2027
q = Polynomial(rand(deg))

# ╔═╡ 36be0072-a89e-46e8-84bc-57e208a9e70a
p * q

# ╔═╡ c312e94d-a461-4d8a-bf56-99699e82fe04
Polynomials.poly_multiplication_fft(p, q)

# ╔═╡ 4580177c-9654-4950-bbf5-a599adb44db6
function poly_time(n, with_fft::Bool)
	p = Polynomial(rand(n))
	q = Polynomial(rand(n))
	@elapsed if with_fft
		Polynomials.poly_multiplication_fft(p, q)
	else
		p * q
	end
end

# ╔═╡ 9500f893-7eb0-4f24-83eb-5a2240279062
let
	n = 2 .^ (1:13)
	std_time(n) = poly_time(n, false)
	fft_time(n) = poly_time(n, true)
    plot(n, std_time ∘ Int, label = "Standard multiplication")
	plot!(n, fft_time ∘ Int, label = "With FFT")
end

# ╔═╡ 44f73082-7cdb-47b9-bb86-f952ac7c4498
md"""
### 1.1.2 Définition de la transformée de Fourier

La transformée de Fourier continue et son inverse:
```math
\begin{align}
  F(\xi) & = \int_{-\infty}^\infty f(t) e^{-i 2\pi\xi t} \text{d}t\\
  f(t) & = \int_{-\infty}^\infty F(\xi) e^{i 2\pi\xi t} \text{d}\xi
\end{align}
```
La transformée de Fourier discrète et son inverse:
```math
\begin{align}
  X_k & = \sum_{n=0}^{N-1} x_n e^{-i 2\pi \frac{k}{N} n}\\
  x_n & = \frac{1}{N}\sum_{k=0}^{N-1} X_k e^{i 2\pi \frac{k}{N} n}
\end{align}
```
"""

# ╔═╡ a2359c44-3150-4b0a-8675-90d1627b5c07
md"""
**Intuition**
Le signal pur "la" est
``f_{\text{la}}(t) = \exp(2\pi 440 t i) = \cos(2\pi 440 t) + i \sin(2\pi 440 t)``.
Si j'évalue en ``t = 1/440 s``, j'obtiens la fin d'une période:
``\exp(2\pi i) = \cos(2\pi) + i \sin(2\pi)``.
Ces exponentielles pures sont *orthogonales*, c'est à dire que ces intégrales entre fréquences différentes sont nulles.
```math
\begin{align}
  F(440) & = \int_{-\infty}^\infty f(t) e^{-i 2\pi 440 t} \text{d}t = \int_{-\infty}^\infty 1 \text{d}t = \infty\\
  F(293.7) & = \int_{-\infty}^\infty f(t) e^{-i 2\pi 293.7 t} \text{d}t = 0
\end{align}
```
Cet intégrale permet donc de détecter la quantité d'une certaine fréquence dans un signal, en ignorant toutes les autres fréquences!
"""

# ╔═╡ 86368d24-c206-42dc-afdb-1e762b3ac14c
md"Rajouter linearité : TODO"

# ╔═╡ f0eda0e0-0a52-403a-a73d-3a604cbe557c
md"TODO : même dimension en temporel et fréquentiel"

# ╔═╡ d29773c7-a413-46a9-8c62-2f54d4de1730
md"""
### 1.1.3 Le replis spectral

On vient de voir que le signal est extrapolé périodiquement par la transformée de fourier discrète.
Si on prend un signal de `0` à `t_max` et qu'on divise la fréquence d'échantillonnage de sa transformée de fourier par 2.
Ça signifie que le signal sera 2 fois moins long.
"""

# ╔═╡ 49d078e2-d865-4f2f-9cfc-3c2d5a535088
md"Illustrons cela avec la fonction `exp(-(t-centre)^2)` qui est initiallement échantillonnée avec `512` échantillons de `0` à `t_max` et qu'on ré-échantillonne ensuite à `512 / ratio` en fréquentiel, ce qui divise `t_max` par `ratio` également."

# ╔═╡ b80ff616-8bb7-4132-88c7-ed10c3c6786c
md"t\_max = $(@bind t_max Slider(range(1, stop=10, length=10), default=5, show_value = true))"

# ╔═╡ 56ff3957-0295-41fe-839f-12dae0fdd14f
md"centre = $(@bind centre Slider(range(0.1, stop=t_max, length=32), default=t_max / 2, show_value = true))"

# ╔═╡ 4980c4ff-7da0-4747-b779-891569d604fd
md"ratio = $(@bind ratio Slider(2 .^ (0:4), default=2, show_value = true))"

# ╔═╡ 5032f1f1-47a1-4ff3-b3fd-49d60e5e28a5
let
	N = 512
	d = ratio
	x = range(0, stop = t_max - t_max / N, length = N)
	gauss(x) = exp(-16(x-centre)^2)/4 + exp(-16(x-centre/2)^2)
	plot(x, gauss, label = "signal original")
	X = fft(gauss.(x))
	Xs = X[1:d:N]
	half = div(N, 2d)
	xs = real.(ifft(Xs))
	plot!(x[1:length(xs)], xs, label = "rééchantilloné")
end

# ╔═╡ 78026d88-7051-11ef-29f0-67f85176a548
md"## 1.2 Fast Fourier Transform"

# ╔═╡ dc4c4d53-feb2-40ce-b20e-3386aab2a45f
md"""
Attardons-nous à présent sur l'algorithme utilisé pour calculé la DFT (transformée de fourier discrète).
On a vu précédemment que cet algorithme a une complexité de ``\Omega(n \log(n))``.
C'est parmis le [top 10 des meilleurs algorithmes du 20ᵉ siècle](https://archive.siam.org/pdf/news/637.pdf) !
Son Implémentation la plus rapide est dans la libraire [FFTW](https://www.fftw.org/) : the Fastest Fourier Transform in the West ! Comment est-elle la plus rapide ? Voir [ici](https://youtu.be/mSgXWpvQEHE?t=1800).
L'algorithme est initialement découvert part Gauss en 1805 puis redécouvert par Cooley et Tukey en 1965.
"""

# ╔═╡ 5f9eeb76-0a9d-40ad-858e-4ab67af46427
md"""
La DFT est en fait une évaluation d'un polynômes aux différentes racines `N`ièmes de l'unité:
```math
\begin{align}
  X_k & = \sum_{n=0}^{N-1} x_n e^{-i 2\pi \frac{k}{N} n}\\
  X_k & = \sum_{n=0}^{N-1} x_n z_{k,N}^{n} \qquad z_{k,N} = e^{-i 2\pi \frac{k}{N}}
\end{align}
```
La DFT peut aussi être vue comme un produit matriciel:
```math
  X =
  \begin{bmatrix}
    1 & 1 & \cdots & 1 & 1\\
    1 & z_N & \cdots & z_N^{N-2} & z_N^{N-1}\\
    1 & z_N^2 & \cdots & z_N^{2(N-2)} & z_N^{2(N-1)}\\
	\vdots & \vdots & \ddots & \vdots & \vdots\\
    1 & z_N^{(N-1)2} & \cdots & z_N^{(N-1)(N-2)} & z_N^{(N-1)(N-1)}
  \end{bmatrix}
  x
```
Un produit matriciel a complexité ``\Omega(n^2)``. Ça ne nous donne pas une complexité de ``\Omega(n \log n)``, il va falloir utiliser la structure assez spéciale de cette matrice...

On a acquis une intuition géométrique sur les racines ``z_N``. On va s'en servir pour visualiser cette matrice.
"""

# ╔═╡ a33ebf12-5c62-4cc9-af1e-7bad9026db21
md"Besoin d'un indice ? $(@bind colored_arrows CheckBox())"

# ╔═╡ c64e6062-3793-4556-9001-bf184cd090e9
md"On peut le vérifier numériquement également:"

# ╔═╡ f952619b-194a-4088-a2bb-b8210392417a
isapprox(0.0, 1e-300, atol = 1e-8)

# ╔═╡ 6119e719-4e14-4516-ac7e-a236cfb4fcfa
md"Vérifions le à nouveau numériquement:"

# ╔═╡ 6d79251e-f7d5-4126-a781-7221a4313c9a
cispi(-2/8)

# ╔═╡ 59bc8221-938c-4919-8b94-7a2d61cf2062
cispi(-2/8) .^ collect(0:3)

# ╔═╡ f53db6db-d9f9-4371-a551-e3bfc2d6b162
md"En se rappelant que l'intuition géométrique des `N` racines de l'unité, on comprend pourquoi les `N/2` racines sont un sous-ensembles des `N` racines."

# ╔═╡ bcf5ecbc-ccb6-48ce-9146-55b641f770f2
md"`max_N` = $(@bind max_N Slider(2 .^ (2:5), default=8, show_value = true))"

# ╔═╡ 59130f1a-4fcd-4ae1-9ecf-1818dfc07612
md"## Utilitaires"

# ╔═╡ 1817a169-0afc-4fff-be23-1ea510ff337a
function play_sound(sound, samplerate)
	# For the static version, there is no sound so we just make the buttons
	# do nothing
	if PortAudio.Pa_GetDefaultOutputDevice() != PortAudio.LibPortAudio.paNoDevice
		PortAudioStream(0, 2; samplerate) do stream
    		write(stream, real.(sound))
		end
	end
end

# ╔═╡ 04ac4af9-5789-4d0a-a958-1630d0a54299
# ╠═╡ disabled = true
#=╠═╡
let
	la_clicked

	échantillonnage = 2^13
	Δt = 1 / échantillonnage
	temps = range(Δt, stop=2, length=échantillonnage)
	la = cispi.(2*440*temps)
	play_sound(la, échantillonnage)
	nothing
end
  ╠═╡ =#

# ╔═╡ 5d64d5ac-ff4f-415a-85bb-c673da00e9a1
let
	ré_clicked

	échantillonnage = 2^13
	Δt = 1 / échantillonnage
	temps = range(Δt, stop=2, length=échantillonnage)
	ré = cispi.(2*293.7*temps)
	play_sound(ré, échantillonnage)
	nothing
end

# ╔═╡ 93534f3d-f9da-46c4-a9c3-d055530f3045
let
	la_ré_clicked
	
	échantillonnage = 2^13
	Δt = 1 / échantillonnage
	temps = range(Δt, stop=2, length=échantillonnage)
	la = cispi.(2*440*temps)
	ré = cispi.(2*293.7*temps)
	play_sound(la + ré, échantillonnage)
	nothing
end

# ╔═╡ 0f87464a-1a55-4f7c-8a88-bf514f5177a0
begin
    struct Join
		list
	    Join(a) = new(a)
	    Join(a, b, args...) = Join(tuple(a, b, args...))
    end
	function Base.show(io::IO, mime::MIME"text/html", d::Join)
		for el in d.list
			show(io, mime, el)
		end
	end
end

# ╔═╡ 9e6732a7-d2df-465f-a7fa-c119a64deb9a
begin
	struct HTMLTag
		tag::String
		parent
	end
	function Base.show(io::IO, mime::MIME"text/html", d::HTMLTag)
		write(io, "<", d.tag, ">")
		show(io, mime, d.parent)
		write(io, "</", d.tag, ">")
	end
end

# ╔═╡ 0e9e8c6e-83ff-4c4c-b9d1-19a0593a3cba
f(i, j, n) = cispi(-2 * i * j / n)

# ╔═╡ 5e5ff9b1-0a78-4e9c-9c90-4fb80b0dbe21
function F(n)
	A = [
		f(k, m, n)
		for k in 0:(n-1), m in 0:(n-1)
	]
end

# ╔═╡ 222cfba3-88f7-4824-92ca-f9ab81a000c4
F(8)

# ╔═╡ c9131a2a-8978-4c81-8588-44e4d7071c0b
partie_mauve = F(8)[1:4, 1:2:8]

# ╔═╡ f24f49d9-10e6-4ff5-ace7-8767c69f8594
partie_bleue = F(8)[5:8, 1:2:8]

# ╔═╡ ff3ec271-5e08-4ab8-8c13-7f5e1996e383
partie_mauve == partie_bleue

# ╔═╡ aa4ff136-6739-46dc-aae6-cebcb3b31263
partie_bleue

# ╔═╡ 8676e10b-6c41-4e8f-a56b-ca548b9690dd
partie_verte = F(8)[1:4, 2:2:8]

# ╔═╡ 4ead6b3a-af8d-4dc6-905e-bfc4b6a00e64
partie_verte

# ╔═╡ faf71e8e-7ba0-44fa-800d-949a42fbacec
partie_rouge = F(8)[5:8, 2:2:8]

# ╔═╡ fd5ef1e7-bc1e-4382-b5c6-e960ecb3c893
partie_verte == -partie_rouge

# ╔═╡ 3567a3cb-a5be-4d7c-a917-d3947b5549f6
F(4)

# ╔═╡ 0c74b1e8-a67f-4e41-9638-bb10babf7138
partie_bleue == F(4)

# ╔═╡ 96c678be-8528-462f-b455-446b9a367fab
cispi(-2/8) .^ (0:3) .* F(4)

# ╔═╡ 9a7a1c48-8168-4a21-b57e-0aaf2dbb9de8
isapprox.(real.(partie_verte), real.(cispi(-2/8) .^ (0:3) .* F(4)), atol=1e-8)

# ╔═╡ 0d9e0a31-b229-4d06-92c9-0136404399da
F(8)

# ╔═╡ f80595c3-6bbb-4b62-b03b-1c80cecaf7f9
function arrows(n, use_color = false)
	plot(showaxis = false)
	for i in 0:(n-1)
		for j in 0:(n-1)
			if use_color
				color = Colors.JULIA_LOGO_COLORS[1 + iseven(j) * 2 + (2i < n)]
			else
				color = :black
			end
			scatter!([j], [-i], markersize = 8 / sqrt(n), markerstrokewidth = 0, legend = nothing; color)
			c = f(i, j, n) / 3 # Divide by 3 so that arrow is not too long in plot
			plot!([j, j + real(c)], [-i, -i + imag(c)], arrow = true, legend = nothing; color)
		end
	end
	plot!()
end

# ╔═╡ cd3b2930-2c33-4028-9b27-d25238c79bac
arrows(8, colored_arrows)

# ╔═╡ 3739ddf4-ecb8-49d5-a70c-aa347c6fb17d
arrows(4, false)

# ╔═╡ ecf827b8-6a01-4e3b-9856-451239a7b175
arrows(8, true)

# ╔═╡ 792723a3-0011-422e-bbf4-798624913780
function draw_roots!(N; kws...)
	for i in 1:N
		b, a = sincospi(2 * i / N)
		plot!([0, a], [0, b]; arrow = true, label = nothing, kws...)
	end
end

# ╔═╡ 43c13895-a651-4e61-8fda-abdde42718da
let
	p = plot(ratio = :equal, ticks = nothing, showaxis = false)
	max_pow = round(Int, log2(max_N))
	colors = []
	for pow in max_pow:-1:2
		draw_roots!(2^pow, linewidth = 2 * (pow - 1), color = Colors.JULIA_LOGO_COLORS[pow - 1])
	end
	p
end

# ╔═╡ b0a265b8-e3fd-4a6c-b813-904d18ea83e7
function qa(question, answer)
	return HTMLTag("details", Join(HTMLTag("summary", question), answer))
end

# ╔═╡ 1dfad36a-05a5-4108-8c67-d4c866be89c0
qa(
	html"Quel est la complexité de la convolution vs produit classique s'ils sont discrets de longueur n ?",
	md"Convolution: ``\Theta(n^2)``, produit classique: ``\Theta(n)``. On verra que la complexité de la FFT est ``\Theta(n \log(n))`` donc ``f \ast g`` a une complexité de ``\Theta(n^2)`` et ``\mathcal{F}^{-1}(\mathcal{F}(f) \cdot \mathcal{F}(g))`` a une complexité de ``\Theta(n + 3n \log(n)) = \Theta(n\log(n))`` ce qui est plus rapide pour un nombre ``n`` suffisamment grand.",
)

# ╔═╡ 13d00d91-dae3-451c-924a-6b120c0a605d
qa(
	html"Pourquoi obtient-on le même résultat avec la convolution et la cross-corrélation ?",
	md"Parce que la gaussienne est symmétrique !",
)

# ╔═╡ 6b9a064c-1097-4754-af85-5dcaa435b893
qa(md"Quelle est la complexité de FIR ? et de de FFT ?", md"""
**FIR**: ``(f \ast g)_{k,l} = \sum_{i=1}^n\sum_{j=1}^n f_{i,j} g_{k - i, l -j}``
Pour une pair ``k,l``, c'est ``n^2`` opérations.
Comme il y a ``n`` valeurs de ``k`` et idem pour ``l``, je dois faire cette double somme ``n^2`` fois.
Donc la complexité totale c'est ``n^4``.

**FFT**: ``\mathcal{F}^{-1}(\mathcal{F}(f) \cdot \mathcal{F}(g))``
* La complexité de ``\mathcal{F}(f)`` est ``n^2\log(n^2)``
* La complexité de ``\mathcal{F}(g)`` est ``n^2\log(n^2)``
* La complexité du produit ``\mathcal{F}(f) \cdot \mathcal{F}(g)`` est ``n^2``.
* La complexité de la FT inverse ``\mathcal{F}^{-1}(\ldots)`` est ``n^2\log(n^2)``
Au total, ``3n^2\log(n^2) + n^2``. Notons que ``\log(n^2) = 2\log(n)``.
En ignorant les constantes, on a donc ``n^2\log(n) + n^2``.
""")

# ╔═╡ 976ac813-8f46-4912-a04a-ec586155867d
qa(md"Quelle est la complexité du produit de polynômes avec ou sans la FFT ?", md"""
**Sans FFT**: ``(f \ast g)_{k} = \sum_{i=1}^n f_{i} g_{l - i}``
Pour un ``k``, c'est ``n`` opérations.
Comme il y a ``n`` valeurs de ``k``, je dois faire cette somme ``n`` fois.
Donc la complexité totale c'est ``n^2``.

**FFT**: ``\mathcal{F}^{-1}(\mathcal{F}(f) \cdot \mathcal{F}(g))``
* La complexité de ``\mathcal{F}(f)`` est ``n\log(n)``
* La complexité de ``\mathcal{F}(g)`` est ``n\log(n)``
* La complexité du produit ``\mathcal{F}(f) \cdot \mathcal{F}(g)`` est ``n``.
* La complexité de la FT inverse ``\mathcal{F}^{-1}(\ldots)`` est ``n\log(n``
Au total, ``3n\log(n) + n``.
En ignorant les constantes, on a donc ``n\log(n) + n``.
""")

# ╔═╡ 87dbf272-1d03-4f9e-b7e1-9b76197b0a83
qa(
	html"Comment utiliser cela pour calculer le produit de grand nombres ?",
	md"Un nombre représenté par les chiffres ``a_na_{n-1}\cdots a_1a_0`` est égal à ``a_n 10^n + a_{n-1} 10^{n-1} + \cdots + a_1 10 + a_0`` ou encore ``p(10)`` où ``p(x) = a_n x^n + a_{n-1} x^{n-1} + \cdots + a_1 x + a_0``. Le produit de nombre peut donc être calculé via un produit de polynôme!",
)

# ╔═╡ 3524f1bf-8e39-4d04-9dc3-e434e65db409
qa(
	html"Que vaut la distance entre deux valeurs successive en temporel et en fréquentiel pour les signaux discrets ?",
	md"""
L'échelle temporelle ou fréquentielle n'est pas contenue dans ``X_k`` ou ``x_n``.
Par contre, elles sont étroitement liées. On a
```math
\frac{k}{N}n = \frac{k}{N\Delta t} n\Delta t
```
Si on décide que l'écart temporel entre ``x_n`` et ``x_{n+1}`` est ``\Delta t``, on a ``t = n\Delta t``. Il faut alors que ``\xi = \frac{k}{N\Delta t}`` et donc l'écart fréquentiel entre ``X_k`` et ``X_{k+1}`` est ``\frac{1}{N\Delta t}``.
""",
)

# ╔═╡ 57373c5b-a34a-4e48-af42-27a0bd05a173
qa(
	html"Que vaut le signal en dehors de ses N points ?",
	md"""On a
```math
\begin{align}
  X_{k + N} & = \sum_{n=0}^{N-1} x_n e^{-i 2\pi \frac{(k + N)}{N} n}\\
  & = \sum_{n=0}^{N-1} x_n e^{-i 2\pi \frac{k}{N} n} e^{-i2\pi n}\\
  & = \sum_{n=0}^{N-1} x_n e^{-i 2\pi \frac{k}{N} n}\\
  & = X_k
\end{align}
```
Le signal est donc périodique !
Plus précisément, la discrétisation en temporelle rend le signal fréquentiel périodique et la discrétisation en fréquentiel rend le signal temporel périodique.
""",
)

# ╔═╡ b915e8c7-bc11-4f3f-88f0-a82dfe447643
qa(
	html"Que va-t-il advenir de la deuxième moitié du signal ?",
	md"""
Le nouveau signal de `0` à `t_max/2` sera la somme du signal de `0` à `t_max/2` et du signal de `t_max/2` à `t_max`.
""",
)

# ╔═╡ ae2adccf-edfc-4cba-97d8-b440946c9b69
qa(
	html"Qu'observe-t-on dans la matrice ci-dessous ?",
	md"""
La partie mauve et bleue sont égales. La partie rouge est l'opposée de la partie verte.
""",
)

# ╔═╡ 7c96fc08-da4a-48e6-a43e-cd24c9aaef42
qa(
	html"Qu'observe-t-on en comparant cela à la matrice 2 fois plus petite ?",
	md"""
Les parties mauve et bleue sont égale à `F(N ÷ 2)`.
La ligne `k` (en commençant à `k = 1`) est une rotation de ``2\pi(k-1)/N`` radians (donc une multiplication par ``z^{-(k-1)}``) par rapport à `F(N ÷ 2)`.
""",
)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
PortAudio = "80ea8bcb-4634-5cb3-8ee8-a132660d1d2d"
TestImages = "5e47fb64-e119-507b-a336-dd2b206d9990"

[compat]
FFTW = "~1.10.0"
ImageFiltering = "~0.7.10"
Images = "~0.26.2"
Plots = "~1.41.1"
PlutoUI = "~0.7.73"
Polynomials = "~4.1.0"
PortAudio = "~1.3.0"
TestImages = "~1.8.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.1"
manifest_format = "2.0"
project_hash = "3745d20831e3d92d304bf21edf0de1697debf4a7"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BerkeleyDB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll"]
git-tree-sha1 = "77a1bd0eed92aae78fa1bb1318ac53d3c617e9d3"
uuid = "cd00e070-8fe2-570d-8212-aefc8f89bd06"
version = "18.1.41+0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BlueZ_jll]]
deps = ["Artifacts", "Dbus_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libical_jll", "Pkg", "Readline_jll", "eudev_jll"]
git-tree-sha1 = "d4c413db1759fa113135800ff2993ee01206126b"
uuid = "471b5b61-da80-5748-8755-67d5084d21f2"
version = "5.54.0+1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "0df00546373af8eee1598fb4b2ba480b1ebe895c"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.10"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Elfutils_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "XZ_jll", "Zlib_jll", "argp_standalone_jll", "fts_jll", "obstack_jll"]
git-tree-sha1 = "ab92028799ddede63b16af075f8a053a2af04339"
uuid = "ab5a07f8-06af-567f-a878-e8bb879eba5a"
version = "0.189.0+1"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FLAC_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "8476481230247b3671a98f8b3072053bb001102a"
uuid = "1d38b3a6-207b-531b-80e8-c83f48dafa73"
version = "1.3.4+2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "d60eb76f37d7e5a40cc2e7c36974d864b82dc802"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.1"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.3.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.GSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "56f1e2c9e083e0bb7cf9a7055c280beb08a924c0"
uuid = "1b77fbbe-d8ee-58f0-85f9-836ddc23a7a4"
version = "2.7.2+0"

[[deps.GStreamer_jll]]
deps = ["Artifacts", "Elfutils_jll", "GMP_jll", "GSL_jll", "Glib_jll", "JLLWrappers", "LibUnwind_jll", "Libdl", "Pkg", "libcap_jll"]
git-tree-sha1 = "455c99eb5cd91f12943d48f54e34b26765867dc0"
uuid = "aaaaf01e-2457-52c6-9fe8-886f7267d736"
version = "1.20.3+0"

[[deps.Gdbm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Readline_jll"]
git-tree-sha1 = "64929c4ee6b015679b8fc9f2dc36f1b738f13abd"
uuid = "54ca2031-c8dd-5cab-9ed4-295edde1660f"
version = "1.19.0+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "50c11ffab2a3d50192a228c313f05b5b5dc5acb2"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7a98c6502f4632dbe9fb1973a4244eaa3324e84d"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

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

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6b1e49820922eca7bfc862442da6e54173a075b4"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "68.2.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "33485b4e40d1df46c806498c73ea32dc17475c59"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.1"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "eea3a5095c0c5f143e62773164ab11f67e43c4bb"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.10"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "8e64ab2f0da7b928c8ae889c514a52741debc1c2"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.4.2"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "afde851466407a99d48829051c36ac80749d8d7c"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "7.1.1048+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "cffa21df12f00ca1a365eb8ed107614b40e8c6da"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.6"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "7196039573b6f312864547eb7a74360d6c0ab8e6"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.9.0"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "dfde81fafbe5d6516fb864dc79362c5c6b973c82"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.2"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "a49b96fd4a8d1a9a718dfd9cde34c154fc84fcd5"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "b842cbff3f44804a84fda409745cc8f04c029a20"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.6"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues", "TranscodingStreams"]
git-tree-sha1 = "d97791feefda45729613fafeccc4fbef3f539151"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.15"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

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

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

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

[[deps.LibUnwind_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "745a5e78-f969-53e9-954f-d19f2f74f4e3"
version = "1.8.1+2"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libical_jll]]
deps = ["Artifacts", "BerkeleyDB_jll", "Glib_jll", "ICU_jll", "JLLWrappers", "Libdl", "Pkg", "XML2_jll"]
git-tree-sha1 = "c61ffd9e8faf24c19a88f369f1966d53967824d1"
uuid = "bce108ef-3f60-5dd0-bcd6-e13a096cb796"
version = "3.0.9+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libtool_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "88c91b79c1d71166340fb7554bc274876cb3d98e"
uuid = "a76c16ae-fb8f-5ff0-8826-da3b7a640f0b"
version = "2.4.7+4"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "fa7fd067dca76cadd880f1ca937b4f387975a9f5"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.16.0+0"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "a9fc7883eb9b5f04f46efb9a540833d1fad974b3"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.173"

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    ForwardDiffNNlibExt = ["ForwardDiff", "NNlib"]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.LoopVectorization.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "e9650bea7f91c3397eb9ae6377343963a22bf5b8"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.8.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Ncurses_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b5e7e7ad16adfe5f68530f9f641955b5b0f12bbb"
uuid = "68e3532b-a499-55ff-9963-d1c0c0748b3a"
version = "6.5.1+0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "7dc7028a10d1408e9103c0a77da19fdedce4de6c"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "386b47442468acfb1add94bf2d85365dea10cbab"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1f7f9bbd5f7a2e5a9f7d96e51c9754454ea7f60b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.4+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "12ce661880f8e309569074a61d3767e5756a199f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.1"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PortAudio]]
deps = ["LinearAlgebra", "SampledSignals", "Suppressor", "alsa_plugins_jll", "libportaudio_jll"]
git-tree-sha1 = "1c485addb6c281f039d406137a71394afdcb3585"
uuid = "80ea8bcb-4634-5cb3-8ee8-a132660d1d2d"
version = "1.3.0"

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

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.PulseAudio_jll]]
deps = ["Artifacts", "BlueZ_jll", "Dbus_jll", "FFTW_jll", "GStreamer_jll", "Gdbm_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Libtool_jll", "OpenSSL_jll", "SBC_jll", "SoXResampler_jll", "SpeexDSP_jll", "alsa_jll", "eudev_jll", "libasyncns_jll", "libcap_jll", "libsndfile_jll"]
git-tree-sha1 = "df6d51f380df6e16fdae052e15bd2c02a17fe98f"
uuid = "02771fc1-bdb7-5db5-8d11-300768e00fbd"
version = "15.0.1+0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Readline_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ncurses_jll"]
git-tree-sha1 = "6044f482a91c7aa2b82ab614aedd726be633ad05"
uuid = "05236dd9-4125-5232-aa7c-9ec0c9b2c25a"
version = "8.2.13+0"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.SBC_jll]]
deps = ["Libdl", "Pkg", "libsndfile_jll"]
git-tree-sha1 = "34755bff50b6b08988cdfe5fee69c1c1b24ff810"
uuid = "da37f231-8920-5702-a09a-bdd970cb6ddc"
version = "1.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SampledSignals]]
deps = ["Base64", "Compat", "DSP", "FFTW", "FixedPointNumbers", "IntervalSets", "LinearAlgebra", "Random", "TreeViews", "Unitful"]
git-tree-sha1 = "0eaf25f56d43267dc58f6989fc79e2043a649ab6"
uuid = "bd7594eb-a658-542f-9e75-4c4d8908c167"
version = "2.1.4"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "3e5f165e58b18204aed03158664c4982d691f454"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.5.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.SoXResampler_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a95ff1842456719a727e23fe28712eb26f7818b8"
uuid = "fbe68eb6-6641-54c6-99e3-f7c7c4d73a57"
version = "0.1.3+0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SpeexDSP_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "ecc65cb4a4e77f624deae8d881787c789af6deaf"
uuid = "f2f9631b-9a4e-5b48-9975-88f638ec36a7"
version = "1.2.0+0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "a136f98cefaf3e2924a66bd75173d1c891ab7453"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.7"

[[deps.StringDistances]]
deps = ["Distances", "StatsAPI"]
git-tree-sha1 = "5b2ca70b099f91e54d98064d5caf5cc9b541ad06"
uuid = "88034a9c-02f8-509d-84a9-84ec65e18404"
version = "0.11.3"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.Suppressor]]
deps = ["Logging"]
git-tree-sha1 = "6dbb5b635c5437c68c28c2ac9e39b87138f37c0a"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.8"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestImages]]
deps = ["AxisArrays", "ColorTypes", "FileIO", "ImageIO", "ImageMagick", "OffsetArrays", "Pkg", "StringDistances"]
git-tree-sha1 = "0567860ec35a94c087bd98f35de1dddf482d7c67"
uuid = "5e47fb64-e119-507b-a336-dd2b206d9990"
version = "1.8.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "02aca429c9885d1109e58f400c333521c13d48a0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.4"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

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

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "83360bda12f61c250835830cc40b64f487cc2230"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.1"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "d1d9a935a26c475ebffd54e9c7ad11627c43ea85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.72"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.alsa_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "10a5bb558f9aa896cfd87daa9856660cfe884cd6"
uuid = "45378030-f8ea-5b20-a7c7-1a9d95efb90e"
version = "1.2.12+0"

[[deps.alsa_plugins_jll]]
deps = ["Artifacts", "FFMPEG_jll", "JLLWrappers", "Libdl", "Pkg", "PulseAudio_jll", "alsa_jll", "libsamplerate_jll"]
git-tree-sha1 = "a43b5bcdfadfbe06c42cd6b007572c4806f2c0f7"
uuid = "5ac2f6bb-493e-5871-9171-112d4c21a6e7"
version = "1.2.2+0"

[[deps.argp_standalone_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "feaf9f6293003c2bf53056fd6930d677ed340b34"
uuid = "c53206cc-00f7-50bf-ad1e-3ae1f6e49bc3"
version = "1.3.1+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fts_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa21810b841ae26d2fc7f780cb1596b4170a4c49"
uuid = "d65627f6-89bd-53e8-8ab5-8b75ff535eee"
version = "1.2.8+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libasyncns_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "38a54b0ebad9bc225a38106ff66b7827fac5bd9e"
uuid = "ed080073-db63-57db-a029-74e11ae80737"
version = "0.8.0+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libcap_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "500e32aac3dda306137e6dc09675a3fdbef91c41"
uuid = "eef66a8b-8d7a-5724-a8d2-7c31ae1e29ed"
version = "2.51.0+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libportaudio_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "alsa_jll"]
git-tree-sha1 = "fbce8030d68816899cd5f068670feaad67e84e4a"
uuid = "2d7b7beb-0762-5160-978e-1ab83a1e8a31"
version = "19.7.0+0"

[[deps.libsamplerate_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "45ba80d9b0a208fd5165d159d93a3725fab0d76b"
uuid = "9427e74d-4e05-59c1-8ff3-7d74b6e52ac8"
version = "0.1.9+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libsndfile_jll]]
deps = ["Artifacts", "FLAC_jll", "JLLWrappers", "Libdl", "Ogg_jll", "Opus_jll", "Pkg", "alsa_jll", "libvorbis_jll"]
git-tree-sha1 = "f35a5fbfb2b18ff837dec4594c7e096ac6604154"
uuid = "5bf562c0-5a39-5b4f-b979-f64ac885830c"
version = "1.1.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "ccbb625a89ec6195856a50aa2b668a5c08712c94"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.4.0+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.obstack_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5353d2b8d19b8ed8d972a4bed38fff85d27f7f73"
uuid = "c88a4935-d25e-5644-aacc-5db6f1b8ef79"
version = "1.2.3+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╟─b8949a96-8d97-4048-8022-ab8fabf326a1
# ╟─4309cf1d-2a7a-4a84-996b-ecc99d62dad1
# ╟─df7fa000-70ce-462b-af19-fb9c4b022012
# ╟─58d25489-0d99-4c33-ad22-5ffb2751c3c6
# ╟─5010b40e-86c1-4ed1-b7ae-25086ad45205
# ╟─27d19b85-8f5b-4d53-93c2-033cca3b9006
# ╟─04ac4af9-5789-4d0a-a958-1630d0a54299
# ╟─41333814-27a7-4ecf-99ce-4f972c0a37da
# ╟─5d64d5ac-ff4f-415a-85bb-c673da00e9a1
# ╟─308757bd-bfbc-43d5-8e72-5903390c1458
# ╟─93534f3d-f9da-46c4-a9c3-d055530f3045
# ╟─8f77bda0-25b0-44fb-b235-afb6c2e725ce
# ╟─83f8f109-50fa-4174-8fcd-6bd9d8ef2d65
# ╟─54e54a58-3fdf-43d8-a87a-22ff047a392e
# ╟─83c1cbf7-e24d-4e9a-a309-c7b60e730344
# ╠═44dbdbff-18a2-43ab-b8e5-399bf8fe9639
# ╟─590c895f-d7d4-467a-90cd-90bfba59c2b5
# ╠═40490a0f-9a06-4283-9d3a-d2c96b6b14df
# ╟─41155d3b-ba20-49ee-af03-babd4ce10c68
# ╠═f18256be-d172-4f10-9c6b-662015f535e4
# ╟─b790ddf3-37ac-45d6-87aa-ebd06c829d82
# ╠═02084ab9-305a-4317-9fad-a4acf8cc5cd3
# ╟─4a11717b-4086-44f3-9a50-b531d51c6944
# ╠═0ba2c1bc-aa82-4fbd-a704-034274349c4d
# ╟─ef156c25-cb7d-4cd3-8412-d9f5f2bc7a3b
# ╟─4265ccd9-3bf2-4655-a0da-2d774baf32ac
# ╟─c26a3bd3-cf79-426a-b439-2b29501cee4b
# ╟─a315cbaf-a4b5-42ea-a0e0-f0b064b9995e
# ╟─bb5eebce-e385-4c4f-b827-0ffcdead0799
# ╟─2488c83e-5e54-4d47-80d2-6ef41ca7ac7f
# ╟─279e5f5e-4fc6-475f-8dab-c5740a470b0e
# ╟─1dfad36a-05a5-4108-8c67-d4c866be89c0
# ╟─7cdf9587-af9f-4346-b619-ba8bdf6d6799
# ╠═d4ac67de-9442-464c-be8b-e8ff468dd996
# ╟─54ced290-6da4-4f24-a18a-e46ac116448c
# ╟─93cfa0b9-dbb6-4956-bd94-575fe233c9f4
# ╠═3855c981-206d-46e6-8a9b-eda011769208
# ╟─5dd93a6b-f355-4927-8f3d-a3522f43934f
# ╟─53045a43-9941-4b3f-a0be-d6de969ee124
# ╠═8f9cccf1-ec8d-448d-a86f-3955f6c238d8
# ╟─5edc270e-a7c6-4e44-97b0-7697adca16a2
# ╠═e1c25e06-3259-4a56-b26a-60903e57f92b
# ╟─13d00d91-dae3-451c-924a-6b120c0a605d
# ╟─802cbd16-2126-4d4f-aca2-df29024d5fdf
# ╠═9fae43f1-61cd-43c6-b383-d9bbd4c2a1ac
# ╠═31aca097-4a85-4d63-b114-6d1828648e87
# ╠═e7b74b94-bf5f-4455-9fc7-922ce7d4c284
# ╟─6b9a064c-1097-4754-af85-5dcaa435b893
# ╟─101663c2-acca-4eb7-858c-20a70c2103e9
# ╠═ae1e22a7-b4af-4f83-8f28-97534d7e0c20
# ╟─eeafd828-82a7-449e-b083-f22eb0d07079
# ╟─66fd9f77-6c2a-418d-9650-39033396b584
# ╠═0230ce41-bf49-4363-b344-3c475ef6474d
# ╠═680c44ba-02b9-41d1-b85f-86f49a85ca7a
# ╠═14a60621-f76e-48f6-8dbe-cf4c871b2027
# ╠═36be0072-a89e-46e8-84bc-57e208a9e70a
# ╠═c312e94d-a461-4d8a-bf56-99699e82fe04
# ╠═4580177c-9654-4950-bbf5-a599adb44db6
# ╠═9500f893-7eb0-4f24-83eb-5a2240279062
# ╟─976ac813-8f46-4912-a04a-ec586155867d
# ╟─87dbf272-1d03-4f9e-b7e1-9b76197b0a83
# ╟─44f73082-7cdb-47b9-bb86-f952ac7c4498
# ╟─a2359c44-3150-4b0a-8675-90d1627b5c07
# ╟─86368d24-c206-42dc-afdb-1e762b3ac14c
# ╟─3524f1bf-8e39-4d04-9dc3-e434e65db409
# ╟─f0eda0e0-0a52-403a-a73d-3a604cbe557c
# ╟─57373c5b-a34a-4e48-af42-27a0bd05a173
# ╟─d29773c7-a413-46a9-8c62-2f54d4de1730
# ╟─b915e8c7-bc11-4f3f-88f0-a82dfe447643
# ╟─49d078e2-d865-4f2f-9cfc-3c2d5a535088
# ╟─b80ff616-8bb7-4132-88c7-ed10c3c6786c
# ╟─56ff3957-0295-41fe-839f-12dae0fdd14f
# ╟─4980c4ff-7da0-4747-b779-891569d604fd
# ╟─5032f1f1-47a1-4ff3-b3fd-49d60e5e28a5
# ╟─78026d88-7051-11ef-29f0-67f85176a548
# ╟─dc4c4d53-feb2-40ce-b20e-3386aab2a45f
# ╟─5f9eeb76-0a9d-40ad-858e-4ab67af46427
# ╠═222cfba3-88f7-4824-92ca-f9ab81a000c4
# ╟─ae2adccf-edfc-4cba-97d8-b440946c9b69
# ╟─a33ebf12-5c62-4cc9-af1e-7bad9026db21
# ╠═cd3b2930-2c33-4028-9b27-d25238c79bac
# ╟─c64e6062-3793-4556-9001-bf184cd090e9
# ╠═c9131a2a-8978-4c81-8588-44e4d7071c0b
# ╠═f24f49d9-10e6-4ff5-ace7-8767c69f8594
# ╠═ff3ec271-5e08-4ab8-8c13-7f5e1996e383
# ╠═8676e10b-6c41-4e8f-a56b-ca548b9690dd
# ╠═faf71e8e-7ba0-44fa-800d-949a42fbacec
# ╠═fd5ef1e7-bc1e-4382-b5c6-e960ecb3c893
# ╠═f952619b-194a-4088-a2bb-b8210392417a
# ╟─7c96fc08-da4a-48e6-a43e-cd24c9aaef42
# ╠═3739ddf4-ecb8-49d5-a70c-aa347c6fb17d
# ╟─ecf827b8-6a01-4e3b-9856-451239a7b175
# ╟─6119e719-4e14-4516-ac7e-a236cfb4fcfa
# ╠═aa4ff136-6739-46dc-aae6-cebcb3b31263
# ╠═3567a3cb-a5be-4d7c-a917-d3947b5549f6
# ╠═0c74b1e8-a67f-4e41-9638-bb10babf7138
# ╠═4ead6b3a-af8d-4dc6-905e-bfc4b6a00e64
# ╠═6d79251e-f7d5-4126-a781-7221a4313c9a
# ╠═59bc8221-938c-4919-8b94-7a2d61cf2062
# ╠═96c678be-8528-462f-b455-446b9a367fab
# ╠═9a7a1c48-8168-4a21-b57e-0aaf2dbb9de8
# ╟─f53db6db-d9f9-4371-a551-e3bfc2d6b162
# ╟─bcf5ecbc-ccb6-48ce-9146-55b641f770f2
# ╟─43c13895-a651-4e61-8fda-abdde42718da
# ╟─59130f1a-4fcd-4ae1-9ecf-1818dfc07612
# ╠═583c3b75-71cf-48f6-b93d-2201c6c7d520
# ╠═34c779d2-5a9d-4c98-ac28-e4be1a85da40
# ╠═e3a7dd72-f106-4865-8b24-d86056782160
# ╠═1817a169-0afc-4fff-be23-1ea510ff337a
# ╠═0f87464a-1a55-4f7c-8a88-bf514f5177a0
# ╠═9e6732a7-d2df-465f-a7fa-c119a64deb9a
# ╠═0e9e8c6e-83ff-4c4c-b9d1-19a0593a3cba
# ╠═5e5ff9b1-0a78-4e9c-9c90-4fb80b0dbe21
# ╠═0d9e0a31-b229-4d06-92c9-0136404399da
# ╠═f80595c3-6bbb-4b62-b03b-1c80cecaf7f9
# ╠═792723a3-0011-422e-bbf4-798624913780
# ╠═b0a265b8-e3fd-4a6c-b813-904d18ea83e7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
