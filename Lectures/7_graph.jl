### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f83fdf0a-a2f5-4d72-9297-a409cf4e2bbb
using PlutoUI, Graphs, GraphPlot, Colors, DataStructures, SimpleWeightedGraphs

# ╔═╡ bd9ab2c4-b098-11ef-3c1e-7fe1458a8556
include("utils.jl")

# ╔═╡ c9595f48-bfec-443e-8f8f-94ed0c5a9530
section("Graph Theory")

# ╔═╡ 15f7cc14-7663-4ad2-9cfe-9ee8aada601c
frametitle("Seven Bridges of Königsberg")

# ╔═╡ b65dd471-10da-4b3a-a14b-cfbb70e70050
frametitle("Définition")

# ╔═╡ d7bbc61f-5cc6-4090-8cc4-313187f28bab
let
	g = Graph(4)
	add_edge!(g, 1, 2)
	add_edge!(g, 1, 4)
	add_edge!(g, 2, 3)
	add_edge!(g, 2, 4)
	add_edge!(g, 3, 4)
	gplot(g, edgelabel = [1, 2, 1, 1, 2])
end

# ╔═╡ b52b2192-df63-4095-aa15-6d3ccafc92ba
frametitle("Graphes et polyhèdres")

# ╔═╡ e41498fa-8fed-4354-bb3b-a2fc8afeacd9
frametitle("Terminologie des parcours")

# ╔═╡ c41901fb-6570-487d-bd77-4708d1f0b627
frametitle("Degré d'un noeud")

# ╔═╡ fe642397-9ae1-42b5-9d0e-d2b96cfb178b
let
	g = Graph(4)
	add_edge!(g, 1, 2)
	add_edge!(g, 1, 4)
	add_edge!(g, 2, 3)
	add_edge!(g, 2, 4)
	add_edge!(g, 3, 4)
	gplot(g, nodelabel = [3, 3, 3, 5], edgelabel = [1, 2, 1, 1, 2])
end

# ╔═╡ 043ff87b-41e4-49c2-aaee-bee39f088b58
frametitle("Composante connexe")

# ╔═╡ 9e065cf0-7483-4ec2-809c-cdd7d22a5564
let
	g = Graph(7)
	add_edge!(g, 1, 2)
	add_edge!(g, 1, 4)
	add_edge!(g, 2, 3)
	add_edge!(g, 2, 4)
	add_edge!(g, 3, 4)
	add_edge!(g, 5, 6)
	add_edge!(g, 5, 7)
	add_edge!(g, 6, 7)
	gplot(g)
end

# ╔═╡ cbc8da76-9adf-4670-b07c-0f5366b55547
frametitle("Calcul de composantes connexes")

# ╔═╡ b184882e-f2ed-4e9c-b763-0173711d6fac
frametitle("Disjoint-Set datastructure")

# ╔═╡ 8e7b6a51-1f8e-43d7-846b-a9e7596052b9
frametitle("Piste Eulérienne")

# ╔═╡ d884c51e-9ddd-4dc4-9f9c-e461f0455a17
frametitle("Arbres et forêts")

# ╔═╡ ceb5aea4-5a69-47cf-bedd-17d47a6e07a8
gplot(uniform_tree(10))

# ╔═╡ b94745d7-ffda-4190-82af-8ee2fe9ea165
frametitle("Spanning tree")

# ╔═╡ 944b6ad9-731b-4479-94bd-e61876987917
function spanning_forest(g, edges_it = edges(g))
	component = IntDisjointSets(nv(g))
	function create_loop!(edge)
		loop = in_same_set(component, src(edge), dst(edge))
		union!(component, src(edge), dst(edge))
		return loop
	end
	return [edge for edge in edges_it if !create_loop!(edge)]
end

# ╔═╡ 3c09769f-13b2-4713-8d2e-06abd7735c33
frametitle("Minimum spanning tree")

# ╔═╡ 3eb8d016-d1a1-4e4f-a307-e2c439c502e4
kruskal(g) = spanning_forest(g, sort(collect(edges(g)), by = weight))

# ╔═╡ 8c384e86-18d4-419b-977b-c5d5bfd1318f
frametitle("Graphes dirigés")

# ╔═╡ 49092fe9-7c24-4564-9d2a-aa313d54ed8a
gplot(random_regular_digraph(16, 2))

# ╔═╡ 264eb1b6-0e3f-4d67-b8b0-5f0831d18b95
frametitle("Strongly connected components")

# ╔═╡ 66c7051b-c377-4376-a472-b2bd5c788267
frametitle("Directed Acyclic Graph (DAG)")

# ╔═╡ 71566570-7d26-4fba-8f9a-8726fb97f2d8
frametitle("Suite de Fibonacci")

# ╔═╡ 2f905018-1d7b-43dc-b576-e4080dc946ee
function fib(n)
	if n <= 2
		return 1
	else
		return fib(n - 1) + fib(n - 2)
	end
end

# ╔═╡ a95116d6-614a-4a8a-b3a6-5a96e9b3c439
frametitle("Visualisation des appels récursifs")

# ╔═╡ 3c82aba6-0539-4578-ad04-648eb3bcf400
frametitle("Nombre d'appels")

# ╔═╡ 30861bc6-4af9-43dc-a1f5-fed9bb157564
frametitle("Bottom-Up / Dynamic Programming")

# ╔═╡ 9b024507-48e2-4247-aa4e-ed0958af7b2d
function fib_bottom_up(n)
	fib = zeros(Int, n)
	fib[1] = fib[2] = 1
	for k in 3:n
		fib[k] = fib[k - 1] + fib[k - 2]
	end
	return fib[n]
end

# ╔═╡ a9379cc3-6ca8-4aed-81f7-992f3664ceb3
frametitle("Top-Down / Mémoïzation")

# ╔═╡ c395343c-3f7c-4dfe-bc0f-8950afc115bf
function fib!(cache, n)
	if cache[n] == 0
		cache[n] = fib!(cache, n - 1) + fib!(cache, n - 2)
	end
	return cache[n]
end

# ╔═╡ 73c0c135-f070-422e-ab1d-20c6faa160d8
function cached_fib(n)
	cache = zeros(Int, n)
	cache[1] = cache[2] = 1
	return fib!(cache, n)
end

# ╔═╡ 16c5f723-77d4-4bee-86c6-90a1b597e7d9
frametitle("Une IA pour le jeu Puissance 4")

# ╔═╡ 54d0cbaf-4bf7-4508-b706-b5042dda1644
frametitle("Puissance 2 de 2 colonnes de hauteur 2")

# ╔═╡ 73d104ac-7aff-46cc-95b5-a67b6811f64d
frametitle("Puissance 2 de 2 colonnes de hauteur 3")

# ╔═╡ 3de604e8-20b6-4ef5-9952-b8f27c28fef1
begin
	slider_longueur = @bind longueur Slider(2:4, default = 2, show_value = true);
	slider_nrows = @bind nrows Slider(2:6, default = 3, show_value = true);
	slider_ncols = @bind ncols Slider(2:6, default = 3, show_value = true);
	checkbox_patient = @bind patient CheckBox()
	nothing
end

# ╔═╡ 1308d56a-f91a-4fe2-8f5d-7d7a15b8f040
frametitle("Puissance $longueur : sans mémoïzation")

# ╔═╡ a7530e0a-a3b5-46c8-9666-fb0219d7ba04
frametitle("Mémoïsation pour Puissance 4")

# ╔═╡ 0a7e65f5-da36-4331-8463-221051286728
begin
	slider_wooclap_nrows = @bind wooclap_nrows Slider(2:6, default = 3, show_value = true);
	slider_wooclap_ncols = @bind wooclap_ncols Slider(2:6, default = 3, show_value = true);
	nothing
end;

# ╔═╡ 3b82fa0a-65b1-46dd-b81a-d1fa117893d5
begin
	slider_mem_longueur = @bind mem_longueur Slider(2:4, default = 2, show_value = true);
	slider_mem_nrows = @bind mem_nrows Slider(2:6, default = 3, show_value = true);
	slider_mem_ncols = @bind mem_ncols Slider(2:7, default = 3, show_value = true);
	checkbox_mem_patient = @bind mem_patient CheckBox()
	nothing
end

# ╔═╡ 3f813480-726b-42e0-a648-e3ad76abcf9d
frametitle("Mémoïsation pour Puissance $mem_longueur")

# ╔═╡ a436f564-8418-477c-9b0d-9d791c893bb1
frametitle("Pièces de monnaie")

# ╔═╡ 9db78023-e29c-41e5-9102-870692111f39
frametitle("Depth-First Search (DFS)")

# ╔═╡ 6dc95039-9e2a-47e0-8250-9c5ad8a2af0c
function dfs!(path, g::DiGraph, current, target)
	if current in path
		return false
	end
	push!(path, current)
	if current == target
		return true
	end
	for next in Graphs.outneighbors(g, current)
		if dfs!(path, g, next, target)
			return true
		end
	end
	pop!(path)
	return false
end

# ╔═╡ 2829213a-8533-4e92-a130-57d53166b376
function dfs(g, current, target)
	path = Int[]
	dfs!(path, g, current, target)
	return path
end

# ╔═╡ 2c7c9b52-0d27-4dc8-ae69-29a6edc87041
frametitle("Breadth-First Search (BFS)")

# ╔═╡ fb29e166-b361-49e3-bc5c-042d409023f1
function bfs(g::DiGraph, source)
	dist = fill(-1, nv(g))
	queue = Queue{Tuple{Int,Int}}()
	dist[source] = 0
	enqueue!(queue, (source, 0))
	while !isempty(queue)
		current, d = dequeue!(queue)
		for next in outneighbors(g, current)
			if dist[next] == -1
				dist[next] = d + 1
				enqueue!(queue, (next, d + 1))
			end
		end
	end
	return dist
end

# ╔═╡ 670aaab7-a101-4dea-ba97-80f1a5346122
frametitle("Dijkstra")

# ╔═╡ ab3997f3-83ef-4435-bb90-7cef1c7bcdcf
function dijkstra(g::SimpleWeightedGraph, source)
	dist = fill(-1, nv(g))
	dist[source] = 0
	queue = PriorityQueue{Int,Int}()
	push!(queue, source => 0)
	while !isempty(queue)
		current = dequeue!(queue)
		for next in outneighbors(g, current)
			d = dist[current] + g[current, next, Val(:weight)]
			if dist[next] == -1
				dist[next] = d
				enqueue!(queue, next => d)
			elseif d < dist[next]
				dist[next] = d
				queue[next] = d
			end
		end
	end
	return dist
end

# ╔═╡ 634dea01-3832-49e7-90af-4e1b834cead7
frametitle("Max Flow / Min Cut")

# ╔═╡ 75bb5e30-76b8-486b-b3fa-9979aaf830a5
section("Utils")

# ╔═╡ accdf77c-fc4b-4777-bebf-9d69ebb36526
import DocumenterCitations, Polyhedra, Makie, WGLMakie, StaticArrays, Bonito, Luxor, Formatting, Random

# ╔═╡ 730e5e09-e818-49c5-8fa1-495b39b4267e
let
	Random.seed!(32)
	g = random_regular_digraph(32, 2)
	c = Set.(strongly_connected_components(g))
	cols = distinguishable_colors(length(c))
	tree = Set(spanning_forest(g))
	gplot(g, nodefillc = [cols[findfirst(comp -> v in comp, c)] for v in vertices(g)])
end

# ╔═╡ c7663ad1-16b4-4b29-8f4a-8993d094a7e1
let
	Random.seed!(32)
	n = 32
	g = random_regular_digraph(n, 2)
	dist = bfs(g, 1)
	gplot(g, nodelabel = string.(dist))
end

# ╔═╡ b01e3f5a-ad9f-4557-ba62-30e9c57b8532
begin
	slider_n = @bind n Slider(1:42, default = 10, show_value = true);
	slider_k = @bind k Slider(1:42, default = 10, show_value = true);
	nothing
end

# ╔═╡ b53e5719-cc0e-4169-9c8d-3eb5a50570b0
@time fib(n)

# ╔═╡ 6183e415-898a-4d1f-a5db-9c09503affff
called = let
	called = zeros(Int, n)
	function fib(k)
		called[k] += 1
		if n > 42
			error("Don't use `n` higher than 42")
		elseif k <= 2
			return 1
		else
			return fib(k - 1) + fib(k - 2)
		end
	end
	fib(n)
	called
end;

# ╔═╡ 46977b09-a43a-427c-8974-a1fd5f2d7da2
Makie.barplot(1:n, called)

# ╔═╡ 4c76167b-a8ac-436f-87a8-d493d489e61f
@time fib_bottom_up(n)

# ╔═╡ 9eeaca05-9e6e-4881-85aa-d81ad5460e8c
@time cached_fib(n)

# ╔═╡ fdd9318e-ee09-4eab-9322-ad97b03da002
function draw_fib(n, depth)
	if n < 3
		return 0
	end
	dy = 140
	δ = round(Int, dy * 0.9^depth) * 1.2
	w = fib(n - 2)
	Δy = 30
	if depth != 1
		Luxor.translate(δ * w/2, 0)
	end
	if n == 3
		Luxor.setdash("dashed")
		Luxor.circle(Luxor.Point(6, -6), 20, :stroke)
	end
	Luxor.sethue("black")
	Luxor.text(string(n))
	Luxor.translate(-δ * w/2, 0)
	cur = 0
	for m in [n - 1, n - 2]
		if m < 3
			continue
		end
		dx = cur * δ
		Luxor.translate(dx, dy)
		dw = draw_fib(m, depth + 1)
		Luxor.sethue("black")
		Luxor.setdash("solid")
		Luxor.setline(δ/40)
		Luxor.line(Luxor.Point(-dx+δ * w/2, -dy+Δy), Luxor.Point(δ*dw/2, -Δy), :stroke)
		Luxor.translate(-dx, -dy)
		cur += dw
	end
	return w
end

# ╔═╡ 039d5065-4ab2-4787-86e2-5612cfcd01b4
function draw_fibo(n)
	Luxor.@draw begin
		Luxor.translate(0, -350)
		Luxor.fontsize(30)
    	Luxor.fontface("Alegreya Sans SC")
		draw_fib(n, 1)
	end 1200 800
end

# ╔═╡ d35df21e-2dd5-478f-b0c1-17cbd96a9cc6
draw_fibo(8)

# ╔═╡ 3c9dd530-8a25-4c1b-b1ff-3b5a3748b74a
function random_weighted(n, k, weights)
	g = random_regular_graph(n, k)
	wg = SimpleWeightedGraph(n)
	for e in edges(g)
		add_edge!(wg, src(e), dst(e), rand(weights))
	end
	return wg
end

# ╔═╡ 0ee692b3-fcf6-418a-b3ba-39febcf6eed2
let
	Random.seed!(0)
	g = random_weighted(10, 3, 1:9)
	dist = dijkstra(g, 1)
	gplot(g, nodelabel=string.(dist), edgelabel=string.(Int.(weight.(edges(g)))))
end

# ╔═╡ 33ac9568-a83d-47b0-a17c-6e30df0dbaa4
colors = Colors.JULIA_LOGO_COLORS

# ╔═╡ d1329a48-5bed-428f-8554-854a5d51488e
let
	g = random_regular_graph(10, 5)
	cols = [colors.red, colors.green]
	tree = Set(spanning_forest(g))
	gplot(g, edgestrokec = [cols[(edge in tree) + 1] for edge in edges(g)], nodefillc = colors.blue)
end

# ╔═╡ ecf6a7d0-3279-4bff-b0f4-73327f60b702
let
	Random.seed!(0)
	g = random_weighted(10, 3, 1:9)
	cols = [colors.red, colors.green]
	tree = Set(kruskal(g))
	gplot(g, edgelabel=string.(Int.(weight.(edges(g)))), edgestrokec = [cols[(edge in tree) + 1] for edge in edges(g)], nodefillc = colors.blue)
end

# ╔═╡ b54d9e8b-83ab-434e-9e6c-f6774558abad
function plot_path(g, path)
	cols = [colors.red, colors.green]
	gplot(g, nodefillc = [cols[v in path ? 1 : 2] for v in vertices(g)], nodelabel = [something(findfirst(isequal(v), path), "") for v in vertices(g)])
end

# ╔═╡ f8d160f8-7e41-4bdf-99f4-d28b344a46a4
let
	Random.seed!(32)
	n = 32
	g = random_regular_digraph(n, 2)
	path = dfs(g, 1, n)
	plot_path(g, path)
end

# ╔═╡ 1b748b69-5262-46c9-aeff-01696c585e6f
biblio = load_biblio!()

# ╔═╡ d57a6a11-0387-4693-a50e-9eba62605a49
polyhedron = Polyhedra.polyhedron(Polyhedra.vrep(Float64[
	0 1 2 3;0 2 1 3; 1 0 2 3; 1 2 0 3; 2 0 1 3; 2 1 0 3;
    0 1 3 2;0 3 1 2; 1 0 3 2; 1 3 0 2; 3 0 1 2; 3 1 0 2;
    0 3 2 1;0 2 3 1; 3 0 2 1; 3 2 0 1; 2 0 3 1; 2 3 0 1;
    3 1 2 0;3 2 1 0; 1 3 2 0; 1 2 3 0; 2 3 1 0; 2 1 3 0
]))

# ╔═╡ 2fedb811-6f34-424f-a46b-806bb0c6b9c0
projected = Polyhedra.project(polyhedron, [1 1 1; -1 1 1; 0 -2 1; 0 0 -3])

# ╔═╡ 12fa3c11-c5d3-4c81-a92f-174f09ec10ff
struct Config{N}
	libre::StaticArrays.SVector{N,UInt8}
	player::StaticArrays.SVector{N,UInt8}
	max::UInt8
end

# ╔═╡ b107f529-f47f-4252-b06a-6e53302be860
begin
_no_winner(::AbstractMatrix) = 0
_no_winner(::Config) = (false, false)
_has_winner(x::Int) = !iszero(x)
_has_winner(x::Tuple{Bool,Bool}) = x[1]
end

# ╔═╡ e2f5854a-c6ce-4600-8e96-62ac8ffc6b6e
const Δs = [(-1, 0), (0, 1), (0, -1), (1, 1), (-1, 1), (1, -1), (-1, -1)]

# ╔═╡ 938e0359-0a86-4755-a92b-f1a7a892961d
function winner(config, longueur, pos)
	for Δ in Δs
		win = winner(config, pos, Δ, longueur)
		if _has_winner(win)
			return win, Δ
		end
	end
	return _no_winner(config), (0, 0)
end

# ╔═╡ 1c762a4a-4ee1-46ee-93ee-7b566376d5a4
begin
_first_row(::AbstractMatrix) = 1
_next_row(row::Int) = row + 1
_done(m::AbstractMatrix, row) = row > size(m, 1)
_first_row(::Config) = 0b1
_next_row(row::UInt8) = row << 1
_done(c::Config, row) = row == c.max
end

# ╔═╡ 8dea1779-e432-4fe3-baf0-5b0d8d3dd97b
begin
_shift(x::Int, y::Int) = x + y
_shift(x::UInt8, y::Int) = x << y
_inrows(m::AbstractMatrix, x::Int, y::Int) = in(x + y, axes(m, 1))
function _inrows(config::Config, x::UInt8, y::Int)
	if y == 0
		return true
	elseif y < 0
		return x >= (0b1 << (UInt8(-y)))
	else
		return (config.max >> UInt8(y)) > x
	end
end
end

# ╔═╡ efbdb43d-3c6b-48d6-be40-51369f3158de
begin
	_shift_bit(x::Int) = x
	function _shift_bit(x::UInt8)
		y = findfirst(i -> (0b1 << (i-1)) == x, 1:7)
		return y
	end
end

# ╔═╡ 3ca057ef-1367-4fef-8e18-bdef855bc35e
begin
_num_rows(mat::AbstractMatrix) = size(mat, 1)
_num_rows(config::Config) = _shift_bit(config.max) - 1
_reverse(config::Config) = Config(reverse(config.libre), reverse(config.player), config.max)
Base.getindex(c::Config, i::UInt8, j::Int) = (!iszero(c.libre[j] & i), !iszero(c.player[j] & i))
_compare(player::Int, a, b) = player * a > player * b
_compare(player::Bool, a, b) = player ? a < b : a > b
_sign(player::Int, a) = player * a
_sign(w::Tuple{Bool,Bool}, a) = w[2] ? -a : a
_first_player(::AbstractMatrix) = 1
_first_player(::Config) = false
_next_player(player::Int) = -player
_next_player(player::Bool) = !player
_columns(m::AbstractMatrix) = axes(m, 2)
_columns(config::Config) = eachindex(config.libre)
function _next(mat::AbstractMatrix, col)
	for r in axes(config, 1)
		if iszero(config[r, col])
			return r
		end
	end
	return 0
end
function _next(config::Config, col)
	mask = config.libre[col] + 0b1
	if mask == config.max
		return 0b0
	end
	return mask
end
function _assign(mat::StaticArrays.SMatrix{I,J}, r, c, value) where {I,J}
	target = r + (c - 1) * I
	return StaticArrays.SMatrix{I,J}(ntuple(i -> i == target ? value : mat[i], Val(I * J)))
	#return @SMatrix [(row, col) == (r, c) ? value : mat[row, col] for row in 1:nrows, col in 1:ncols]
end
function _assign(v::StaticArrays.SVector{N}, mask, col) where {N}
	StaticArrays.SVector{N}(ntuple(c -> c == col ? (v[c] | mask) : v[c], Val(N)))
end
function _assign(config::Config, mask, col, value::Bool)
	libre = _assign(config.libre, mask, col)
	if value
		player = _assign(config.player, mask, col)
	else
		player = config.player
	end
	return Config(libre, player, config.max)
end
function _assign(mat::Matrix, row, col, value)
	mat[row, col] = value
	return mat
end
_unassign(::StaticArrays.SMatrix, _, _) = nothing
_unassign(::Config, _, _) = nothing
function _unassign(mat::Matrix, row, col)
	mat[row, col] = 0
	return
end
_copy(m::StaticArrays.SMatrix) = m
_copy(m::Config) = m
_copy(m::Matrix) = copy(m)
end

# ╔═╡ faca26cb-14a3-4a45-ba0b-ec3cd57d8d29
HAlign(md"""
Est-il possible de prendre tous les ponts de la ville une et une seule fois (le point de départ est au choix) ?
""",
img("Konigsberg_bridges"))

# ╔═╡ 8bec0324-274c-4ae6-99f9-08b613460064
md"""
Un graph a un ensemble ``V`` de noeuds (nodes) / sommets (vertices) et ``E`` d'arêtes (edges) / arcs.

| Math | Informatique | Non-dirigé | Dirigé |
|------|--------------|------------|--------|
| ``V`` | Nodes | Sommet / vertex | Noeud / node |
| ``E`` | Edges | Arête / edge | Arc |

Le graphe du problème de Königsberg test:
"""

# ╔═╡ c87e1a66-62e6-41aa-bcc0-fb1704d7e259
md"""
Prenons un polyhèdre, où se cache le graphe ?
"""

# ╔═╡ 55283461-b2e6-4e81-884f-48ef81e605d4
md"""
On associe un sommet du graphe à chaque sommet du polyhèdre et une arêtes à chaque arêtes (c'est la même terminologie, coincidence ?).
"""

# ╔═╡ 73642ba1-918f-450c-a555-7bcd393ec54e
md"""
**Relation d'Euler**: le nombre de faces est égal à
``
2 - |V| + |E|
``.
"""

# ╔═╡ 27ee39a0-0eb2-42dc-bf90-eb2fa0d7497e
md"""
Un parcours / walk de longueur ``k`` d'un graphe ``G = (V, E)`` est une suite de ``k+1`` noeuds ``v_0, v_1, \ldots, v_k \in V`` et ``k`` arêtes ``(v_0, v_1), (v_1, v_2), \ldots, (v_{k-1}, v_k) \in E``.
Si ``v_0 = v_k``, le parcours est fermé.

|    | Noeuds distincts | Arêtes distinctes |
|----|------------------|-------------------|
| Ouvert | Chemin | Piste / Trail  |
| Fermé  | Cycle  | Circuit |

Piste (resp. circuit) *Eulérienne* : Une piste (resp. circuit) qui visite **toutes** les arêtes.
"""

# ╔═╡ ff87b1bb-7903-4954-8e7a-f377aef97b2d
md"""
Le degré d'un noeud ``v`` est le nombre d'arête qui lui sont incidentes.
"""

# ╔═╡ 2ebc96d9-e526-408a-a264-a87771b96258
md"""
Une composante connexe ``C`` est un ensemble de noeuds tels que toutes paires de noeuds est connectée par un chemin.
"""

# ╔═╡ 338527b3-ec32-4794-8b92-63cce6ce4285
md"""
On démarre en assignant une composant connexe différente pour chaque noeud.
Pour chaque arête, on fusionne les composantes connexes des deux noeuds liés par l'arête. Comment calculer cette fusion efficacement ?

Supposons qu'on ait 4 noeuds. On démarre avec le tableau ``(1, 2, 3, 4)`` signifiant que chaque noeud est dans une composante connexe différente.

* En commençant avec l'arête ``(1, 2)``, on fusionne et on arrive au tableau ``(1, 1, 3, 4)``.
* Si on voit ensuite l'arête ``(3, 4)``, on arrive alors au tableau ``(1, 1, 3, 3)``.
* Si on voit ensuite l'arête ``(1, 3)``, on update le tableau à ``(1, 1, 1, 1)``.
"""

# ╔═╡ 269320fe-a575-4fe4-86b6-67b0dd4484a5
qa(md"Quelle est la complexité de cet algorithm ?", md"""
À la 3ième étape, on doit mettre à jours deux éléments du tableau. Dans le pire cas, si on doit fusionner deux composante connexes dont la taille est ``O(|V|)``, la complexité d'une fusion peut être ``O(|V|)``.
Comme on a ``|E|`` fusion, la complexité de l'algorithm est ``O(|V|\cdot |E|)``.
""")

# ╔═╡ caa3a9dd-e7d6-493b-9586-8d08d9c8b9e5
md"""
Pour améliorer l'algorithme, pour la fusion de l'algorithm ``(1, 3)``, on peut updater le tableau à ``(1, 1, 1, 3)``.
Pour le noeud 4, on encode donc qu'il est dans la même connected component que le noeud 3 qui qui est dans la même connected component que le noeud 1. Le noeud 4 est donc lié indirectement au noeud 1 qui est appelé sa racine.
"""

# ╔═╡ fb0dc703-d81e-40e3-9cee-59744a2359c1
qa(md"Quelle est complexité après cette amélioration ?", md"""
Pour fusionner, il ne faut qu'updater une seule entrée du tableau, la racine d'un des deux noeuds de l'arête. Seulement, il pourrait faloir ``O(|V|)`` étapes pour trouver la racine. On a donc toujours ``O(|V|\cdot|E|)``.
""")

# ╔═╡ 9cd44e27-1c13-4831-a17b-bff99514194e
md"""
De façon surprenante, si on met à jour le valeur dans le tableau pour mettre directement le root après l'avoir calculé, la complexité passe à ``O(|V| \alpha(|V|))`` où ``\alpha`` est la réciproque de la [fonction d'Ackermann](https://fr.wikipedia.org/wiki/Fonction_d%27Ackermann), une fonction qui augmente extrèmement lentement donc en pratique, c'est presque ``O(|V|)``.
"""

# ╔═╡ 0d914458-7ab1-4d0e-afc6-1bd0db52076a
qa(md"Proof", md"Pour chaque noeud de degré pair, quand on arrive à ce noeud, on saura le quitter par une arête non-utilisé. Son degré d'arête non-utilisé va être diminué par deux et donc rester pair. On va alors s'arêter sur le noeud de départ si tous les degrés sont pairs. Si 2 degrés sont impairs, on doit partir d'un noeud à degré impair et on arrivera à l'autre noeud à degré impair. Ceci construit un circuit ou piste et laisse le graphe avec un degré pair d'arêtes non-utilisées à chaque noeud. On peut alors trouver un autre circuit avec la même approche et continuer à construire de nouveaux circuit jusqu'à ce que tous les degrés soient zero. Si le graphe a une seule composante connecté, on sait ensuite fusionner tous ces circuits en un seul circuit.")

# ╔═╡ 4253daed-b737-49bf-9927-6b9db39d43c1
md"""
**Définition**

* Une *forêt* (forest) est un graphe qui n'a **pas** de **cycle**.
* Un *arbre* (tree) est une forêt avec **une** seule **composante** connexe.

Par abus de language, on parle souvent d'arbre sans se soucier de si le graphe est connecté ou pas.
"""

# ╔═╡ 60e42a23-9acc-4df0-8a9d-4c5f87a74a90
md"""
Étant donné un graphe ``G(V, E)``, sa *forêt sous-tendante* (spanning forest) est une forêt ``G'(V, E')`` où ``E' \subseteq E``. Comment la trouver ?
"""

# ╔═╡ b5088341-c171-45bc-ab93-c49059da6c6a
qa(md"Quelle est la complexité de Kruskal ?", md"??")

# ╔═╡ dc4766e9-e4a7-4b8b-a6db-360df927cc5d
md"""
Dans un graphe dirigé, une arête ``(i, j)`` a un **sens**. Intuitivement, on peut aller de ``i`` vers ``j`` mais on ne peut qu'aller de ``j`` vers ``i`` si il y a aussi une autre arête de ``(j, i)``.

Que deviennent les notions de **connected components**, **spanning tree**, etc...
"""

# ╔═╡ ea6fabad-d23a-46b7-ba8c-e4822216cd4e
md"""
**Définition** Une **composante fortement connexe** d'un graphe dirigé ``G(V, E)`` est un ensemble ``V' \subseteq V`` tel qu'il existe un chemin de ``u`` vers ``v`` pour tous ``u, v \in V'``.
"""

# ╔═╡ d819bd0c-966d-46b4-88fc-8246989b070e
md"""
Si un graph contient un cycle, tous les noeuds de ce cycle sont contenus dans la même composante fortement connexe.
Si chaque composante fortement connexe est fusionnée en un noeud, il ne reste plus de cycle.
Le graphe résultat n'a donc plus de cycle, il est dit *acyclique* (DAG).

On va voir que beaucoup de problème peuvent se voir comme un problème de calcul de chaque noeud d'un graphe où les arêtes représentent les dépendences de calcul.
"""

# ╔═╡ d4cb08b6-7829-45fe-ab61-b45b39b6c3ec
md"""
```math
F_n = F_{n-1} + F_{n-2} \qquad F_2 = F_1 = 1
```
"""

# ╔═╡ 8eca3dcf-e69f-40b6-b1cd-65785af1d2c4
md"n = $slider_n"

# ╔═╡ 0a2f66ad-bff3-42e5-9f91-2f8c46ab41a7
md"Combien de **fois** `fib(3)` est appelé quand l'utilisateur appelle `fib(8)` ?"

# ╔═╡ ae970602-a2f0-472f-8e77-35b6bfd4cb39
md"n = $slider_n"

# ╔═╡ cf1fd0a7-0512-4b0f-800a-912be0c6138a
md"Nombre d'appels de `fib(k)` venant de `fib(n)`"

# ╔═╡ 76ca63ec-10cb-485c-a0f4-fc4209201032
md"""
**Définition** *Ordre topologique*: Ordre des noeuds dans lequel chaque noeud ``u`` apparait après tous les noeuds ``v`` tels qu'il y a un chemin de ``u`` vers ``v``.
"""

# ╔═╡ f9e49361-cdbc-4b77-baf6-1cdb8ed93bce
md"n = $slider_n"

# ╔═╡ 09192bc0-ac3e-47fe-8307-67f9c66d4f8e
md"n = $slider_n"

# ╔═╡ 28cfa917-cf00-4790-b71d-e629a18e4e10
md"""
* Le Puissance 4 de 7 colonnes de hauteur 6 est résolu par ordinateur en 1988
* Utilisation de l'algorithme **minimax** vu au cours S4...
* ...mais avec l'aide de la programmation dynamique!
"""

# ╔═╡ 877e5e67-b8a2-4d42-8043-47057a6df1ad
HAlign(md"Colonnes = $slider_ncols", md"Hauteur = $slider_nrows")

# ╔═╡ a252ae16-783c-4e8f-a467-1bcbf1e12565
md"Longueur = $slider_longueur" # Cols(md"Longueur = $slider_longueur", md"Gagnant patient ? $checkbox_patient")

# ╔═╡ 9e45076a-a26d-4964-a31f-da58f40f9d50
HAlign(md"Colonnes = $slider_wooclap_ncols", md"Hauteur = $slider_wooclap_nrows")

# ╔═╡ 9e3a4cb5-5ed9-42ac-b0d1-11aa5189089c
md"Quelle **structure de donnée** utiliser pour la mémoïsation : état du jeu → action optimale ?"

# ╔═╡ 7613cf26-b429-4483-bf55-bb110e49f45e
HAlign(md"Colonnes = $slider_mem_ncols", md"Hauteur = $slider_mem_nrows")

# ╔═╡ ff30250c-8844-45c5-b5fb-2f09bf05ac4e
md"Longueur = $slider_mem_longueur" # Cols(md"Longueur = $slider_mem_longueur", md"Gagnant patient ? $checkbox_mem_patient")

# ╔═╡ dc54e4b3-2cc6-4528-9490-6638a566a529
md"""
Comme on a vu avec Fibonacci et le puissance 4, beaucoup de problème de calcul peuvent se représenter par un graphe où le but est de calculer la valeur d'un noeud du graphe et chaque arête ``(u, v)`` signifie que pour calculer la valeur du noeud ``u``, il faut d'abord calculer la valeur du noeud ``v``.

Considérons le problème de calculer de combien de manières différentes ``x_n``, on sait faire une certaine somme avec des pièces de monnaie.
Par exemple, pour faire 3 € avec les pièces de 1 € et 2 €, on peut soit utiliser 3 pièces de 1 € ou 1 pièce de 2 € et une pièce de 1 € donc il y a 2 possibilitées.

Quel graphe peut-on utiliser pour modéliser ce problème ? On peut utiliser un noeud par somme en €, représentant la solution du problème pour cette somme là. Si on ne travaille pas avec les cents, et que avec les billets jusque 10 €, on a
```math
x_n = x_{n - 10} + x_{n - 5} + x_{n - 2} + x_{n - 1}
```
On a donc les arêtes ``(n, n - 10)``, ``(n, n - 5)``, ``(n, n - 5)`` et ``(n, n - 1)``.
"""

# ╔═╡ 1e6c4a4a-6c76-4f55-ac9d-bc93ccedc190
qa(md"Est-ce que cette formule donne le bon résultat ?",
md"""
Non. Ça considère plusieurs fois la même combinaison. Par exemple, ``x_3`` serait calculé comme ``x_2 + x_1 = 2 + 1 = 3`` alors que ``x_3 = 2``. Le problème c'est que ça considère que prendre une pièce de 2 € puis une pièce de 1 € est différent d'une pièce de 1 € puis une pièce de 2 €.

Il y a une symmétrie dans le problème due au fait qu'on ne se souce pas de l'ordre des pièces et on aimerait casser cette symmétrie.
Pour cela, on empêche de considérer les solutions où une pièce/billet arrive après une de moins grande valeur. Donc on ne peut pas prendre une pièce de 1 € puis une de 2 €.

Au noeud 2, on ne peut que prendre l'arête correspondant à la pièce de 2 € si on n'a pas encore pris de pièce de 1 € avant.
Seulement, la possibilité de prendre une arête ne dépend que du noeud, pas du passé des arêtes déjà prises.
On a donc besoin de considérer un graphe différent.
""")

# ╔═╡ f4562517-53ec-429c-a398-3f146a82df22
md"""
Considérons la valeur ``y_{n,k}`` qui correspond à la façon de faire la somme ``n`` avec des billets/pièces de valeur max ``k``.
"""

# ╔═╡ 174b7c95-a9f7-4cfd-a679-8000a26b39c9
qa(md"Quelle est la formule pour ``y_{n,10}`` ? Quel est le graphe ? Est-ce que la bonne valeur est calculée ?",
md"""
```math
y_{n, 10} = y_{n - 10, 10} + y_{n - 5, 5} + y_{n - 2, 2} + y_{n - 1, 1}.
```
Le graphe a donc les arêtes ``((n, 10), (n - 10, 10))``, ... La bonne valeur est maintenant calculée, on a bien ``y_{3, 2} = y_{1, 2} + y_{2, 1} = 1 + 1 = 2``.
""")

# ╔═╡ 783532d8-910e-434d-bf59-0c8dc75cfef3
md"""
**Exercice**: Implémenter le calcul avec l'approche top-down et bottom-up. Laquelle est la plus rapide ? Pourquoi ?

**Supplémentaire**: On peut aussi utiliser la formule
```math
y_{n, 10} = y_{n - 10, 5} + y_{n - 20, 5} + ... + y_{n - 5, 2} + y_{n - 10, 2} + ...
```
Modifier l'implémentation pour utiliser cette formule. Est-ce plus rapide ? Pourquoi ?
"""

# ╔═╡ d8d21aa6-36b8-4498-82e5-5fcbbd9dbe46
qa(md"Quelle est la complexité de cet algorithme ?", md"???")

# ╔═╡ 46192b0c-e856-4e96-94bc-d691d201b391
qa(md"Trouve-t-il le plus court chemin ?", md"Non, l'algorithme DFS a d'autres utilité comme le calcul de composantes fortement connexes ou tri topologique mais pas le calcul de plus court chemin. Un exception est le calcul de plus cours chemin dans un DAG. On peut alors le calculer en top-down avec memoïsation ce qui ressemble à un DFS.")

# ╔═╡ 03986809-5a0d-4029-af2c-55eaf68ff5e2
md"L'algorithme BFS résoud le problème de shortest path dans un graphe dirigé dont tous les poids sont 1."

# ╔═╡ 18e81794-8275-4f90-a18a-8052aaf34bc6
qa(md"Quelle est la complexité de l'algorithme BFS ?", md"???")

# ╔═╡ 2b1634c2-71fd-4fc3-afda-04721bef791d
md"L'algorithme de Dijkstra résoud le problème de shortest path dans un graph dirigé quand tous les poids sont positifs."

# ╔═╡ f34ffdae-4c4f-4358-9b7d-4125e82e781d
qa(md"Quelle est la complexité de l'algorithme Dijkstra ?", md"???")

# ╔═╡ e37ed884-1fd7-4fbe-ac7f-9931d595b3d3
qa(md"Un fois qu'un noeud `v` est retourné par `dequeue!`, peut-il encore satisfaire `d < dist[next]` pour `v == next` par la suite ?", md"Non, les noeuds suivant retournés par `dequeue!` auront une valeur dans `dist` plus élevée donc `d` sera encore plus élevé donc on aura `d > d[v]` durant le reste de l'algorithme.")

# ╔═╡ 61dc575d-44cb-4afb-be2b-c9a1c15e86de
qa(md"Que se peut-il se passer si une des arêtes a un poids négatif ?", md"L'algorithme peut ne jamais terminer... Il faut alors utiliser l'algorithme de Bellman-Ford (qui ne fait pas partie du cours 😢)")

# ╔═╡ 0d96847c-6ff2-4811-baec-50705344ce0e
md"""
* Le problème de Max Flow détermine le flot le plus important d'un noeud source à un noeud target étant donné une capacité maximale de chaque arête.
* Le problème de Min Cut détermine la somme des capacité minimale d'une coupe du graphe séparant la source de la target.

**Théorème** Pour tout graphe, le flot maximal est égal à la cut minimale !
"""

# ╔═╡ f0ba7d64-a360-450b-a9bb-ee6d5751e9ec
function winner(config, longueur)
	for col in _columns(config)
		row = _first_row(config)
		while true
			pos = (row, col)
			win, Δ = winner(config, longueur, pos)
			if _has_winner(win)
				return win, pos, Δ
			end
			row = _next_row(row)
			if _done(config, row)
				break
			end
		end
	end
	return _no_winner(config), (0, 0), (0, 0)
end

# ╔═╡ 92956d94-f128-46fd-9bc5-f96d6ad265fd
function winner(config, x, Δ::Tuple{Int,Int}, longueur)
	if !_inrows(config, x[1], (longueur::Int - 1) * Δ[1]) ||
	    !in(x[2] + (longueur::Int - 1) * Δ[2], _columns(config))
	#if !all(in.(x .+ (longueur - 1) .* Δ, axes(config))) 
		return _no_winner(config)
	end
	cur = config[x[1], x[2]]
	if iszero(cur[1])
		return _no_winner(config)
	end
	for i in 1:(longueur::Int - 1)
		if config[_shift(x[1], i * Δ[1]), _shift(x[2], i * Δ[2])] != cur
		#if config[(x .+ i .* Δ)...] != cur
			return _no_winner(config)
		end
	end
	return cur
end

# ╔═╡ 3b494bea-74b1-43c8-9ab8-6ad744294c2f
function play!(config, col::Int, player)
	row = _next(config, col)
	#row = @time findfirst(row -> iszero(config[row, col]), axes(config, 1))
	#if !isnothing(row)
	if !iszero(row)
		config = _assign(config, row, col, player)
	end
	return row, config
end

# ╔═╡ 2ef98325-0497-4a76-b2b4-161e663575c4
function _draw(config, longueur, δ = 50)
	colors = Dict(
		(false, false) => "white", (true, true) => "red", (true, false) => "palegreen",
		0 => "white", 1 => "palegreen", -1 => "red",
	)
	sz = (length(_columns(config)), _num_rows(config))
	c = (sz .+ 1) ./ 2 .* δ
	screen = sz .* δ
	pos(row, col) = Luxor.Point(δ * col, δ * (sz[2] - _shift_bit(row) + 1))
	Luxor.translate(-Luxor.Point(c...))
	Luxor.sethue("royalblue")
	Luxor.polysmooth(Luxor.box(Luxor.Point(c...), screen..., vertices=true), δ * 0.2, :fill)
	for col in _columns(config)
		row = _first_row(config)
		while true
			Luxor.sethue(colors[config[row, col]])
			Luxor.circle(pos(row, col), δ / 3, :fill)
			row = _next_row(row)
			if _done(config, row)
				break
			end
		end
	end
	win, x, Δ = winner(config, longueur)
	if _has_winner(win)
		Luxor.sethue(colors[win])
		#polysmooth(Luxor.box(pos((x .+ Δ .* (longueur - 1) ./ 2)...), reverse((Δ .* δ .* (longueur - 0.4)) .+ δ/5)...), 4, :fill)
		Luxor.setline(δ/10)
		Luxor.line(pos(x...), pos((_shift_bit.(x) .+ Δ .* (longueur-1))...), :stroke)
	end
	Luxor.translate(Luxor.Point(c...))
end

# ╔═╡ 74e93e11-f07d-4f8d-a5aa-7b6e38825270
function best_move_rec(config, player, longueur, depth, cache; patient, maxturn)
	if depth > maxturn + 1
		error("$depth > $maxturn + 1")
	end
	if cache !== nothing
		if haskey(cache, config)
			return convert(Tuple{Int,Int}, cache[config])
		end
		reversed = _reverse(config)
		if haskey(cache, reversed)
			return convert(Tuple{Int,Int}, cache[reversed])
		end
	end
	config = _copy(config)
	best = 0
	best_col = 0
	for col in _columns(config)
		row, new_config = play!(config, col, player)
		if !iszero(row)
			win, _ = winner(new_config, longueur, (row, col))
			if _has_winner(win)
				best = _sign(win, 2maxturn - depth)
				best_col = col
			end
		end
	end
	if iszero(best)
	for col in _columns(config)
		row, new_config = play!(config, col, player)
		if !iszero(row)
			cur, _ = best_move_rec(new_config, _next_player(player), longueur, depth + 1, cache; patient, maxturn)
			if best_col == 0 || _compare(player, cur, best)
				best = cur
				best_col = col
				if patient && _compare(player, best, 0)
					break
				end
			end
			_unassign(config, row, col)
		end
	end
	end
	if cache !== nothing
		cache[config] = (best, best_col)
	end
	return best, best_col
end

# ╔═╡ 2429f085-20f2-4941-820e-c49499331874
function game!(config, cols, longueur, memoization; patient, maxturn)
	if memoization
		cache = Dict{typeof(config),Tuple{Int8,Int8}}()
	else
		cache = nothing
	end
	player = _first_player(config)
	for i in 1:maxturn
		_, col = best_move_rec(config, player, longueur, 1, cache; patient, maxturn)
		if col == 0
			break
		end
		push!(cols, col)
		_, config = play!(config, col, player)
		w = winner(config, longueur)
		if _has_winner(w[1])
			break
		end
		player = _next_player(player)
	end
	return cols
end

# ╔═╡ 2766c3d6-6a48-4b30-bf4e-56ac9a341145
function draw(config, longueur, δ = 50)
	sz = (length(_columns(config)), _num_rows(config))
	screen = sz .* δ
	Luxor.@draw begin
		_draw(config, longueur, δ)
	end screen[1] screen[2]
end

# ╔═╡ 103f4378-34f6-46ae-9ff0-1ffd1fe4d8dd
md"""
Observations clés:
* Arbre de recherche:
  - ``\text{colonnes} = `` $wooclap_ncols sous-branches par noeuds
  - profondeur : ``\text{colonnes} \times \text{hauteur}`` donc ...
  - ... $wooclap_ncols``{}^{\text{colonnes} \times \text{hauteur}} = `` $(replace(Formatting.format(wooclap_ncols^(wooclap_nrows * wooclap_ncols), commas=true), "," => " ")) feuilles !
* Stratégie optimale dépend de la **grille** mais **pas du passé**
* Chaque case a 3 états possibles: $(draw(zeros(Int, 1, 1), 1, 20)), $(draw(ones(Int, 1, 1), 1, 20)) et $(draw(-ones(Int, 1, 1), 1, 20))
* ``3^{\text{colonnes} \times \text{hauteur}} = `` $(replace(Formatting.format(3^(wooclap_nrows * wooclap_ncols), commas=true), "," => " ")) grilles différentes...
* ... mais beaucoup de ces grilles sont invalides ou **pas explorées**
"""

# ╔═╡ 8729cfca-e8f2-4eb8-adaf-301dc0d9f061
function _init(nrows, ncols)
	z = StaticArrays.SVector{ncols}(ntuple(_ -> zero(UInt8), Val(ncols)))
	Config(z, z, (0b1 << nrows))
end

# ╔═╡ fa7956f3-822a-4d34-bb64-f94104b3d03e
game(nrows, ncols, longueur; memoization = false, patient = true) = game!(_init(nrows, ncols), Int[], longueur, memoization; patient, maxturn = nrows * ncols)

# ╔═╡ 17dc0577-35dd-4580-b2d1-958f2d3336e5
cols = @time game(nrows, ncols, longueur)

# ╔═╡ 797a4c59-6978-482e-b085-9924ea5f63bc
slider_tour = @bind tour Slider(1:length(cols), default = length(cols), show_value = true);

# ╔═╡ a2596069-63c1-4950-b9be-bf10a402b5bb
mem_cols = @time game(mem_nrows, mem_ncols, mem_longueur, memoization = true)

# ╔═╡ f0f1539a-0a11-4c81-95d8-b6375fbc9160
slider_mem_tour = @bind mem_tour Slider(1:length(mem_cols), default = length(mem_cols), show_value = true);

# ╔═╡ 5b9132b5-6638-413f-83f3-2bf1f2f9e3e6
function play(nrows, ncols, cols)
	config = _init(nrows, ncols)
	player = _first_player(config)
	for col in cols
		_, config = play!(config, col, player)
		player = _next_player(player)
	end
	return config
end

# ╔═╡ 9e114f3b-77f2-4b06-bff2-35e36c806b56
HAlign(
	md"tour = $slider_tour",
	HTML(html(draw(play(nrows, ncols, cols[1:tour]), longueur)))
)

# ╔═╡ 35cffe19-61e3-4eb1-9ec1-aeb36f91c93b
HAlign(
	md"tour = $slider_mem_tour",
	HTML(html(draw(play(mem_nrows, mem_ncols, mem_cols[1:mem_tour]), mem_longueur)))
)

# ╔═╡ 1f057543-d832-42e3-9d5d-32ba4c825e63
function minimax(nrows::Integer, ncols, longueur; vert = -200, w = 1200, h = 600, dy = 140, damping = 0.9, Δy = 55, scaling)
	config = _init(nrows, ncols)
	Luxor.@draw minimax(config, _first_player(config), longueur, 1; dy, damping, vert, Δy, scaling) w h
end

# ╔═╡ b97cff7f-26c4-4343-986c-193b5980e1c6
function width(config, player, longueur, damping, depth)
	if _has_winner(winner(config, longueur)[1])
		return damping^depth
	end
	w = 0
	for col in _columns(config)
		row, new_config = play!(config, col, player)
		if !iszero(row)
			w += width(new_config, _next_player(player), longueur, damping, depth + 1)
		end
	end
	return max(w, 1)
end;

# ╔═╡ 3098f82f-cd44-4fd8-92ae-6ee4550846fb
function minimax(config, player, longueur, depth; dy, damping, vert, Δy, scaling)
	δ = round(Int, dy * damping^depth) * scaling
	w = width(config, player, longueur, damping, depth)
	if all(iszero(config.libre))
		Luxor.translate(0, vert)
	else
		Luxor.translate(δ * w/2, 0)
	end
	_draw(config, longueur)
	Luxor.translate(-δ * w/2, 0)
	if _has_winner(winner(config, longueur)[1])
		return w
	end
	cur = 0
	for col in _columns(config)
		row, new_config = play!(config, col, player)
		if !iszero(row)
			dx = cur * δ
			Luxor.translate(dx, dy)
			win, _, _ = winner(new_config, longueur)
			score, _ = best_move_rec(new_config, _next_player(player), longueur, depth + 1, nothing, patient = true, maxturn = 30)
			dw = minimax(new_config, _next_player(player), longueur, depth + 1; dy, damping, vert, Δy, scaling)
			if _has_winner(win)
				if player
					Luxor.sethue("red")
				else
					Luxor.sethue("palegreen")
				end
			elseif score == 0
				Luxor.sethue("black")
			elseif score > 0
				Luxor.sethue("palegreen")
			else
				Luxor.sethue("red")
			end
			Luxor.setline(δ/20)
			Luxor.line(Luxor.Point(-dx+δ * w/2, -dy+Δy), Luxor.Point(δ*dw/2, -Δy), :stroke)
			Luxor.translate(-dx, -dy)
			cur += dw
		end
	end
	return w
end;

# ╔═╡ 005b2830-0510-4ef2-bb06-5a63ea5c372d
minimax(2, 2, 2, scaling = 1.7)

# ╔═╡ 9e9f07dd-edcb-471c-9fb4-613553e10353
minimax(3, 2, 2; vert = -350, w = 1400, h = 900, dy = 170, Δy = 70, scaling = 0.8, damping = 1)

# ╔═╡ f1af8644-4daf-464f-aba1-b306541508bc
mesh = Polyhedra.Mesh(projected)

# ╔═╡ 064d8f08-13b3-43f4-8492-6b9b4c7297d4
Makie.mesh(mesh)

# ╔═╡ 2d0c6a4f-f55d-41d2-89cf-501a97bcd89d
cite(args...) = bibcite(biblio, args...)

# ╔═╡ 7080a31b-82da-4ca2-9ddb-1a48c0685d68
md"""
* Voir $(cite("west2022Introduction", "Chapter 1-2"))
* ou $(cite("cormen2022Introduction", "Chapter 22-24"))
"""

# ╔═╡ 84926c9f-5ed7-4720-8153-01700575d3c5
md"""
**Théorème** $(cite("west2022Introduction", "Theorem 1.2.26"))
Considérons un graphe avec une seule composante connectée.
Il existe un circuit Eulérienne si et seulement si tous les noeuds ont un degré pair.
Si tous les noeuds ont un degré pair sauf 2 alors il n'existe pas de circuit mais une piste Eulérienne.
"""

# ╔═╡ 37b6ef7f-2b07-4189-a3fe-ae3ab1661eb9
refs(args...) = bibrefs(biblio, args...)

# ╔═╡ 912276f9-96d4-407d-bbde-c643c81a15ec
refs(["west2022Introduction", "cormen2022Introduction"])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bonito = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
DocumenterCitations = "daee34ce-89f3-4625-b898-19384cb65244"
Formatting = "59287772-0a20-5a39-b81b-1366585eb4c0"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
Bonito = "~4.0.0"
Colors = "~0.12.11"
DataStructures = "~0.18.20"
DocumenterCitations = "~1.3.5"
Formatting = "~0.4.3"
GraphPlot = "~0.6.0"
Graphs = "~1.9.0"
Luxor = "~4.1.0"
Makie = "~0.21.17"
PlutoUI = "~0.7.60"
Polyhedra = "~0.7.8"
SimpleWeightedGraphs = "~1.4.0"
StaticArrays = "~1.9.8"
WGLMakie = "~0.10.17"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "5bf96e00a2e13ab0a5de00e8daf8c8bd9034ba1f"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

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

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

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

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.BibInternal]]
deps = ["TestItems"]
git-tree-sha1 = "b3107800faf461eca3281f89f8d768f4b3e99969"
uuid = "2027ae74-3657-4b95-ae00-e2f7d55c3e64"
version = "0.3.7"

[[deps.BibParser]]
deps = ["BibInternal", "DataStructures", "Dates", "JSONSchema", "TestItems", "YAML"]
git-tree-sha1 = "33478bed83bd124ea8ecd9161b3918fb4c70e529"
uuid = "13533e5b-e1c2-4e57-8cef-cac5e52f6474"
version = "0.2.2"

[[deps.Bibliography]]
deps = ["BibInternal", "BibParser", "DataStructures", "Dates", "FileIO", "YAML"]
git-tree-sha1 = "520c679daed011ce835d9efa7778863aad6687ed"
uuid = "f1be7e48-bf82-45af-a471-ae754a193061"
version = "0.2.20"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bonito]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "LinearAlgebra", "Markdown", "MsgPack", "Observables", "RelocatableFolders", "SHA", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "262f58917d5d9644d16ec6f53480e11a6e128db2"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
version = "4.0.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+2"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "TranscodingStreams"]
git-tree-sha1 = "e7c529cc31bb85b97631b922fa2e6baf246f5905"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "c785dfb1b3bfddd1da557e861b919819b82bbe5b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.1"

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
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "bf6570a34c850f99407b494757f5d7ad233a7257"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.5"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "e1371a23fd9816080c828d0ce04373857fe73d33"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.3"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3101c32aab536e7a27b1763c0797dba151b899ad"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.113"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Documenter]]
deps = ["ANSIColoredPrinters", "AbstractTrees", "Base64", "CodecZlib", "Dates", "DocStringExtensions", "Downloads", "Git", "IOCapture", "InteractiveUtils", "JSON", "LibGit2", "Logging", "Markdown", "MarkdownAST", "Pkg", "PrecompileTools", "REPL", "RegistryInstances", "SHA", "TOML", "Test", "Unicode"]
git-tree-sha1 = "d0ea2c044963ed6f37703cead7e29f70cba13d7e"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "1.8.0"

[[deps.DocumenterCitations]]
deps = ["AbstractTrees", "Bibliography", "Dates", "Documenter", "Logging", "Markdown", "MarkdownAST", "OrderedCollections", "Unicode"]
git-tree-sha1 = "5a72f3f804deb1431cb79f5636163e4fdf8ed8ed"
uuid = "daee34ce-89f3-4625-b898-19384cb65244"
version = "1.3.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e51db81749b0777b2147fbe7b783ee79045b8e99"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.4+1"

[[deps.Extents]]
git-tree-sha1 = "81023caa0021a41712685887db1fc03db26f41f5"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.4"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2dd20384bf8c6d411b5c7370865b1e9b26cb2ea3"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.6"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "21fac3c77d7b5a9fc03b0ec503aa1a6392c34d2b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.15.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "786e968a8d2fb167f2e4880baba62e0e26bd8e4e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.3+1"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "d52e255138ac21be31fa633200b65e4e71d26802"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.6"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random", "libblastrampoline_jll"]
git-tree-sha1 = "c4f9c87b74aedf20920034bd4db81d0bffc527d2"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.3.14"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "59107c179a586f0fe667024c5eb7033e81333271"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.2"

[[deps.GeoInterface]]
deps = ["Extents", "GeoFormatTypes"]
git-tree-sha1 = "826b4fd69438d9ce4d2b19de6bc2f970f45f0f88"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.8"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0224cce99284d997f6880a42ef715a37c99338d1"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.2+0"

[[deps.Git]]
deps = ["Git_jll"]
git-tree-sha1 = "04eff47b1354d702c3a85e8ab23d539bb7d5957e"
uuid = "d7ba0133-e1db-5d97-8f8c-041e4b3a1eb2"
version = "1.3.1"

[[deps.Git_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libiconv_jll", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "399f4a308c804b446ae4c91eeafadb2fe2c54ff9"
uuid = "f8c6e375-362e-5223-8a59-34ff63f689eb"
version = "2.47.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "48b5d4c75b2c9078ead62e345966fa51a25c05ad"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.2+1"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "f76a7a0f10af6ce7f227b7a921bfe351f628ed45"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.6.0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "01979f9b37367603e2848ea225918a3b3861b606"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+1"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "dc6bed05c15523624909b3953686c5f5ffa10adc"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "6c22309e9a356ac1ebc5c8a217045f9bae6f8d9a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.13"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "b1c2585431c382e3fe5805874bda6aea90a95de9"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.25"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

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

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "LinearAlgebra", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "24c095b1ec7ee58b936985d31d5df92f9b9cfebb"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.19"
weakdeps = ["DiffRules", "ForwardDiff", "IntervalSets", "RecipesBase"]

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "1d322381ef7b087548321d3f878cb4c9bd8f8f9b"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.1"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JSONSchema]]
deps = ["Downloads", "JSON", "JSON3", "URIs"]
git-tree-sha1 = "243f1cdb476835d7c249deb9f29ad6b7827da7d3"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "1.4.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MacroTools", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays"]
git-tree-sha1 = "866dd0bf0474f0d5527c2765c71889762ba90a27"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.23.5"

    [deps.JuMP.extensions]
    JuMPDimensionalDataExt = "DimensionalData"

    [deps.JuMP.weakdeps]
    DimensionalData = "0703355e-b756-11e9-17c0-8b28908087d0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "0f2da712350b020bc3957f269c9caad516383ee0"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.3.0"

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
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6ce1e19f3aec9b59186bdf06cdf3c4fc5f5f3e6"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.50.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "61dfdba58e585066d8bce214c5a51eaa0539f269"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "84eef7acd508ee5b3e956a2ae51b05024181dee0"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.2+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "edbf5309f9ddf1cab25afc344b1e8150b7c832f9"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "PolygonAlgorithms", "PrecompileTools", "Random", "Rsvg"]
git-tree-sha1 = "134570038473304d709de27384621bd0810d23fa"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "4.1.0"
weakdeps = ["LaTeXStrings", "MathTeXEngine"]

    [deps.Luxor.extensions]
    LuxorExtLatex = ["LaTeXStrings", "MathTeXEngine"]

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "260d6e1ac8abcebd939029e6eedeba4e3870f13a"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.21.17"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "b774d0563bc332f64d136d50d0420a195d9bdcc6"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.8.11"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MarkdownAST]]
deps = ["AbstractTrees", "Markdown"]
git-tree-sha1 = "465a70f0fc7d443a00dcdc3267a497397b8a3899"
uuid = "d0879d2d-cac2-40c8-9cee-1863dc0c7391"
version = "0.1.2"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "e065ca5234f53fd6f920efaee4940627ad991fb4"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.34.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "f45c8916e8385976e1ccd055c9874560c257ab13"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.2"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

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
version = "2023.12.12"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "a2710df6b0931f987530f59427441b21245d8f5e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "39d000d9c33706b8364817d8894fae1548f40295"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.2"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

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

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PolygonAlgorithms]]
git-tree-sha1 = "a5ded6396172cff3bacdd1354d190b93cb667c4b"
uuid = "32a0d02f-32d9-4438-b5ed-3a2932b48f96"
version = "0.2.0"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polyhedra]]
deps = ["GenericLinearAlgebra", "GeometryBasics", "JuMP", "LinearAlgebra", "MutableArithmetics", "RecipesBase", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "e2e81cdf047d2921d99f7ff0b6dece45b3195185"
uuid = "67491407-f73d-577b-9b50-8179a7c68029"
version = "0.7.8"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "52af86e35dd1b177d051b12681e1c581f53c281b"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "4b33e0e081a825dbfaf314decf58fa47e53d6acb"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "777657803913ffc7e8cc20f0fd04b634f871af8f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.8"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "9537ef82c42cdd8c5d443cbc359110cbb36bae10"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.21"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

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

[[deps.TestItems]]
git-tree-sha1 = "42fd9023fef18b9b78c8343a4e2f3813ffbcefcb"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.0.0"

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "0248b1b2210285652fbc67fd6ced9bf0394bcfec"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.1"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

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
git-tree-sha1 = "01915bfcd62be15329c9a07235447a89d588327c"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.1"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.WGLMakie]]
deps = ["Bonito", "Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "db71caa2e1ac6b3f806333c9de32393ed75d60e6"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.10.17"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "a2fccc6559132927d4c5dc183e3e01048c6dcbd6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "15e637a697345f6743674f1322beefbc5dcd5cfc"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+1"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2b0e27d52ec9d8d483e2ca0b72b3cb1a8df5c27a"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+1"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "02054ee01980c90297412e4c809c8694d7323af3"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+1"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d7155fea91a4123ef59f42c4afb5ab3b4ca95058"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+1"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee57a273563e273f0f53275101cd41a8153517a"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+1"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "1a74296303b6524a0472a8cb12d3d87a78eb3612"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+1"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b9ead2d2bdb27330545eb14234a2e300da61232e"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+1"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "dea63ff72079443240fbd013ba006bcbc8a9ac00"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.12"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "86e7731be08b12fa5e741f719603ae740e16b666"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.10+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "7dfa0fd9c783d3d0cc43ea1af53d69ba45c447df"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "ccbb625a89ec6195856a50aa2b668a5c08712c94"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.4.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

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
"""

# ╔═╡ Cell order:
# ╟─c9595f48-bfec-443e-8f8f-94ed0c5a9530
# ╟─7080a31b-82da-4ca2-9ddb-1a48c0685d68
# ╟─912276f9-96d4-407d-bbde-c643c81a15ec
# ╟─15f7cc14-7663-4ad2-9cfe-9ee8aada601c
# ╟─faca26cb-14a3-4a45-ba0b-ec3cd57d8d29
# ╟─b65dd471-10da-4b3a-a14b-cfbb70e70050
# ╟─8bec0324-274c-4ae6-99f9-08b613460064
# ╟─d7bbc61f-5cc6-4090-8cc4-313187f28bab
# ╟─b52b2192-df63-4095-aa15-6d3ccafc92ba
# ╟─c87e1a66-62e6-41aa-bcc0-fb1704d7e259
# ╟─064d8f08-13b3-43f4-8492-6b9b4c7297d4
# ╟─55283461-b2e6-4e81-884f-48ef81e605d4
# ╟─73642ba1-918f-450c-a555-7bcd393ec54e
# ╟─e41498fa-8fed-4354-bb3b-a2fc8afeacd9
# ╟─27ee39a0-0eb2-42dc-bf90-eb2fa0d7497e
# ╟─c41901fb-6570-487d-bd77-4708d1f0b627
# ╟─ff87b1bb-7903-4954-8e7a-f377aef97b2d
# ╠═fe642397-9ae1-42b5-9d0e-d2b96cfb178b
# ╟─043ff87b-41e4-49c2-aaee-bee39f088b58
# ╟─2ebc96d9-e526-408a-a264-a87771b96258
# ╟─9e065cf0-7483-4ec2-809c-cdd7d22a5564
# ╟─cbc8da76-9adf-4670-b07c-0f5366b55547
# ╟─338527b3-ec32-4794-8b92-63cce6ce4285
# ╟─269320fe-a575-4fe4-86b6-67b0dd4484a5
# ╟─b184882e-f2ed-4e9c-b763-0173711d6fac
# ╟─caa3a9dd-e7d6-493b-9586-8d08d9c8b9e5
# ╟─fb0dc703-d81e-40e3-9cee-59744a2359c1
# ╟─9cd44e27-1c13-4831-a17b-bff99514194e
# ╟─8e7b6a51-1f8e-43d7-846b-a9e7596052b9
# ╟─84926c9f-5ed7-4720-8153-01700575d3c5
# ╟─0d914458-7ab1-4d0e-afc6-1bd0db52076a
# ╟─d884c51e-9ddd-4dc4-9f9c-e461f0455a17
# ╟─4253daed-b737-49bf-9927-6b9db39d43c1
# ╠═ceb5aea4-5a69-47cf-bedd-17d47a6e07a8
# ╟─b94745d7-ffda-4190-82af-8ee2fe9ea165
# ╟─60e42a23-9acc-4df0-8a9d-4c5f87a74a90
# ╠═944b6ad9-731b-4479-94bd-e61876987917
# ╠═d1329a48-5bed-428f-8554-854a5d51488e
# ╟─3c09769f-13b2-4713-8d2e-06abd7735c33
# ╠═3eb8d016-d1a1-4e4f-a307-e2c439c502e4
# ╟─ecf6a7d0-3279-4bff-b0f4-73327f60b702
# ╟─b5088341-c171-45bc-ab93-c49059da6c6a
# ╟─8c384e86-18d4-419b-977b-c5d5bfd1318f
# ╟─dc4766e9-e4a7-4b8b-a6db-360df927cc5d
# ╠═49092fe9-7c24-4564-9d2a-aa313d54ed8a
# ╟─264eb1b6-0e3f-4d67-b8b0-5f0831d18b95
# ╟─ea6fabad-d23a-46b7-ba8c-e4822216cd4e
# ╠═730e5e09-e818-49c5-8fa1-495b39b4267e
# ╟─66c7051b-c377-4376-a472-b2bd5c788267
# ╟─d819bd0c-966d-46b4-88fc-8246989b070e
# ╟─71566570-7d26-4fba-8f9a-8726fb97f2d8
# ╟─d4cb08b6-7829-45fe-ab61-b45b39b6c3ec
# ╠═2f905018-1d7b-43dc-b576-e4080dc946ee
# ╟─8eca3dcf-e69f-40b6-b1cd-65785af1d2c4
# ╠═b53e5719-cc0e-4169-9c8d-3eb5a50570b0
# ╟─0a2f66ad-bff3-42e5-9f91-2f8c46ab41a7
# ╟─a95116d6-614a-4a8a-b3a6-5a96e9b3c439
# ╠═d35df21e-2dd5-478f-b0c1-17cbd96a9cc6
# ╟─3c82aba6-0539-4578-ad04-648eb3bcf400
# ╟─ae970602-a2f0-472f-8e77-35b6bfd4cb39
# ╟─cf1fd0a7-0512-4b0f-800a-912be0c6138a
# ╠═46977b09-a43a-427c-8974-a1fd5f2d7da2
# ╟─6183e415-898a-4d1f-a5db-9c09503affff
# ╟─30861bc6-4af9-43dc-a1f5-fed9bb157564
# ╟─76ca63ec-10cb-485c-a0f4-fc4209201032
# ╟─f9e49361-cdbc-4b77-baf6-1cdb8ed93bce
# ╠═9b024507-48e2-4247-aa4e-ed0958af7b2d
# ╠═4c76167b-a8ac-436f-87a8-d493d489e61f
# ╟─a9379cc3-6ca8-4aed-81f7-992f3664ceb3
# ╟─09192bc0-ac3e-47fe-8307-67f9c66d4f8e
# ╠═c395343c-3f7c-4dfe-bc0f-8950afc115bf
# ╠═73c0c135-f070-422e-ab1d-20c6faa160d8
# ╠═9eeaca05-9e6e-4881-85aa-d81ad5460e8c
# ╟─16c5f723-77d4-4bee-86c6-90a1b597e7d9
# ╟─28cfa917-cf00-4790-b71d-e629a18e4e10
# ╟─54d0cbaf-4bf7-4508-b706-b5042dda1644
# ╠═005b2830-0510-4ef2-bb06-5a63ea5c372d
# ╟─73d104ac-7aff-46cc-95b5-a67b6811f64d
# ╟─9e9f07dd-edcb-471c-9fb4-613553e10353
# ╟─1308d56a-f91a-4fe2-8f5d-7d7a15b8f040
# ╟─877e5e67-b8a2-4d42-8043-47057a6df1ad
# ╟─a252ae16-783c-4e8f-a467-1bcbf1e12565
# ╠═17dc0577-35dd-4580-b2d1-958f2d3336e5
# ╟─9e114f3b-77f2-4b06-bff2-35e36c806b56
# ╟─797a4c59-6978-482e-b085-9924ea5f63bc
# ╟─3de604e8-20b6-4ef5-9952-b8f27c28fef1
# ╟─a7530e0a-a3b5-46c8-9666-fb0219d7ba04
# ╟─9e45076a-a26d-4964-a31f-da58f40f9d50
# ╟─103f4378-34f6-46ae-9ff0-1ffd1fe4d8dd
# ╟─9e3a4cb5-5ed9-42ac-b0d1-11aa5189089c
# ╟─0a7e65f5-da36-4331-8463-221051286728
# ╟─3f813480-726b-42e0-a648-e3ad76abcf9d
# ╟─7613cf26-b429-4483-bf55-bb110e49f45e
# ╠═a2596069-63c1-4950-b9be-bf10a402b5bb
# ╟─ff30250c-8844-45c5-b5fb-2f09bf05ac4e
# ╟─35cffe19-61e3-4eb1-9ec1-aeb36f91c93b
# ╟─3b82fa0a-65b1-46dd-b81a-d1fa117893d5
# ╟─f0f1539a-0a11-4c81-95d8-b6375fbc9160
# ╟─a436f564-8418-477c-9b0d-9d791c893bb1
# ╟─dc54e4b3-2cc6-4528-9490-6638a566a529
# ╟─1e6c4a4a-6c76-4f55-ac9d-bc93ccedc190
# ╟─f4562517-53ec-429c-a398-3f146a82df22
# ╟─174b7c95-a9f7-4cfd-a679-8000a26b39c9
# ╟─783532d8-910e-434d-bf59-0c8dc75cfef3
# ╟─9db78023-e29c-41e5-9102-870692111f39
# ╠═6dc95039-9e2a-47e0-8250-9c5ad8a2af0c
# ╠═2829213a-8533-4e92-a130-57d53166b376
# ╠═f8d160f8-7e41-4bdf-99f4-d28b344a46a4
# ╠═b54d9e8b-83ab-434e-9e6c-f6774558abad
# ╟─d8d21aa6-36b8-4498-82e5-5fcbbd9dbe46
# ╟─46192b0c-e856-4e96-94bc-d691d201b391
# ╟─2c7c9b52-0d27-4dc8-ae69-29a6edc87041
# ╟─03986809-5a0d-4029-af2c-55eaf68ff5e2
# ╠═fb29e166-b361-49e3-bc5c-042d409023f1
# ╠═c7663ad1-16b4-4b29-8f4a-8993d094a7e1
# ╟─18e81794-8275-4f90-a18a-8052aaf34bc6
# ╟─670aaab7-a101-4dea-ba97-80f1a5346122
# ╟─2b1634c2-71fd-4fc3-afda-04721bef791d
# ╠═ab3997f3-83ef-4435-bb90-7cef1c7bcdcf
# ╠═0ee692b3-fcf6-418a-b3ba-39febcf6eed2
# ╟─f34ffdae-4c4f-4358-9b7d-4125e82e781d
# ╟─e37ed884-1fd7-4fbe-ac7f-9931d595b3d3
# ╟─61dc575d-44cb-4afb-be2b-c9a1c15e86de
# ╟─634dea01-3832-49e7-90af-4e1b834cead7
# ╟─0d96847c-6ff2-4811-baec-50705344ce0e
# ╟─75bb5e30-76b8-486b-b3fa-9979aaf830a5
# ╠═f83fdf0a-a2f5-4d72-9297-a409cf4e2bbb
# ╠═accdf77c-fc4b-4777-bebf-9d69ebb36526
# ╠═bd9ab2c4-b098-11ef-3c1e-7fe1458a8556
# ╠═b01e3f5a-ad9f-4557-ba62-30e9c57b8532
# ╠═fdd9318e-ee09-4eab-9322-ad97b03da002
# ╠═039d5065-4ab2-4787-86e2-5612cfcd01b4
# ╠═3c9dd530-8a25-4c1b-b1ff-3b5a3748b74a
# ╠═33ac9568-a83d-47b0-a17c-6e30df0dbaa4
# ╠═1b748b69-5262-46c9-aeff-01696c585e6f
# ╠═d57a6a11-0387-4693-a50e-9eba62605a49
# ╠═2fedb811-6f34-424f-a46b-806bb0c6b9c0
# ╠═12fa3c11-c5d3-4c81-a92f-174f09ec10ff
# ╠═3ca057ef-1367-4fef-8e18-bdef855bc35e
# ╠═b107f529-f47f-4252-b06a-6e53302be860
# ╠═e2f5854a-c6ce-4600-8e96-62ac8ffc6b6e
# ╠═938e0359-0a86-4755-a92b-f1a7a892961d
# ╠═f0ba7d64-a360-450b-a9bb-ee6d5751e9ec
# ╠═1c762a4a-4ee1-46ee-93ee-7b566376d5a4
# ╠═92956d94-f128-46fd-9bc5-f96d6ad265fd
# ╠═8dea1779-e432-4fe3-baf0-5b0d8d3dd97b
# ╠═3b494bea-74b1-43c8-9ab8-6ad744294c2f
# ╠═2ef98325-0497-4a76-b2b4-161e663575c4
# ╠═efbdb43d-3c6b-48d6-be40-51369f3158de
# ╠═74e93e11-f07d-4f8d-a5aa-7b6e38825270
# ╠═fa7956f3-822a-4d34-bb64-f94104b3d03e
# ╠═2429f085-20f2-4941-820e-c49499331874
# ╠═2766c3d6-6a48-4b30-bf4e-56ac9a341145
# ╠═5b9132b5-6638-413f-83f3-2bf1f2f9e3e6
# ╠═8729cfca-e8f2-4eb8-adaf-301dc0d9f061
# ╠═1f057543-d832-42e3-9d5d-32ba4c825e63
# ╠═b97cff7f-26c4-4343-986c-193b5980e1c6
# ╠═3098f82f-cd44-4fd8-92ae-6ee4550846fb
# ╠═f1af8644-4daf-464f-aba1-b306541508bc
# ╠═2d0c6a4f-f55d-41d2-89cf-501a97bcd89d
# ╠═37b6ef7f-2b07-4189-a3fe-ae3ab1661eb9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
