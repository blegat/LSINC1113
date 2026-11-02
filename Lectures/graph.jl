### A Pluto.jl notebook ###
# v1.0.4

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

# ╔═╡ 8e31c397-10cd-4ef3-85ce-30b1bb888135
md"""
Kruskal est un algorithm glouton (ou greedy). Il prend les arêtes de poids minimales successivement tant qu'elles ne créent pas de cycle.
Bien qu'un algorithme greedy n'arrive souvent pas à la solution optimale, on peut prouver que Kruskal trouve systématiquement l'optimal.
"""

# ╔═╡ b5088341-c171-45bc-ab93-c49059da6c6a
qa(md"Quelle est la complexité de Kruskal ?", md"``O(|E|\log(|E|))`` to sort the edges then ``O(|E|\alpha(|V|))`` for the selection of edge with the Disjoint-Set datastructure.")

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
  - ... $wooclap_ncols``{}^{\text{colonnes} \times \text{hauteur}} = `` $(replace(Formatting.format(big(wooclap_ncols)^(wooclap_nrows * wooclap_ncols), commas=true), "," => " ")) feuilles !
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
Bonito = "~5.2.0"
Colors = "~0.13.4"
DataStructures = "~0.19.6"
DocumenterCitations = "~1.5.0"
Formatting = "~0.4.3"
GraphPlot = "~0.6.2"
Graphs = "~1.15.0"
Luxor = "~4.5.0"
Makie = "~0.24.17"
PlutoUI = "~0.7.84"
Polyhedra = "~0.9.0"
SimpleWeightedGraphs = "~1.6.0"
StaticArrays = "~1.9.24"
WGLMakie = "~0.13.17"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.3"
manifest_format = "2.1"
project_hash = "6b071a3a66b745b2a94d8cc1d71f45ce60523cbc"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
registries = "General"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "8c0bc12397dabf5d2dfc4b5763ec463854cc44d3"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.46"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "dee8388a9bec89c8713ad176f96960ae2f98598a"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.4"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptiveArrayPools]]
deps = ["Preferences", "Printf", "Random"]
git-tree-sha1 = "86d0ae34bc46373170064ef33c12a303994ab736"
registries = "General"
uuid = "4f381ef7-9af0-4cbe-99d4-cf36d7b0f233"
version = "0.3.12"

    [deps.AdaptiveArrayPools.extensions]
    AdaptiveArrayPoolsCUDAExt = "CUDA"
    AdaptiveArrayPoolsMetalExt = "Metal"

    [deps.AdaptiveArrayPools.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
registries = "General"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
registries = "General"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
registries = "General"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
registries = "General"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
registries = "General"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.BibInternal]]
deps = ["TestItems"]
git-tree-sha1 = "ec1f1c78c1a6c917125abdee98518685fe653e64"
registries = "General"
uuid = "2027ae74-3657-4b95-ae00-e2f7d55c3e64"
version = "0.4.0"

[[deps.BibParser]]
deps = ["BibInternal", "DataStructures", "TestItems"]
git-tree-sha1 = "5a6aeb593475157a5911cda3e0297f1552d45102"
registries = "General"
uuid = "13533e5b-e1c2-4e57-8cef-cac5e52f6474"
version = "0.3.0"

    [deps.BibParser.extensions]
    BibParserCFFExt = ["Dates", "JSONSchema", "YAML"]
    BibParserCSLExt = "JSON3"
    BibParserXMLExt = "EzXML"

    [deps.BibParser.weakdeps]
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    EzXML = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
    JSONSchema = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
    YAML = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"

[[deps.Bibliography]]
deps = ["BibInternal", "BibParser", "DataStructures", "Dates", "FileIO", "JSONSchema", "TestItems", "YAML"]
git-tree-sha1 = "2246611f8da01251b73665fe0d57ed9d0e3b880b"
registries = "General"
uuid = "f1be7e48-bf82-45af-a471-ae754a193061"
version = "0.4.0"

    [deps.Bibliography.extensions]
    BibliographyCSLExt = "JSON3"
    BibliographyXMLExt = "EzXML"

    [deps.Bibliography.weakdeps]
    EzXML = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
registries = "General"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.Bonito]]
deps = ["ANSIColoredPrinters", "Base64", "CodecZlib", "Colors", "CommonMark", "Dates", "Deno_jll", "HTTP", "Hyperscript", "JSON", "LinearAlgebra", "Markdown", "MbedTLS", "MsgPack", "Observables", "OrderedCollections", "PrecompileTools", "Random", "RelocatableFolders", "SHA", "Scratch", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "f947a8b5967714bf5b7d82070b647631de1f10d9"
registries = "General"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
version = "5.2.0"
weakdeps = ["Documenter", "MarkdownAST"]

    [deps.Bonito.extensions]
    BonitoDocumenterExt = ["Documenter", "MarkdownAST"]

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "0db159493e83ecb2d122f9c40138d89b70d43963"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.1.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
registries = "General"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "0cfb779e24fdced77a05df1dbe99200821a9f563"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.8+1"

[[deps.ChainRulesCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "295fa9b898db109bc56df62c0fc20bd105cec00e"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.8"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "TranscodingStreams"]
git-tree-sha1 = "84990fa864b7f2b4901901ca12736e45ee79068c"
registries = "General"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "970758a3d591a2a5c2a907c53f2e2f8c1b1d3537"
registries = "General"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.9"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
registries = "General"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
registries = "General"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "060adaf7a5c2e14bca025edf0eea0206ffdebd87"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.32.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "af523e59de94de05b9b367a2ca4d63233aacd222"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.5"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "96fce5c80e0834e96ab6d24fb384fd945c9cb2d4"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.2"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "0737a3e59f9d94957f5b177fa2ca40901636b85b"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.4"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "7c8fe02c7eb6fe22e89d3990123a2a931ca8fcfb"
registries = "General"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "1.0.4"
weakdeps = ["Markdown", "MarkdownAST"]

    [deps.CommonMark.extensions]
    CommonMarkMarkdownASTExt = "MarkdownAST"
    CommonMarkMarkdownExt = "Markdown"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
registries = "General"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "d75431a71f82758e218779ccb3369f3243fd2bc1"
registries = "General"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.7"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "b56df0c58a9d37b4a926f84f3f4f283a5987680e"
registries = "General"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.9"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
registries = "General"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "e1e28f37cbdca4b594ad86ecfe2ee8f25796ebff"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.7.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
registries = "General"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
registries = "General"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
registries = "General"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
registries = "General"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "2ddbc341251ed713d39e503b504593a5d6789551"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.132"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Documenter]]
deps = ["ANSIColoredPrinters", "AbstractTrees", "Base64", "CodecZlib", "Dates", "DocStringExtensions", "Downloads", "Git", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "MarkdownAST", "Pkg", "PrecompileTools", "REPL", "RegistryInstances", "SHA", "TOML", "Test", "Unicode"]
git-tree-sha1 = "191e6bef0cf32cac3a3913cd4787d851715254c2"
registries = "General"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "1.19.0"

[[deps.DocumenterCitations]]
deps = ["AbstractTrees", "Bibliography", "Bijections", "Dates", "Documenter", "Logging", "Markdown", "MarkdownAST", "OrderedCollections", "SHA", "Unicode"]
git-tree-sha1 = "5699464c66535cb99ad6ab075513cb215217b9f2"
registries = "General"
uuid = "daee34ce-89f3-4625-b898-19384cb65244"
version = "1.5.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
registries = "General"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
registries = "General"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0f666deaf8099a7d7ac229a74402ada88898be1"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.9.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "7bd13840b4148949e290748071bbe6826b3eabcd"
registries = "General"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d9d3cd382f2c999f684e6bfa4039fcc16c153b6f"
registries = "General"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "9.0.2+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "16aa62b6b2c5c65a642155bd1bce90633b8d3130"
registries = "General"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.2"

[[deps.FastInterpolations]]
deps = ["AdaptiveArrayPools", "HelpPlots", "LinearAlgebra", "Preferences", "Printf"]
git-tree-sha1 = "f7bcd496f0a2e4e1174e7e1f667f27fb671fd7df"
registries = "General"
uuid = "9ea80cae-fc13-4c00-8066-6eaedb12f34b"
version = "0.4.20"

    [deps.FastInterpolations.extensions]
    FastInterpolationsChainRulesCoreExt = "ChainRulesCore"
    FastInterpolationsColorVectorSpaceExt = ["ColorTypes", "ColorVectorSpace"]
    FastInterpolationsEnzymeExt = "Enzyme"
    FastInterpolationsForwardDiffExt = "ForwardDiff"
    FastInterpolationsRecipesBaseExt = "RecipesBase"
    FastInterpolationsSymbolicsExt = ["Symbolics", "SymbolicUtils"]

    [deps.FastInterpolations.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ColorTypes = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
    ColorVectorSpace = "c3611d14-8923-5661-9e6a-0046d554d3a4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
registries = "General"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
registries = "General"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
registries = "General"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "831f86fde78ab7f6470d6ba4d1161a11fb6b4ec0"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.18.2"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "e1e9ba3a8f5d71e6ae698f138a1362aeabb42e9f"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.18.3+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
registries = "General"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "3b0f72e2ffef1a139ac450725c9ecd01c0a3c050"
registries = "General"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.4.6"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
registries = "General"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d95a7a3a7e0f2fffc6291bd33035830aab01310b"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.15.0+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
registries = "General"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e894be022a2eba58463939f1bf1f7f7006f27f52"
registries = "General"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.19+0"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random", "libblastrampoline_jll"]
git-tree-sha1 = "3b0e97a0de86e54fa897c6e431c179ad4cd1dea0"
registries = "General"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.4.2"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "ec46c5825710fa1a15d468acb2d93cc939a7a5fe"
registries = "General"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.13"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
registries = "General"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3efbbc027441271444dcd0c0a46f2d119dc4329"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "6.1.3+0"

[[deps.Git]]
deps = ["Git_LFS_jll", "Git_jll", "JLLWrappers", "OpenSSH_jll"]
git-tree-sha1 = "824a1890086880696fc908fe12a17bcf61738bd8"
registries = "General"
uuid = "d7ba0133-e1db-5d97-8f8c-041e4b3a1eb2"
version = "1.5.0"

[[deps.Git_LFS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d310b22511a51569401e975194cc76077bcd3786"
registries = "General"
uuid = "020c3dae-16b3-5ae5-87b3-4cb189e250b2"
version = "3.8.1+0"

[[deps.Git_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libiconv_jll", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "3842d761162baf8efce819cdfb5ea070fe44ac40"
registries = "General"
uuid = "f8c6e375-362e-5223-8a59-34ff63f689eb"
version = "2.56.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "ea9f1fb37ff10597353b56f3ba0ebd525acef6e0"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.90.1+0"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "066c87e33a8fcc3518c9e9970a1cbf85aa79fd6c"
registries = "General"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.6.2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "dbb2b976f605b220dfe139d6db9f3859680da09e"
registries = "General"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.15.0"
weakdeps = ["Distributed", "SharedArrays"]

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "ef70da5e123a06a29e2d6ddff0f09985bc226491"
registries = "General"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.3"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "EnumX", "Logging", "PrecompileTools", "Random", "Reseau", "SHA", "TranscodingStreams", "URIs", "UUIDs", "Zlib_jll"]
git-tree-sha1 = "779aaa01f905311a6a2c21b04517da84efb0be3c"
registries = "General"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "2.14.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f700c32c1e9e11b576592043f270e6b33868a56f"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.15000.0+0"

[[deps.HelpPlots]]
git-tree-sha1 = "07d389d81ca31e42d68bd88cfa1e900a63b7d378"
registries = "General"
uuid = "5880b2e8-ad0a-48c9-bacd-30dc9c1357e0"
version = "1.1.0"

    [deps.HelpPlots.extensions]
    HelpPlotsExt = "Plots"

    [deps.HelpPlots.weakdeps]
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "c28dceeadc1c9c699656a3ac9d63f79cfaa17472"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.31"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "f0f005f997dfb8c5fe23920d99458a9619873893"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.10"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "1c531bf0f8a5c60a340926e058fd3f209b5eef5d"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.12"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticMakieExt = "Makie"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "0db6aea5b64caa1c47e9fdd27394f0296f05e8bf"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.15"
weakdeps = ["Makie", "Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsMakieExt = "Makie"
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "5cc6ad3e2bfc75a11eed92632ee341d37f4a4536"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.18"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "d22e8bdb3f03c300a103ff23c0dd5a654c6701de"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.7"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
registries = "General"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "123311cf9834c6ee6a8347bba9af68f255615391"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.11.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSONSchema]]
deps = ["Downloads", "JSON", "URIs"]
git-tree-sha1 = "5864ebcd5d9c6966c35bc29e536c8b40d713c00c"
registries = "General"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "1.6.0"

    [deps.JSONSchema.extensions]
    JSONSchemaJSON3Ext = "JSON3"

    [deps.JSONSchema.weakdeps]
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e06b451acfb1c6fa5ea753b5042207226271f455"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.1+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
registries = "General"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
registries = "General"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d7c65411dda4bfa0beb7490476122b4f805cf31"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.3.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "0f2da712350b020bc3957f269c9caad516383ee0"
registries = "General"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.3.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.22.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.7+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b8a91037379e1fd661e9b5ff3d7a770424b1f9f2"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.4+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "Cairo_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pango_jll", "XML2_jll", "gdk_pixbuf_jll"]
git-tree-sha1 = "e6ab5dda9916d7041356371c53cdc00b39841c31"
registries = "General"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.7+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9e887184f306b28dbfcff88cefbb0ec80f1aa9e0"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.4+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "28d620bcdb37a3fb02de9f521928ee453690f414"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.3"

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

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "PolygonAlgorithms", "PrecompileTools", "Random", "Rsvg"]
git-tree-sha1 = "fe8060b3d693f682e14f1019b058c64effb62b43"
registries = "General"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "4.5.0"

    [deps.Luxor.extensions]
    LuxorExtLatex = ["LaTeXStrings", "MathTeXEngine"]
    LuxorExtTypstry = ["Typstry"]

    [deps.Luxor.weakdeps]
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    MathTeXEngine = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FastInterpolations", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "8ff1cd11162828e9f517e32d83b49d2af63391c6"
registries = "General"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.17"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MarkdownAST]]
deps = ["AbstractTrees", "Markdown"]
git-tree-sha1 = "93c718d892e73931841089cdc0e982d6dd9cc87b"
registries = "General"
uuid = "d0879d2d-cac2-40c8-9cee-1863dc0c7391"
version = "0.1.3"

[[deps.MathOptInterface]]
deps = ["CodecBzip2", "CodecZlib", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test"]
git-tree-sha1 = "637aa0115265323850b125287833341a6f7b9b54"
registries = "General"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.54.0"

    [deps.MathOptInterface.extensions]
    MathOptInterfaceBenchmarkToolsExt = "BenchmarkTools"
    MathOptInterfaceCliqueTreesExt = "CliqueTrees"

    [deps.MathOptInterface.weakdeps]
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
registries = "General"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "0e4b139a1f0f0b8c054fcdfd3be5f610cfb8c787"
registries = "General"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.2.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
registries = "General"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
registries = "General"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "d9d20d6eb2a3290ba5c70a1c013b45305cc39c9e"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.10.20"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
registries = "General"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ffd1161361638e5f0bb796a095b0fdc230d6f019"
registries = "General"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.9.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "676c0f1c1c34224d6869675047a54b42fdc7b8c5"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.5"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
registries = "General"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
registries = "General"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "OpenJPH_jll", "Zlib_jll", "Zstd_jll", "libdeflate_jll"]
git-tree-sha1 = "afa66dc1246b28b68600d053715e00a379268980"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.5.3+0"

[[deps.OpenJPH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "1b1f3cabb1130670812b848f4623d69d578ea805"
registries = "General"
uuid = "95196869-6d61-5868-8eee-c38c33804e0c"
version = "0.32.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll"]
git-tree-sha1 = "1f640b82a627df3dbbaaff756a201a9ceb10dd1f"
registries = "General"
uuid = "9bd350c2-7e96-507f-8002-3f2e150b4e1b"
version = "10.6.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.8+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
registries = "General"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "b079a6affebde167dad8ffe0d670b53917fbcb68"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.47.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
registries = "General"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
registries = "General"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "2facd7b5863092567d54590997caf8f6db8e0b1e"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.11"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "ddcc36fd83fb4cd459a53c891981697bde0d3de6"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "9869376035f47756b8b8e711fd0d06a06d51da6f"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.84"

[[deps.PolygonAlgorithms]]
git-tree-sha1 = "c1092ada65e6d59d6361d5086ddb0a5ea63ae204"
registries = "General"
uuid = "32a0d02f-32d9-4438-b5ed-3a2932b48f96"
version = "0.4.0"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
registries = "General"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polyhedra]]
deps = ["GenericLinearAlgebra", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "7d4df6c9c039abc28c1e9ab4e2e20f638b9eb918"
registries = "General"
uuid = "67491407-f73d-577b-9b50-8179a7c68029"
version = "0.9.0"

    [deps.Polyhedra.extensions]
    PolyhedraGeometryBasicsExt = "GeometryBasics"
    PolyhedraJuMPExt = "JuMP"
    PolyhedraRecipesBaseExt = "RecipesBase"

    [deps.Polyhedra.weakdeps]
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "aef1987d90dcc2e2cbd919ad10b2d25bcce4f418"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Printf"]
git-tree-sha1 = "ce5ba3edfdeba107a721748d53b5f5887f415853"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.13.0"

    [deps.ProgressMeter.extensions]
    ProgressMeterDistributedExt = "Distributed"
    ProgressMeterIJuliaExt = "IJulia"

    [deps.ProgressMeter.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0c0dcffc464cda2ff33fc8ffc1dd1bdd53cd98be"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.4.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
registries = "General"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
registries = "General"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Reseau]]
deps = ["NetworkOptions", "OpenSSL_jll", "PrecompileTools", "Random", "SHA"]
git-tree-sha1 = "13d0465bbd89b45183c0363b635dc82d5843593f"
registries = "General"
uuid = "802f3686-a58f-41ce-bb0c-3c43c75bba36"
version = "1.6.6"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "c5e8a0a52c0d5d043214de04f11f4d9c2c6ccdb8"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.13"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "e53dad0507631c0b8d5d946d93458cbabd0f05d7"
registries = "General"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.1.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "37ac518b3c924f110d82d50173eaf4ebeaf27ad7"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.8.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "57aa595158717ef165e6f5ab639fe2e3178c0a2b"
registries = "General"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
registries = "General"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "93f0ac4ed92a34bf859f1e64813a9f41925a3447"
registries = "General"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.6.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2c79185e5261b159474903598c571ddd9f12ff7b"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "88ef0d98204515bc53dd4aaea5f7788882d1b66c"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "8a8a52f887f38d44014e24ee46b8a7311c625f6d"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.24"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
registries = "General"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "4daaddbe37aa80b4551bf7c56adf53e7cfcdf35a"
registries = "General"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.4"

    [deps.StructArrays.extensions]
    AdaptExt = "Adapt"
    GPUArraysCoreExt = "GPUArraysCore"
    KernelAbstractionsExt = "KernelAbstractions"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"
    StaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestItems]]
git-tree-sha1 = "a6dd904babc04670c784f81b67957a3545271610"
registries = "General"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.1.0"

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
registries = "General"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
registries = "General"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
deps = ["Serialization"]
git-tree-sha1 = "eb3b85ca904a1d1647b15acc091b113b8874b648"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.3"

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
registries = "General"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WGLMakie]]
deps = ["Bonito", "Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "798486379d2b1ef502c5f68e18885b8274183b66"
registries = "General"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.13.17"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
registries = "General"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
registries = "General"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "ce91b54a5f90c4f27df6950706234b34e49cf79e"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.14+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
registries = "General"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a54cf54677c9f6eac0d97e9d757c233b56bc21b9"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.13+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
registries = "General"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "9c685a6e94e255019e2bf9486318147cd42fe34a"
registries = "General"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.17"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.2+0"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "54fdde8cc7d9350412f75aaa98b1d57f9d4efca2"
registries = "General"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.44.8+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
registries = "General"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "21e05cef1411bfe6b22a1ca8cbdc64a297d4b447"
registries = "General"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.15.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cb007192783c56d8249db4cf0e3495001edfe414"
registries = "General"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.5+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdeflate_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "11c07205bc3b3b7833b8fb36beb4e4764d0126a7"
registries = "General"
uuid = "46979653-d7f6-5232-b59e-dd310c4598de"
version = "1.26.0+1"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "28e57478e8a160d346a19c28b3fffb9273bcc9c2"
registries = "General"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.134+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
registries = "General"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1fcd6cf05d23adab0b5ffc496d020537215751dc"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "e067c8bae65bb40866552296a2d5b4dc65b8928f"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.80.702+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "21cbc8fc2ac0bb182ed461e556d9c062ae6a3164"
registries = "General"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.24.1+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
registries = "General"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "52d3b9475133c3bc8c0a7f90f18bfc8cdb5a443a"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.70.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.3+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
registries = "General"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
registries = "General"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
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
# ╟─8e31c397-10cd-4ef3-85ce-30b1bb888135
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
