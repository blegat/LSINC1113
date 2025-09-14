using Graphs, GraphPlot, SparseArrays, Colors

function myplot(graph, components)
    cols = Colors.distinguishable_colors(length(components))
    node_cols = Vector{eltype(cols)}(undef, nv(graph))
    for (comp, col) in zip(components, cols)
        node_cols[comp] .= col
    end
    gplot(graph, nodefillc = node_cols)
end

n = 40
adjacency = sprand(Bool, n, n, 0.02)

graph = Graph(max.(adjacency, adjacency'))

connected_components(graph)

myplot(graph, connected_components(graph))

function my_connected_component(graph)
    # TODO
    # Vous pouvez commencer par utiliser https://juliacollections.github.io/DataStructures.jl/stable/disjoint_sets/
    # pour vous aider
    # Ensuite, essayez d'implémenter DisjointSet vous-même !
    for edge in edges
        # Arête liant les noeuds `edge.src` et `edge.dst`
    end
end
