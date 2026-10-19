# Construire le graphe de flot pour le problème de matching
function build_flow_graph(scanners::Vector{Scanner}, employees::Vector{Employee})
    n_scanners = length(scanners)
    n_employees = length(employees)
    
    # Nœuds: 0=source, 1:n_scanners=scanners, n_scanners+1:n_scanners+n_employees=employés, n_scanners+n_employees+1=sink
    source = 0
    sink = n_scanners + n_employees + 1
    
    edges = Edge[]
    
    # Arêtes de la source vers les scanners (capacité 1)
    for i in 1:n_scanners
        push!(edges, Edge(source, i, 1, 0))
    end
    
    # Arêtes des scanners vers les employés (capacité 1)
    for (emp_idx, employee) in enumerate(employees)
        employee_node = n_scanners + emp_idx
        for scanner_id in employee.can_use_scanners
            push!(edges, Edge(scanner_id, employee_node, 1, 0))
        end
    end
    
    # Arêtes des employés vers le puits (capacité 1)
    for i in 1:n_employees
        employee_node = n_scanners + i
        push!(edges, Edge(employee_node, sink, 1, 0))
    end
    
    return edges, source, sink, n_scanners, n_employees
end

# Algorithme de Ford-Fulkerson pour trouver le flot maximum
function ford_fulkerson(edges::Vector{Edge}, source::Int, sink::Int, n_nodes::Int)
    # Créer la matrice d'adjacence
    capacity = zeros(Int, n_nodes + 1, n_nodes + 1)
    flow = zeros(Int, n_nodes + 1, n_nodes + 1)
    
    for edge in edges
        capacity[edge.from + 1, edge.to + 1] = edge.capacity
    end
    
    # Algorithme de Ford-Fulkerson
    max_flow = 0
    
    while true
        # Trouver un chemin augmentant avec BFS
        parent = fill(-1, n_nodes + 1)
        queue = [source]
        parent[source + 1] = source
        
        found_path = false
        while !isempty(queue) && !found_path
            u = popfirst!(queue)
            for v in 0:n_nodes
                if parent[v + 1] == -1 && capacity[u + 1, v + 1] > flow[u + 1, v + 1]
                    parent[v + 1] = u
                    if v == sink
                        found_path = true
                        break
                    end
                    push!(queue, v)
                end
            end
        end
        
        if !found_path
            break
        end
        
        # Trouver la capacité minimale le long du chemin
        path_flow = Inf
        v = sink
        while v != source
            u = parent[v + 1]
            path_flow = min(path_flow, capacity[u + 1, v + 1] - flow[u + 1, v + 1])
            v = u
        end
        
        # Mettre à jour les flots
        v = sink
        while v != source
            u = parent[v + 1]
            flow[u + 1, v + 1] += path_flow
            flow[v + 1, u + 1] -= path_flow
            v = u
        end
        
        max_flow += path_flow
    end
    
    return max_flow, flow
end

# Algorithme de Max-Flow pour trouver le nombre maximum de scans simultanés
function find_maximum_scans(scanners::Vector{Scanner}, employees::Vector{Employee})
    # Construire le graphe de flot
    edges, source, sink, n_scanners, n_employees = build_flow_graph(scanners, employees)
    n_nodes = n_scanners + n_employees + 1
    
    # Appliquer l'algorithme de Ford-Fulkerson
    max_flow, flow_matrix = ford_fulkerson(edges, source, sink, n_nodes)
    
    # Extraire les assignations optimales
    assignments = Tuple{Int, Int}[]
    for i in 1:n_scanners
        for j in 1:n_employees
            employee_node = n_scanners + j
            if flow_matrix[i + 1, employee_node + 1] > 0
                push!(assignments, (i, j))
            end
        end
    end
    
    return max_flow, assignments
end

# Fonction pour trouver les formations optimales
function find_optimal_training(scanners::Vector{Scanner}, employees::Vector{Employee})
    # Trouver la capacité actuelle
    current_max, _ = find_maximum_scans(scanners, employees)
    
    best_improvement = 0
    best_training = Tuple{Int, Int}[]
    best_new_capacity = current_max
    
    # Essayer de former chaque employé sur chaque scanner qu'il ne connaît pas encore
    for (emp_idx, employee) in enumerate(employees)
        for scanner_id in 1:length(scanners)
            if !(scanner_id in employee.can_use_scanners)
                # Créer une copie des employés avec cette formation
                new_employees = deepcopy(employees)
                push!(new_employees[emp_idx].can_use_scanners, scanner_id)
                
                # Calculer la nouvelle capacité
                new_max, _ = find_maximum_scans(scanners, new_employees)
                improvement = new_max - current_max
                
                if improvement > best_improvement
                    best_improvement = improvement
                    best_training = [(emp_idx, scanner_id)]
                    best_new_capacity = new_max
                elseif improvement == best_improvement && improvement > 0
                    push!(best_training, (emp_idx, scanner_id))
                end
            end
        end
    end
    
    return best_training, best_new_capacity
end
