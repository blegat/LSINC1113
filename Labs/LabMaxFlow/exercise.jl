# Exercice : Optimisation de Scanners Hospitaliers avec l'Algorithme de Max-Flow
# Contexte : Un hôpital possède des scanners et des employés. Chaque employé peut utiliser
# certains scanners. L'objectif est de maximiser le nombre de scans simultanés possibles
# en utilisant l'algorithme de Max-Flow.

using Random
using DataStructures

# Structure pour représenter un employé
struct Employee
    name::String
    can_use_scanners::Set{Int}  # Ensemble des scanners que cet employé peut utiliser
end

# Structure pour représenter un scanner
struct Scanner
    id::Int
    name::String
end

# Structure pour représenter une arête dans le graphe de flot
struct Edge
    from::Int
    to::Int
    capacity::Int
    flow::Int
end

# Générer des données de test pour l'exercice
function generate_hospital_data(n_scanners::Int, n_employees::Int)
    Random.seed!(42)  # Pour des résultats reproductibles
    
    scanners = [Scanner(i, "Scanner_$i") for i in 1:n_scanners]
    employees = Employee[]
    
    for i in 1:n_employees
        employee_name = Char(64 + i)  # A, B, C, D, E, ...
        
        # Chaque employé peut utiliser 1-3 scanners aléatoires
        n_scanners_can_use = min(rand(1:3), n_scanners)
        can_use = Set(rand(1:n_scanners, n_scanners_can_use))
        
        push!(employees, Employee(employee_name, can_use))
    end
    
    return scanners, employees
end

# Fonction utilitaire pour afficher les assignations
function print_assignments(assignments::Vector{Tuple{Int, Int}}, scanners::Vector{Scanner}, employees::Vector{Employee})
    println("Assignations optimales :")
    for (scanner_id, emp_idx) in assignments
        scanner = scanners[scanner_id]
        employee = employees[emp_idx]
        println("  $(scanner.name) → $(employee.name)")
    end
end

# Données du problème original (selon le diagramme du README)
function create_original_problem_data()
    scanners = [Scanner(i, "Scanner_$i") for i in 1:5]
    
    # Selon le diagramme : A peut utiliser 1,2,5 ; B peut utiliser 1,2,3,4 ; etc.
    employees = [
        Employee("A", Set([1, 2, 5])),
        Employee("B", Set([1, 2, 3, 4])),
        Employee("C", Set([5])),
        Employee("D", Set([4])),
        Employee("E", Set([4, 5]))
    ]
    
    return scanners, employees
end

# Fonction principale pour démontrer l'exercice
function main()
    println("=== Exercice d'Optimisation de Scanners Hospitaliers ===\n")
    
    # Utiliser les données du problème original
    scanners, employees = create_original_problem_data()
    
    println("Scanners disponibles :")
    for scanner in scanners
        println("  $(scanner.name)")
    end
    
    println("\nEmployés et leurs compétences :")
    for employee in employees
        scanner_names = [scanners[s].name for s in employee.can_use_scanners]
        println("  $(employee.name): peut utiliser $(join(sort(scanner_names), ", "))")
    end
    
    println("\nRésolution avec l'algorithme de Max-Flow...")
    
    max_scans, assignments = find_maximum_scans(scanners, employees)
    
    println("\nRésultat de l'optimisation :")
    println("Nombre maximum de scans simultanés : $max_scans")
    print_assignments(assignments, scanners, employees)
    
    println("\nAnalyse des formations optimales...")
    training_recommendations, new_capacity = find_optimal_training(scanners, employees)
    
    println("Formations recommandées :")
    for (emp_idx, scanner_id) in training_recommendations
        employee = employees[emp_idx]
        scanner = scanners[scanner_id]
        println("  Former $(employee.name) sur $(scanner.name)")
    end
    println("Nouvelle capacité maximale après formation : $new_capacity")
end
