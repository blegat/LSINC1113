# Exercice : Optimisation de Traitement Médical avec l'Algorithme de Dijkstra
# Contexte : Un patient présente plusieurs pathologies qui nécessitent un traitement. 
# Chaque médicament peut guérir certaines pathologies mais peut également en causer 
# d'autres comme effets secondaires. Trouvez le plan de traitement de coût minimal 
# pour guérir toutes les pathologies.

using Random
using DataStructures
using StatsBase

# Structure pour représenter un médicament
struct Medication
    name::String
    cures::Int          # Bitset des pathologies que ce médicament guérit
    side_effects::Int   # Bitset des pathologies que ce médicament cause
    price::Float64      # Prix en euros
end

# Générer des données de test pour l'exercice
function generate_medical_data(n_pathologies::Int, n_medications::Int)
    Random.seed!(42)  # Pour des résultats reproductibles
    
    pathology_names = ["Pathology_$i" for i in 1:n_pathologies]
    medications = Medication[]
    
    for i in 1:n_medications
        # Chaque médicament guérit 1-3 pathologies aléatoires
        n_cures = min(rand(1:3), n_pathologies)
        cures_bits = 0
        cured_pathologies = sample(1:n_pathologies, n_cures, replace=false)
        for p in cured_pathologies
            cures_bits |= (1 << (p-1))
        end
        
        # Chaque médicament peut causer 0-2 effets secondaires (différents de ce qu'il guérit)
        n_side_effects = rand(0:min(2, n_pathologies - n_cures))
        side_effects_bits = 0
        if n_side_effects > 0
            remaining_pathologies = setdiff(1:n_pathologies, cured_pathologies)
            if !isempty(remaining_pathologies)
                side_effect_pathologies = sample(remaining_pathologies, 
                                                min(n_side_effects, length(remaining_pathologies)), 
                                                replace=false)
                for p in side_effect_pathologies
                    side_effects_bits |= (1 << (p-1))
                end
            end
        end
        
        # Prix aléatoire entre 10 et 100 euros
        price = round(rand() * 90 + 10, digits=2)
        
        push!(medications, Medication("Med_$i", cures_bits, side_effects_bits, price))
    end
    
    return pathology_names, medications
end

# Fonction utilitaire pour convertir un bitset en liste lisible de pathologies
function bitset_to_pathologies(bitset::Int, pathology_names::Vector{String})
    pathologies = String[]
    for i in 1:length(pathology_names)
        if (bitset & (1 << (i-1))) != 0
            push!(pathologies, pathology_names[i])
        end
    end
    return pathologies
end

# Appliquer un médicament à l'état actuel des pathologies
function apply_medication(current_state::Int, medication::Medication)
    # Retirer les pathologies guéries et ajouter les effets secondaires
    new_state = current_state & (~medication.cures)  # Retirer les pathologies guéries
    new_state |= medication.side_effects             # Ajouter les effets secondaires
    return new_state
end

# TODO: Implémenter l'algorithme de Dijkstra pour trouver le traitement optimal
function find_optimal_treatment(initial_pathologies::Int, medications::Vector{Medication})
    # Votre implémentation ici
    # Utilisez une file de priorité pour explorer les états par ordre de coût croissant
    # Retournez (coût_total, séquence_de_médicaments)
    
    error("À implémenter : utilisez l'algorithme de Dijkstra pour trouver le traitement optimal")
end

# Fonction principale pour démontrer l'exercice
function main()
    n_pathologies = 5
    n_medications = 8
    
    println("=== Exercice d'Optimisation de Traitement Médical ===\n")
    
    # Générer les données
    pathology_names, medications = generate_medical_data(n_pathologies, n_medications)
    
    println("Pathologies disponibles :")
    for (i, name) in enumerate(pathology_names)
        println("  Bit $i: $name")
    end
    
    println("\nMédicaments disponibles :")
    for (i, med) in enumerate(medications)
        cures = bitset_to_pathologies(med.cures, pathology_names)
        side_effects = bitset_to_pathologies(med.side_effects, pathology_names)
        println("  $(med.name): €$(med.price)")
        println("    Guérit: $(isempty(cures) ? "Aucune" : join(cures, ", "))")
        println("    Effets secondaires: $(isempty(side_effects) ? "Aucun" : join(side_effects, ", "))")
    end
    
    # Exemple de patient avec les pathologies 1, 3, et 4
    initial_pathologies = (1 << 0) | (1 << 2) | (1 << 3)  # Bits 1, 3, 4
    
    println("\nPathologies initiales du patient :")
    initial_names = bitset_to_pathologies(initial_pathologies, pathology_names)
    println("  $(join(initial_names, ", ")) (bitset: $initial_pathologies)")
    
    println("\nRésolution avec l'algorithme de Dijkstra...")
    
    cost, treatment_path = find_optimal_treatment(initial_pathologies, medications)
    
    if cost == Inf
        println("Aucune solution trouvée ! Impossible de guérir toutes les pathologies.")
    else
        println("\nTraitement optimal trouvé !")
        println("Coût total : €$cost")
        println("Séquence de traitement :")
        
        current_state = initial_pathologies
        for (step, med_idx) in enumerate(treatment_path)
            med = medications[med_idx]
            new_state = apply_medication(current_state, med)
            
            println("  Étape $step: Prendre $(med.name) (€$(med.price))")
            
            if med.cures != 0
                cured = bitset_to_pathologies(current_state & med.cures, pathology_names)
                println("    Guérit: $(join(cured, ", "))")
            end
            
            if med.side_effects != 0
                new_side_effects = bitset_to_pathologies(med.side_effects & (~current_state), pathology_names)
                if !isempty(new_side_effects)
                    println("    Nouveaux effets secondaires: $(join(new_side_effects, ", "))")
                end
            end
            
            remaining = bitset_to_pathologies(new_state, pathology_names)
            println("    Pathologies restantes: $(isempty(remaining) ? "Aucune" : join(remaining, ", "))")
            
            current_state = new_state
        end
    end
end

# Exécuter l'exercice
main()
