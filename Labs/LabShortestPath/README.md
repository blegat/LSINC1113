# Exercice : Optimisation de Traitement Médical avec l'Algorithme de Dijkstra

## Contexte

Un patient présente plusieurs pathologies qui nécessitent un traitement. Chaque médicament peut guérir certaines pathologies mais peut également en causer d'autres comme effets secondaires. L'objectif est de trouver le plan de traitement de coût minimal pour guérir toutes les pathologies.

## Problème

Étant donné :
- Un ensemble de pathologies initiales du patient
- Une liste de médicaments disponibles, chacun avec :
  - Un nom
  - Les pathologies qu'il guérit
  - Les effets secondaires qu'il peut causer
  - Un prix en euros

Trouver la séquence de médicaments qui :
1. Guérit toutes les pathologies initiales
2. Minimise le coût total du traitement
3. Prend en compte les effets secondaires qui peuvent apparaître

## Structure du Code

### Structure `Medication`
```julia
struct Medication
    name::String           # Nom du médicament
    cures::Int            # Bitset des pathologies guéries
    side_effects::Int     # Bitset des effets secondaires
    price::Float64        # Prix en euros
end
```

### Fonctions Fournies

1. **`generate_medical_data(n_pathologies, n_medications)`**
   - Génère des données de test aléatoires
   - Crée des pathologies et médicaments avec des propriétés réalistes

2. **`bitset_to_pathologies(bitset, pathology_names)`**
   - Convertit un bitset en liste lisible de pathologies

3. **`apply_medication(current_state, medication)`**
   - Applique un médicament à l'état actuel des pathologies
   - Retire les pathologies guéries et ajoute les effets secondaires

### Fonction à Implémenter

**`find_optimal_treatment(initial_pathologies, medications)`**
- Implémentez l'algorithme de Dijkstra
- Trouvez le chemin de coût minimal vers l'état "toutes pathologies guéries"
- Retournez `(coût_total, séquence_de_médicaments)`

## Algorithme de Dijkstra

L'algorithme utilise :
- **États** : Combinaison de pathologies actives (représentées par un bitset)
- **Transitions** : Application d'un médicament
- **Coût** : Prix du médicament
- **Objectif** : Atteindre l'état 0 (aucune pathologie)

### Pseudo-code
```
1. Initialiser une file de priorité avec l'état initial
2. Tant que la file n'est pas vide :
   a. Extraire l'état de coût minimal
   b. Si l'état est 0, retourner le chemin
   c. Pour chaque médicament :
      - Calculer le nouvel état
      - Ajouter à la file si meilleur coût
3. Retourner "pas de solution" si aucun chemin trouvé
```

## Utilisation

### Exécution
```bash
julia --project=. main.jl
```

### Exemple de Sortie Attendu
```
=== Exercice d'Optimisation de Traitement Médical ===

Pathologies disponibles :
  Bit 1: Pathology_1
  Bit 2: Pathology_2
  ...

Médicaments disponibles :
  Med_1: €45.67
    Guérit: Pathology_1, Pathology_3
    Effets secondaires: Pathology_2

Pathologies initiales du patient :
  Pathology_1, Pathology_3, Pathology_4 (bitset: 13)

Résolution avec l'algorithme de Dijkstra...

Traitement optimal trouvé !
Coût total : €78.45
Séquence de traitement :
  Étape 1: Prendre Med_2 (€32.10)
    Guérit: Pathology_1
    Nouveaux effets secondaires: Pathology_5
    Pathologies restantes: Pathology_3, Pathology_4, Pathology_5
  ...
```

## Conseils d'Implémentation

1. **File de priorité** : Utilisez `PriorityQueue` de `DataStructures`
2. **États visités** : Gardez une trace des états déjà explorés avec leur coût minimal
3. **Bitsets** : Manipulez efficacement les ensembles de pathologies avec les opérations bit à bit
4. **Optimisation** : Évitez de revisiter les mêmes états avec un coût supérieur

## Objectifs d'Apprentissage

1. **Comprendre l'algorithme de Dijkstra** dans un contexte non-standard
2. **Manipuler les bitsets** pour représenter efficacement les états
3. **Modéliser un problème réel** avec des contraintes complexes
4. **Optimiser les performances** avec des structures de données appropriées

## Dépendances

- `Random` : Génération de données aléatoires
- `DataStructures` : File de priorité pour Dijkstra
- `StatsBase` : Fonction `sample` pour l'échantillonnage

## Fichiers

- `exercise.jl` : Code de l'exercice (sans solution)
- `solution.jl` : Solution complète (pour les enseignants)
- `main.jl` : Point d'entrée pour exécuter l'exercice
- `README.md` : Ce fichier d'instructions
