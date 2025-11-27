# Lab Max Flow - Optimisation de Scanners Hospitaliers

Un hôpital possède 5 scanners différents et 5 employés. Seulement, tous les employés ne savent pas utiliser tous les scanners.
Les scanners sont numérotés de 1 à 5 et les employés de A à E. Les employés sachant utiliser chaque scanner sont
représentés par le diagramme suivant :

```mermaid
graph LR
    S1[Scanner 1] --> E1[Employé A]
    S1 --> E2[Employé B]
    S2[Scanner 2] --> E1
    S2 --> E2
    S3[Scanner 3] --> E2
    S4[Scanner 4] --> E2
    S4 --> E4[Employé D]
    S4 --> E5[Employé E]
    S5[Scanner 5] --> E1
    S5 --> E3[Employé C]
    S5 --> E5
```

*Une flèche de scanner i vers employé j représente que j sait utiliser le scanner i.*

## Questions

### Question 1
En supposant que les 5 employés sont disponibles. Combien de scans peuvent être effectués en même temps ?

<details>
<summary>Solution</summary>

Les scanners 1, 2, 3 ne peuvent qu'être utilisés par A et B donc on ne sait pas utiliser ces 3 scanners en même temps.
On ne sait donc utiliser que 4 scanners à la fois.

</details>

### Question 2
Proposer une formation à l'utilisation d'un scan pour un employé qui augmenterait le nombre de scans pouvant être effectués simultanément.

<details>
<summary>Solution</summary>

Il faut former C, D ou E à l'utilisation d'un des scanners 1, 2 ou 3.

</details>

Suite à l'apparition d'une pandémie, l'état a investi massivement dans l'hôpital.
Il y a maintenant 100 scanners et 100 employés pouvant effectuer les scans.

### Question 3
Quel algorithme pouvez-vous utiliser pour résoudre ce problème efficacement ? (*Indice: il faut peut-être ajouter des nœuds fictifs au graphe ci-dessus pour que ça corresponde à un des problèmes vu en cours...*)

### Question 4
Comment trouver quelle formation proposer pour augmenter la capacité de scan de l'hôpital à partir de la solution du problème ?

<details>
<summary>Solution</summary>

On ajoute une source qui relie à tous les scanners et on relie tous les employés à une source.
Ça donne un problème de Max-Flow avec une capacité de 1 pour toutes les arêtes.

De façon équivalente, le Max-Flow est égal à Min-Cut. La solution du Min-Cut est de grouper les nœuds
bleus et verts comme ci-dessous.
La coupe est alors formée par les arêtes en rouge qui relient un nœud bleu
à un nœud vert.

Il faut former un employé dont la solution du Max-Flow donne une valeur de 0
à son arête vers T à utiliser un scanner dont la solution du Max-Flow donne
une valeur de 0 à l'arête le reliant à S.

```mermaid
graph LR
    S[Source S] --> S1[Scanner 1]
    S --> S2[Scanner 2]
    S --> S3[Scanner 3]
    S -.->|"Coupe (rouge)"| S4[Scanner 4]
    S -.->|"Coupe (rouge)"| S5[Scanner 5]
    
    S1 --> E1[Employé A]
    S1 --> E2[Employé B]
    S2 --> E1
    S2 --> E2
    S3 --> E2
    S4 --> E2
    S4 --> E4[Employé D]
    S4 --> E5[Employé E]
    S5 --> E1
    S5 --> E3[Employé C]
    S5 --> E5
    
    E1 -.->|"Coupe (rouge)"| T[Puits T]
    E2 -.->|"Coupe (rouge)"| T
    E3 --> T
    E4 --> T
    E5 --> T
    
    classDef blue fill:#e1f5fe
    classDef green fill:#e8f5e8
    classDef red stroke:#f44336,stroke-width:3px
    
    class S,S1,S2,S3,E1,E2 blue
    class S4,S5,E3,E4,E5,T green
```

</details>
