/-
Copyright (c) 2025 Quantum Play. All rights reserved.
Spin Networks - graph structures with SU(2) irrep labels
-/

import SU2Rep

/-- A finite oriented graph with vertex and edge sets -/
structure Graph where
  vertices : Finset Nat
  edges : Finset (Nat × Nat)
  edge_valid : ∀ (v, w) ∈ edges, v ∈ vertices ∧ w ∈ vertices

/-- SU(2) irrep label on an edge (stored as twice_j) -/
structure EdgeLabel where
  edge : Nat × Nat
  twice_j : Int
  non_neg : twice_j ≥ 0

/-- Intertwiner at a vertex (placeholder for now) -/
structure Intertwiner where
  vertex : Nat
  incident_spins : List Int
  gauge_invariant : Bool  -- Will be refined later

/-- Definition 1: Spin network -/
structure SpinNetwork where
  graph : Graph
  edge_labels : List EdgeLabel
  intertwiners : List Intertwiner
  labels_match_edges : ∀ label ∈ edge_labels, label.edge ∈ graph.edges
  intertwiner_gauge : ∀ i ∈ intertwiners, i.gauge_invariant = true

namespace SpinNetwork

/-- Partition of vertices into regions A and B -/
structure Partition (sn : SpinNetwork) where
  region_A : Finset Nat
  region_B : Finset Nat
  covers_all : region_A ∪ region_B = sn.graph.vertices
  disjoint : region_A ∩ region_B = ∅

/-- Definition 2: Cut - edges crossing the partition -/
def cut (sn : SpinNetwork) (p : Partition sn) : Finset (Nat × Nat) :=
  sn.graph.edges.filter (fun (v, w) =>
    (v ∈ p.region_A ∧ w ∈ p.region_B) ∨ (v ∈ p.region_B ∧ w ∈ p.region_A))

/-- Boundary multiset of spins on the cut -/
def boundary_spins (sn : SpinNetwork) (p : Partition sn) : List Int :=
  (cut sn p).toList.map (fun edge =>
    match sn.edge_labels.find? (fun label => label.edge = edge) with
    | some label => label.twice_j
    | none => 0)  -- Default case (shouldn't happen with well-formed networks)

/-- Definition 3: Relational entropy -/
def relational_entropy (sn : SpinNetwork) (p : Partition sn) : Nat :=
  singlet_dimension (boundary_spins sn p)

/-- Helper: Check if boundary has even parity -/
def boundary_has_even_parity (sn : SpinNetwork) (p : Partition sn) : Bool :=
  boundary_parity (boundary_spins sn p)

end SpinNetwork
