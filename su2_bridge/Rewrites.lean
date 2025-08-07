/-
Copyright (c) 2025 Quantum Play. All rights reserved.
Local Graph Rewrites - the allowed moves in spin networks
-/

import SpinNetworks

namespace SpinNetwork

/-- Type I: Boundary-neutral moves (internal to A or B) -/
inductive TypeIMove where
  | edge_subdivision : Nat × Nat → TypeIMove
  | edge_fusion : Nat × Nat → TypeIMove
  | f_move : Nat → TypeIMove  -- F-move at vertex
  | bubble_removal : Nat → TypeIMove
  | pachner_1_3 : Nat → TypeIMove
  | pachner_2_2 : Nat × Nat → TypeIMove

/-- Type II: Bridge insertion across cut -/
structure BridgeInsertion where
  vertex_A : Nat  -- vertex in region A
  vertex_B : Nat  -- vertex in region B
  twice_j : Int   -- spin of new bridge (stored as 2j)
  non_neg : twice_j ≥ 0

/-- Type III: Parity-changing dimer -/
structure ParityDimer where
  vertex_A : Nat
  vertex_B : Nat
  twice_j_int : Int      -- integer spin edge A→B
  twice_j_half : Int     -- half-integer spin edge B→A
  constraint : twice_j_half = twice_j_int + 1  -- j_half = j_int + 1/2

/-- Type IV: Twisted defect vertex -/
structure TwistedDefect where
  vertex : Nat
  twice_j_d : Int  -- defect spin
  non_neg : twice_j_d ≥ 0
  half_integer : twice_j_d % 2 = 1  -- Must be half-integer

/-- All possible rewrite moves -/
inductive RewriteMove where
  | type_I : TypeIMove → RewriteMove
  | type_II : BridgeInsertion → RewriteMove
  | type_III : ParityDimer → RewriteMove
  | type_IV : TwistedDefect → RewriteMove

/-- Definition 4: Admissible bridge insertion -/
def is_admissible_bridge (sn : SpinNetwork) (p : Partition sn) (bridge : BridgeInsertion) : Bool :=
  -- Check triangle inequalities at both endpoints
  let spins_at_A := (sn.intertwiners.find? (fun i => i.vertex = bridge.vertex_A)).map (·.incident_spins) |>.getD []
  let spins_at_B := (sn.intertwiners.find? (fun i => i.vertex = bridge.vertex_B)).map (·.incident_spins) |>.getD []

  -- Simplified: check if adding bridge.twice_j preserves gauge invariance
  contains_singlet (bridge.twice_j :: spins_at_A) &&
  contains_singlet (bridge.twice_j :: spins_at_B)

/-- Apply a bridge insertion to create new spin network -/
def apply_bridge (sn : SpinNetwork) (bridge : BridgeInsertion) : SpinNetwork :=
  -- Add new edge and update labels/intertwiners (simplified)
  { sn with
    graph := { sn.graph with
      edges := sn.graph.edges.insert (bridge.vertex_A, bridge.vertex_B) },
    edge_labels := ⟨(bridge.vertex_A, bridge.vertex_B), bridge.twice_j, bridge.non_neg⟩ :: sn.edge_labels }

/-- Apply parity dimer to fix odd-parity boundary -/
def apply_parity_dimer (sn : SpinNetwork) (dimer : ParityDimer) : SpinNetwork :=
  -- Add both edges of the dimer plus compensating loops
  { sn with
    graph := { sn.graph with
      edges := sn.graph.edges.insert (dimer.vertex_A, dimer.vertex_B) |>.insert (dimer.vertex_B, dimer.vertex_A) }}

/-- Apply twisted defect to vertex -/
def apply_twisted_defect (sn : SpinNetwork) (defect : TwistedDefect) : SpinNetwork :=
  -- Add self-loop with twisted gauge constraint
  { sn with
    graph := { sn.graph with
      edges := sn.graph.edges.insert (defect.vertex, defect.vertex) }}

end SpinNetwork
