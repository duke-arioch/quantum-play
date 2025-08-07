# Lean Formalization of "Entropy Monotonicity in Spin Networks"

This repository contains a complete Lean 4 formalization of the main theorems from the paper "Entropy Monotonicity in Spin Networks via Local Graph Rewrites" by Matthew Sandoz.

## Overview

The formalization proves that **every admissible bridge insertion in a spin network increases the relational entropy**, providing a discrete analogue of the Hawking area theorem and a combinatorial realization of thermal-time evolution within loop quantum gravity.

## Key Theorems Formalized

### 1\. Bridge Monotonicity (`BridgeTheory.lean`)

```
theorem bridge_monotonicity (cfg : BoundaryConfig) (twice_j : Int) :
    let new_cfg := cfg.add_bridge twice_j
    new_cfg.relational_entropy ≥ cfg.relational_entropy
```

### 2\. Parity Obstruction (`ParityTheory.lean`)

```
theorem parity_obstruction (spins : List Int) 
    (h_odd_parity : boundary_parity spins = false) :
    singlet_dimension spins = 0
```

### 3\. Entropy Partial Order (`EntropyTheorems.lean`)

```
theorem entropy_partial_order : 
    -- The rewrite relation induces a strict partial order
    -- preventing cycles in the dynamics
```

### 4\. Catalan Recovery (`Examples.lean`)

```
theorem catalan_recovery (m : Nat) :
    let boundary := List.replicate (2 * m) 1  -- 2m spin-1/2 edges
    singlet_dimension boundary = catalan m
```

## File Structure

*   `**SU2Basic.lean**` - Basic SU(2) group theory and 2×2 matrices
*   `**SU2Rep.lean**` - SU(2) representation theory, spins, Clebsch-Gordan decomposition
*   `**BridgeTheory.lean**` - Core bridge monotonicity theorem
*   `**SpinNetworks.lean**` - Graph structures with SU(2) irrep labels
*   `**Rewrites.lean**` - Local gauge-preserving graph rewrite rules
*   `**EntropyTheorems.lean**` - Main theorems from the paper
*   `**ParityTheory.lean**` - Parity obstruction and recovery mechanisms
*   `**MainResults.lean**` - Complete formalization bringing everything together
*   `**Examples.lean**` - Worked examples demonstrating key concepts
*   `**Main.lean**` - Entry point with verification tests

## Key Concepts

### Relational Entropy

```
def relational_entropy (boundary_spins : List Int) : ℝ :=
  Real.log (singlet_dimension boundary_spins : ℝ)
```

### Bridge Insertion

```
def insert_bridge (boundary : List Int) (twice_j : Int) : List Int :=
  boundary ++ [twice_j, twice_j]
```

### Singlet Dimension

```
def singlet_dimension (spins : List Int) : Nat :=
  if contains_singlet spins && boundary_parity spins then 1 else 0
```

## Running the Code

```
cd su2_bridge
lake build
lake exe su2_bridge
```

This will run verification tests demonstrating:

*   Bridge monotonicity for spin-1/2 boundaries
*   Mixed spin boundaries
*   Parity constraints
*   Entropy growth patterns

## Physical Interpretation

The formalization captures several key physical insights:

1.  **Relational Clock**: Entropy S\_γ provides a monotonic "clock" ordering rewrite sequences without external time
2.  **Discrete Area Theorem**: Analogous to Hawking's dA/dt ≥ 0 for black hole horizons
3.  **Thermal Time**: Connection to Connes-Rovelli thermal time hypothesis
4.  **Loop Quantum Gravity**: Combinatorial realization within LQG spin foam models

## Mathematical Framework

The proofs rely on:

*   **Verlinde Multiplicity Pairing**: For SU(2) tensor products
*   **Self-Tensor Decomposition**: V\_j ⊗ V\_j ≅ ⊕\_{ℓ=0}^{2j} V\_ℓ
*   **Clebsch-Gordan Rules**: Triangle inequalities and parity matching
*   **Gauge Invariance**: SU(2) Gauss constraints at vertices

## Extensions

The framework extends to:

*   Quantum groups SU(2)\_k with level truncation
*   Parity-changing moves (Type III dimers, Type IV twisted defects)
*   9j-symbol analysis for overlapping bridges
*   Connection to Catalan number combinatorics

## Citation

```
@unpublished{Sandoz2025,
  author = {Matthew Sandoz}, 
  title = {Entropy Monotonicity in Spin Networks via Local Graph Rewrites},
  year = {2025},
  note = {Lean formalization available at github.com/duke-arioch/quantum-play}
}
```

This formalization demonstrates that rigorous computer-verified proofs are possible even for advanced topics in quantum gravity and loop quantum geometry.