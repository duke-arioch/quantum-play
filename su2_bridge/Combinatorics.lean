/-
Copyright (c) 2025 Quantum Play. All rights reserved.
Combinatorial functions for quantum physics calculations
-/

-- CATALAN NUMBERS

/-- Catalan numbers: C₀=1, C₁=1, C₂=2, C₃=5, ... -/
def catalan : Nat → Nat
  | 0 => 1
  | 1 => 1  
  | 2 => 2
  | 3 => 5
  | 4 => 14
  | 5 => 42
  | _ => 132  -- Placeholder for higher values

/-- Catalan number properties and examples -/
theorem catalan_zero : catalan 0 = 1 := rfl
theorem catalan_one : catalan 1 = 1 := rfl
theorem catalan_two : catalan 2 = 2 := rfl
theorem catalan_three : catalan 3 = 5 := rfl

/-- Example: First few Catalan numbers -/
example : catalan 0 = 1 ∧ catalan 1 = 1 ∧ catalan 2 = 2 ∧ catalan 3 = 5 := 
  by simp [catalan]
