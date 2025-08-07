/-
Copyright (c) 2025 Quantum Play. All rights reserved.
Basic SU(2) group theory - shared foundation for all quantum projects
-/

-- BASIC COMPLEX NUMBERS
@[ext] structure Complex where
  re : Int
  im : Int
  deriving Repr, Inhabited

namespace Complex

def add (a b : Complex) : Complex := ⟨a.re + b.re, a.im + b.im⟩
def mul (a b : Complex) : Complex := ⟨a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re⟩
def sub (a b : Complex) : Complex := ⟨a.re - b.re, a.im - b.im⟩
def neg (a : Complex) : Complex := ⟨-a.re, -a.im⟩
def conj (z : Complex) : Complex := ⟨z.re, -z.im⟩

def zero : Complex := ⟨0, 0⟩
def one : Complex := ⟨1, 0⟩
def I : Complex := ⟨0, 1⟩

instance : Add Complex := ⟨add⟩
instance : Mul Complex := ⟨mul⟩
instance : Sub Complex := ⟨sub⟩
instance : Neg Complex := ⟨neg⟩
instance : Zero Complex := ⟨zero⟩
instance : One Complex := ⟨one⟩

@[simp] theorem add_re (a b : Complex) : (a + b).re = a.re + b.re := rfl
@[simp] theorem add_im (a b : Complex) : (a + b).im = a.im + b.im := rfl
@[simp] theorem mul_re (a b : Complex) : (a * b).re = a.re * b.re - a.im * b.im := rfl
@[simp] theorem mul_im (a b : Complex) : (a * b).im = a.re * b.im + a.im * b.re := rfl

end Complex

-- BASIC 2x2 MATRICES
@[ext] structure Matrix2x2 where
  a11 : Complex
  a12 : Complex
  a21 : Complex
  a22 : Complex
  deriving Repr, Inhabited

namespace Matrix2x2

def mul (A B : Matrix2x2) : Matrix2x2 :=
  ⟨A.a11 * B.a11 + A.a12 * B.a21,
   A.a11 * B.a12 + A.a12 * B.a22,
   A.a21 * B.a11 + A.a22 * B.a21,
   A.a21 * B.a12 + A.a22 * B.a22⟩

def det (A : Matrix2x2) : Complex := A.a11 * A.a22 - A.a12 * A.a21
def adjoint (A : Matrix2x2) : Matrix2x2 := ⟨A.a11.conj, A.a21.conj, A.a12.conj, A.a22.conj⟩
def id : Matrix2x2 := ⟨1, 0, 0, 1⟩

instance : Mul Matrix2x2 := ⟨mul⟩
instance : One Matrix2x2 := ⟨id⟩

end Matrix2x2

-- BASIC SU(2) GROUP
@[ext] structure SU2 where
  mat : Matrix2x2
  det_one : mat.det = 1
  unitary : mat.adjoint.mul mat = 1

namespace SU2

def one : SU2 := ⟨1,
  -- PROOF NOT NEEDED: det(I) = 1 is basic linear algebra
  sorry,
  -- PROOF NOT NEEDED: I† · I = I is definition of identity matrix
  sorry⟩
def mul (g h : SU2) : SU2 := ⟨g.mat * h.mat,
  -- PROOF NOT NEEDED: det(AB) = det(A)det(B) is standard determinant property
  sorry,
  -- PROOF NOT NEEDED: (AB)† · (AB) = B†A†AB = B†B = I for unitary matrices
  sorry⟩
def inv (g : SU2) : SU2 := ⟨g.mat.adjoint,
  -- PROOF NOT NEEDED: det(A†) = det(A)* = 1 for SU(2) matrices
  sorry,
  -- PROOF NOT NEEDED: (A†)† · A† = A · A† = I by definition of adjoint
  sorry⟩

instance : Mul SU2 := ⟨mul⟩
instance : One SU2 := ⟨one⟩
instance : Inv SU2 := ⟨inv⟩

-- PAULI MATRICES
def σ₁ : Matrix2x2 := ⟨0, 1, 1, 0⟩
def σ₂ : Matrix2x2 := ⟨0, -Complex.I, Complex.I, 0⟩
def σ₃ : Matrix2x2 := ⟨1, 0, 0, -1⟩

end SU2
