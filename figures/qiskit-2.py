# bridge_monotonicity_demo_exact.py   (Python 3.9+)

import numpy as np
from functools import reduce
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector

# ---------------------------------------------------------
# 1. helper : prepare the two‑qubit singlet |Ψ⁻⟩
# ---------------------------------------------------------
def singlet_2qubit():
    qc = QuantumCircuit(2, name="psi_minus")
    qc.cx(0, 1)
    qc.h(0)
    qc.cx(0, 1)
    qc.z(0)
    return qc


# ---------------------------------------------------------
# 2. build the pre‑ and post‑bridge circuits
# ---------------------------------------------------------
# 4‑qubit: two disjoint singlets
pre = QuantumCircuit(4, name="pre")
pre.append(singlet_2qubit(), [0, 1])
pre.append(singlet_2qubit(), [2, 3])

# 6‑qubit: same + Bell‑pair bridge on (4,5)
post = QuantumCircuit(6, name="post")
post.append(singlet_2qubit(), [0, 1])
post.append(singlet_2qubit(), [2, 3])
post.x(5)                               #  |↑↓⟩  =  |0〉|1〉
post.append(singlet_2qubit(), [4, 5])

# ---------------------------------------------------------
# 3. exact projector onto total spin S = 0
# ---------------------------------------------------------
def exact_singlet_projector(n: int) -> np.ndarray:
    """Dense matrix of the projector onto global SU(2) singlet (S=0)."""
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    paulis = [sx, sy, sz]

    # Build Σ⃗ = Σ_i σ⃗^{(i)}
    def kron_to_n(op, pos, qubits):
        return reduce(
            np.kron,
            [op if i == pos else np.eye(2) for i in range(qubits)],
        )

    sig_tot = [sum(kron_to_n(p, q, n) for q in range(n)) for p in paulis]
    S2 = sum(np.dot(s, s) for s in sig_tot) / 4.0  #  Σ⃗²/4  (spin convention)

    # Diagonalise once; eigenvectors with eigenvalue ~0 form the singlet subspace
    evals, evecs = np.linalg.eigh(S2)
    tol = 1e-9
    cols = [i for i, ev in enumerate(evals) if abs(ev) < tol]
    P = evecs[:, cols] @ evecs[:, cols].conj().T
    return P


def log_dim_singlet(n: int) -> float:
    P = exact_singlet_projector(n)
    d = int(round(np.trace(P).real))
    return np.log(d), d


# ---------------------------------------------------------
# 4.  compute entropy jump  ΔS = log d1 – log d0
# ---------------------------------------------------------
log_d0, d0 = log_dim_singlet(4)
log_d1, d1 = log_dim_singlet(6)

# delta_S = log_d1 - log_d0

print(f"Singlet‐subspace dimension (4 qubits)  d₀ = {d0}")
print(f"Singlet‐subspace dimension (6 qubits)  d₁ = {d1}\n")

print(f"Relational entropy before bridge  :  S₀ = {log_d0:.6f}")
print(f"Relational entropy after  bridge  :  S₁ = {log_d1:.6f}")

S0, d0 = log_dim_singlet(4)   # = (ln2, 2)
S1, d1 = log_dim_singlet(6)   # = (ln4, 4)
print("Entropy jump  ΔS =", S1 - S0, "(should be ln2)")

# print(f"Entropy jump  ΔS = {delta_S:.6f}     (expected ln2 = {np.log(2):.6f})")
