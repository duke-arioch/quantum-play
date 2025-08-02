# -----------------------------------------------
#  SU(2) invariant–subspace counter and benchmark
# -----------------------------------------------
from functools import lru_cache
from sympy import S          # exact rationals for half‑integers
import numpy  as np
import pandas as pd

# ---------- 1. the raw invariant counter ----------
@lru_cache(None)
def dim_inv(spins):
    """
    Dimension of Inv( ⊗_e V_{j_e} )
    `spins` is a tuple of sympy.Rational entries (e.g.  S(1)/2, 1, 3/2, …).
    """
    if not spins:                 # empty tensor product → scalar rep
        return 1
    if len(spins) == 1:           # single irrep is singlet only for j = 0
        return 1 if spins[0] == 0 else 0

    j1, j2, *rest = spins
    total = 0
    j_min, j_max = abs(j1 - j2), j1 + j2
    # iterate in half‑integer steps: j = j_min, j_min+½, …, j_max
    step_count = int(2*(j_max - j_min)) + 1
    for n in range(step_count):
        j = j_min + S(n, 2)
        total += dim_inv((j, *rest))
    return total


# ---------- 2. helper wrappers ----------
def boundary_dim(m):
    """d₀ for a boundary with 2m spin‑½ edges (even parity)"""
    return dim_inv(tuple([S(1)/2] * (2*m)))


def bridged_dim(m, jb):
    """d₁ after inserting ONE antiparallel spin‑j_b bridge"""
    # two extra edges, both spin‑j_b
    spins = tuple([S(1)/2]*(2*m) + [S(jb), S(jb)])
    return dim_inv(spins)


# ---------- 3. analytic benchmark WITHOUT brute force ----------
@lru_cache(None)
def analytic_dim(m, jb):
    """
    Proven formula:
      d₁(m,j_b) = Σ_{k=0}^{2j_b}  dim Inv( (V_{½})^{⊗2m} ⊗ V_k )
    Each summand is still Catalan‑like, so we reuse `dim_inv`.
    """
    total = 0
    for k in range(int(2*jb)+1):                # k = 0,1,…,2j_b   (all integers)
        total += dim_inv(tuple([S(1)/2]*(2*m) + [S(k)]))
    return total


# ---------- 4. sweep the parameter range and print a nice table ----------
rows = []
for m in range(1, 5):          # m = 1 … 4  (feel free to extend)
    for jb in [S(1)/2, 1, S(3)/2, 2]:
        d0   = boundary_dim(m)
        dnum = bridged_dim(m, jb)
        dana = analytic_dim(m, jb)
        err  = dnum - dana
        rel  = err / dana if dana else float('nan')
        rows.append([int(m), float(jb), int(dnum), int(dana), int(err), rel])

df = pd.DataFrame(rows,
    columns=['m', 'jb', 'd_num', 'd_ana', 'err', 'rel'])

print("\nNumerical vs analytic (even parity, single bridge)")
print(df.to_string(index=False, formatters={'rel': '{:.2e}'.format}))

# ---------- 5. entropy growth simulation and plot ----------
import matplotlib.pyplot as plt

T = 50          # number of time steps
lam = 0.2       # bridge insertion rate per step
jb  = S(1)/2    # bridge spin
m = 2           # use m=2 for the boundary dimension calculation
d0 = boundary_dim(m)

# Generate cumulative bridge insertions
traj = np.random.poisson(lam, size=T).cumsum()

# Calculate entropy at each time step
S_t = np.log(float(d0)) + traj * np.log(float(2*jb + 1))

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(range(T), S_t, 'b-', linewidth=2, label=f'Entropy growth (λ={lam}, j_b={float(jb)})')
plt.xlabel('Time step t')
plt.ylabel('Entropy S_t')
plt.title('Entropy Growth Due to Bridge Insertions')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# Save the plot
plt.savefig('growth.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.show()

print(f"\nEntropy growth plot saved as 'growth.pdf'")
print(f"Initial entropy S_0 = ln({d0}) ≈ {np.log(float(d0)):.3f}")
print(f"Final entropy S_{T-1} ≈ {S_t[-1]:.3f}")

