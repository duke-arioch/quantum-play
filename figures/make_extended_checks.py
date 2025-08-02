#!/usr/bin/env python3
# -------------------------------------------
# Bridge‑monotonicity numerical verification
# -------------------------------------------

from collections import Counter
from fractions import Fraction as Q    # exact rationals, eg Q(3,2) == 3/2
from itertools import product
from math import comb, prod, log
import pandas as pd

# --------------------------------------------------------------------
# Singlet‑dimension routine (multiplicity‑aware Clebsch–Gordan fusion)
# --------------------------------------------------------------------
def singlet_dim(spins):
    """
    Return dim Inv(⊗_e V_{j_e}) for an arbitrary iterable *spins*.
    The algorithm keeps a Counter of multiplicities for each total
    spin J encountered so far, so it sums over *all* binary fusion trees.
    Complexity: O(n · J_max^2), fast for the cases relevant here.
    """
    counter = Counter({0: 1})        # start with the trivial rep
    for j in spins:
        new = Counter()
        for J, mult in counter.items():
            # Clebsch–Gordan range |J−j| .. J+j in integer steps
            low, high = abs(J - j), J + j
            # step = 1 when J+j and |J−j| differ by an integer,
            # but for half‑integers that still enumerates every allowed value
            k = 0
            while True:
                Jtot = low + k
                if Jtot > high: break
                new[Jtot] += mult
                k += 1
        counter = new
    return counter[0]                # multiplicity of total spin‑0


# -------------------------------
# Analytic formulae from the paper
# -------------------------------
def catalan(m):                  # C_m = (2m choose m)/(m+1)
    return comb(2 * m, m) // (m + 1)

def analytic_dim_even(m, bridges):
    """
    Analytic dim Inv for an even‑parity cut carrying 2m spin‑½ legs
    plus a list 'bridges' of additional antiparallel (j_b,j_b) pairs.
    """
    return catalan(m) * prod(2 * jb + 1 for jb in bridges)

# ---------------
# Test‑bench grid
# ---------------
max_m   = 4
bridge_js = [Q(1,1), Q(3,2), Q(2,1)]      # 1, 3/2, 2
max_k   = 4                               # up to four identical bridges each
rows = []

for m, jb, k in product(range(1, max_m + 1), bridge_js, range(1, max_k + 1)):
    base_legs   = [Q(1,2)] * (2 * m)      # the 2m spin‑½ edges on the cut
    bridge_legs = [jb, jb] * k            # k antiparallel pairs
    spins       = base_legs + bridge_legs

    d_num = singlet_dim(spins)
    d_ana = analytic_dim_even(m, [jb] * k)

    rows.append(dict(m=m,
                     jb=float(jb),    # for nicer printing
                     k=k,
                     d_num=d_num,
                     d_ana=d_ana,
                     err=abs(d_num - d_ana) / d_ana))

# --------------------------
# Pretty print & export
# --------------------------
df = pd.DataFrame(rows)
df.sort_values(['m', 'jb', 'k'], inplace=True)

# Convert error column to float for proper formatting
df['err'] = df['err'].astype(float)

print("\nExtended‑checks summary")
print("=".ljust(70, '='), "\n")
print(df.to_string(index=False,
                   formatters={'err': '{:8.2e}'.format}))

perfect  = (df.err == 0).sum()
good     = (df.err < 1e-10).sum()           # exact integers -> tiny float eps
print(f"\nPerfect matches : {perfect}/{len(df)}")
print(f"err < 1e‑10      : {good}/{len(df)}")

# --- Optional export -------------------------------------------------
df.to_csv("extended_checks.csv", index=False)
with open("extended_checks.tex","w") as f:
    f.write(df.to_latex(index=False, float_format="%.0f",
                        columns=['m','jb','k','d_num','d_ana']))

