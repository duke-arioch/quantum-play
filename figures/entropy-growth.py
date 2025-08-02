import math
import numpy as np
import matplotlib.pyplot as plt

def catalan(n):
    return math.comb(2*n, n)//(n+1)

k_vals = list(range(0, 7))

# Even parity case: start with 2 spin-1/2 edges (m=1)
S_even = [math.log(catalan(1 + k)) for k in k_vals]

# Odd parity case: start with 1 spin-1/2 edge (odd); any bridge keeps parity odd => d0 = 0 => entropy undefined.
# Represent undefined as NaN so it shows as missing points.
S_odd = [np.nan for _ in k_vals]

plt.figure(figsize=(6,4))
plt.plot(k_vals, S_even, marker='o', label='Even parity (2m edges)')
plt.plot(k_vals, S_odd, marker='x', label='Odd parity (2m+1 edges)')
plt.xlabel("Number of bridges inserted $k$")
plt.ylabel("Entropy $S_\\gamma$")
plt.title("Entropy growth vs. parity of half‑integer edge count")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("entropy_growth_even_vs_odd.pdf")
plt.show()
