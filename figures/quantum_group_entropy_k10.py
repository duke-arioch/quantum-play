import numpy as np, matplotlib.pyplot as plt
k = 10
j = np.linspace(0, k/2, 200)
ΔS = np.log(np.minimum(2*j, k-2*j) + 1)
print(f"ΔS for k={k}:", ΔS, "...")
plt.plot(j, ΔS, lw=2)
plt.axvline(k/4, ls='--'), plt.axhline(np.log(6), ls=':')
plt.plot(k/4, np.log(6), 'ko', markersize=8)  # solid dot at peak (2.5, ln 6)
plt.xlabel(r'spin $j_b$'); plt.ylabel(r'$\Delta S_\gamma$ (nats)')
plt.tight_layout();# plt.show()
plt.savefig('quantum_group_entropy_k10.pdf')