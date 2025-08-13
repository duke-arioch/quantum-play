# Re-execute the cell (previous kernel state reset).

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.linalg import eigh

# Simple replacement for caas_jupyter_tools.display_dataframe_to_user
def display_dataframe_to_user(title, df):
    print(f"\n{title}")
    print("=" * len(title))
    print(df.to_string(index=False))
    print()

N = 1500
Rmax = 100.0
lam = 1.0
ell_list = [0,1,2]
num_eigs = 6

r = np.linspace(1e-6, Rmax, N)
dr = r[1]-r[0]

def hamiltonian_matrix(ell, lam):
    main = np.full(N, 1.0/dr**2)
    off  = np.full(N-1, -0.5/dr**2)
    T = np.diag(-0.5*( -2.0*main )) + np.diag(-0.5*off, k=1) + np.diag(-0.5*off, k=-1)
    V = -lam / r + 0.5 * ell*(ell+1) / (r**2)
    H = T + np.diag(V)
    return H

def solve_spectrum(ell, lam, num_eigs):
    H = hamiltonian_matrix(ell, lam)
    w, v = eigh(H)
    mask = w < 0
    w = w[mask]; v = v[:,mask]
    idx = np.argsort(w)[:num_eigs]
    return w[idx], v[:,idx]

def analytic_energy(n, lam):
    return - (lam**2) / (2.0 * n**2)

levels = {
    (1,0): "1s",
    (2,0): "2s", (2,1): "2p",
    (3,0): "3s", (3,1): "3p", (3,2): "3d",
    (4,0): "4s", (4,1): "4p", (4,2): "4d"
}

rows = []
eigvecs = {}
for ell in ell_list:
    w, v = solve_spectrum(ell, lam, num_eigs=num_eigs)
    for k, E in enumerate(w):
        rows.append({"ell": ell, "k_index": k, "E_numeric": E})
    eigvecs[ell] = (w, v)

df = pd.DataFrame(rows).sort_values(["ell","E_numeric"]).reset_index(drop=True)
display_dataframe_to_user("Discrete radial energies (simulation units)", df.round(8))

def normalize(P, dr):
    norm = np.trapezoid(P, dx=dr)
    return P / norm

# Extract bound states safely, checking if they exist
bound_states = {}
for ell in ell_list:
    w, v = eigvecs[ell]
    if len(w) > 0:
        bound_states[ell] = {
            'energy': w[0],
            'wavefunction': v[:,0],
            'probability': normalize(np.abs(v[:,0])**2, dr)
        }
        print(f"Found bound state for ℓ={ell}: E = {w[0]:.6f}")
    else:
        print(f"No bound states found for ℓ={ell}")

# Use only available bound states for plotting
available_states = list(bound_states.keys())

# Plots
n_vals = np.arange(1,6)
E_analytic = [analytic_energy(n, lam) for n in n_vals]

plt.figure()
plt.plot(n_vals, E_analytic, marker='o')
plt.xlabel("n"); plt.ylabel("Analytic E_n (a.u.)"); plt.title("Hydrogenic analytic energies E_n = -1/(2 n^2)")
plt.savefig("analytic_energies.png", bbox_inches="tight"); plt.close()

for ell in ell_list:
    w, _ = eigvecs[ell]
    if len(w) > 0:  # Only plot if bound states exist
        plt.figure()
        plt.plot(np.arange(len(w)), w, marker='o')
        plt.xlabel("State index (per ℓ)"); plt.ylabel("Numeric E (sim units)"); plt.title(f"Numeric bound energies for ℓ={ell}")
        plt.savefig(f"numeric_energies_ell{ell}.png", bbox_inches="tight"); plt.close()

# Plot radial probability densities for available bound states
state_labels = {0: "1s", 1: "2p", 2: "3d"}
for ell in available_states:
    if ell in bound_states:
        plt.figure()
        plt.plot(r, bound_states[ell]['probability'])
        plt.xlabel("r (sim units)")
        plt.ylabel(f"|u_{state_labels.get(ell, f'ell{ell}')}(r)|^2")
        plt.title(f"Radial probability density: {state_labels.get(ell, f'ℓ={ell}')}")
        
        # Save plots with consistent naming
        if ell == 0:
            plt.savefig("radial_1s.png", bbox_inches="tight")
        else:
            plt.savefig(f"radial_{state_labels.get(ell, f'ell{ell}')}.png", bbox_inches="tight")
        plt.close()

E1s_analytic = analytic_energy(1, lam)
E1s_phys_eV = -13.605693
s_E = E1s_phys_eV / E1s_analytic  # eV per sim-unit

summary = {
    "grid": {"N": N, "Rmax": Rmax, "dr": float(dr)},
    "lam": lam,
    "E_numeric_first_per_l": {int(ell): float(eigvecs[ell][0][0]) for ell in ell_list if len(eigvecs[ell][0]) > 0},
    "E1s_analytic_sim": float(E1s_analytic),
    "E1s_anchor_eV": E1s_phys_eV,
    "energy_scale_eV_per_simunit": float(s_E),
    "plots": {
        "analytic_energies": "/mnt/data/analytic_energies.png",
        "numeric_energies_s": "/mnt/data/numeric_energies_ell0.png",
        "numeric_energies_p": "/mnt/data/numeric_energies_ell1.png",
        "numeric_energies_d": "/mnt/data/numeric_energies_ell2.png",
        "radial_1s": "/mnt/data/radial_1s.png",
        "radial_2p": "/mnt/data/radial_2p.png",
        "radial_3d": "/mnt/data/radial_3d.png"
    }
}
summary
