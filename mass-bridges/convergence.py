# Re-run after state reset

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.linalg import eigh
from caas_jupyter_tools import display_dataframe_to_user

def solve_ground_state(N=600, Rmax=40.0, lam=1.0):
    r = np.linspace(1e-6, Rmax, N)
    dr = r[1]-r[0]
    main = np.full(N, 1.0/dr**2)
    off  = np.full(N-1, -1.0/(2.0*dr**2))
    T = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    V = -lam / r
    H = T + np.diag(V)
    w, v = eigh(H)
    neg_mask = w < 0
    w = w[neg_mask]; v = v[:, neg_mask]
    idx0 = np.argmin(w)
    E0 = w[idx0]
    u0 = v[:, idx0]
    norm = np.trapz(np.abs(u0)**2, x=r)
    u0 /= np.sqrt(norm)
    return r, u0, E0

def fit_tail(r, u, r_min=6.0, r_max=20.0):
    mask = (r >= r_min) & (r <= r_max)
    x = r[mask]
    y = np.log(np.maximum(np.abs(u[mask])**2, 1e-300))
    A = np.vstack([np.ones_like(x), x]).T
    coeff, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    a, b = coeff
    return a, b

configs = [
    {"N": 400, "Rmax": 30.0},
    {"N": 600, "Rmax": 40.0},
    {"N": 900, "Rmax": 60.0},
]

rows = []
paths = {}
for cfg in configs:
    r, u, E0 = solve_ground_state(N=cfg["N"], Rmax=cfg["Rmax"], lam=1.0)
    a, b = fit_tail(r, u, r_min=6.0, r_max=min(20.0, 0.6*cfg["Rmax"]))
    plt.figure()
    mask = (r >= 0) & (r <= cfg["Rmax"])
    plt.plot(r[mask], np.log(np.maximum(np.abs(u[mask])**2, 1e-300)), label="log |u|^2")
    xfit = np.linspace(6.0, min(20.0, 0.6*cfg["Rmax"]), 100)
    yfit = a + b*xfit
    plt.plot(xfit, yfit, '--', label=f"fit slope ≈ {b:.3f}")
    plt.xlabel("r (a.u.)"); plt.ylabel("log |u(r)|^2"); plt.title(f"1s tail fit (N={cfg['N']}, Rmax={cfg['Rmax']})\nE0 ≈ {E0:.6f} a.u. (theory -0.5)")
    plt.legend()
    outpath = f"./mass-bridges/tailfit_N{cfg['N']}_R{int(cfg['Rmax'])}.png"
    plt.savefig(outpath, bbox_inches="tight"); plt.close()
    paths[f"N{cfg['N']}_R{int(cfg['Rmax'])}"] = outpath
    rows.append({
        "N": cfg["N"], "Rmax": cfg["Rmax"], "E0_numeric": E0, "E0_theory": -0.5, 
        "tail_slope_numeric": b, "tail_slope_theory": -2.0
    })

df = pd.DataFrame(rows)
display_dataframe_to_user("1s tail fit results (numeric vs theory)", df.round(6))

csv_path = "./mass-bridges/1s_tail_fit_results.csv"
df.to_csv(csv_path, index=False)

{"results_csv": csv_path, "plots": paths}
