# Appendix: Parameter Estimates & Example Curves for I_B(t)
# This script produces a simple ΛCDM-based example evolution for the index budget I_B,
# using the ODE from the paper. It also computes an example instability timescale τ_inst.
#
# Outputs in current directory:
#   - ibudget_example_curves.png / .pdf
#   - ibudget_example_data.json
#   - inst_timescale_sweep.png / .pdf
#
# Notes:
# - Values chosen are illustrative (order-of-magnitude), not fits. You can tune alphas later.
# - Cosmology: Planck-like parameters (H0, Ω_m, Ω_r, Ω_Λ).


import numpy as np, json, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- Cosmology & units ----------
# H0 in s^-1: 67.4 km/s/Mpc
H0_km_s_Mpc = 67.4
Mpc_in_km = 3.0856775814913673e19
H0 = H0_km_s_Mpc * 1000.0 / Mpc_in_km  # s^-1
Omega_m = 0.315
Omega_r = 9.0e-5
Omega_L = 1.0 - Omega_m - Omega_r

# CMB temperature
T0 = 2.7255  # K
k_B = 1.380649e-23  # J/K
hbar = 1.054571817e-34  # J*s
c = 299792458.0  # m/s

# Critical density today (SI): ρ_c0 = 3 H0^2 / (8πG)
G = 6.67430e-11  # m^3 / (kg s^2)
rho_c0 = 3 * H0**2 / (8 * math.pi * G)  # kg/m^3

# Entropy density s ~ (2π^2/45) g_* T^3 (natural units). We'll keep as an arbitrary normalized unit.
# For late times, g_* ≈ 3.36 effective relativistic dof for photons+neutrinos.
g_star = 3.36
pref_s = (2 * math.pi**2 / 45.0) * g_star  # in natural units; we will keep "arb" units

def H_of_a(a):
    return H0 * math.sqrt(Omega_r / a**4 + Omega_m / a**3 + Omega_L)

def s_of_a(a):
    # T(a) = T0/a; scale s ∝ T^3 ∝ a^-3; return in arbitrary units normalized to s(a=1)=pref_s*T0^3
    return pref_s * (T0**3) / (a**3)

def rho_m_of_a(a):
    return Omega_m * rho_c0 / a**3  # kg/m^3

def rho_r_of_a(a):
    return Omega_r * rho_c0 / a**4  # kg/m^3

def rho_L_of_a(a):
    return Omega_L * rho_c0  # kg/m^3

# Weighted matter functional Θ = Σ ν_i ρ_i. Choose illustrative ν weights (dimensionful absorbed into α_χ).
nu_m, nu_r, nu_L = 1.0, 0.2, 0.0

def Theta_of_a(a):
    return nu_m * rho_m_of_a(a) + nu_r * rho_r_of_a(a) + nu_L * rho_L_of_a(a)

# ---------- Budget ODE parameters ----------
# Coefficients (illustrative; order-of-magnitude). You will tune to your algebra later.
alpha_H = 0.8     # dimensionless
alpha_s = 1.0e-3  # sets strength of entropy screening
alpha_chi = 1.0e-3  # sets strength of matter drag
# No micro-source in baseline:
def S_micro(a): return 0.0

# ---------- Integrate dI_B/dt via scale factor a ----------
# Use chain rule: dI/dt = (dI/da) * aH  =>  dI/da = [ -α_H H I  - α_s s  - α_χ Θ + S ] / (a H)
# Integrate from recombination (z≈1100, a≈1/1101) to today (a=1).
a_start = 1.0 / 1101.0
a_end = 1.0
N = 1800
agrid = np.linspace(a_start, a_end, N)
I_init = 10.0  # "nats per unit area"—illustrative starting budget at recombination

def integrate_ibudget(I0, alpha_H, alpha_s, alpha_chi):
    I = I0
    I_hist = [I0]
    for i in range(1, len(agrid)):
        a0, a1 = agrid[i-1], agrid[i]
        H0a = H_of_a(a0)
        s0 = s_of_a(a0)
        Th0 = Theta_of_a(a0)
        rhs0 = -alpha_H * H0a * I - alpha_s * s0 - alpha_chi * Th0 + S_micro(a0)
        dIda0 = rhs0 / (a0 * H0a)

        # simple RK2 (midpoint)
        amid = 0.5*(a0+a1)
        Hm = H_of_a(amid)
        sm = s_of_a(amid)
        Thm = Theta_of_a(amid)
        Imid = I + dIda0 * (a1 - a0)/2.0
        rhsm = -alpha_H * Hm * Imid - alpha_s * sm - alpha_chi * Thm + S_micro(amid)
        dIdam = rhsm / (amid * Hm)

        I = I + dIdam * (a1 - a0)
        I_hist.append(I)
    return np.array(I_hist)

I_hist = integrate_ibudget(I_init, alpha_H, alpha_s, alpha_chi)

# ---------- Produce comparison curves for different α's ----------
curves = {}
curves['baseline'] = I_hist.tolist()

I_softH = integrate_ibudget(I_init, 0.4, alpha_s, alpha_chi)
curves['low_alpha_H'] = I_softH.tolist()

I_strong_matter = integrate_ibudget(I_init, alpha_H, alpha_s, 4.0e-3)
curves['strong_matter_drag'] = I_strong_matter.tolist()

# ---------- Instability timescale example ----------
# τ_inst ≈ χ / (D_I k_*^2) with k_*^2 ~ χ^{-1} / σ  → τ_inst ~ χ / (D_I * χ^{-1}/σ) = σ / D_I
# So in this simple estimate τ_inst is set by σ (stiffness) and D_I (mobility).
# We'll sweep D_I and σ to show a band of possibilities.
D_I_values = np.logspace(0, 8, 30)   # s^-1 (from 1 Hz to 10^8 Hz attempt rate)
sigma_values = [1e-4, 1e-3, 1e-2]    # arbitrary stiffness units (set by rewrite cost scale)
tau_inst_curves = {f"sigma={s}": (1.0*np.array([s/Di for Di in D_I_values])).tolist() for s in sigma_values}

# ---------- Save data ----------
out = {
    "cosmology": {
        "H0_s^-1": H0,
        "Omega_m": Omega_m,
        "Omega_r": Omega_r,
        "Omega_L": Omega_L,
        "rho_c0_kg_m3": rho_c0
    },
    "alphas": {
        "alpha_H": alpha_H,
        "alpha_s": alpha_s,
        "alpha_chi": alpha_chi
    },
    "initial_conditions": {
        "a_start": a_start,
        "z_start": 1.0/a_start - 1.0,
        "I_init": I_init
    },
    "agrid": agrid.tolist(),
    "I_curves": curves,
    "tau_inst": {
        "D_I_values_s^-1": D_I_values.tolist(),
        "sigma_values": sigma_values,
        "curves": tau_inst_curves
    }
}

with open("ibudget_example_data.json","w") as f:
    json.dump(out, f, indent=2)

# ---------- Plots ----------
# 1) I_B(a) vs redshift z
zgrid = 1.0/agrid - 1.0

plt.figure(figsize=(8,5))
plt.plot(zgrid, curves['baseline'], label=r'baseline $\alpha_H=0.8,\ \alpha_\chi=10^{-3}$')
plt.plot(zgrid, curves['low_alpha_H'], '--', label=r'weaker expansion $\alpha_H=0.4$')
plt.plot(zgrid, curves['strong_matter_drag'], ':', label=r'stronger matter drag $\alpha_\chi=4\!\times\!10^{-3}$')
plt.gca().invert_xaxis()
plt.xlabel("redshift z")
plt.ylabel(r"$I_B$ (nats / unit area, illustrative)")
plt.title(r"Example $I_B$ evolution from recombination ($z\!\sim\!1100$) to today")
plt.legend()
plt.tight_layout()
plt.savefig("ibudget_example_curves.png", dpi=160)
plt.savefig("ibudget_example_curves.pdf")
plt.close()

# 2) τ_inst vs D_I for different σ
plt.figure(figsize=(8,5))
for s, y in tau_inst_curves.items():
    plt.loglog(D_I_values, y, label=s)
plt.xlabel(r"index mobility $D_I$ (s$^{-1}$)")
plt.ylabel(r"instability timescale $\tau_{\rm inst}$ (s)")
plt.title(r"Example $\tau_{\rm inst}\approx \sigma/D_I$ bands")
plt.legend()
plt.tight_layout()
plt.savefig("inst_timescale_sweep.png", dpi=160)
plt.savefig("inst_timescale_sweep.pdf")
plt.close()

"ibudget_example_curves.png|ibudget_example_curves.pdf|inst_timescale_sweep.png|inst_timescale_sweep.pdf|ibudget_example_data.json"
