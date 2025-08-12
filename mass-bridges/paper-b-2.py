# make_ibudget_figs.py
# Generates: ibudget_example_curves.(png|pdf), inst_timescale_sweep.(png|pdf),
#            ibudget_vs_temperature.(png|pdf), ibudget_example_data.json

import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

# -------------------- Cosmology (ΛCDM-ish) --------------------
H0_km_s_Mpc = 67.4
Mpc_in_km = 3.0856775814913673e19
H0 = H0_km_s_Mpc * 1000.0 / Mpc_in_km  # s^-1
Omega_m = 0.315
Omega_r = 9.0e-5
Omega_L = 1.0 - Omega_m - Omega_r

T0_K = 2.7255
k_B_eV_per_K = 8.617333262145e-5  # eV/K
T0_eV = T0_K * k_B_eV_per_K
T0_GeV = T0_eV * 1e-9

G = 6.67430e-11  # m^3/(kg*s^2)
rho_c0 = 3 * H0**2 / (8 * math.pi * G)  # kg/m^3

def H_of_a(a):
    return H0 * math.sqrt(Omega_r / a**4 + Omega_m / a**3 + Omega_L)

# Entropy density (arb. units): s ∝ T^3 ∝ a^-3
g_star = 3.36
pref_s = (2 * math.pi**2 / 45.0) * g_star
def s_of_a(a): return pref_s * (T0_GeV**3) / (a**3)

# Matter functional Θ = ν_m ρ_m + ν_r ρ_r (arb. units absorbed in α_χ)
def rho_m_of_a(a): return Omega_m * rho_c0 / a**3
def rho_r_of_a(a): return Omega_r * rho_c0 / a**4
nu_m, nu_r, nu_L = 1.0, 0.2, 0.0
def Theta_of_a(a): return nu_m * rho_m_of_a(a) + nu_r * rho_r_of_a(a) + nu_L * (1 - Omega_m - Omega_r) * rho_c0

# -------------------- Budget ODE coefficients --------------------
# TUNE these if you want different behavior/crossings:
alpha_H  = 0.8    # expansion coupling
alpha_s  = 1.0e-3 # entropy screening
alpha_ch = 1.0e-3 # matter drag
def S_micro(a): return 0.0  # source term (set to 0 for baseline)

def dI_da(a, I):
    H = H_of_a(a)
    rhs = -alpha_H * H * I - alpha_s * s_of_a(a) - alpha_ch * Theta_of_a(a) + S_micro(a)
    return rhs / (a * H)

# -------------------- Helper: integrate I_B(a) --------------------
def integrate_ibudget(I0, agrid):
    I = I0
    I_hist = [I0]
    for i in range(1, len(agrid)):
        a0, a1 = agrid[i-1], agrid[i]
        # RK2 midpoint
        k1 = dI_da(a0, I)
        k2 = dI_da(0.5*(a0+a1), I + 0.5*k1*(a1-a0))
        I = I + k2 * (a1 - a0)
        I_hist.append(I)
    return np.array(I_hist)

# ================================================================
# Figure 1: I_B vs redshift (recombination -> today), plus variations
# ================================================================
a_start = 1.0 / 1101.0
a_end   = 1.0
N       = 1800
agrid1  = np.linspace(a_start, a_end, N)
zgrid1  = 1.0/agrid1 - 1.0
I_init1 = 10.0  # initial budget at recombination (nats/area, illustrative)

# Baseline
I_hist_baseline = integrate_ibudget(I_init1, agrid1)

# Variation: weaker expansion coupling
alpha_H_save = alpha_H
alpha_H = 0.4
I_hist_lowH = integrate_ibudget(I_init1, agrid1)
alpha_H = alpha_H_save

# Variation: stronger matter drag
alpha_ch_save = alpha_ch
alpha_ch = 4.0e-3
I_hist_strongM = integrate_ibudget(I_init1, agrid1)
alpha_ch = alpha_ch_save

# Save data bundle
curves = {
    "baseline": I_hist_baseline.tolist(),
    "low_alpha_H": I_hist_lowH.tolist(),
    "strong_matter_drag": I_hist_strongM.tolist(),
}

out = {
    "cosmology": {
        "H0_s^-1": H0, "Omega_m": Omega_m, "Omega_r": Omega_r, "Omega_L": Omega_L, "rho_c0_kg_m3": rho_c0
    },
    "alphas": {"alpha_H": alpha_H, "alpha_s": alpha_s, "alpha_ch": alpha_ch},
    "initial_conditions": {"a_start": a_start, "z_start": 1.0/a_start - 1.0, "I_init": I_init1},
    "agrid": agrid1.tolist(),
    "I_curves": curves
}
with open("ibudget_example_data.json","w") as f:
    json.dump(out, f, indent=2)

# Plot fig 1
plt.figure(figsize=(8,5))
plt.plot(zgrid1, I_hist_baseline, label=r'baseline $\alpha_H=0.8,\ \alpha_\chi=10^{-3}$')
plt.plot(zgrid1, I_hist_lowH, '--', label=r'weaker expansion $\alpha_H=0.4$')
plt.plot(zgrid1, I_hist_strongM, ':', label=r'stronger matter drag $\alpha_\chi=4\times10^{-3}$')
plt.gca().invert_xaxis()
plt.xlabel("redshift z"); plt.ylabel(r"$I_B$ (nats / unit area, illustrative)")
plt.title(r"$I_B$ from recombination ($z\!\sim\!1100$) to today")
plt.legend()
plt.tight_layout()
plt.savefig("ibudget_example_curves.png", dpi=160)
plt.savefig("ibudget_example_curves.pdf")
plt.close()

# ================================================================
# Figure 2: Instability timescale τ_inst ~ σ / D_I
# ================================================================
D_I_values = np.logspace(0, 8, 30)   # s^-1
sigma_values = [1e-4, 1e-3, 1e-2]    # arbitrary units
plt.figure(figsize=(8,5))
for s in sigma_values:
    tau = np.array([s/Di for Di in D_I_values])
    plt.loglog(D_I_values, tau, label=fr"$\sigma={s}$")
plt.xlabel(r"index mobility $D_I$ (s$^{-1}$)"); plt.ylabel(r"instability timescale $\tau_{\rm inst}$ (s)")
plt.title(r"Example $\tau_{\rm inst}\approx \sigma/D_I$ bands")
plt.legend()
plt.tight_layout()
plt.savefig("inst_timescale_sweep.png", dpi=160)
plt.savefig("inst_timescale_sweep.pdf")
plt.close()

# ================================================================
# Figure 3: I_B vs Temperature with EW/QCD thresholds
# ================================================================
# Temperature range (GeV)
T_max, T_min = 1e2, 1e-4  # 100 GeV -> 100 keV
a_min2 = T0_GeV / T_max
a_max2 = T0_GeV / T_min
agrid2 = np.linspace(a_min2, a_max2, 2500)
Tgrid  = T0_GeV / agrid2

# Choose initial budget at high T; tune this + alphas to move the curve
I_init2 = 6.0
I_hist_T = integrate_ibudget(I_init2, agrid2)

# Thresholds from your paper
I_th_EW  = 1.10  # single spin-1 bridge
I_th_QCD = 2.08  # three spin-1/2 bridges

def find_cross(Tgrid, I_hist, Ith):
    for i in range(1, len(I_hist)):
        if I_hist[i-1] >= Ith and I_hist[i] <= Ith:
            t = (Ith - I_hist[i]) / (I_hist[i-1] - I_hist[i] + 1e-30)
            return Tgrid[i] * (1-t) + Tgrid[i-1] * t
    return None

T_cross_QCD = find_cross(Tgrid, I_hist_T, I_th_QCD)
T_cross_EW  = find_cross(Tgrid, I_hist_T, I_th_EW)

plt.figure(figsize=(8.2,5.2))
plt.semilogx(Tgrid, I_hist_T, label=r'$I_B(T)$ (illustrative)')
plt.axhline(I_th_QCD, color='tab:red',   linestyle='--', label=r'$I_{B,\mathrm{QCD}} \approx 2.08$')
plt.axhline(I_th_EW,  color='tab:green', linestyle='--', label=r'$I_{B,\mathrm{EW}} \approx 1.10$')

if T_cross_QCD is not None:
    plt.scatter([T_cross_QCD], [I_th_QCD], color='tab:red')
    plt.text(T_cross_QCD*1.05, I_th_QCD+0.05, f"cross @ {T_cross_QCD:.2e} GeV", fontsize=8)
if T_cross_EW is not None:
    plt.scatter([T_cross_EW], [I_th_EW], color='tab:green')
    plt.text(T_cross_EW*1.05, I_th_EW+0.05, f"cross @ {T_cross_EW:.2e} GeV", fontsize=8)

plt.gca().invert_xaxis()
plt.xlabel("Temperature T (GeV)"); plt.ylabel(r"$I_B$ (nats / unit area, illustrative)")
plt.title(r"Illustrative $I_B(T)$ with EW/QCD thresholds")
plt.legend()
plt.tight_layout()
plt.savefig("ibudget_vs_temperature.png", dpi=160)
plt.savefig("ibudget_vs_temperature.pdf")
plt.close()

print("Wrote:",
      "ibudget_example_curves.png, ibudget_example_curves.pdf,",
      "inst_timescale_sweep.png, inst_timescale_sweep.pdf,",
      "ibudget_vs_temperature.png, ibudget_vs_temperature.pdf,",
      "ibudget_example_data.json")

