#!/usr/bin/env python3
"""--------------------
Estimate a Green–Kubo transport coefficient D from synthetic bridge-flux time series.

Features
- Synthetic flux generator (Ornstein–Uhlenbeck) to emulate stationary bridge current J(t)
- Multiprocessing over configurations for linear speedup
- FFT-based autocorrelation (O(N log N)) for efficiency
- Ensemble average over configurations
- Running Green–Kubo integral I(t) = ∫_0^t C(τ) dτ, plateau detection + window suggestion
- CSV + JSON summaries and PNG plots
- Reproducible via random seed, configurable via CLI
"""

import numpy as np
import argparse, json, time, math
from pathlib import Path
import multiprocessing as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------
# Synthetic OU generator (per-config)
# -----------------------------
def ou_process(nstep, dt, tau, sigma, rng):
    """Ornstein–Uhlenbeck J(t) with correlation time tau and stationary std sigma.
    Exact AR(1): J_t = a J_{t-1} + s * N(0,1), a=exp(-dt/tau), s=sigma*sqrt(1-a^2).
    """
    a = math.exp(-dt/tau)
    s = math.sqrt(max(0.0, sigma**2 * (1 - a*a)))
    J = np.empty(nstep, dtype=float)
    J[0] = s * rng.standard_normal()
    for t in range(1, nstep):
        J[t] = a * J[t-1] + s * rng.standard_normal()
    return J

def gen_flux_traj(args):
    """Worker function to generate one trajectory and return it."""
    idx, nstep, dt, tau, sigma, seed = args
    rng = np.random.default_rng(seed + idx*7919)
    return ou_process(nstep, dt, tau, sigma, rng)

# -----------------------------
# Correct FFT-based autocorrelation (unbiased)
# -----------------------------
def autocorr_fft_unbiased(x):
    """
    Unbiased autocorrelation via FFT:
    - zero-mean
    - zero-pad to nfft >= 2N-1 (next power of 2)
    - use rfft/irfft
    - divide by (N-k) so C[0] = variance
    """
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = x.size
    nfft = 1 << int(np.ceil(np.log2(2*n - 1)))  # next pow2 >= 2N-1
    f = np.fft.rfft(x, n=nfft)
    s = f * np.conjugate(f)
    acf_full = np.fft.irfft(s, n=nfft)
    acf = acf_full[:n].real
    # unbiased normalization
    denom = np.arange(n, 0, -1, dtype=float)
    acf = acf / denom
    return acf  # C[0] ≈ Var(x) (stationary)

# -----------------------------
# Robust plateau finder
# -----------------------------
def estimate_tau_from_acf(C, t):
    """Estimate correlation time as first index where C drops below C(0)/e."""
    c0 = C[0]
    if c0 <= 0:
        return max(t[-1] * 0.1, t[1]-t[0])
    target = c0 / math.e
    idx = np.where(C <= target)[0]
    if len(idx) == 0:
        return max(t[-1] * 0.1, t[1]-t[0])
    return t[idx[0]]

def plateau_window(I, t, C, rel_tol=0.02, width_factor=5.0):
    """
    Choose the earliest window of width ~ width_factor * tau_est where the
    running integral is stable: std(window)/max(|mean|,1e-3) < rel_tol.
    Returns (t_start, t_end).
    """
    tau_est = estimate_tau_from_acf(C, t)
    dt = t[1] - t[0]
    W = max(10, int(round((width_factor * tau_est) / dt)))
    n = len(I)
    best = None
    for i in range(0, n - W):
        w = I[i:i+W]
        m = float(np.mean(w))
        s = float(np.std(w))
        if s / max(abs(m), 1e-3) < rel_tol:
            best = (t[i], t[i+W-1])
            break
    if best is None:
        # fallback: a window centered at ~5 tau_est
        center = min(t[-1], 5.0 * tau_est)
        i0 = max(0, int(center/dt) - W//2)
        i1 = min(n-1, i0 + W - 1)
        return t[i0], t[i1]
    return best

# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nconf", type=int, default=200, help="number of independent configurations")
    ap.add_argument("--nstep", type=int, default=20000, help="time steps per configuration")
    ap.add_argument("--dt", type=float, default=1.0, help="time step")
    ap.add_argument("--tau", type=float, default=20.0, help="OU correlation time")
    ap.add_argument("--sigma", type=float, default=1.0, help="OU stationary std")
    ap.add_argument("--jobs", type=int, default=0, help="parallel workers (0=cpu_count)")
    ap.add_argument("--seed", type=int, default=12345, help="base random seed")
    ap.add_argument("--outdir", type=str, default="gk_out", help="output directory")
    ap.add_argument("--slope_tol", type=float, default=0.015, help="(legacy) plateau slope tolerance")
    ap.add_argument("--min_plateau", type=float, default=0.25, help="(legacy) min fraction of timeline for plateau")
    ap.add_argument("--rel_tol", type=float, default=0.02, help="plateau relative-variance tolerance")
    ap.add_argument("--width_factor", type=float, default=5.0, help="plateau width ~ factor * tau_est")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    jobs = args.jobs if args.jobs>0 else mp.cpu_count()

    # Generate flux trajectories
    worker_args = [(i, args.nstep, args.dt, args.tau, args.sigma, args.seed) for i in range(args.nconf)]
    if jobs == 1:
        series = [gen_flux_traj(a) for a in worker_args]
    else:
        with mp.Pool(processes=jobs) as pool:
            series = list(pool.map(gen_flux_traj, worker_args))
    t1 = time.time()

    # Autocorrelation per config (FFT, unbiased), then ensemble average
    acfs = [autocorr_fft_unbiased(J) for J in series]
    acfs = np.asarray(acfs)
    C = acfs.mean(axis=0)  # ensemble-averaged autocorrelation

    # Time axis and running Green–Kubo integral
    t = np.arange(len(C)) * args.dt
    I = np.cumsum(C) * args.dt

    # Plateau detection (new robust method; legacy args kept for API compatibility)
    t_start, t_end = plateau_window(I, t, C, rel_tol=args.rel_tol, width_factor=args.width_factor)

    # Estimate D from plateau mean
    mask = (t >= t_start) & (t <= t_end)
    D_est = float(np.mean(I[mask]))

    # Analytic OU check: D_true = sigma^2 * tau
    D_true = args.sigma**2 * args.tau
    ratio = D_est / D_true if D_true != 0 else float("nan")

    t2 = time.time()

    # Save CSVs
    import csv
    with open(outdir/"gk_autocorr.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t","C(t)"]); w.writerows(zip(t, C))
    with open(outdir/"gk_running_integral.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t","I(t)"]); w.writerows(zip(t, I))

    # Plots
    plt.figure()
    plt.plot(t, C)
    plt.xlabel("t"); plt.ylabel("C(t)"); plt.title("Autocorrelation (ensemble-averaged)")
    plt.savefig(outdir/"gk_autocorr.png", bbox_inches="tight"); plt.close()

    plt.figure()
    plt.plot(t, I, label="I(t)")
    plt.axvspan(t_start, t_end, alpha=0.2, label="plateau")
    plt.xlabel("t"); plt.ylabel("I(t) = ∫_0^t C(τ) dτ")
    plt.title("Green–Kubo Running Integral (plateau shaded)")
    plt.legend()
    plt.savefig(outdir/"gk_running_integral.png", bbox_inches="tight"); plt.close()

    # Summary + timings
    summary = {
        "nconf": args.nconf, "nstep": args.nstep, "dt": args.dt,
        "tau": args.tau, "sigma": args.sigma, "jobs": jobs,
        "D_estimate": D_est,
        "D_true_OU": D_true,
        "relative_error": (D_est - D_true) / D_true if D_true != 0 else None,
        "plateau_window": [float(t_start), float(t_end)],
        "timings_seconds": {
            "generate_series": t1 - t0,
            "autocorr_and_integral": t2 - t1,
            "total": t2 - t0
        }
    }
    with open(outdir/"gk_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Scaling log (very simple; users can sweep externally)
    with open(outdir/"gk_scaling_log.json", "w") as f:
        json.dump({
            "O_cost": "O(nconf * nstep log nstep) dominated by FFT",
            "nconf": args.nconf, "nstep": args.nstep, "jobs": jobs,
            "timings_seconds": summary["timings_seconds"]
        }, f, indent=2)

    # Also print summary to stdout (CLI compatibility)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
