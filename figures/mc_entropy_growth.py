import math, random, collections, numpy as np
import matplotlib.pyplot as plt

# ------- helper: dim Inv -------
def invariant_dim(spins):
    """Return dim Inv(⊗ V_j) for a list of spins."""
    current = {0: 1}                         # maps total spin×2 → mult
    for j in spins:
        j2 = int(round(j*2))
        nxt = collections.defaultdict(int)
        for t2, mult in current.items():
            for s2 in range(abs(t2-j2), t2+j2+1, 2):
                nxt[s2] += mult
        current = nxt
    return current.get(0, 0)

# ------- helper: admissibility -------
def is_admissible(v_spins, new_spin):
    """Bridge admissible iff invariant subspace survives."""
    return invariant_dim(v_spins + [new_spin]) > 0

# ------- one Monte-Carlo history -------
def run_history(steps=10, rng=random.Random()):
    u_spins, v_spins = [0.5], [0.5]          # initial even-parity cut
    S_vals  = [math.log(invariant_dim(u_spins + v_spins))]
    acc = rej = 0
    for _ in range(steps):
        while True:                          # rejection sampling
            j_b = rng.choice([0.5, 1.0, 1.5])
            if is_admissible(u_spins, j_b) and is_admissible(v_spins, j_b):
                # accept the bridge
                u_spins.append(j_b); v_spins.append(j_b)
                acc += 1
                S_vals.append(math.log(invariant_dim(u_spins + v_spins)))
                break
            else:
                rej += 1                     # reject and resample
    return np.array(S_vals), acc, rej

# ------- aggregate many histories -------
runs, steps = 300, 10
all_S = np.empty((runs, steps + 1))
acc_tot = rej_tot = 0
rng = random.Random(0)

for r in range(runs):
    S, acc, rej = run_history(steps, rng)
    all_S[r] = S
    acc_tot += acc
    rej_tot += rej

mean_S = np.mean(all_S, axis=0)
print(f"accepted {acc_tot}, rejected {rej_tot}, acceptance {acc_tot/(acc_tot+rej_tot):.3f}")

# ------- plot & save -------
plt.figure(figsize=(6,4))
plt.plot(range(steps + 1), mean_S, marker='o')
plt.xlabel("Monte-Carlo step")
plt.ylabel(r"Mean $S_\gamma$")
plt.title("Monte-Carlo growth of boundary entropy")
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/mc_entropy_growth.pdf", dpi=300)
