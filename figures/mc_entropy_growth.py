import math, random, collections
import numpy as np
import matplotlib.pyplot as plt

def invariant_dim(spins):
    current = {0:1}
    for j in spins:
        j2 = int(round(j*2))
        new = collections.defaultdict(int)
        for t2,mult in current.items():
            for s2 in range(abs(t2-j2), t2+j2+1, 2):
                new[s2]+=mult
        current=new
    return current.get(0,0)

def simulate(steps=10, runs=200):
    rng = random.Random(0)
    all_S = np.empty((runs, steps+1))
    for r in range(runs):
        spins=[0.5,0.5]  # even half-integer start
        all_S[r,0]=math.log(invariant_dim(spins))
        for t in range(1, steps+1):
            j_b = rng.choice([0.5,1.0,1.5])
            spins.extend([j_b,j_b])
            d0=invariant_dim(spins)
            all_S[r,t]=math.log(d0) if d0 else np.nan
    return all_S

steps=10
all_S=simulate(steps=steps, runs=300)
mean_S=np.nanmean(all_S, axis=0)

plt.figure(figsize=(6,4))
plt.plot(range(steps+1), mean_S, marker='o')
plt.xlabel("Monte‑Carlo step")
plt.ylabel("Mean $S_\\gamma$")
plt.title("Monte‑Carlo growth of boundary entropy")
plt.grid(True)
plt.tight_layout()
plt.savefig("mc_entropy_growth.pdf")
plt.show()
