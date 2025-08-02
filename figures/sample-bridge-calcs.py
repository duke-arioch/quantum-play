# response to:
comment = """
Intersecting Bridges: Even a simple example calculation for two intersecting spin-1/2 bridges would strengthen the paper significantly.
Numerical Comparison: If possible, include a comparison with existing spin-foam numerical results to demonstrate the utility of the analytical formula.

can we run that with a simple python program and grab the results?
"""


import itertools

def total_spin_range(j1, j2):
    return range(abs(j1 - j2), j1 + j2 + 1, 2)

def multiplicity_singlet(spins):
    """
    Return the number of singlet channels in the tensor product
    of SU(2) reps labelled by `spins` (list of half‑integers).
    """
    twice = [int(2*j) for j in spins]         # work with integers
    partial = {0: 1}                          # multiplicity dictionary

    for j in twice:
        new = {}
        for t, mult in partial.items():
            for k in range(abs(t - j), t + j + 1, 2):
                new[k] = new.get(k, 0) + mult
        partial = new
    return partial.get(0, 0)

# simple progression: start with 2 half‑integer edges (even parity)
configs = {
    "no bridge (2×½)": [0.5, 0.5],
    "+1 bridge":       [0.5, 0.5, 0.5, 0.5],
    "+2 bridges":      [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
}

results = [(name, len(spins), multiplicity_singlet(spins))
           for name, spins in configs.items()]

import pandas as pd, ace_tools as tools
tools.display_dataframe_to_user("Singlet multiplicities for successive bridges", 
                                pd.DataFrame(results, columns=["configuration",
                                                               "edges", 
                                                               "dim Inv"]))
