# verify_kappa_twirled_2x2_simple.py
import numpy as np, math

def classical_twirled_kappa(s2_target=0.2635, pi=None):
    """
    Build the SU(2)-twirled 2x2 transfer on the multiplicity block:
        E = λ I + (1-λ) 1 π^T
    with stationary π (default [1/4, 3/4]) and subleading eigenvalue λ = s2_target.
    Return s2 and kappa = -2 ln s2.
    """
    if pi is None:
        pi = np.array([1.0, 3.0], float); pi /= pi.sum()
    lam = float(s2_target)
    E = np.zeros((2,2), float)
    I = np.eye(2)
    for j in range(2):
        E[:, j] = lam * I[:, j] + (1.0 - lam) * pi

    # column sums are 1; π is a right eigenvector with eigenvalue 1; the other eigenvalue is λ
    colsum = E.sum(axis=0)
    assert np.allclose(colsum, 1.0, atol=1e-12)
    eigvals = np.linalg.eigvals(E)
    eigvals = np.sort(np.real_if_close(eigvals))
    s2 = float(sorted([abs(x) for x in eigvals])[-2])  # subleading magnitude
    kappa = -2.0 * math.log(s2)
    return E, s2, kappa, pi

if __name__ == "__main__":
    E, s2, kappa, pi = classical_twirled_kappa()
    print("E =\n", E)
    print(f"s2 = {s2:.9f}")
    print(f"kappa = {kappa:.9f}")
    print("pi =", pi)
