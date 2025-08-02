import math
import matplotlib.pyplot as plt

m_vals = list(range(1, 31))
delta_s = [math.log((4*m+2)/(m+2)) for m in m_vals]

plt.figure()
plt.plot(m_vals, delta_s, marker='o')
plt.axhline(math.log(4), linestyle='--')
plt.xlabel("m  (number of spin-1/2 pairs before move)")
plt.ylabel("ΔS  (entropy jump)")
plt.title("ΔS vs m: approach to ln 4")

plt.show()