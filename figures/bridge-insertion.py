import matplotlib.pyplot as plt
import networkx as nx

# Set up figure and axis
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Define common node positions for both panels
positions = {
    "A": (-1, 0), "B": (1, 0), "C": (0, 1.5), "O": (0, 0.5),
    "A'": (-1, 0), "B'": (1, 0), "C'": (0, 1.5), "O'": (0, 0.5)
}

# Panel 1: GHZ state (no bridge)
G1 = nx.Graph()
G1.add_edges_from([("A", "O"), ("B", "O"), ("C", "O")])
nx.draw(G1, pos=positions, ax=ax[0], with_labels=True, node_size=800, node_color='lightblue', edge_color='gray')
ax[0].set_title("GHZ State (No Entanglement in AB)")
ax[0].axis('off')

# Panel 2: Bridge inserted (AB entangled)
G2 = nx.Graph()
G2.add_edges_from([("A'", "O'"), ("B'", "O'"), ("C'", "O'"), ("A'", "B'")])
nx.draw(G2, pos=positions, ax=ax[1], with_labels=True, node_size=800, node_color='lightgreen', edge_color='gray')
ax[1].set_title("After $U_C$ / Bridge Insertion")
ax[1].axis('off')

# Add dashed boundary for AB region in both panels
for i in range(2):
    ax[i].plot([-1.4, 1.4], [-0.5, -0.5], 'k--', linewidth=1)
    ax[i].text(0, -0.8, "Boundary $\gamma$", ha='center', fontsize=10)

plt.tight_layout()
plt.savefig("contextual_activation_bridge_diagram.png")
plt.show()
