import matplotlib.pyplot as plt
import networkx as nx

# --- Step 1: Define example entropy-equivalence layers ---
# We'll model this as a layered directed graph: horizontal edges ΔS=0, vertical edges ΔS>0

G = nx.DiGraph()

# Example: 3 layers, each with 4 configurations
layers = {
    0: ['A1', 'A2', 'A3', 'A4'],
    1: ['B1', 'B2', 'B3', 'B4'],
    2: ['C1', 'C2', 'C3', 'C4']
}

# Add nodes with layer attribute
for layer, nodes in layers.items():
    for n in nodes:
        G.add_node(n, layer=layer)

# Add horizontal entropy-neutral edges
for layer, nodes in layers.items():
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            G.add_edge(nodes[i], nodes[j], deltaS=0)
            G.add_edge(nodes[j], nodes[i], deltaS=0)

# Add vertical ΔS>0 edges
for l in range(len(layers)-1):
    for src in layers[l]:
        for dst in layers[l+1]:
            G.add_edge(src, dst, deltaS=1)

# --- Step 2: Layout ---
pos = {}
layer_gap = 5
node_gap = 2
for layer, nodes in layers.items():
    for i, n in enumerate(nodes):
        pos[n] = (i * node_gap, -layer * layer_gap)

# --- Step 3: Draw ---
plt.figure(figsize=(10, 6))
# Draw ΔS=0 edges in blue
zero_edges = [(u,v) for u,v,d in G.edges(data=True) if d['deltaS']==0]
nx.draw_networkx_edges(G, pos, edgelist=zero_edges, edge_color='blue', alpha=0.5, arrows=False)

# Draw ΔS>0 edges in red (directed)
positive_edges = [(u,v) for u,v,d in G.edges(data=True) if d['deltaS']>0]
nx.draw_networkx_edges(G, pos, edgelist=positive_edges, edge_color='red', alpha=0.6, arrows=True)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_color='lightgray', node_size=800, edgecolors='black')
nx.draw_networkx_labels(G, pos, font_size=10)

# Labels for layers
for layer in layers:
    plt.text(-2.5, -layer*layer_gap, f"Layer {layer}", fontsize=12, verticalalignment='center')

plt.title("Entropy-Equivalence Layers: ΔS=0 (blue), ΔS>0 (red)", fontsize=14)
plt.axis('off')

# Save PDF for paper integration
pdf_path = "entropy_layers_diagram.pdf"
plt.savefig(pdf_path, bbox_inches='tight')

pdf_path
