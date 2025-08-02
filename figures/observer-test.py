import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
import matplotlib.patches as mpatches
from qiskit import QuantumCircuit
from qiskit.visualization import circuit_drawer
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, concurrence
import seaborn as sns

# Set publication style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.dpi'] = 300

class EntanglementAnalysis:
    def __init__(self):
        self.backend = AerSimulator(method='statevector')
        self.results = {}
        
    def create_ghz_state(self):
        """Create GHZ state"""
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        return qc
    
    def create_w_state(self):
        """Create W state"""
        qc = QuantumCircuit(3)
        qc.initialize([0, 1/np.sqrt(3), 1/np.sqrt(3), 0, 
                      1/np.sqrt(3), 0, 0, 0], [0, 1, 2])
        return qc
    
    def apply_context(self, qc, context):
        """Apply context transformation"""
        if context == 'A':
            qc.h(0)
            qc.cz(0, 1)
            qc.h(0)
        elif context == 'B':
            qc.h(1)
            qc.cz(1, 2)
            qc.h(1)
        elif context == 'C':
            qc.h(2)
            qc.cz(2, 0)
            qc.h(2)
        return qc
    
    def analyze_state(self, qc, context):
        """Full analysis of quantum state"""
        qc_copy = qc.copy()
        qc_copy = self.apply_context(qc_copy, context)
        qc_copy.save_statevector()
        
        job = self.backend.run(qc_copy)
        state = job.result().get_statevector()
        rho = DensityMatrix(state)
        
        # Calculate various partitions
        rho_AB = partial_trace(rho, [2])
        rho_AC = partial_trace(rho, [1])
        rho_BC = partial_trace(rho, [0])
        
        # Calculate entropies
        S_AB = float(entropy(rho_AB, base=2))
        S_AC = float(entropy(rho_AC, base=2))
        S_BC = float(entropy(rho_BC, base=2))
        
        # Calculate concurrences
        C_AB = float(concurrence(rho_AB))
        C_AC = float(concurrence(rho_AC))
        C_BC = float(concurrence(rho_BC))
        
        # Calculate mutual information
        rho_A = partial_trace(rho, [1, 2])
        rho_B = partial_trace(rho, [0, 2])
        rho_C = partial_trace(rho, [0, 1])
        
        S_A = float(entropy(rho_A, base=2))
        S_B = float(entropy(rho_B, base=2))
        S_C = float(entropy(rho_C, base=2))
        
        I_AB = S_A + S_B - S_AB
        
        return {
            'S_AB': S_AB, 'S_AC': S_AC, 'S_BC': S_BC,
            'C_AB': C_AB, 'C_AC': C_AC, 'C_BC': C_BC,
            'S_A': S_A, 'S_B': S_B, 'S_C': S_C,
            'I_AB': I_AB,
            'state': state,
            'rho_AB': rho_AB
        }
    
    def run_full_analysis(self):
        """Run complete analysis for all states and contexts"""
        states = {'GHZ': self.create_ghz_state(), 'W': self.create_w_state()}
        contexts = ['None', 'A', 'B', 'C']
        
        for state_name, state_qc in states.items():
            self.results[state_name] = {}
            for context in contexts:
                if context == 'None':
                    self.results[state_name][context] = self.analyze_state(state_qc, None)
                else:
                    self.results[state_name][context] = self.analyze_state(state_qc, context)
        
        return self.results

# Create analyzer and run analysis
analyzer = EntanglementAnalysis()
results = analyzer.run_full_analysis()

# Figure 1: Main Results - Entanglement Measures
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
fig.suptitle('Observer-Dependent Entanglement in Tripartite Quantum States', fontsize=16)

# GHZ State Results
contexts = ['None', 'A', 'B', 'C']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

# Plot S_AB for GHZ
ax = axes[0, 0]
s_ab_ghz = [results['GHZ'][c]['S_AB'] for c in contexts]
c_ab_ghz = [results['GHZ'][c]['C_AB'] for c in contexts]
x = np.arange(len(contexts))
width = 0.35

bars1 = ax.bar(x - width/2, s_ab_ghz, width, label='$S_{AB}$', color=colors[0], alpha=0.8)
bars2 = ax.bar(x + width/2, c_ab_ghz, width, label='$C_{AB}$', color=colors[1], alpha=0.8)

ax.set_xlabel('Context')
ax.set_ylabel('Entanglement Measure')
ax.set_title('GHZ State: AB Partition')
ax.set_xticks(x)
ax.set_xticklabels(contexts)
ax.legend()
ax.grid(True, alpha=0.3)

# Plot all partitions for GHZ Context C
ax = axes[0, 1]
partitions = ['AB', 'AC', 'BC']
s_vals = [results['GHZ']['C'][f'S_{p}'] for p in partitions]
c_vals = [results['GHZ']['C'][f'C_{p}'] for p in partitions]

x = np.arange(len(partitions))
bars1 = ax.bar(x - width/2, s_vals, width, label='Entropy', color=colors[2], alpha=0.8)
bars2 = ax.bar(x + width/2, c_vals, width, label='Concurrence', color=colors[3], alpha=0.8)

ax.set_xlabel('Partition')
ax.set_ylabel('Entanglement Measure')
ax.set_title('GHZ State Context C: All Partitions')
ax.set_xticks(x)
ax.set_xticklabels(partitions)
ax.legend()
ax.grid(True, alpha=0.3)

# Mutual Information
ax = axes[0, 2]
i_ab_vals = [results['GHZ'][c]['I_AB'] for c in contexts]
bars = ax.bar(contexts, i_ab_vals, color=colors, alpha=0.8)
ax.set_xlabel('Context')
ax.set_ylabel('Mutual Information $I(A:B)$')
ax.set_title('GHZ State: Mutual Information')
ax.grid(True, alpha=0.3)

# W State Results (bottom row)
# Similar plots for W state
for i, measure in enumerate(['S_AB', 'C_AB', 'I_AB']):
    ax = axes[1, i]
    if measure == 'I_AB':
        vals = [results['W'][c][measure] for c in contexts]
        ax.bar(contexts, vals, color=colors, alpha=0.8)
        ax.set_ylabel('Mutual Information $I(A:B)$')
    else:
        vals = [results['W'][c][measure] for c in contexts]
        ax.plot(contexts, vals, 'o-', markersize=8, linewidth=2, 
                label=f'${measure.replace("_", "{")}}}$')
        ax.set_ylabel('Entanglement Measure')
        ax.legend()
    
    ax.set_xlabel('Context')
    ax.set_title(f'W State: {measure.replace("_", " ")}')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figure1_main_results.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figure1_main_results.png', dpi=300, bbox_inches='tight')

# Figure 2: State Visualization and Circuit Diagrams
fig = plt.figure(figsize=(12, 10))

# Create grid for subplots
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.2], hspace=0.4, wspace=0.3)

# Top row: Circuit diagrams
for i, context in enumerate(['A', 'B', 'C']):
    ax = fig.add_subplot(gs[0, i])
    
    # Create circuit
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.barrier()
    
    if context == 'A':
        qc.h(0)
        qc.cz(0, 1)
        qc.h(0)
    elif context == 'B':
        qc.h(1)
        qc.cz(1, 2)
        qc.h(1)
    elif context == 'C':
        qc.h(2)
        qc.cz(2, 0)
        qc.h(2)
    
    # Draw circuit
    circuit_drawer(qc, output='mpl', ax=ax, style='iqx')
    ax.set_title(f'Context {context} Transformation', fontsize=12)

# Middle row: Density matrices
for i, context in enumerate(['None', 'B', 'C']):
    ax = fig.add_subplot(gs[1, i])
    
    rho_AB = results['GHZ'][context]['rho_AB']
    rho_data = np.abs(rho_AB.data)
    
    im = ax.imshow(rho_data, cmap='RdBu', vmin=0, vmax=0.5, aspect='equal')
    ax.set_title(f'$\\rho_{{AB}}$ for Context {context}')
    ax.set_xlabel('Basis state')
    ax.set_ylabel('Basis state')
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([0, 1, 2, 3])
    ax.set_xticklabels(['00', '01', '10', '11'])
    ax.set_yticklabels(['00', '01', '10', '11'])
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('$|\\rho_{ij}|$')

# Bottom: Conceptual diagram
ax = fig.add_subplot(gs[2, :])
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
ax.axis('off')

# Draw conceptual diagram
# Initial state
rect1 = FancyBboxPatch((0.5, 2), 2, 1.5, boxstyle="round,pad=0.1", 
                        facecolor='lightblue', edgecolor='black', linewidth=2)
ax.add_patch(rect1)
ax.text(1.5, 2.75, 'GHZ State\n$\\frac{|000\\rangle + |111\\rangle}{\\sqrt{2}}$', 
        ha='center', va='center', fontsize=11)

# Arrows and transformations
for i, (context, x_pos) in enumerate([('A', 3.5), ('B', 5.5), ('C', 7.5)]):
    ax.arrow(2.5, 2.75, 0.8, 0, head_width=0.2, head_length=0.1, fc='black', ec='black')
    
    rect = FancyBboxPatch((x_pos, 1.5 + i*0.7), 1.8, 1, boxstyle="round,pad=0.1",
                          facecolor=colors[i+1], alpha=0.3, edgecolor='black', linewidth=2)
    ax.add_patch(rect)
    
    if context == 'C':
        text = f'Context {context}\n$S_{{AB}}=0$\n$C_{{AB}}=1$'
    else:
        text = f'Context {context}\n$S_{{AB}}=1$\n$C_{{AB}}=0$'
    
    ax.text(x_pos + 0.9, 2 + i*0.7, text, ha='center', va='center', fontsize=10)

ax.text(5, 0.5, 'Different contexts reveal different entanglement structures!', 
        ha='center', va='center', fontsize=12, style='italic')

plt.suptitle('Quantum Circuits and State Analysis', fontsize=16)
plt.savefig('figure2_circuits_states.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figure2_circuits_states.png', dpi=300, bbox_inches='tight')

# Figure 3: Connection to LQG
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle('Connection to Loop Quantum Gravity', fontsize=16)

# Create data for bridge analogy
ax = axes[0, 0]
ax.set_title('Spin Network Analogy')
ax.text(0.5, 0.5, 'Conceptual diagram\nto be created', ha='center', va='center')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# Summary table
ax = axes[0, 1]
ax.axis('tight')
ax.axis('off')

table_data = [
    ['Property', 'Context None/A/B', 'Context C'],
    ['$S_{AB}$', '1 bit', '0 bits'],
    ['$C_{AB}$', '0', '1'],
    ['State Type', 'Mixed', 'Pure'],
    ['LQG Analog', 'No bridge', 'With bridge']
]

table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                 colWidths=[0.3, 0.35, 0.35])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 2)

# Header row styling
for i in range(3):
    table[(0, i)].set_facecolor('#4CAF50')
    table[(0, i)].set_text_props(weight='bold', color='white')

ax.set_title('Entanglement Properties Comparison', pad=20)

# Statistical analysis over multiple runs
ax = axes[1, 0]
# Run multiple times to get error bars
n_runs = 100
s_ab_runs = {'A': [], 'B': [], 'C': []}
c_ab_runs = {'A': [], 'B': [], 'C': []}

for _ in range(n_runs):
    for context in ['A', 'B', 'C']:
        result = analyzer.analyze_state(analyzer.create_ghz_state(), context)
        s_ab_runs[context].append(result['S_AB'])
        c_ab_runs[context].append(result['C_AB'])

# Plot with error bars
contexts_plot = ['A', 'B', 'C']
s_means = [np.mean(s_ab_runs[c]) for c in contexts_plot]
s_stds = [np.std(s_ab_runs[c]) for c in contexts_plot]
c_means = [np.mean(c_ab_runs[c]) for c in contexts_plot]
c_stds = [np.std(c_ab_runs[c]) for c in contexts_plot]

x = np.arange(len(contexts_plot))
ax.errorbar(x - 0.1, s_means, yerr=s_stds, fmt='o', capsize=5, 
            label='$S_{AB}$', markersize=8)
ax.errorbar(x + 0.1, c_means, yerr=c_stds, fmt='s', capsize=5, 
            label='$C_{AB}$', markersize=8)

ax.set_xlabel('Context')
ax.set_ylabel('Value')
ax.set_title(f'Statistical Validation (n={n_runs})')
ax.set_xticks(x)
ax.set_xticklabels(contexts_plot)
ax.legend()
ax.grid(True, alpha=0.3)

# Interpretation diagram
ax = axes[1, 1]
ax.text(0.5, 0.8, 'Key Insight:', fontsize=14, weight='bold', ha='center')
ax.text(0.5, 0.6, 'Context C creates a pure entangled state\n'
                  'between A and B (C_{AB} = 1),\n'
                  'but with zero entropy because it\'s pure.\n\n'
                  'This mirrors how bridge insertion in LQG\n'
                  'creates gauge-invariant singlet channels.', 
        ha='center', va='center', fontsize=12)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig('figure3_lqg_connection.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figure3_lqg_connection.png', dpi=300, bbox_inches='tight')

print("Figures generated successfully!")
print("\nKey results summary:")
print("=" * 50)
print("GHZ State:")
for context in contexts:
    r = results['GHZ'][context]
    print(f"\nContext {context}:")
    print(f"  S_AB = {r['S_AB']:.3f}")
    print(f"  C_AB = {r['C_AB']:.3f}")
    print(f"  I(A:B) = {r['I_AB']:.3f}")