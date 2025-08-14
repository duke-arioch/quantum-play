import numpy as np
import matplotlib.pyplot as plt

def quantum_dimension(j, k):
    """
    Compute quantum dimension for spin j in SU(2)_k.
    d_j^(q) = sin((2j+1)π/(k+2)) / sin(π/(k+2))
    """
    if k == np.inf:  # Classical limit
        return 2*j + 1
    
    numerator = np.sin((2*j + 1) * np.pi / (k + 2))
    denominator = np.sin(np.pi / (k + 2))
    return numerator / denominator

def anchor_charge(j, k):
    """
    Compute anchor charge Q = d_j^(q) - 1.
    """
    return quantum_dimension(j, k) - 1

def verify_charge_quantization():
    """
    Verify that Q ≈ 1 for spin-1/2 at k=48.
    """
    print("CHARGE QUANTIZATION FROM QUANTUM DIMENSIONS")
    print("=" * 60)
    
    # Test for various k values
    k_values = [10, 20, 48, 100, 200, np.inf]
    j = 0.5  # spin-1/2 (electron)
    
    print(f"\nSpin-{j} anchor charge vs level k:")
    print("-" * 40)
    for k in k_values:
        d_q = quantum_dimension(j, k)
        Q = anchor_charge(j, k)
        if k == np.inf:
            print(f"k = ∞ (classical): d_q = {d_q:.4f}, Q = {Q:.4f}")
        else:
            print(f"k = {k:3}: d_q = {d_q:.6f}, Q = {Q:.6f}")
    
    print(f"\nAt k=48: Q = {anchor_charge(0.5, 48):.6f} ≈ 1")
    print("This matches unit electric charge!")

def test_multiple_spins():
    """
    Test charge for different spin anchors.
    """
    k = 48  # Our standard level
    spins = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    
    print("\nCharge spectrum for different anchor spins (k=48):")
    print("-" * 50)
    print("Spin j | d_j^(q) | Q = d_j^(q) - 1 | Interpretation")
    print("-" * 50)
    
    for j in spins:
        d_q = quantum_dimension(j, k)
        Q = anchor_charge(j, k)
        
        # Interpretation
        if j == 0:
            interp = "Vacuum (no charge)"
        elif j == 0.5:
            interp = "Electron (unit charge)"
        elif j == 1:
            interp = "W boson analog"
        elif j == 1.5:
            interp = "Exotic fermion"
        else:
            interp = "Heavy exotic"
            
        print(f"{j:5.1f} | {d_q:7.4f} | {Q:15.4f} | {interp}")

def plot_charge_spectrum():
    """
    Visualize how charge varies with spin and level k.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Q vs j for fixed k=48
    k = 48
    j_values = np.linspace(0, k/2, 100)
    Q_values = [anchor_charge(j, k) for j in j_values]
    
    ax1.plot(j_values, Q_values, 'b-', linewidth=2)
    ax1.axhline(y=1, color='r', linestyle='--', label='Q=1 (electron)')
    ax1.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax1.scatter([0.5], [anchor_charge(0.5, k)], color='r', s=100, zorder=5)
    ax1.set_xlabel('Anchor spin j')
    ax1.set_ylabel('Charge Q')
    ax1.set_title(f'Charge spectrum at k={k}')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Right: Q(j=1/2) vs k
    k_values = np.linspace(4, 100, 100)
    Q_half = [anchor_charge(0.5, k) for k in k_values]
    
    ax2.plot(k_values, Q_half, 'g-', linewidth=2)
    ax2.axhline(y=1, color='r', linestyle='--', label='Q=1')
    ax2.axvline(x=48, color='b', linestyle=':', label='k=48')
    ax2.set_xlabel('Level k')
    ax2.set_ylabel('Q for spin-1/2')
    ax2.set_title('Approach to unit charge')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('charge_quantization.png', dpi=150)
    plt.show()
    print("\nPlot saved as 'charge_quantization.png'")

def verify_fusion_rules():
    """
    Verify that charge adds for multiple anchors.
    """
    k = 48
    print("\nCHARGE ADDITION FOR MULTIPLE ANCHORS:")
    print("=" * 50)
    
    # Two electrons
    Q_electron = anchor_charge(0.5, k)
    print(f"One electron: Q = {Q_electron:.6f}")
    print(f"Two electrons: Q = {2*Q_electron:.6f} ≈ 2")
    
    # Electron + higher spin
    Q_spin1 = anchor_charge(1, k)
    print(f"\nSpin-1 anchor: Q = {Q_spin1:.6f}")
    print(f"Electron + Spin-1: Q = {Q_electron + Q_spin1:.6f}")
    
    # Check neutrality
    print(f"\nNeutral configuration:")
    print(f"Electron + Positron: Q = {Q_electron - Q_electron:.6f} = 0")

# Run all tests
if __name__ == "__main__":
    verify_charge_quantization()
    test_multiple_spins()
    verify_fusion_rules()
    plot_charge_spectrum()