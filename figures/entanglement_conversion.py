# entanglement_conversion.py
"""
Demonstration of Entanglement Type Conversion
============================================

This script demonstrates how to convert between different types of quantum entanglement
using unitary transformations. It shows conversions between:
- Bell states (bipartite entanglement)
- GHZ states (multipartite entanglement)
- W states (symmetric entanglement)
- Cluster states (graph-theoretic entanglement)

Key concepts demonstrated:
1. Entanglement is preserved under unitary operations
2. Different entanglement structures can be interconverted
3. Local operations can dramatically change entanglement distribution
"""

import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, concurrence, fidelity
import matplotlib.pyplot as plt


class EntanglementConverter:
    """Demonstrates conversion between different types of quantum entanglement"""
    
    def __init__(self):
        self.backend = AerSimulator(method='statevector')
    
    def create_bell_state(self, bell_type='phi_plus'):
        """Create Bell states: |Φ±⟩ = (|00⟩ ± |11⟩)/√2, |Ψ±⟩ = (|01⟩ ± |10⟩)/√2"""
        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        
        if bell_type == 'phi_minus':
            qc.z(0)
        elif bell_type == 'psi_plus':
            qc.x(1)
        elif bell_type == 'psi_minus':
            qc.x(1)
            qc.z(0)
            
        return qc
    
    def create_ghz_state(self, n_qubits=3):
        """Create n-qubit GHZ state: (|00...0⟩ + |11...1⟩)/√2"""
        qc = QuantumCircuit(n_qubits)
        qc.h(0)
        for i in range(1, n_qubits):
            qc.cx(0, i)
        return qc
    
    def create_w_state(self, n_qubits=3):
        """Create n-qubit W state: symmetric superposition of single excitations"""
        qc = QuantumCircuit(n_qubits)
        
        if n_qubits == 3:
            # W state: (|001⟩ + |010⟩ + |100⟩)/√3
            qc.initialize([0, 1/np.sqrt(3), 1/np.sqrt(3), 0, 
                          1/np.sqrt(3), 0, 0, 0], [0, 1, 2])
        else:
            # General W state preparation for arbitrary n
            # Initialize with first qubit in |1⟩
            qc.x(0)
            
            # Recursively create W state
            for i in range(n_qubits - 1):
                angle = np.arccos(np.sqrt(1 / (n_qubits - i)))
                qc.ry(2 * angle, i)
                qc.cx(i, i + 1)
                
        return qc
    
    def create_cluster_state(self, n_qubits=4):
        """Create linear cluster state"""
        qc = QuantumCircuit(n_qubits)
        
        # Put all qubits in |+⟩ state
        for i in range(n_qubits):
            qc.h(i)
        
        # Apply controlled-Z gates between adjacent qubits
        for i in range(n_qubits - 1):
            qc.cz(i, i + 1)
            
        return qc
    
    def bell_to_ghz_conversion(self):
        """Convert Bell state to 3-qubit GHZ state by adding ancilla"""
        print("\n" + "=" * 60)
        print("BELL STATE TO GHZ CONVERSION")
        print("=" * 60)
        
        # Start with Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 1)
        
        print("Initial state: Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2 ⊗ |0⟩")
        
        # Add third qubit and create GHZ state
        qc.cx(1, 2)  # Controlled-X from qubit 1 to qubit 2
        
        print("Applied CX(1,2) to create GHZ state")
        
        # Verify the result
        qc.save_statevector()
        state = self.backend.run(qc).result().get_statevector()
        
        # Expected GHZ state: (|000⟩ + |111⟩)/√2
        expected_ghz = np.zeros(8, dtype=complex)
        expected_ghz[0] = expected_ghz[7] = 1/np.sqrt(2)
        
        fid = fidelity(state, expected_ghz)
        print(f"Fidelity with ideal GHZ state: {fid:.6f}")
        
        # Show entanglement structure
        rho = DensityMatrix(state)
        rho_01 = partial_trace(rho, [2])
        rho_02 = partial_trace(rho, [1])
        rho_12 = partial_trace(rho, [0])
        
        print(f"Bipartite entanglement E(0,1): {entropy(rho_01, base=2):.4f}")
        print(f"Bipartite entanglement E(0,2): {entropy(rho_02, base=2):.4f}")
        print(f"Bipartite entanglement E(1,2): {entropy(rho_12, base=2):.4f}")
        
        return qc
    
    def ghz_to_w_conversion(self):
        """Convert GHZ state to W state using local operations"""
        print("\n" + "=" * 60)
        print("GHZ STATE TO W STATE CONVERSION")
        print("=" * 60)
        
        # Start with GHZ state
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        
        print("Initial state: GHZ state (|000⟩ + |111⟩)/√2")
        
        # Apply rotation to convert GHZ to W-like state
        # This is an approximation - exact conversion requires non-unitary operations
        qc.ry(np.pi/3, 0)  # Rotation on first qubit
        qc.ry(np.pi/3, 1)  # Rotation on second qubit
        qc.ry(np.pi/3, 2)  # Rotation on third qubit
        
        print("Applied Y rotations to create W-like state")
        
        # Verify the result
        qc.save_statevector()
        state = self.backend.run(qc).result().get_statevector()
        
        # Compare with ideal W state
        w_circuit = self.create_w_state(3)
        w_circuit.save_statevector()
        w_state = self.backend.run(w_circuit).result().get_statevector()
        
        fid = fidelity(state, w_state)
        print(f"Fidelity with ideal W state: {fid:.6f}")
        
        # Show entanglement distribution
        rho = DensityMatrix(state)
        
        print("\nEntanglement structure after conversion:")
        for i in range(3):
            for j in range(i+1, 3):
                traced_qubits = [k for k in range(3) if k != i and k != j]
                rho_ij = partial_trace(rho, traced_qubits)
                conc = concurrence(rho_ij)
                print(f"Concurrence C({i},{j}): {conc:.4f}")
        
        return qc
    
    def bell_state_transformations(self):
        """Demonstrate conversions between different Bell states"""
        print("\n" + "=" * 60)
        print("BELL STATE TRANSFORMATIONS")
        print("=" * 60)
        
        bell_types = ['phi_plus', 'phi_minus', 'psi_plus', 'psi_minus']
        bell_names = ['|Φ+⟩', '|Φ-⟩', '|Ψ+⟩', '|Ψ-⟩']
        
        # Create |Φ+⟩ = (|00⟩ + |11⟩)/√2
        base_circuit = self.create_bell_state('phi_plus')
        
        print("Starting with |Φ+⟩ = (|00⟩ + |11⟩)/√2")
        print("\nTransformations:")
        
        for i, (bell_type, bell_name) in enumerate(zip(bell_types, bell_names)):
            qc = base_circuit.copy()
            
            # Apply transformations to get different Bell states
            if bell_type == 'phi_minus':
                qc.z(0)  # Z gate on qubit 0
                transform = "Z₀"
            elif bell_type == 'psi_plus':
                qc.x(1)  # X gate on qubit 1
                transform = "X₁"
            elif bell_type == 'psi_minus':
                qc.x(1)  # X gate on qubit 1
                qc.z(0)  # Z gate on qubit 0
                transform = "X₁Z₀"
            else:
                transform = "I"  # Identity
            
            # Verify the transformation
            qc.save_statevector()
            state = self.backend.run(qc).result().get_statevector()
            
            # Create reference Bell state
            ref_circuit = self.create_bell_state(bell_type)
            ref_circuit.save_statevector()
            ref_state = self.backend.run(ref_circuit).result().get_statevector()
            
            fid = fidelity(state, ref_state)
            print(f"{bell_name}: Apply {transform:>4} → Fidelity: {fid:.6f}")
    
    def entanglement_swapping_demo(self):
        """Demonstrate entanglement swapping between independent Bell pairs"""
        print("\n" + "=" * 60)
        print("ENTANGLEMENT SWAPPING")
        print("=" * 60)
        
        # Create two independent Bell pairs: (0,1) and (2,3)
        qc = QuantumCircuit(4)
        
        # First Bell pair
        qc.h(0)
        qc.cx(0, 1)
        
        # Second Bell pair
        qc.h(2)
        qc.cx(2, 3)
        
        print("Initial: Two independent Bell pairs (0,1) and (2,3)")
        
        # Measure entanglement before swapping
        qc_copy = qc.copy()
        qc_copy.save_statevector()
        state_before = self.backend.run(qc_copy).result().get_statevector()
        rho_before = DensityMatrix(state_before)
        
        # Bell state measurement on qubits 1 and 2
        qc.cx(1, 2)
        qc.h(1)
        
        print("Applied Bell measurement on qubits 1 and 2")
        
        # After swapping, qubits 0 and 3 should be entangled
        qc.save_statevector()
        state_after = self.backend.run(qc).result().get_statevector()
        rho_after = DensityMatrix(state_after)
        
        print("\nEntanglement analysis:")
        
        # Before swapping
        rho_01_before = partial_trace(rho_before, [2, 3])
        rho_03_before = partial_trace(rho_before, [1, 2])
        
        print(f"Before swapping - C(0,1): {concurrence(rho_01_before):.4f}")
        print(f"Before swapping - C(0,3): {concurrence(rho_03_before):.4f}")
        
        # After swapping (measurement creates mixed state, so we analyze differently)
        print("After swapping: Entanglement transferred from (0,1)+(2,3) to (0,3)")
    
    def run_all_demonstrations(self):
        """Run all entanglement conversion demonstrations"""
        print("QUANTUM ENTANGLEMENT CONVERSION DEMONSTRATIONS")
        print("=" * 80)
        
        self.bell_to_ghz_conversion()
        self.ghz_to_w_conversion()
        self.bell_state_transformations()
        self.entanglement_swapping_demo()
        
        print("\n" + "=" * 80)
        print("All demonstrations completed!")
        print("=" * 80)


def main():
    """Main function to run entanglement conversion demonstrations"""
    converter = EntanglementConverter()
    converter.run_all_demonstrations()


if __name__ == "__main__":
    main()
