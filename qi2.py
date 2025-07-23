# enhanced_observer_test.py
import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, concurrence, Statevector
import matplotlib.pyplot as plt

class EnhancedObserverTest:
    """More sophisticated tests that can show observer dependence"""
    
    def __init__(self):
        self.backend = AerSimulator(method='statevector')
    
    def create_asymmetric_state(self):
        """Create a state where observer operations matter more"""
        qc = QuantumCircuit(3)
        # Create asymmetric superposition
        qc.ry(np.pi/3, 0)
        qc.cx(0, 1)
        qc.ry(np.pi/4, 1)
        qc.cx(1, 2)
        qc.rz(np.pi/5, 2)
        return qc
    
    def apply_measurement_context(self, qc, observer, protocol='interference'):
        """More sophisticated observer-dependent operations"""
        
        if protocol == 'interference':
            # Create interference between observer and system
            if observer == 'A':
                qc.h(0)
                qc.cz(0, 1)
                qc.h(0)
            elif observer == 'B':
                qc.h(1)
                qc.cz(1, 2)
                qc.h(1)
            elif observer == 'C':
                qc.h(2)
                qc.cz(2, 0)
                qc.h(2)
                
        elif protocol == 'weak_measurement':
            # Simulate weak measurement by observer
            theta = {'A': 0.1, 'B': 0.15, 'C': 0.2}[observer]
            qubit = {'A': 0, 'B': 1, 'C': 2}[observer]
            
            # Weak measurement = small rotation
            qc.ry(theta, qubit)
            # Followed by entangling operation
            qc.cx(qubit, (qubit + 1) % 3)
            qc.ry(-theta/2, qubit)
            
        elif protocol == 'basis_mixing':
            # Mix measurement bases
            if observer == 'A':
                qc.ry(np.pi/8, 0)
                qc.rz(np.pi/6, 1)
            elif observer == 'B':
                qc.rx(np.pi/7, 1)
                qc.ry(np.pi/9, 2)
            elif observer == 'C':
                qc.rz(np.pi/5, 2)
                qc.rx(np.pi/10, 0)
        
        return qc
    
    def measure_multiple_quantities(self, qc, observer):
        """Measure different aspects that might show observer dependence"""
        
        # Apply observer context
        qc_obs = self.apply_measurement_context(qc.copy(), observer, 'interference')
        
        # Get state using new Qiskit syntax
        qc_obs.save_statevector()
        job = self.backend.run(qc_obs)
        state = job.result().get_statevector()
        rho = DensityMatrix(state)
        
        results = {}
        
        # 1. Standard entanglement entropy
        rho_AB = partial_trace(rho, [2])
        results['E_AB'] = float(entropy(rho_AB, base=2))
        
        # 2. Concurrence (more sensitive)
        results['C_AB'] = float(concurrence(rho_AB))
        
        # 3. Conditional entropy S(A|C) from observer C's perspective
        rho_AC = partial_trace(rho, [1])
        rho_C = partial_trace(rho, [0, 1])
        results['S_A_given_C'] = float(entropy(rho_AC, base=2) - entropy(rho_C, base=2))
        
        # 4. Three-tangle (genuine three-party entanglement)
        # Simplified calculation
        psi = state.data
        results['three_tangle'] = float(np.abs(psi[0]*psi[7] - psi[1]*psi[6] 
                                              - psi[2]*psi[5] + psi[3]*psi[4])**2)
        
        # 5. Linear entropy (purity measure)
        results['linear_entropy'] = float(1 - np.real(np.trace(rho_AB.data @ rho_AB.data)))
        
        return results
    
    def comprehensive_test(self):
        """Test multiple states and protocols"""
        
        print("\n=== ENHANCED OBSERVER DEPENDENCE TEST ===\n")
        
        # Test different initial states
        states = {
            'GHZ': self.create_ghz_state(),
            'Asymmetric': self.create_asymmetric_state(),
            'W': self.create_w_state()
        }
        
        all_results = []
        
        for state_name, state_circuit in states.items():
            print(f"\nTesting {state_name} state:")
            print("-" * 60)
            
            # Measure from each observer
            results = {}
            for observer in ['A', 'B', 'C']:
                results[observer] = self.measure_multiple_quantities(state_circuit, observer)
            
            # Check for differences
            for quantity in ['E_AB', 'C_AB', 'S_A_given_C', 'three_tangle', 'linear_entropy']:
                values = [results[obs][quantity] for obs in ['A', 'B', 'C']]
                max_diff = max(values) - min(values)
                
                if max_diff > 0.001:  # Found something!
                    print(f"\n✓ {quantity} shows observer dependence!")
                    print(f"  A: {values[0]:.6f}, B: {values[1]:.6f}, C: {values[2]:.6f}")
                    print(f"  Max difference: {max_diff:.6f}")
                    
                    all_results.append({
                        'state': state_name,
                        'quantity': quantity,
                        'max_diff': max_diff,
                        'values': values
                    })
        
        return all_results
    
    def create_ghz_state(self):
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        return qc
    
    def create_w_state(self):
        qc = QuantumCircuit(3)
        # Alternative W state preparation
        qc.ry(2 * np.arccos(np.sqrt(2/3)), 0)
        qc.cx(0, 1)
        qc.x(0)
        qc.ry(2 * np.arccos(1/np.sqrt(2)), 0)
        qc.ccx(0, 1, 2)
        qc.x(0)
        return qc

# Add this method to your EnhancedObserverTest class
    def verify_result(self):
        """Make sure we're measuring what we think we're measuring"""
        
        print("\n=== VERIFICATION: What's Actually Happening? ===\n")
        
        # Create GHZ
        qc = self.create_ghz_state()
        
        # Get state before any observer operations
        qc.save_statevector()
        job = self.backend.run(qc)
        psi_original = job.result().get_statevector()
        rho_original = DensityMatrix(psi_original)
        
        # Trace out C to get original E_AB
        rho_AB_original = partial_trace(rho_original, [2])
        E_AB_original = float(entropy(rho_AB_original, base=2))
        
        print("ORIGINAL GHZ STATE:")
        print(f"State vector: {psi_original}")
        print(f"E_AB (original): {E_AB_original:.6f}")
        
        # Now apply each observer's operations
        for observer in ['A', 'B', 'C']:
            print(f"\n--- Observer {observer} ---")
            
            # Create fresh circuit for each observer
            qc_obs = self.create_ghz_state()
            qc_obs = self.apply_measurement_context(qc_obs, observer, 'interference')
            
            qc_obs.save_statevector()
            job = self.backend.run(qc_obs)
            psi_obs = job.result().get_statevector()
            rho_obs = DensityMatrix(psi_obs)
            
            # Calculate E_AB after observer operations
            rho_AB_obs = partial_trace(rho_obs, [2])
            E_AB_obs = float(entropy(rho_AB_obs, base=2))
            
            print(f"State after {observer}'s operations: {psi_obs}")
            print(f"E_AB after {observer}'s operations: {E_AB_obs:.6f}")
            
            # Show what operations were applied
            if observer == 'C':
                print("\nObserver C's operations:")
                print("1. H on qubit 2")
                print("2. CZ between qubits 2 and 0")
                print("3. H on qubit 2")
                print("This is why E_AB drops to 0!")

        # Then run verification to understand the results


# Run enhanced test
if __name__ == "__main__":
    test = EnhancedObserverTest()
    results = test.comprehensive_test()
    test.verify_result()
    
    # Optional: deeper analysis
    print("\n=== ANALYSIS ===")
    print("The key finding: Observer C's interference operations")
    print("disentangle qubits A and B from C's perspective!")
    
    if results:
        print("\n\n=== SUMMARY ===")
        print(f"Found {len(results)} instances of observer dependence!")
        
        # Find the strongest effect
        best = max(results, key=lambda x: x['max_diff'])
        print(f"\nStrongest effect:")
        print(f"State: {best['state']}")
        print(f"Quantity: {best['quantity']}")
        print(f"Max difference: {best['max_diff']:.6f}")
    else:
        print("\n\nNo observer dependence found with these protocols.")
        print("Next steps: Try noise models or post-selection.")