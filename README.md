# Loop Quantum Gravity and Entropy Study

This repository began as a holding area for a fun python program I thought my be relavent in QC around entanglement. It changed into a bit more.

The current direction aims at creating some proofs that can lead from a purely combinatoric and graph theoretic view of observer independence in LQG, to:

1.  derivation of known physical constants
2.  predictions of observable phenomena
3.  operator theory work
4.  a greater link between LQG and QI/QC
5.  insights into the "why" behind some of how nature acts
6.  a deeper understanding of time

I am working through some minor biblio issues, as my source quoting was not as rigorous as it should be. If anyone wants to volunteer to help and has some better handle on what they are doing - please reach out!

My current WIP:

*   bridge-monotocity.pdf/tex - describes a way of calculating relational entropy
*   entropy-spin-networks.pdf/tex - extends bridge-monotocity and shows an arrow of time and a chain leading to derivation of a known value
*   operator-theory.pdf/tex

---

Here's the below part of the readme that contained just the first experiment.

# Quantum Contextuality Demonstration

This project demonstrates quantum contextuality in three-party quantum systems using Qiskit. It shows how different measurement contexts (reference frames) reveal different entanglement structures in the same quantum state, based on the Plávala-Gühne theorem connecting contextuality and entanglement.

## Features

*   **Quantum Contextuality**: Demonstrates how measurement context affects observed quantum correlations
*   **Multiple Quantum States**: Tests with GHZ and W-states
*   **Reference Frame Transformations**: Uses unitary transformations representing different measurement contexts
*   **Comprehensive Analysis**: Measures entanglement entropy, concurrence, and provides mathematical verification
*   **Educational Focus**: Clear explanations of what the results mean and what they don't mean

## Setup and Installation

### Prerequisites

*   Python 3.8 or higher
*   pip (Python package installer)

### Installation

**Clone the repository:**

**Create and activate a virtual environment:**

On Windows:

On macOS/Linux:

**Install dependencies:**

## Usage

### Main Demonstration

Run the main contextuality demonstration:

```
python qi2.py
```

This will show:

*   How different measurement contexts (A, B, C) reveal different entanglement structures
*   Context A & B: E\_AB = 1.0, Concurrence = 0.0 (mixed state)
*   Context C: E\_AB = 0.0, Concurrence = 1.0 (pure Bell state)
*   W-state verification with perfect fidelity
*   Mathematical verification of the quantum transformations

### Analysis Scripts

For deeper technical analysis:

```
# Debug the actual quantum measurements
python debug_analysis.py

# Get corrected mathematical interpretation  
python final_check.py
```

## Understanding the Results

### Key Insight

The demonstration shows that **Context C** transforms the GHZ state `(|000⟩ + |111⟩)/√2` into `(|000⟩ + |011⟩)/√2`, which creates a Bell state between qubits A and B when we trace out qubit C.

### What the Metrics Mean

*   **E\_AB = 0**: The A-B subsystem is in a pure state (not mixed)
*   **Concurrence = 1**: The A-B subsystem is maximally entangled (Bell state)
*   **Different contexts reveal different entanglement structures** - this is quantum contextuality!

### Important Notes

*   This is **NOT** about measurement collapse or wavefunction reduction
*   These are **unitary transformations** representing different reference frames
*   The effect demonstrates the **Plávala-Gühne theorem** connecting contextuality and entanglement

## Files Overview

*   `**qi2.py**`: Main quantum contextuality demonstration script
*   `**debug_analysis.py**`: Technical analysis showing actual quantum state transformations
*   `**final_check.py**`: Corrected mathematical interpretation of the results
*   `**requirements.txt**`: Python dependencies
*   `**setup_and_run.bat/.sh**`: Setup scripts for different platforms

## Scientific Background

Based on: Plávala & Gühne, "Contextuality as a Precondition for Quantum Entanglement" (arXiv:2209.09942)

This demonstrates that quantum entanglement is fundamentally contextual - the observed correlations depend on the measurement framework used to observe them. Different unitary reference frame transformations can reveal completely different entanglement structures within the same quantum state.

## Dependencies

*   `qiskit>=1.0.0`: Main quantum computing framework
*   `qiskit-aer>=0.13.0`: High-performance quantum circuit simulator
*   `numpy>=1.24.0`: Numerical computations
*   `matplotlib>=3.7.0`: Plotting and visualization (optional)

## Contributing

Feel free to fork this repository and submit pull requests for improvements or additional quantum contextuality demonstrations.

## License

This project is open source. Please include appropriate attribution if you use this code in your research or projects. If using for academic work, please cite the relevant papers on quantum contextuality.

```
pip install -r requirements.txt
```

```
python -m venv venv
source venv/bin/activate
```

```
python -m venv venv
venv\Scripts\activate
```

```
git clone <your-repo-url>
cd qplay
```