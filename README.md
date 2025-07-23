# Quantum Observer Dependence Test

This project demonstrates enhanced quantum observer dependence tests using Qiskit. It explores how different measurement contexts and observer operations can affect quantum entanglement measures.

## Features

*   **Enhanced Observer Tests**: Sophisticated tests showing observer dependence in quantum measurements
*   **Multiple Quantum States**: Tests with GHZ, W-state, and custom asymmetric quantum states
*   **Comprehensive Analysis**: Measures entanglement entropy, concurrence, three-tangle, and other quantum information metrics
*   **Observer Protocols**: Implements interference, weak measurement, and basis mixing protocols

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

Run the main script to see the enhanced observer dependence tests:

```
python qi2.py
```

The script will:

*   Test multiple quantum states (GHZ, Asymmetric, W-state)
*   Apply different observer protocols
*   Measure various quantum information quantities
*   Report any observer-dependent effects found
*   Provide verification and analysis of the results

## Understanding the Output

The program tests for observer dependence by:

1.  Creating entangled quantum states
2.  Applying observer-specific operations
3.  Measuring entanglement quantities from each observer's perspective
4.  Comparing results to identify observer-dependent effects

Key metrics measured:

*   **E\_AB**: Entanglement entropy between qubits A and B
*   **C\_AB**: Concurrence (a measure of entanglement)
*   **S\_A\_given\_C**: Conditional entropy
*   **three\_tangle**: Three-party entanglement measure
*   **linear\_entropy**: Purity measure

## Dependencies

*   `qiskit>=1.0.0`: Main quantum computing framework
*   `qiskit-aer>=0.13.0`: High-performance quantum circuit simulator
*   `numpy>=1.24.0`: Numerical computations
*   `matplotlib>=3.7.0`: Plotting and visualization
*   `scipy>=1.10.0`: Additional scientific computing functions

## Deactivating the Virtual Environment

When you're done working with the project, deactivate the virtual environment:

```
deactivate
```

## Contributing

Feel free to fork this repository and submit pull requests for improvements or additional quantum observer tests.

## License

This project is open source. Please include appropriate attribution if you use this code in your research or projects.

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