# Development Guide

## Setting up the Development Environment

### Manual Setup

**Clone and navigate to the repository:**

**Create virtual environment:**

**Install dependencies:**

### Automated Setup

**Windows:**

```
setup_and_run.bat
```

**macOS/Linux:**

```
./setup_and_run.sh
```

## Project Structure

```
qplay/
├── qi2.py                 # Main quantum observer test script
├── requirements.txt       # Python dependencies
├── README.md             # Project documentation
├── .gitignore            # Git ignore rules
├── setup_and_run.bat     # Windows setup script
├── setup_and_run.sh      # Unix setup script
├── DEVELOPMENT.md        # This file
└── venv/                 # Virtual environment (not in git)
```

## Running Tests

```
# Activate virtual environment first
source venv/bin/activate  # Unix
# or
venv\Scripts\activate     # Windows

# Run the main script
python qi2.py
```

## Key Components

### EnhancedObserverTest Class

*   `**create_asymmetric_state()**`: Creates quantum states with intentional asymmetry
*   `**apply_measurement_context()**`: Applies observer-specific operations
*   `**measure_multiple_quantities()**`: Measures various quantum information metrics
*   `**comprehensive_test()**`: Runs full test suite across multiple states

### Observer Protocols

1.  **Interference**: Creates interference patterns between observer and system
2.  **Weak Measurement**: Simulates weak measurement interactions
3.  **Basis Mixing**: Mixes measurement bases for different observers

### Quantum States Tested

*   **GHZ State**: Maximally entangled three-qubit state
*   **W State**: Symmetric entangled state
*   **Asymmetric State**: Custom state with built-in asymmetries

## Adding New Tests

To add new observer dependence tests:

1.  Create new state preparation methods
2.  Add new observer protocols to `apply_measurement_context()`
3.  Implement new quantum measures in `measure_multiple_quantities()`
4.  Update the test loop in `comprehensive_test()`

## Dependencies Management

When adding new dependencies:

1.  Install in your activated virtual environment
2.  Update `requirements.txt`:
3.  Test that others can install from the updated requirements

## Deactivating Environment

```
deactivate
```

## Publishing Changes

1.  Ensure virtual environment is in `.gitignore`
2.  Test setup scripts on clean checkout
3.  Update README.md with any new features
4.  Commit only source files, not virtual environment

```
pip freeze > requirements.txt
```

```
pip install -r requirements.txt
```

```
# Windows
python -m venv venv
venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

```
git clone <your-repo-url>
cd qplay
```