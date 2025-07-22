# Quantum Algorithm Demo

This project implements quantum algorithms for solving molecular Hamiltonians and performing clustering analysis.

## Project Overview
- Uses Qiskit for constructing and solving molecular Hamiltonians
- Supports models such as Heisenberg, Spin Glass, and molecules (e.g., LiH)
- Integrates VITE (Variational Imaginary Time Evolution) and clustering analysis (KMeans, Hopkins, etc.)

## Directory Structure
```
.
├── molecular/           # Core algorithm modules
│   ├── heisenberg.py    # Heisenberg model
│   ├── spin_glass.py    # Spin glass model
│   ├── molecular.py     # Molecular Hamiltonian and related tools
│   ├── vqe.py           # Quantum circuits and VQE, VITE
│   ├── utils.py         # Utility functions
│   └── clustering.py    # Clustering and analysis
├── scripts/
│   └── run_experiment.py # Main experiment entry script
├── README.md            # Project documentation
├── requirements.txt     # Dependencies 
└── .gitignore           # Git ignore file (to be completed)
```

## Installation
Python 3.8+ is recommended. Install dependencies with:

```sh
pip install -r requirements.txt
```

Main dependencies include:
- qiskit
- qiskit-nature
- numpy
- scipy
- pandas
- matplotlib
- seaborn
- scikit-learn

## Quick Start
1. Clone this repository:
   ```sh
   git clone https://github.com/renmengzhen/molecular-demo.git
   cd molecular-demo
   ```
2. Install dependencies:
   ```sh
   pip install -r requirements.txt
   ```
3. Run the main experiment script:
   ```sh
   python scripts/run_experiment.py
   ```

## Contribution & Feedback
Contributions are welcome! Please open an issue or submit a pull request for suggestions and improvements.
