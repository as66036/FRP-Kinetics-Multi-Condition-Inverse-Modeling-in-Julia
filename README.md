# 🧪 CFRP Oxidative Decomposition — Multi-Condition Inverse Modeling

This project was developed for **Homework 2** in *Introduction to Scientific Computing* (Instructor: **Dr. Hayri Sezer**).

It presents a complete workflow for modeling and parameter estimation of **Carbon Fiber Reinforced Polymer (CFRP)** thermal degradation using **Julia**.

---

## 📌 Overview

CFRP decomposition under high temperature involves multiple coupled reactions:

| Reaction | Description | Rate Law |
|----------|------------|----------|
| **1. Matrix Pyrolysis** | M → ν_char·C + (1-ν_char)·G₁ | r₁ = A₁·exp(-E₁/RT)·[M]^n₁ |
| **2. Char Oxidation** | C + O₂ → G₂ | r₂ = A₂·exp(-E₂/RT)·[C]^n₂·PO₂^m₂ |
| **3. Fiber Oxidation** | F + O₂ → G₂ | r₃ = A₃·exp(-E₃/RT)·[F]^n₃·PO₂^m₃ |

The objective is to:
- Simulate the forward reaction model  
- Generate synthetic TGA data  
- Solve the inverse problem to recover kinetic parameters  
- Validate the model using multi-condition experiments  

---

## 📂 Project Structure

```
cfrp_kinetics/
├── Project.toml
├── README.md
├── run.jl
├── src/
│   ├── model.jl
│   ├── data_generation.jl
│   ├── inverse_problem.jl
│   └── visualization.jl
└── results/
```

---

## ⚙️ How to Run

```julia
cd("path/to/cfrp_kinetics")

using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("run.jl")
```

---

## 🧪 Experimental Setup

| Exp | β (K/min) | PO₂ | Purpose |
|-----|-----------|------|---------|
| 1 | 2.5 | 0.21 | Baseline kinetics |
| 2 | 5.0 | 0.21 | Baseline kinetics |
| 3 | 10.0 | 0.21 | Activation energy |
| 4 | 5.0 | 0.05 | Oxygen order |

---

## 📊 Outputs

- results/measurements.csv  
- results/parameter_comparison.csv  
- results/all_experiments.png  
- results/validation_plot.png  

---

## 📦 Dependencies

- Catalyst.jl  
- PEtab.jl  
- OrdinaryDiffEq.jl  
- Plots.jl  
- DataFrames.jl  
- CSV.jl  

---

## 📝 Author

Abdullah Al Sayeed  
Georgia Southern University
