# Enzyme Kinetics and Thermo Model

## Project Description

This project implements a comparison of enzyme kinetics models, primarily including calculations for three-step and four-step models. It is implemented in the Julia language and utilizes PyCall to invoke Python modules for comparison. The project employs symbolic computation and data visualization to help researchers understand the performance of different models in enzyme kinetics.

## File Structure

- `model_cal.jl`: Defines the equations for the enzyme kinetics model and calculates the net reaction rate.
- `model_comparison_pycall.jl`: Uses PyCall to call Python modules and compare the results of different models.
- `thermo_kinetics_modified.jl`: Defines thermodynamic parameters and calculation functions, supporting the calculations for both three-step and four-step models.

## Environment Setup

### 1. Install Julia

Download and install the latest stable version of Julia from the official website: [https://julialang.org/downloads/](https://julialang.org/downloads/)

Follow the installation instructions specific to your operating system (Windows, macOS, or Linux).

### 2. Configure Python Environment for PyCall

This project uses the `PyCall` package to interact with Python. It's recommended to use `Conda.jl` to manage the Python environment used by `PyCall`.

Open a Julia REPL (type `julia` in your terminal) and run the following commands:

```julia
using Pkg
Pkg.add("Conda")
using Conda
Conda.add("python") # Or specify a specific Python version, e.g., Conda.add("python=3.9")
# Add any required Python packages here if needed, e.g., Conda.add("numpy")
```

This will create an isolated Conda environment for Julia to use. `PyCall` will automatically use this environment.

## Installation

Ensure that you have Julia and the necessary dependencies installed. You can install the required Julia packages using the following commands:

```julia
using Pkg
Pkg.add("Symbolics")
Pkg.add("Plots")
Pkg.add("Measures")
Pkg.add("DataFrames")
Pkg.add("PyCall")
```

## Usage Instructions

1. **Run Model Calculations**:
   The equations defined in `model_cal.jl` will calculate the net reaction rate and output the results.
2. **Compare Models**:
   In `model_comparison_pycall.jl`, you can set different parameters and run the comparison function, with results output as a DataFrame.
3. **Generate Reports**:
   Running the `generate_markdown_report` function will generate a Markdown report containing the results of the model comparison.

## Example

In `model_comparison_pycall.jl`, you can find an example usage of the `main` function, where parameters are set and model comparisons are run.

```julia
# Set parameters
ΔGT = -40.0  # kJ/mol
ΔG1 = -15.0  # kJ/mol
S_values = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]  # mM

# Compare models
results_df = compare_models(ΔGT, ΔG1, S_values)
```

## Results

After running the models, you will obtain comparison results for Km and reaction rates at different substrate concentrations, along with corresponding charts and reports.
