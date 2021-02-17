Here is the source-code for the dynamic biological community parameter estimation 
case study where we compare the utility of a formulation that uses the discrete 
experimental data directly and our proposed formulation that uses continuous fitted 
empirical functions of the experimental data. The file structure is as follows:
- `synthetic_data_generator.jl`: Simulates the system and generates the experimental data files.
- `sigmoidal_function_fitter.jl`: Fits a sigmoidal function to each experiment and stores fitted parameters via a data file.
- `continuous_estimation_problem.jl`: Implements the proposed continuous formulation and saves the results.
- `discrete_estimation_problem.jl`: Implements the traditional formulation and saves the results.
- `result_analyzer.jl`: Takes the above results and analyzes the accuracies.

To configure the required packages, we recommend creating a Julia environment 
using the included `Project.toml` file. Thus, we can configure the environment and 
run the case study via:
```julia
julia> cd("[INSERT_PATH_TO_FILES]")

julia> ]

(@v1.5) pkg> activate .

julia> include("[DESIRED_SCRIPT_FILE]")
```
