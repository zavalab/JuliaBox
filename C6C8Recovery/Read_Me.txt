This folder provides supporting codes for the paper "A Supply Chain Framework for the Analysis of the Recovery of Biogas and Fatty Acids from Organic Waste".

* All original data, data treatment and case study results are summarized in the file "Data_results.xlsx".
* In each folder, "c8_mix.jl" is the model script, "c8_Run.jl" is the running script, "c8_Flow.jl" is used for printing flow(logistic) results, "Technologies_Sited.jl" is used for printing technology location results.
* We recommend use Julia 0.5.2 and Gurobi 7.0.2 to run all code files.
* The model takes "demand_matrix.csv", "node_matrix.csv", "supply_matrix.csv" and "temperatures.csv" as input data. Other data are set in the script. Other .csv files are output files.

All case studies have corresponding files in this folder:

1. Four base cases
* do nothing case: see folder "Do_nothing"
* biogas only case: see folder "MIX_biogas"
* C6/C8 only case: see folder "MIX_C6C8"
* hybrid case: see folder "MIX_Combo"

2. Mixing effect cases
* see folder "Mixing_effect"
* We decompose the whole network into three sub-networks to study the case when mixing is not allowed.

3. Sensitivity analysis of biogas price
* see folder "Biogas_price"
* numbers in each subfolder name represent the biogas price

4. Sensitivity analysis of social cost of carbon
* see folder "SCC_sensitivity"
* numbers in each subfolder name represent the value of social cost of carbon

5. Sensitivity analysis of yield factors of C6 and C8
* see folder "yield_sensitivity"
* numbers in each subfolder name represent the COD conversion of C6 and C8 repectively
