# SchwarzOCP Examples

To reproduce the results, run the following script
```julia
julia --project=. -p 20 # make sure that your hardware has enough cores

] add https://github.com/sshin23/SimpleNLModels.jl#2a0023d802ba8211f60a8c920a1a5bbd18eb04ed
] add https://github.com/sshin23/SimpleSchwarz.jl#4b6f7c04817b98ede8362f621e9b02924a99c376
] add https://github.com/sshin23/SimpleADMM.jl#408cb08d86fb76c8a08e07d0df9ac7d3241e1268
] instantiate

include("run_case_study.jl")
```
