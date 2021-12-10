# SchwarzOCP Examples

To reproduce the results, run the following script
```julia
julia --project=. -t 20 # make sure that your hardware has enough cores

] add https://github.com/sshin23/SimpleNLModels.jl#2a0023d802ba8211f60a8c920a1a5bbd18eb04ed
] add https://github.com/sshin23/SimpleSchwarz.jl#6ef555336f30ee2e6cedaf874e89f60a136cb2b2
] add https://github.com/sshin23/SimpleADMM.jl#408cb08d86fb76c8a08e07d0df9ac7d3241e1268
] instantiate

include("run_case_study.jl")
```
