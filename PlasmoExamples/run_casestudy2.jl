println("Activating Julia Environment")
using Pkg
Pkg.activate(@__DIR__)

println("\nRunning DC OPF Case Study\n")
include("case_studies/dcpowergrid/partition_and_solve.jl")
