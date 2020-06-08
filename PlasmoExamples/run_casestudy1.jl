println("Activating Julia Environment")
using Pkg
Pkg.activate(@__DIR__)

println("\nRunning Gas Network Case Study\n")
include("case_studies/gas_network/partition_and_solve.jl")
