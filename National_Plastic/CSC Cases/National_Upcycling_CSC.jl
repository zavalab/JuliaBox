using CoordinatedSupplyChains
using DataFrames
using DelimitedFiles
using Gurobi

#########################
# Run Model on Path
#########################

# Employ CoordinatedSupplyChains.jl package

case_path = pwd()   # Choose case study folder

global T, N, P, Q, A, D, G, V, M, L, Subsets, Pars, MOD, ModelStats, SOL, POST = RunCSC(case_path, optimizer=Gurobi.Optimizer, UseArcLengths=true, Output=true)
