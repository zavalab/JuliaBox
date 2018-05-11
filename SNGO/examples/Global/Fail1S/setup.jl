push!(LOAD_PATH, ENV["SINGODIR"])
using PlasmoOld, Ipopt, SCIP
using JuMP
include("../../bb.jl")
NS = 1