push!(LOAD_PATH, ENV["SINGODIR"])
using PlasmoOld, Ipopt, SCIP
using JuMP
NS = 1000