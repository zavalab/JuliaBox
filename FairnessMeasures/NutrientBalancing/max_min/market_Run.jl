## Market coordination case study
## Market analysis of the rock river area

println("----------Reading Model-----------")
tic()
include("market_model.jl"); # Reading model file
toc()
println("----------Begin solving-----------")

tic()
solve(m) 
toc()

# print main results
println("\n----------Print results----------")
include("market_print.jl")
