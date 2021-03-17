@isdefined(SchwarzOCPExamples) || include("SchwarzOCPExamples.jl")
using SimpleNLModels, SimpleSchwarz, SimpleADMM, Ipopt, .SchwarzOCPExamples

include("exponential_decay.jl")
include("demonstrating_convergence.jl")
include("effect_of_overlap.jl")
include("benchmark_quadrotor.jl") 
include("benchmark_quadrotor.jl") # needs to be run twice for exact timing
include("benchmark_thin_plate.jl") 
include("benchmark_thin_plate.jl") # needs to be run twice for exact timing

println("case study ran successfully")
