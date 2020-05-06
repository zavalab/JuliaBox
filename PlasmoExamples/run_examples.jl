#This script runs all of the examples in the examples folder and plots the results in examples/figures

println("Activating Julia Environment")
using Pkg
Pkg.activate(@__DIR__)

#Model Building Examples
#Example 1 - simple modelgraph
println("\nRunning Example 1\n")
include("examples/example1.jl")
savefig(plt_graph1,"examples/figures/example1_layout.pdf")
savefig(plt_matrix1,"examples/figures/example1_matrix.pdf")

#Example 2 - distributed modelgraph
println("\nRunning Example 2\n")
include("examples/example2.jl")
savefig(plt_graph2,"examples/figures/example2_layout.pdf")
savefig(plt_matrix2,"examples/figures/example2_matrix.pdf")

#Example 3 - hierarchical modelgraph
println("\nRunning Example 3\n")
include("examples/example3.jl")
savefig(plt_graph3,"examples/figures/example3_layout.pdf")
savefig(plt_matrix3,"examples/figures/example3_matrix.pdf")

#Partitioning Examples - Optimal Control Problem
#Example 4
println("\nRunning Example 4\n")
include("examples/example4.jl")
# ocp structure
savefig(plt_graph4,"examples/figures/example4_layout1.pdf")
savefig(plt_matrix4,"examples/figures/example4_matrix1.pdf")

# partition structure
savefig(plt_graph5,"examples/figures/example4_layout2.pdf")
savefig(plt_matrix5,"examples/figures/example4_matrix2.pdf")

# combined structure
savefig(plt_graph6,"examples/figures/example4_layout3.pdf")
savefig(plt_matrix6,"examples/figures/example4_matrix3.pdf")

#Example 5 - overlap structure
println("\nRunning Example 5\n")
include("examples/example5.jl")
savefig(plt_graph7,"examples/figures/example5_layout.pdf")
savefig(plt_matrix7,"examples/figures/example5_matrix.pdf")
