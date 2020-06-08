using KaHyPar
using SchwarzSolver

include((@__DIR__)*"/powergrid.jl")

overlap = 10
n_parts = 4
max_imbalance = 0.1

hypergraph,hyper_map = gethypergraph(powergrid) #create hypergraph object based on graph
n_vertices = length(vertices(hypergraph))
node_sizes = [num_variables(node) for node in all_nodes(powergrid)]
edge_weights = [num_linkconstraints(edge) for edge in all_edges(powergrid)]

node_vector = KaHyPar.partition(hypergraph,n_parts,configuration = :edge_cut,
imbalance = max_imbalance, node_sizes = node_sizes, edge_weights = edge_weights)

#Create a model partition
# partition = Partition(hypergraph,node_vector,hyper_map)
partition = Partition(powergrid,node_vector,hyper_map)

#Setup subgraphs based on partition
make_subgraphs!(powergrid,partition)


subgraphs = getsubgraphs(powergrid)
expanded_subs = expand.(Ref(powergrid),subgraphs,Ref(overlap))


schwarz_solve(powergrid,expanded_subs;sub_optimizer = optimizer_with_attributes(Ipopt.Optimizer,"tol" => 1e-12,"print_level" => 0),max_iterations = 100,tolerance = 1e-10,
primal_links = primal_links,dual_links = dual_links)


#TODO
# Plot generator control policy
# Plot schwarz residuals
