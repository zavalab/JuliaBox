#
# Partition with KaHyPar
#
using KaHyPar

#NOTE: define n_parts and max_imbalance

hypergraph,hyper_map = gethypergraph(gas_network)

#Setup node and edge weights
n_vertices = length(vertices(hypergraph))
node_weights = [num_variables(node) for node in all_nodes(gas_network)]
edge_weights = [num_link_constraints(edge) for edge in all_edges(gas_network)]

#Use KaHyPar to partition the hypergraph
node_vector = KaHyPar.partition(hypergraph,n_parts,configuration = :edge_cut,
imbalance = max_imbalance, node_weights = node_weights, edge_weights = edge_weights)

#Create a model partition
partition = Partition(node_vector,ref_map)

#Setup subgraphs based on partition
make_subgraphs!(gas_network,partition)

#Combine the subgraphs into model-nodes
combined_graph , combine_map  = combine(gas_network,0)


#
# Solve with PIPS-NLP
#
using Distributed
using MPIClusterManagers

manager=MPIManager(np=2) # specify, number of mpi workers
addprocs(manager)        # start mpi workers and map them to julia workers

#Setup the worker environments
@everywhere using Plasmo
@everywhere using PipsSolver
@everywhere using MPI

#get the julia ids of the mpi workers
julia_workers = collect(values(manager.mpi2j))

#Use PipsSolver.jl to distribute the modelgraph among the workers
#Here, we create the variable `pipsgraph` on each worker
remote_references = PipsSolver.distribute(combined_graph,julia_workers,remote_name = :pipsgraph)

#Solve with PIPS-NLP on each mpi rank
@mpi_do manager begin
    PipsSolver.pipsnlp_solve(pipsgraph)
end
