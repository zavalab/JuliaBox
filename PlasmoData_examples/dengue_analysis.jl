using Revise
using PlasmoData, Graphs, DelimitedFiles, Statistics
using DataGraphPlots

# The data and methods in this example are primarily from the paper of
# Barros de Souze at al., 2022. https://doi.org/10.1088/1742-5468/aca0e5
# "The Euler characteristic as a topological marker for outbreaks in
# vector borne disease"
@time begin
full_data = readdlm((@__DIR__)*"/Recife_data/Dengue_Recife_new_cases_jan_2014_to_dec_2021.csv", ',', skipstart = 1)[:, 2:94]

# Set the matrix type to be Int
data = Matrix{Int}(undef, size(full_data))
data .= full_data

# Define a function to get the smallest filtered graph that only contains one connected component
function find_smallest_filter_graph(graph, k = 5)
    # Define variables as local
    local f_val, EC_val, max_cl, len_max_cl, k_clique, num_k_comms

    iter_values = sort(get_edge_data(graph)[:], rev = true)

    # iterate through filtration levels

    for i in iter_values
        filtered_graph = filter_edges(graph, i; fn = Base.isgreater)

        if length(connected_components(filtered_graph)) == 1
            f_val = i
            EC_val = get_EC(filtered_graph) # Euler characteristic
            max_cl = maximal_cliques(filtered_graph)
            len_max_cl  = length(max_cl) # Number of maximal cliques
            k_clique    = length(max_cl[1]) # Length of maximal cliques (k value)
            num_k_comms = length(clique_percolation(filtered_graph, k = k)) # Number of communities for k-clique percolation
            break
        end
    end

    return f_val, EC_val, len_max_cl, k_clique, num_k_comms
end

# Create array for storing solutions
sols = zeros((2916, 5))

# Run through all 7 day windows in the data
for i in 1:2916
    # Form a correlation matrix based on 7 days of data
    cor_mat = cor(data[i:(i + 6), :], dims = 1)

    # Remove the matrix entries that are NaNs
    bit_vec = (!).(isnan.(cor_mat[:, 1]))
    sym_mat = cor_mat[bit_vec, bit_vec]

    # If there are not more than 2 nodes in the graph, skip this iteration
    node_count = size(sym_mat, 1)

    if node_count <= 2
        continue
    end

    # Build the edge-weighted graph from the correlation matrix
    sym_graph = symmetric_matrix_to_graph(sym_mat)

    # Find the smallest filtration level possible and get the TDA metrics
    f_val, EC_val, len_max_cl, k_clique, num_k_comms = find_smallest_filter_graph(sym_graph, 25)

    sols[i, 1] = f_val
    sols[i, 2] = EC_val
    sols[i, 3] = len_max_cl
    sols[i, 4] = k_clique
    sols[i, 5] = num_k_comms


    if i%5 == 0
        println("DONE WITH $i")
    end
end
end
using DelimitedFiles

#writedlm((@__DIR__)*"/sols_k25_new.csv", sols, ',')
