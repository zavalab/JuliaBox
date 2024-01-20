using Revise
using PlasmoData, Random, DataGraphPlots, Plots

function get_cols_from_cgrad(cgrad_name, vec, max_val = maximum(vec); rev = false)
    node_cgrad = cgrad(cgrad_name, rev = rev)
    len_cgrad = length(node_cgrad)

    min_val = minimum(vec)
    max_val = max_val
    span = max_val - min_val

    col_list = []
    for i in 1:length(vec)
        index = Int(floor(((vec[i] - min_val) / span) * (len_cgrad - 1) + 1))

        node_col = node_cgrad[index]
        push!(col_list, node_col)
    end

    return col_list
end

Random.seed!(15)
random_matrix = rand(12, 12)

matrix_graph = matrix_to_graph(random_matrix; diagonal = true)

set_matrix_node_positions!(matrix_graph, random_matrix)

nodes_for_aggregation = [(3, 7), (3, 8), (3, 9), (4, 7), (4, 8)]

plot_graph(
    matrix_graph,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :binary,
    node_z = get_node_data(matrix_graph, "weight"),
    linecolor = "gray",
    rev = false,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_diags_agg.pdf"
)

aggregated_graph = aggregate(matrix_graph, nodes_for_aggregation, "agg_node")

plot_graph(
    aggregated_graph,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :binary,
    node_z = get_node_data(matrix_graph, "weight"),
    linecolor = "gray",
    rev = false,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_diags_post_agg.pdf"
)
