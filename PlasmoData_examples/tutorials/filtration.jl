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
    fig_name = (@__DIR__)*"/images/mat_graph_diags_filter.pdf"
)

filtered_graph = filter_nodes(matrix_graph, .7, "weight", fn = Base.isless)

node_cols = get_cols_from_cgrad(:binary, get_node_data(matrix_graph)[:, 1][:])

bool_vec = get_node_data(matrix_graph)[:, 1] .< .7

plot_graph(
    filtered_graph,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = node_cols[bool_vec],
    linecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter1.pdf"
)

function comp(a, b)
    return ((a <= .2) || a >= .8)
end

filtered_graph = filter_nodes(matrix_graph, .5, "weight", fn = comp)

node_cols = get_cols_from_cgrad(:binary, get_node_data(matrix_graph)[:, 1][:])

bool_vec = ((get_node_data(matrix_graph)[:, 1] .< .2) .|| (get_node_data(matrix_graph)[:, 1] .> .8))

plot_graph(
    filtered_graph,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = node_cols[bool_vec],
    linecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter2.pdf"
)

Random.seed!(15)
add_edge_dataset!(matrix_graph, rand(506), "weight")

plot_graph(
    matrix_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linecolor = :binary,
    linewidth = 5,
    nodecolor = :gray,
    line_z = get_edge_data(matrix_graph)[:],
    rev = false,
    save_fig = true,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_edge_color.pdf"
)



filtered_graph = filter_edges(matrix_graph, .8, "weight", fn = Base.isgreater)

edge_cols = get_cols_from_cgrad(:binary, get_edge_data(matrix_graph)[:, 1][:])

bool_vec = ((get_edge_data(matrix_graph)[:, 1] .> .8))


plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :gray,
    linecolor = edge_cols[bool_vec],
    save_fig = true,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_edge_filter1.pdf"
)


filtered_graph = filter_edges(matrix_graph, .5, "weight", fn = Base.isless)

edge_cols = get_cols_from_cgrad(:binary, get_edge_data(matrix_graph)[:, 1][:])

bool_vec = ((get_edge_data(matrix_graph)[:, 1] .< .5))


plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :gray,
    linecolor = edge_cols[bool_vec],
    save_fig = true,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_edge_filter2.pdf"
)



Random.seed!(15)
random_matrix = rand(12, 12)

matrix_graph = matrix_to_graph(random_matrix; diagonal = true)


filter_nodes_graph1 = filter_nodes(matrix_graph, .7, "weight", fn = Base.isless)

function extreme_vals(a, b)
    return ((a <= .2) || a >= .8)
end

filter_edges_graph2 = filter_nodes(matrix_graph, .5, "weight", fn = extreme_vals)

Random.seed!(15)
add_edge_dataset!(matrix_graph, rand(506), "weight")

filter_edges_graph1 = filter_edges(matrix_graph, .8, "weight", fn = Base.isgreater)

filter_edges_graph2 = filter_edges(matrix_graph, .5, "weight", fn = Base.isless)
