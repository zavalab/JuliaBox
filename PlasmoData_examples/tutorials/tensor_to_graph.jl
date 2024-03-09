using PlasmoData, Random, PlasmoDataPlots

Random.seed!(15)
random_tensor = rand(4, 5, 6)

tensor_graph_2d = matrix_to_graph(random_tensor; diagonal = true)
set_matrix_node_positions!(tensor_graph_2d, random_tensor[:, :, 1])
plot_graph(
    tensor_graph_2d,
    nodesize = 12,
    xdim = 400,
    ydim = 400,
    linewidth = 5,
    nodecolor = "black",
    linecolor = "gray",
    rev = true,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/tensor_mat_graph_diags.pdf"
)

tensor_graph_3d = tensor_to_graph(random_tensor)
plot_graph(
    tensor_graph_3d,
    nodesize = 8,
    xdim = 500,
    ydim = 500,
    linewidth = 4,
    nodecolor = :grays,
    node_z = get_node_data(tensor_graph_3d, "weight"),
    linecolor = "gray",
    rev = true,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/tensor_graph.pdf"
)

#=
graph_diags = matrix_to_graph(random_matrix; diagonal = true)
set_matrix_node_positions!(graph_diags, random_matrix)
plot_graph(
    graph_diags,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :grays,
    node_z = get_node_data(graph, "weight"),
    linecolor = "gray",
    rev = true,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_diags.pdf"
)

graph_no_diags = matrix_to_graph(random_matrix; diagonal = false)
set_matrix_node_positions!(graph_no_diags, random_matrix)
plot_graph(
    graph_no_diags,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :grays,
    node_z = get_node_data(graph, "weight"),
    linecolor = "gray",
    rev = true,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_no_diags.pdf"
)
=#


#using DelimitedFiles
#writedlm((@__DIR__)*"/rand_tensor.csv", random_tensor, ',')
