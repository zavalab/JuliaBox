using PlasmoData, Random, PlasmoDataPlots

Random.seed!(15)
random_matrix = rand(12, 12)
#using DelimitedFiles
#writedlm((@__DIR__)*"/rand_mat.csv", random_matrix, ',')

matrix_graph_diags = matrix_to_graph(random_matrix; diagonal = true)
set_matrix_node_positions!(matrix_graph_diags, random_matrix)
plot_graph(
    matrix_graph_diags,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :grays,
    node_z = get_node_data(matrix_graph_diags, "weight"),
    linecolor = "gray",
    rev = true,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_diags.pdf"
)

matrix_graph_no_diags = matrix_to_graph(random_matrix; diagonal = false)
set_matrix_node_positions!(matrix_graph_no_diags, random_matrix)
plot_graph(
    matrix_graph_no_diags,
    nodesize = 12,
    xdim = 500,
    ydim = 500,
    linewidth = 5,
    nodecolor = :grays,
    node_z = get_node_data(matrix_graph_no_diags, "weight"),
    linecolor = "gray",
    rev = true,
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_no_diags.pdf"
)
