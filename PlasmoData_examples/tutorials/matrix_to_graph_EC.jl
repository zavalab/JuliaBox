using Revise
using PlasmoData, Random, DataGraphPlots, LinearAlgebra

Random.seed!(5)
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

thresh = 0:.0001:1
EC_vec = PlasmoData.run_EC_on_nodes(matrix_graph, thresh)

using Plots
plot(thresh, EC_vec, legend = :none, grid = false, color = :black, linewidth = 3,
xtickfontsize = 12, ytickfontsize = 12, xlabelfontsize = 16, ylabelfontsize = 16)
xlabel!("Filtration Threshold")
ylabel!("Euler Characteristic")


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


filtered_graph = filter_nodes(matrix_graph, .8, "weight", fn = Base.isless)

node_cols = get_cols_from_cgrad(:binary, get_node_data(matrix_graph)[:, 1][:], rev = true)

bool_vec = get_node_data(matrix_graph)[:, 1] .< .8

plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 3,
    nodecolor = node_cols[bool_vec],
    linecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter1.pdf"
)



filtered_graph = filter_nodes(matrix_graph, .4, "weight", fn = Base.isless)

node_cols = get_cols_from_cgrad(:binary, get_node_data(matrix_graph)[:, 1][:], rev = true)

bool_vec = get_node_data(matrix_graph)[:, 1] .< .4

plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 3,
    nodecolor = node_cols[bool_vec],
    linecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter1.pdf"
)
