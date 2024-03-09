using Revise
using PlasmoData, Random, PlasmoDataPlots, LinearAlgebra

Random.seed!(5)
random_matrix = rand(20, 20)

symmetric_matrix = (random_matrix .+ random_matrix') / 2
symmetric_matrix[diagind(symmetric_matrix)] .= 1

symmetric_matrix_graph = symmetric_matrix_to_graph(symmetric_matrix)
set_circle_node_positions!(symmetric_matrix_graph)

plot_graph(
    symmetric_matrix_graph,
    nodesize = 10,
    nodecolor = "gray",
    xdim = 500,
    ydim = 500,
    linewidth = 3,
    linecolor = :binary,
    rev = true,
    line_z = get_edge_data(symmetric_matrix_graph, "weight"),
    framestyle = :none,
    save_fig = false,
    fig_name = (@__DIR__)*"/images/symmetric_graph.pdf"
)

thresh = 0:.001:1
EC_vec = PlasmoData.run_EC_on_edges(symmetric_matrix_graph, thresh)

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


filtered_graph = filter_edges(symmetric_matrix_graph, .6, "weight", fn = Base.isless)

node_cols = get_cols_from_cgrad(:binary, get_edge_data(symmetric_matrix_graph)[:, 1][:], rev = true)

bool_vec = get_edge_data(symmetric_matrix_graph)[:, 1] .< .6

plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 3,
    linecolor = node_cols[bool_vec],
    nodecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter1.pdf"
)



filtered_graph = filter_edges(symmetric_matrix_graph, .3, "weight", fn = Base.isless)

bool_vec = get_edge_data(symmetric_matrix_graph)[:, 1] .< .3

plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 3,
    linecolor = node_cols[bool_vec],
    nodecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter1.pdf"
)


filtered_graph = filter_edges(symmetric_matrix_graph, 0, "weight", fn = Base.isless)

bool_vec = get_edge_data(symmetric_matrix_graph)[:, 1] .< 0

plot_graph(
    filtered_graph,
    nodesize = 10,
    xdim = 500,
    ydim = 500,
    linewidth = 3,
    linecolor = node_cols[bool_vec],
    nodecolor = "gray",
    save_fig = false,
    framestyle = :none,
    fig_name = (@__DIR__)*"/images/mat_graph_node_filter1.pdf"
)
