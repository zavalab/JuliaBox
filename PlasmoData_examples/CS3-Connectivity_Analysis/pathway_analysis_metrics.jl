using Revise, Plots
using PlasmoData, PlasmoDataPlots, Graphs
using DelimitedFiles, Makie

raw_data  = readdlm((@__DIR__)*"/pathway_example/rawmaterial_data_array.csv", ',')
prod_data = readdlm((@__DIR__)*"/pathway_example/product_data_array.csv", ',')
int_data  = readdlm((@__DIR__)*"/pathway_example/intermediates_data_array.csv", ',')
tech_data = readdlm((@__DIR__)*"/pathway_example/technology_data_array.csv", ',')
edge_data = readdlm((@__DIR__)*"/pathway_example/edge_data.csv", ',')

function get_tech_graph()

    dg = DataDiGraph()

    for i in 1:size(raw_data, 1)
        add_node!(dg, raw_data[i, 1])
        add_node_data!(dg, raw_data[i, 1], 1, "Raw Material")
        add_node_data!(dg, raw_data[i, 1], raw_data[i, 2], "Cost")
        add_node_data!(dg, raw_data[i, 1], raw_data[i, 3], "CO2 Cost")
        add_node_data!(dg, raw_data[i, 1], raw_data[i, 4], "Max Supply")
    end

    for i in 1:size(prod_data, 1)
        add_node!(dg, prod_data[i, 1])
        add_node_data!(dg, prod_data[i, 1], 1, "Product")
        add_node_data!(dg, prod_data[i, 1], prod_data[i, 2], "Demand Limit")
    end

    for i in 1:size(int_data, 1)
        add_node!(dg, int_data[i, 1])
        add_node_data!(dg, int_data[i, 1], 1, "Intermediate")
    end

    for i in 1:size(tech_data, 1)
        add_node!(dg, tech_data[i, 1])
        add_node_data!(dg, tech_data[i, 1], 1, "Technology")
        add_node_data!(dg, tech_data[i, 1], tech_data[i, 2], "Cost")
        add_node_data!(dg, tech_data[i, 1], tech_data[i, 3], "CO2 Cost")
    end

    for i in 1:size(edge_data, 1)
        edge = (edge_data[i, 1], edge_data[i, 2])
        add_edge!(dg, edge)
        PlasmoData.add_edge_data!(dg, edge, edge_data[i, 3], "Optimal Flow")
    end

    return dg
end

dg = get_tech_graph()

function compute_metrics!(dg)
    for r in raw_data[:, 1]
        count = 0
        num_downstream = length(downstream_nodes(dg, r))
        add_node_data!(dg, r, num_downstream, "Number Downstream")
        for p in prod_data[:, 1]
            if PlasmoData.has_path(dg, r, p)
                count += 1
            end
        end
        add_node_data!(dg, r, count, "Connected Products")
    end


    for p in prod_data[:, 1]
        count = 0
        num_upstream = length(upstream_nodes(dg, p))
        add_node_data!(dg, p, num_upstream, "Number Upstream")
        for r in raw_data[:, 1]
            if PlasmoData.has_path(dg, r, p)
                count += 1
            end
        end
        add_node_data!(dg, p, count, "Connected Raw")
    end
end

compute_metrics!(dg)


using SankeyPlots

src = [i for (i, j) in dg.edges]
dst = [j for (i, j) in dg.edges]
weights = get_edge_data(dg, ["Optimal Flow"])

plt_s = sankey(src, dst, weights; size = (1400, 800), node_labels = String.(dg.nodes), edge_color = :black, label_size = 12)
#savefig(plt_s, (@__DIR__)*"/Sankey_plot.pdf")

node_cols = Vector{Any}(undef, 54)

for i in 1:54
    if get_node_data(dg, dg.nodes[i], "Raw Material") == 1
        node_cols[i] = :orange
    elseif get_node_data(dg, dg.nodes[i], "Product") == 1
        node_cols[i] = :yellow
    elseif get_node_data(dg, dg.nodes[i], "Intermediate") == 1
        node_cols[i] = :red
    else
        node_cols[i] = :deepskyblue3
    end
end

# Plot original graph
plt = plot_graph(dg, dag_positions = true, nlabels = dg.nodes, node_color = node_cols, nlabels_fontsize = 28)
#Makie.save((@__DIR__)*"/tech_pathway.png", plt)

# Plot the path between corn stover and nylon
# by default, the path will be black
plt_path = plot_graph_path(dg, "Corn Stover", "Nylon"; nlabels = dg.nodes, node_color = node_cols, nlabels_fontsize = 28)
#Makie.save((@__DIR__)*"/tech_pathway_path.png", plt_path)

# Plot the path between corn stover and nylon that
# goes through the intermediate C.-hexane
plt3 = plot_graph_path(dg, "Corn Stover", "Nylon", nlabels = dg.nodes, node_color = node_cols, intermediate = true, int = "C.-hexane", path_color = :purple)

for i in 1:length(dg.nodes)
    num_downstream = length(downstream_nodes(dg, dg.nodes[i]))
    num_upstream = length(upstream_nodes(dg, dg.nodes[i]))
    add_node_data!(dg, dg.nodes[i], num_downstream, "Number Downstream")
    add_node_data!(dg, dg.nodes[i], num_upstream, "Number Upstream")
end

node_size = get_node_data(dg, "Number Downstream") .+ 1
node_size = node_size ./ maximum(node_size) .* 200

# Plot the original graph with node sizes based on downstream connections
plt = plot_graph(dg, dag_positions = true, nlabels = dg.nodes, node_color = node_cols, node_size = node_size, nlabels_fontsize = 28)
#Makie.save((@__DIR__)*"/tech_pathway_downstream.png", plt)


fg = filter_edges(dg, 0, "Optimal Flow"; fn = Base.isequal)

plt_fg = plot_graph(fg, nlabels = dg.nodes, node_color = node_cols, edge_width = 5, nlabels_fontsize = 28)
#Makie.save((@__DIR__)*"/filtered_pathway.png", plt_fg)


# Compute metrics after removing different nodes in the graph
dg = get_tech_graph()
remove_node!(dg, "Ethylene")
compute_metrics!(dg)
ethylene_metrics_products = get_node_data(dg, ["Connected Products", "Number Downstream"])[1:5, :]
ethylene_metrics_raw = get_node_data(dg, ["Connected Raw", "Number Upstream"])[6:12, :]

dg = get_tech_graph()
remove_node!(dg, "Ter. Acid")
compute_metrics!(dg)
ter_acid_metrics_products = get_node_data(dg, ["Connected Products", "Number Downstream"])[1:5, :]
ter_acid_metrics_raw = get_node_data(dg, ["Connected Raw", "Number Upstream"])[6:12, :]


dg = get_tech_graph()
remove_node!(dg, "C.-hexane")
compute_metrics!(dg)
c_hexane_metrics_products = get_node_data(dg, ["Connected Products", "Number Downstream"])[1:5, :]
c_hexane_metrics_raw = get_node_data(dg, ["Connected Raw", "Number Upstream"])[6:12, :]
