using Revise
using PlasmoData
using DataGraphPlots

dg = DataGraph()

dg = DataGraph{Int, Any, Any, Any, Matrix{Any}, Matrix{Any}}()

add_node!(dg, 1)
add_node!(dg, 2)
add_node!(dg, 3)
add_node!(dg, "node4")
add_node!(dg, :node5)


PlasmoData.add_edge!(dg, 1, 2)
PlasmoData.add_edge!(dg, 2, 3)
PlasmoData.add_edge!(dg, "node4", 1)
PlasmoData.add_edge!(dg, (:node5, 2))
PlasmoData.add_edge!(dg, (3, "node4"))

add_node_data!(dg, 1, [6, 3, 4], "node_weight")
add_node_data!(dg, 2, 3.4, "node_weight")
add_node_data!(dg, 3, "this is on node 3", "node_weight")
add_node_data!(dg, "node4", [1 2; 3 4], "node_weight")
add_node_data!(dg, :node5, DataGraph(), "node_weight")

add_edge_data!(dg, 1, 2, DataGraph(), "edge_weight")
add_edge_data!(dg, 2, 3, [1 2 ; 5 7], "edge_weight")
add_edge_data!(dg, "node4", 1, 1.0, "edge_weight")
add_edge_data!(dg, (:node5, 2), -.00001, "edge_weight")
add_edge_data!(dg, (3, "node4"), Dict(), "edge_weight")

add_graph_data!(dg, 1.0, "graph_weight")

DataGraphPlots.plot_graph(dg; xdim = 400, ydim = 400)
DataGraphPlots.plot_graph(dg; xdim = 400, ydim = 400, save_fig = true, fig_name = (@__DIR__)*"/images/general_overview.pdf")
