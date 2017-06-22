module Plasmo

using Compat, MathProgBase,JuMP

export PlasmoGraph, PlasmoNode, PlasmoEdge, Node, Edge, NodeOrEdge,

#Plasmo graph functions
getnode,getedge,getnodeoredge,getconnectedto,getconnectedfrom,getsupportingnodes,getsupportingedges,

#Graphs package functions
add_vertex!,add_edge!,add_node!,add_subgraph!,create_node,
nv,vertices,edges,src,dst,ne,
getvertexindex,getedgeindex,getnodeindex,getnodes,getedges,out_edges,out_degree,
out_neighbors,in_edges,in_degree,in_neighbors,degree,is_connected,

#Attributes
addattribute!,rmattribute!,getattribute,hasattribute,_copy_subgraphs!,

#Plasmo Model
GraphModel,FlatGraphModel,create_flat_graph_model,getmodel,hasmodel,
getgraph,getnodes,getedges,getnodesandedges,getnodevariables,getnodeobjective,getnodeconstraints,getlinkconstraints,getlinkvarmap,getnodedata,getnode,is_graphmodel,
@linkconstraint,@getconstraintlist,setmodel!,solve,setsolution!,setvalue,

create_flat_graph_model


include("PlasmoGraph.jl")
include("PlasmoModel.jl")
include("PlasmoJuMP.jl")
include("macros.jl")


end

#@nodeconstraint,@nodevariable,@nodeobjective,
