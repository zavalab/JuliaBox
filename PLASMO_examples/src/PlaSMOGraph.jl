import Base.show
import Base.print
import Base.string

########################################################################################
# Graph
########################################################################################
type Graph
    name::AbstractString
    nodes#Node
    edges#Edge
    node_attributes::Dict#{Node,Dict{}}
    edge_attributes::Dict#{Edge,Dict{}}
    # node_callback_map::Dict#{Node,Vector{Function}}
    # edge_callback_map::Dict#{Edge,Vector{Function}}
    node_coupling_functions::Dict
    edge_coupling_functions::Dict
    node_coupling_args::Dict
    edge_coupling_args::Dict
    ext::Dict
    parent
    master
end

Graph() = Graph("graph"*string(gensym()),Node[],Edge[],Dict{Node,Dict}(),Dict{Edge,Dict}(),Dict{Node,Function}(),Dict{Edge,Function}(),Dict(),Dict(),Dict(),nothing,nothing)

#Generate a graph from an incidence matrix
function Graph(A::Matrix)
    graph = Graph()
    nrows, ncols = size(A)
    #Make sure each column has a single -1 and 1
    for i = 1:ncols
        @assert count(x -> x == -1, A[:,i]) == 1
        @assert count(x -> x ==  1, A[:,i]) == 1
        @assert count(x -> x == 0 , A[:,i]) == nrows - 2  #make sure rest of values are zero
    end
    #Add nodes
    for j = 1:nrows
        node = Node(graph)
    end
    #Add edges
    for i = 1:ncols
        from = graph.nodes[find(A[:,i] .== -1)[1]]
        to = graph.nodes[find(A[:,i] .== 1)[1]]
        edge = Edge(graph,from,to)
    end
    return graph
end

getnumnodes(graph::Graph) = length(graph.nodes)
getnumedges(graph::Graph) = length(graph.edges)
getnodes(graph::Graph) = graph.nodes
getedges(graph::Graph) = graph.edges
########################################################################################
# Node
########################################################################################
type Node
    id::Int
    name::AbstractString
    graph::Graph
    edges_in
    edges_out
    subgraph
end

function Node(graph::Graph;name = "node")
    id = getnumnodes(graph) + 1
    edges_in = Edge[]
    edges_out = Edge[]
    node = Node(id,name,graph,edges_in,edges_out,nothing)
    push!(graph.nodes,node)
    graph.node_attributes[node] = Dict()
    #graph.node_callback_map[node] = []
    #graph.node_coupling_fucntions[node] =
    return node
end

addnode!(graph::Graph) = Node(graph)
addnode!(graph::Graph,name) = Node(graph,name = name)

#add a graph as a node to an existing graph
function addgraph!(graph::Graph,subgraph::Graph)
    node = addnode!(graph)
    node.subgraph = subgraph
    node.subgraph.parent = node
    return node
end

function addgraph!(node::Node,subgraph::Graph)
    node.subgraph = subgraph
    node.subgraph.parent = node
    return node
end

getgraph(node::Node) = node.graph

#add a subnode to an existing node
function addnode!(node::Node)
    if node.subgraph == nothing         #if the existing node has no subgraph, create one
        node.subgraph = Graph()
        node.subgraph.parent = node     #set the parent of the subgraph to the existing node
        node = addnode!(node.subgraph)  #run the normal addnode! routine, and add a node to the subgraph
    else
        node = addnode!(node.subgraph)  #if the node already has a subgraph, add a new node to the subgraph
    end
    return node
end

edges_in(node::Node) = node.edges_in
edges_in(graph::Graph,node::Node) = filter(edge -> edge in graph.edges, node.edges_in)
edges_out(node::Node) = node.edges_out
edges_out(graph::Graph,node::Node) = filter(edge -> edge in graph.edges, node.edges_out)
getsupportingedges(node::Node) = Vector([node.edges_in;node.edges_out])

function setmaster!(graph::Graph,node::Node)
    graph.master = node
    #also change up connections here if I need to
end
##############################################################################
# Edge
##############################################################################
type Edge
    id::Int
    name::AbstractString
    graph::Graph
    node_pair::Pair{Node,Node}
end

function Edge(graph::Graph,node1::Node,node2::Node;name = "edge")
    id = getnumedges(graph) + 1
    node_pair = Pair{Node,Node}(node1,node2)
    edge = Edge(id,name,graph,node_pair)
    push!(graph.edges,edge)
    push!(node1.edges_out,edge)
    push!(node2.edges_in,edge)
    graph.edge_attributes[edge] = Dict()
    #graph.edge_callback_map[edge] = []
    return edge
end

function addedge!(node1::Node,node2::Node,name = "edge")
    @assert node1.graph == node2.graph
    graph = node1.graph
    edge = Edge(graph,node1,node2,name = name)
    return edge
end

getedges(graph::Graph) = graph.edges
getconnectedfrom(edge::Edge) = edge.node_pair.first
getconnectedto(edge::Edge) = edge.node_pair.second
getconnectednodes(edge::Edge) = Vector([edge.node_pair.first,edge.node_pair.second])

typealias NodeOrEdge Union{Node,Edge}

##############################################################################
# Querying Functions
##############################################################################
function neighbors(node::Node)
    #If I use an adjacency matrix, I can do fast computations for neighbors and subgraphs
    node_list = Node[]
    for edge in edges_in(node)
        push!(node_list,getconnectedfrom(edge))
    end
    for edge in edges_out(node)
        push!(node_list,getconnectedto(edge))
    end
    return node_list
end

function neighbors_in(node::Node)
    node_list = Node[]
    for edge in edges_in(node)
        push!(node_list,getconnectedfrom(edge))
    end
    return node_list
end

function neighbors_out(node::Node)
    node_list = Node[]
    for edge in edges_out(node)
        push!(node_list,getconnectedto(edge))
    end
    return node_list
end

#There is a more julia-esque way to write this stuff
function hasparentgraph(node::Node)
    if node.graph.parent != nothing
        return true
    else
        return false
    end
end

hassubgraph(node) = node.subgraph == nothing? false : true

neighbors(graph::Graph,node::Node) = filter(node -> node in graph.nodes, neighbors(node))
neighbors_in(graph::Graph,node::Node) = filter(node -> node in graph.nodes, neighbors_in(node))
neighbors_out(graph::Graph,node::Node) = filter(node -> node in graph.nodes, neighbors_out(node))

n_neighbors_in(node::Node) = length(neighbors_in(node))
n_neighbors_out(node::Node) = length(neighbors_out(node))
n_neighbors_in(graph::Graph,node::Node) = length(neighbors_in(graph,node))
n_neighbors_out(graph::Graph,node::Node) = length(neighbors_out(graph,node))

getsupportingedges(node::Node) = Vector([node.edges_in;node.edges_out])

function getedge(node1::Node,node2::Node)
    #Find the edge between two nodes if there is one
    for edge in unique([getsupportingedges(node1);getsupportingedges(node2)])
        nodes = getconnectednodes(edge)
        if all(node -> (node in nodes),[node1,node2])
            return edge
        end
    end
    return []
end

function getnodeinducedsubgraph(nodes)
    graph = nodes[1].graph  #get the graph from one of the nodes
    @assert all(node -> (node in graph.nodes),nodes)  #make sure all the nodes are in the same graph
    edges = []
    for node in nodes
        node_neighbors = neighbors(node)
        for neighbor in node_neighbors
            if neighbor in nodes
                push!(edges,getedge(node,neighbor)) #edge between node and its neighbor
            end
        end
    end
    graph_edges = unique(edges)
    sub_graph = Graph()
    sub_graph.nodes = nodes
    sub_graph.edges = graph_edges
    for node in nodes
        sub_graph.node_attributes[node] = graph.node_attributes[node]
        #delete!(graph.node_attributes,node)
    end
    for edge in graph_edges
        sub_graph.edge_attributes[edge] = graph.edge_attributes[edge]
    end

    #sub_graph.node_attributes = graph
    sub_graph.ext = copy(graph.ext)
    return sub_graph
end

##############################################################################
# Coupling Functions
##############################################################################
# Need to do checking for a lot of things in here

#nodes or edges might not have coupling arguments!
# getcouplingargs(graph,edge::Edge) = graph.edge_coupling_args[edge]
# getcouplingargs(graph,node::Node) = graph.node_coupling_args[node]

function getcouplingargs(graph,node::Node)
    #this can be done using standard julia one-liner constructs
    if haskey(graph.node_coupling_args,node)
        return graph.node_coupling_args[node]
    else
        return nothing
    end
end

function getcouplingargs(graph,edge::Edge)
    #this can be done using standard julia one-liner constructs
    if haskey(graph.edge_coupling_args,edge)
        return graph.edge_coupling_args[edge]
    else
        return nothing
    end
end

function setcouplingfunction(graph::Graph,node::Node,func::Function)
    graph.node_coupling_functions[node] = func
end

function setcouplingfunction(graph::Graph,edge::Edge,func::Function)
    graph.edge_coupling_functions[edge] = func
end

function setcouplingarguments(graph::Graph,edge::Edge,args...)
    graph.edge_coupling_args[edge] = args
end

function setcouplingarguments(graph::Graph,node::Node,args...)
    graph.node_coupling_args[node] = args
end

#Set coupling functions given the higher level node
function setcouplingfunction(master::Node,child::Node,func::Function)
    graph = master.subgraph
    setcouplingfunction(graph,child,func)
end

function setcouplingfunction(master::Node,child::Edge,func::Function)
    graph = master.subgraph
    setcouplingfunction(graph,child,func)
end
##############################################################################
# Attributes
##############################################################################
#Function to add an attribute to a node or edge
function addattribute!(item::NodeOrEdge,attribute::Symbol,value)
    if isa(item,Node)
        item.graph.node_attributes[item][attribute] = value
    elseif isa(item,Edge)
        item.graph.edge_attributes[item][attribute] = value
    else
        error("Can only add attributes to nodes or edges")
    end
end

#add an attribute to a graph
function addattribute!(graph::Graph,attribute::Symbol,value)
    graph.ext[attribute] = value
end


#Function to map attributes to multiple nodes or edges
function mapattribute!(graph::Graph,nodes_or_edges::Vector{NodeOrEdge},attribute::Symbol,value)
    for ne in nodes_or_edges
        addattribute!(graph,ne,attribute,value)
    end
end

function mapattribute!(graph::Graph,nodes::Vector{Node},edges::Vector{Edge},attribute::Symbol,value)
    for node in nodes
        addattribute!(graph,node,attribute,value)
    end
    for edge in edges
        addattribute!(graph,edge,attribute,value)
    end
end

mapattribute!(graph,attribute,value) = mapattribute!(graph,getnodes(graph),getedges(graph),attribute,value)

getattributenames(node::Node) = collect(keys(node.graph.node_attributes[node]))
getattributenames(edge::Edge) = collect(keys(edge.graph.edge_attributes[edge]))
getattributenames(graph::Graph) = collect(keys(graph.ext))

getattribute(node::Node,attr::Symbol) = node.graph.node_attributes[node][attr]
getattribute(edge::Edge,attr::Symbol) = edge.graph.edge_attributes[edge][attr]
getattribute(graph::Graph,attr::Symbol) = graph.ext[attr]
hasattribute(ne::NodeOrEdge,attr::Symbol) = attr in getattributenames(ne)
hasattribute(graph::Graph,attr::Symbol) = attr in getattributenames(graph)

#Get every node in the graph
function getnodelist(graph::Graph)
    nodelist = Node[]
    for node in getnodes(graph)
        if node.subgraph != nothing
            sublist = getnodelist(node.subgraph)
            nodelist = [nodelist;sublist]
        end
    end
    nodelist = [nodelist;getnodes(graph)]
    return nodelist
end
#Get every edge in the graph
function getedgelist(graph::Graph)
    edgelist = Node[]
    for node in getnodes(graph)
        if node.subgraph != nothing
            sublist = getedgelist(node.subgraph)
            edgelist = [edgelist;sublist]
        end
    end
    edgelist = [edgelist;getedges(graph)]
    return edgelist
end

#Get nodes at the bottom level of the graph
function getbottomnodes(graph::Graph)
    nodelist = Node[]
    for node in getnodes(graph)
        if node.subgraph != nothing
            sublist = getbottomnodes(node.subgraph)
            append!(nodelist,sublist)
        else
            push!(nodelist,node)
        end
    end
    return nodelist
end

#get all edges at the bottom level of the graph
function getbottomedges(graph::Graph)
    edgelist = Edge[]
    for node in getnodes(graph)
        if node.subgraph != nothing
            sublist = getbottomedges(node.subgraph)
            append!(edgelist,sublist)
        else
            append!(edgelist,getsupportingedges(node))
        end
    end
    edgelist = unique(edgelist)
    return edgelist
end

#Aggregate nodes.  Combine models into a single model that combines same variable names
function aggregate!(nodes::Array{Node})
    #make sure graph context is not the graph of the given nodes
    #make sure all nodes are in same graph
    g1 = nodes[1].graph  #get graph context for one node
    @assert all(node -> (node in g1.nodes),nodes)
    agraph = getnodeinducedsubgraph(nodes)
    anode = addnode!(g1)  #add new aggragated node the graph

    #Remove nodes and attributes from original graph.  They are now in the aggregated node
    filter!(n -> !(n in agraph.nodes), g1.nodes)
    filter!(e -> !(e in agraph.edges), g1.edges)
    for node in agraph.nodes
        delete!(g1.node_attributes,node)
        node.graph = agraph
    end
    for edge in agraph.edges
        delete!(g1.edge_attributes,edge)
        edge.graph = agraph
    end

    new_coupling_funcs = Dict()
    new_coupling_args = Dict()
    for edge in agraph.edges
        new_coupling_funcs[edge] = g1.edge_coupling_functions[edge]
        new_coupling_args[edge] = g1.edge_coupling_args[edge]
    end
    #Should I copy over node couplings as well?
    agraph.edge_coupling_functions = new_coupling_funcs
    agraph.edge_coupling_args = new_coupling_args
    agraph.parent = anode
    anode.subgraph = agraph
    return anode
end

##########################################################
# Printing overloads
##########################################################
macro display(objecttype)
    expr = quote
      string(object::$objecttype) = object.name * string(object.id)
      print(io::IO, object::$objecttype) = print(io, string(object))
      show(io::IO,object::$objecttype) = print(io,object)
    end
    return esc(expr)
end
@display Node
@display Edge
string(graph::Graph) = graph.name
print(io::IO, graph::Graph) = print(io, string(graph))
show(io::IO,graph::Graph) = print(io,graph)
