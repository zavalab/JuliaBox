#A PlasmoGraph has special attributes for managing structured models
import JuMP:AbstractModel,AbstractConstraint,Variable,ConstraintRef
import MathProgBase

#Specialized Attributes
#Set a coupling function attributes
setcouplingfunction!(graph::PlasmoGraph,ne::NodeOrEdge,func::Function) = push!(ne.attributes[:couplingfunctions][graph],func)
getcouplingfunctions(ne::NodeOrEdge) = ne.attributes[:couplingfunctions]
setcouplingarguments!(graph::PlasmoGraph,ne::NodeOrEdge,func::Function,args...) = ne.attributes[:couplingarguments][graph][func] = args

_assignedtonode(m::AbstractModel) = haskey(m.ext,:node)

#set a variable name when its model gets assigned to a node or edge
_setvarname(nodeoredge::NodeOrEdge,v::Variable) = JuMP.setname(v,string(nodeoredge.label)*"["*getname(v)*"]")

#Component Models (A node or edge can have a single model)
function setmodel!(nodeoredge::NodeOrEdge,m::AbstractModel)
    #make sure the model isn't already assigned to a different node
    @assert !_assignedtonode(m)
    nodeoredge.attributes[:model] = m
    m.ext[:node] = nodeoredge
    #set variable names
    for i = 1:MathProgBase.numvar(m)
        _setvarname(nodeoredge,Variable(m,i))
    end
end

getnode(m::AbstractModel) = m.ext[:node]
getnode(var::Variable) = var.m.ext[:node]

removemodel!(nodeoredge::NodeOrEdge) = nodeoredge.attributes[:model] = nothing
getmodel(nodeoredge::NodeOrEdge) = nodeoredge.attributes[:model]
hasmodel(nodeoredge::NodeOrEdge) = haskey(nodeoredge.attributes,:model)
getindex(nodeoredge::NodeOrEdge,s::Symbol) = getmodel(nodeoredge)[s]  #get a node or edge variable

#run coupling functions for nodes and edges
function _runcouplingfunctions(m::JuMP.Model,graph::Graph)
    for node in getnodes(graph)
        if node.subgraph != nothing
            _runcouplingfunctions(m,node.subgraph)
        end
    end
    #coupling arguments are custom arguments to pass to coupling functions
    for (node, func) in graph.node_coupling_functions
        coupling_args = getcouplingargs(graph,node)
        #new_node = getnode(new_graph,node.id)

        if coupling_args != nothing
            func(m,node,coupling_args...)
        else
            func(m,node)
        end
    end
    for (edge,func) in graph.edge_coupling_functions
        coupling_args = getcouplingargs(graph,edge)
        #new_edge = getedge(new_graph,edge.id)
        #println(coupling_args)
        if coupling_args != nothing
            func(m,edge,coupling_args...)
        else
            func(m,edge)
        end
    end
end

#Create a Plasmo Model from a graph template
function generate_graph_model(graph::Graph)
    m = Model()
    new_graph = copy(graph)      #Copy the topology
    m.ext[:Graph] = new_graph

    #start using iterators for nodes and edges
    for node in getnodelist(new_graph)
        addattribute!(node,:GraphModelData,GraphModelData())
    end

    for edge in getedgelist(new_graph)
        addattribute!(edge,:GraphModelData,GraphModelData())
    end

    nodelist = getnodelist(graph)  #all nodes
    edgelist = getedgelist(graph)  #all edges

    #build each model
    objs = []
    for i = 1:length(edgelist)
        edge = edgelist[i]
        new_edge = getedgelist(new_graph)[i]
        _buildgraphnodemodel!(edge,new_edge,m)
        if hasmodel(edge)
            push!(objs,getnodeobjective(new_edge))
        end
    end
    for i = 1:length(nodelist)
        node = nodelist[i]
        new_node = getnodelist(new_graph)[i]
        _buildgraphnodemodel!(node,new_node,m)
        if hasmodel(node)
            push!(objs,getnodeobjective(new_node))
        end
    end

    # for node in getnodes(graph)
    #     if node.subgraph != nothing
    #         _runcouplingfunctions(m,node.subgraph)
    #     end
    # end
    _runcouplingfunctions(m,getgraph(m))
    #@objective(m,Min,sum(objs))
    #println(objs)
    m.obj = sum(objs)
    return m
end
