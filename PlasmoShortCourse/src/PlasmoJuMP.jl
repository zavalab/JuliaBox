#This file contains all of the constructs to create and manage a JuMP GraphModel.  The idea is that you use PLASMO to create your graph, associate models, and build the
import JuMP:AbstractModel,Model,Variable,ConstraintRef,getvariable,@variable,@constraint,@objective,GenericQuadExpr,GenericAffExpr,solve,setvalue
import MathProgBase

#typealias GenericExpr Union{GenericQuadExpr,GenericAffExpr}
type NodeData
    objective#::GenericExpr                      #Individual objective expression....    #Need the nlp evaluator with the Expr graph to do this?
    variablemap::Dict{Symbol,Any}                #Dictionary of symbols to Model variables
    constraintlist::Vector{ConstraintRef}        #Vector of model constraints (make this a dictionary too)
    indexmap::Dict{Int,Int}                      #linear index of node variable in flat model to the original index of the component model
end
NodeData() = NodeData(0,Dict{Symbol,Any}(),ConstraintRef[],Dict{Int,Int}())

type LinkData
    linkconstraintmap::Dict{PlasmoGraph,Vector{AbstractConstraint}}  #keep track of constraints linking nodes (or edges) to their neighbors
    variablemap::Dict{Variable,NodeOrEdge}       #map which variables belong to which neighbors
end
LinkData() = LinkData(Dict{PlasmoGraph,Vector{AbstractConstraint}}(),Dict{Variable,NodeOrEdge}())

is_nodevar(nodeoredge::NodeOrEdge,var::Variable) = getmodel(nodeoredge) == var.m #checks whether a variable belongs to a node or edge

#Store link constraint references for nodes and edges in their data
function _addlinkconstraint!(graph::PlasmoGraph,nodeoredge::NodeOrEdge,con::AbstractConstraint)
    local_vars = Vector{Variable}() #this should be a unique set
    neighbor_vars = Dict{Variable,NodeOrEdge}()
    vars = con.terms.vars
    for var in vars
        if is_nodevar(nodeoredge,var)  #if the variable is on this node or edge
            push!(local_vars,var)
        elseif is_connected(graph,nodeoredge,getnode(var))
            neighbor_vars[var] = getnode(var)  #should be what it's connecte to!
        else
            error("A link constraint contains variables not connected to $nodeoredge")
        end
    end

    for var in local_vars
        nodeoredge.attributes[:LinkData].variablemap[var] = nodeoredge
    end
    merge!(nodeoredge.attributes[:LinkData].variablemap,neighbor_vars)

    #add the constraint data to the node or edge
    haskey(nodeoredge.attributes[:LinkData].linkconstraintmap,graph)? nothing : nodeoredge.attributes[:LinkData].linkconstraintmap[graph] = Vector{AbstractConstraint}()
    push!(nodeoredge.attributes[:LinkData].linkconstraintmap[graph],con)
end

function _addlinkconstraint!{T}(graph::PlasmoGraph,nodeoredge::NodeOrEdge,cons_refs::Array{AbstractConstraint,T})
    array_type = typeof(cons_refs)
    array_type.parameters.length > 1? cons_refs = vec(cons_refs): nothing
    local_vars = Vector{Variable}() #this should be a unique set
    neighbor_vars = Dict{Variable,NodeOrEdge}()
    vars_checked = Vector{Variable}()
    #Ensure that variables in the constraint are part of nodeoredge or its neighbors)
    for con in cons_refs
        vars = con.terms.vars
        for var in vars
            if !(var in vars_checked)
                if is_nodevar(nodeoredge,var)  #if the variable is on this node or edge
                    push!(local_vars,var)
                elseif is_connected(graph,nodeoredge,getnode(var))
                    neighbor_vars[var] = getnode(var)
                else
                    error("A link constraint contains variables not connected to $nodeoredge")
                end
                push!(vars_checked,var)
            end
        end
    end
    #Add the variable information to the variable map
    for var in local_vars
        nodeoredge.attributes[:LinkData].variablemap[var] = nodeoredge
    end
    merge!(nodeoredge.attributes[:LinkData].variablemap,neighbor_vars)
    #add the constraint data to the node or edge
    haskey(nodeoredge.attributes[:LinkData].linkconstraintmap,graph)? nothing : nodeoredge.attributes[:LinkData].linkconstraintmap[graph] = Vector{AbstractConstraint}()
    append!(nodeoredge.attributes[:LinkData].linkconstraintmap[graph],cons_refs)
end

#Construct a Structured Graph Model from multiple node and edge models (use this for parallel solver interfaces)
function GraphModel()
    m = JuMP.Model()
    m.ext[:Graph] = Plasmo.PlasmoGraph()
    m.ext[:model_type] = :GraphModel
    #m.ext[:GlobalConstraints]  #possibly track this as separate information for global constraints (e.g. sum of all flows in a network)
    return m
end

#Construct a structured model, but roll it all into one JuMP model (this is how we solve with JuMP accessible solvers)
function FlatGraphModel()
    m = JuMP.Model()
    m.ext[:Graph] = Plasmo.PlasmoGraph()
    m.ext[:model_type] = :FlatModel
    return m
end

is_structured_graph_model(m::Model) = m.ext[:model_type] == :GraphModel? true : false
is_flat_graph_model(m::Model) = m.ext[:model_type] == :FlatModel? true : false
is_graphmodel(m::Model) = haskey(m.ext,:Graph)? true : false

#Create a PlasmoModel from a PlasmoGraph (think of the PlasmoGraph as a template)
#TODO
#GraphModel(graph::PlasmoGraph) = generate_graph_model(graph)

function add_node!(m::Model; index = nv(getgraph(m).graph)+1)
    @assert is_graphmodel(m)
    if is_structured_graph_model(m)
        node = PlasmoNode(Dict{PlasmoGraph,Int}(), Symbol("node"),Dict{Any,Any}(:Model => Model(),:LinkData => LinkData()))
    else #it's a flatgraph
        node = PlasmoNode(Dict{PlasmoGraph,Int}(), Symbol("node"),Dict{Any,Any}(:NodeData => NodeData(),:LinkData => LinkData()))
    end
    add_node!(getgraph(m),node,index = index)
    return node
end

function add_edge!(m::Model,node1::PlasmoNode,node2::PlasmoNode)
    @assert is_graphmodel(m)
    if is_structured_graph_model(m)
        edge = PlasmoEdge(Dict{PlasmoGraph,Edge}(), Symbol("edge"),Dict{Any,Any}(:Model => Model(),:LinkData => LinkData()))
    else
        edge = PlasmoEdge(Dict{PlasmoGraph,Edge}(), Symbol("edge"),Dict{Any,Any}(:NodeData => NodeData(),:LinkData => LinkData()))
    end
    add_edge!(getgraph(m),edge,node1,node2)
    return edge
end
#Define all of the PlasmoGraph functions for a PlasmoModel
getgraph(m::Model) = haskey(m.ext, :Graph)? m.ext[:Graph] : error("Model does not have a graph")
getnodes(m::Model) = getnodes(getgraph(m))
getedges(m::Model) = getedges(getgraph(m))
getnode(m::Model,id::Integer) = getnodes(getgraph(m))[id]
getedge(m::Model,id::Edge) = getedges(getgraph(m))[id]
#write node constraints, edge constraints, and coupling constraints
getnodedata(nodeoredge::NodeOrEdge) = getattribute(nodeoredge,:NodeData)
getnodeobjective(nodeoredge::NodeOrEdge) = hasattribute(nodeoredge,:NodeData)? getattribute(nodeoredge,:NodeData).objective : JuMP.getobjective(getmodel(nodeoredge))
getnodevariables(nodeoredge::NodeOrEdge) =  hasattribute(nodeoredge,:NodeData)? getattribute(nodeoredge,:NodeData).variablemap : getmodel(nodeoredge).objDict
getnodeconstraints(nodeoredge::NodeOrEdge) = getattribute(nodeoredge,:NodeData).constraintlist
# getvariable(m::Model,index::Integer,s::Symbol) = getvariable(getnode(getgraph(m),index),s)
# getvariable(nodeoredge::NodeOrEdge,s::Symbol) = getattribute(nodeoredge,:NodeData).variablemap[s]

#is_connected(graph::PlasmoGraph,n1::PlasmoNode,n2::PlasmoNode)
#TODO
#get all of the link constraints from a JuMP model
function getlinkconstraints(m::JuMP.Model)
    @assert is_graphmodel(m)
    cons = Dict()
    for nodeoredge in getnodes(getgraph(m))
        push!(cons,getlinkconstraints(node))
    end
    return cons
    #check that it's a Graph model
    #inspect the graph nodes and edges.  return constraint references for their corresponding link constraints
end

#This will return all of the constraint references that connect to the node or edge in each of its graphs
getlinkconstraints(nodeoredge::NodeOrEdge) = nodeoredge.attributes[:LinkData].linkconstraintmap
getlinkvarmap(nodeoredge::NodeOrEdge) = nodeoredge.attributes[:LinkData].variablemap
#copy the subgraph structure from one graph to another
function _copy_subgraphs!(graph1::PlasmoGraph,graph2::PlasmoGraph)
    for i = 1:length(graph1.subgraphlist)
        subgraph = PlasmoGraph()
        add_subgraph!(graph2,subgraph)
    end
end

#create a copy of the graph.  Make sure indices of everything matches up.
function create_flat_graph_model(m::Model)
    @assert is_structured_graph_model(m)
    graph = getgraph(m)
    flat_model = FlatGraphModel()
    flat_graph = getgraph(flat_model)
    #copy number of subgraphs (might need recursive function here!)
    _copy_subgraphs!(graph,flat_graph)
    #first copy all nodes, then setup all the subgraphs
    for (index,node) in getnodes(graph)
        new_node = add_node!(flat_model,index = index)  #create the node and add a vertex to the top level graph.  We pass the index explicity for this graph
        node_index = getindex(node) #returns dict of {graph => index}
        #setup subgraph references
        for igraph in keys(node_index)
            if igraph.index != 0 #if it's not the top level graph
                graph_index = igraph.index #the index of this subgraph
                subgraph = flat_graph.subgraphlist[graph_index]
                add_node!(subgraph,new_node,index = node.index[igraph.subgraphlist[index]])
                #new_node.index[subgraph] = node.index[igraph.subgraphlist[index]]  #make sure indices match the original graph
            end
        end
        if hasmodel(node)
            node_model = getmodel(node)
            _buildnodemodel!(flat_model,new_node,node_model)
        end
    end

    #copy edges  #edge indices need to work like node indices
    for (index,edge) in getedges(graph)
        pair = getindex(graph,edge)
        new_nodes = getsupportingnodes(flat_graph,pair)
        new_edge = add_edge!(flat_model,new_nodes[1],new_nodes[2])
        #new_edge.index[flat_graph] = index
        for igraph in keys(edge.index)
            if igraph.index != 0
                index = igraph.index
                subgraph = flat_graph.subgraphlist[index]
                pair = getindex(igraph.subgraphlist[index],new_edge)
                add_edge!(subgraph,pair)
            end
        end

        if hasmodel(edge)
            edge_model = getmodel(edge)
            _buildnodemodel!(flat_model,new_edge,edge_model)
        end
    end

    #add the linking constraints
    #inspect the link constraints, and map them to variables within flat model
    for (index,ne) in getnodesandedges(graph)
        #get one list of linkconstraints for each node
        graphlinkconstraints = getlinkconstraints(ne)
        linkvarmap = getlinkvarmap(ne)  #{variable => nodeoredge}
        indexmap = Dict() #{node variable => flat variable index} Need index of node variables to flat model variables
        for (var,nodeoredge) in linkvarmap
            var_index = JuMP.linearindex(var)
            node_index = getindex(graph,nodeoredge)                 #node index in graph.  This could be a problem with edges and subgraphs
            flat_nodeoredge = getnodeoredge(flat_graph,node_index)  #should be the corresponding node
            flat_indexmap = getnodedata(flat_nodeoredge).indexmap
            indexmap[var] = flat_indexmap[var_index]
        end
        link_cons = []
        for con in values(graphlinkconstraints)
            append!(link_cons,con)
        end
        for linkconstraint in link_cons
            t = []
            for terms in linearterms(linkconstraint.terms)
                push!(t,terms)
            end
            con_reference = @constraint(flat_model, linkconstraint.lb <= sum(t[i][1]*JuMP.Variable(flat_model,indexmap[(t[i][2])]) for i = 1:length(t)) + linkconstraint.terms.constant <= linkconstraint.ub)
        end
    end
    #sum the objectives by default
    @objective(flat_model,Min,sum(getnodeobjective(nodeoredge) for nodeoredge in values(getnodesandedges(flat_graph))))
    return flat_model
end

#Function to build a node model for a flat graph model
#THIS IS ALL GOING TO BREAK WHEN JuMP and MathProgBase CHANGE
function _buildnodemodel!(m::Model,nodeoredge::NodeOrEdge,node_model::Model)
    #@assert nodeoredge in [getgraph(m).nodes;getgraph(m).edges]
    num_vars = MathProgBase.numvar(node_model)
    var_map = Dict()              #this dict will map linear index of the node model variables to the new model JuMP variables {index => JuMP.Variable}
    node_map = Dict()             #nodemap. {varkey => [var1,var2,...]}
    index_map = Dict()            #{var index in node => var index in flat model}
    #add the node model variables to the new model
    for i = 1:num_vars
        x = JuMP.@variable(m)            #create an anonymous variable
        setlowerbound(x,node_model.colLower[i])
        setupperbound(x,node_model.colUpper[i])
        var_name = string(Variable(node_model,i))
        new_name = "$(nodeoredge.label)$(nodeoredge.index[getgraph(m)])."*var_name
        setname(x,new_name)       #rename the variable to the node model variable name plus the node or edge name
        setcategory(x,node_model.colCat[i])                                  #set the variable to the same category
        setvalue(x,node_model.colVal[i])                                     #set the variable to the same value
        var_map[i] = x                                                       # map the linear index of the node model variable to the new variable
        index_map[i] = linearindex(x)
        m.objDict[Symbol(new_name)] = x #Update master model variable dictionary
    end
    #setup the node_map dictionary.  This maps the node model's variable keys to variables in the newly constructed model.
    for key in keys(node_model.objDict)  #this contains both variable and constraint references
        if isa(node_model.objDict[key],Union{JuMP.JuMPArray{Variable},Array{Variable}})     #if the JuMP variable is an array or a JuMPArray
            vars = node_model.objDict[key]
            isa(vars,JuMP.JuMPArray)? vars = vars.innerArray : nothing
            # if isa(vars,JuMP.JuMPArray)
            #     vars = vars.innerArray
            # end
            dims = JuMP.size(vars)
            node_map[key] = Array{JuMP.Variable}(dims)
            for j = 1:length(vars)
                var = vars[j]
                node_map[key][j] = var_map[linearindex(var)]
            end
        elseif isa(node_model.objDict[key],Variable) #else it's a single variable
            node_map[key] = var_map[linearindex(node_model.objDict[key])]
            #node_map[key] = var_map[node_model.objDict[key].col]
        end
    end

    getattribute(nodeoredge,:NodeData).variablemap = node_map
    getattribute(nodeoredge,:NodeData).indexmap = index_map
    #getnodevariables(nodeoredge) = node_map
    #copy the linear constraints to the new model
    #cons_refs = ConstraintRef[]
    for i = 1:length(node_model.linconstr)
        con = node_model.linconstr[i]
        #t = collect(linearterms(con.terms))  #This is broken in julia 0.5
        t = []
        for terms in linearterms(con.terms)
            push!(t,terms)
        end
        reference = @constraint(m, con.lb <= sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + con.terms.constant <= con.ub)
        #@constraint(m,reference, con.lb <= sum(t[j][1]*var_map[linearindex(t[j][2])] for j = 1:length(t)) + con.terms.constant <= con.ub)
        #@constraint(m,reference, )
        push!(getattribute(nodeoredge,:NodeData).constraintlist,reference)
    end
    #getattribute(nodeoredge,:NodeData).constraintlist = cons_refs
    #Copy the non-linear constraints to the new model
    d = JuMP.NLPEvaluator(node_model)         #Get the NLP evaluator object.  Initialize the expression graph
    MathProgBase.initialize(d,[:ExprGraph])
    num_cons = MathProgBase.numconstr(node_model)
    for i = 1:num_cons
        if !(MathProgBase.isconstrlinear(d,i))    #if it's not a linear constraint
            expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
            _splicevars!(expr,var_map)              #splice the variables from var_map into the expression
            con = JuMP.addNLconstraint(m,expr)    #raw expression input for non-linear constraint
            push!(getattribute(nodeoredge,:NodeData).constraintlist,con)  #Add the nonlinear constraint reference to the node
        end
    end

    # #If the objective is linear, store it as a node object
    getobjectivesense(node_model) == :Min? sense = 1: sense = -1
    if MathProgBase.isobjlinear(d)                #get the node model objective.  set it to a node attribute
        t = []
        for terms in linearterms(node_model.obj.aff)
            push!(t,terms)
        end
        #t = collect(linearterms(node_model.obj.aff))
        #Make the objective a minimization

        obj = @objective(m,Min,sense*sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + node_model.obj.aff.constant)
        getattribute(nodeoredge,:NodeData).objective = m.obj
    elseif MathProgBase.isobjquadratic(d)
        #Get the linear terms
        t = []
        for terms in linearterms(node_model.obj.aff)
            push!(t,terms)
        end
        #Get the quadratic terms
        qcoeffs = node_model.obj.qcoeffs
        qvars1 = node_model.obj.qvars1
        qvars2 = node_model.obj.qvars2
        obj = @objective(m,Min,sense*(sum(qcoeffs[i]*var_map[linearindex(qvars1[i])]*var_map[linearindex(qvars2[i])] for i = 1:length(qcoeffs)) + sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + node_model.obj.aff.constant))
        getattribute(nodeoredge,:NodeData).objective = m.obj
    end
    return m
    #I don't have the non-linear objective yet, but it shouldn't be any different than the constraints.
end

#splice variables into a constraint expression
function _splicevars!(expr::Expr,var_map::Dict)
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splicevars!(expr.args[i],var_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]     #this is the actual index in x[1], x[2], etc...
                new_var = :($(var_map[var_index]))   #get the JuMP variable from var_map using the index
                expr.args[i] = new_var               #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end

# function getnumvars(nodeoredge::NodeOrEdge)
# end



#copy the solution from one graph to another where nodes and variables match
function set_solution!(graph1::PlasmoGraph,graph2::PlasmoGraph)
    for (index,nodeoredge) in getnodesandedges(graph1)
        nodeoredge2 = getnodeoredge(graph2,index)       #get the corresponding node or edge in graph2
        # for i = 1:MathProgBase.getnumvars(nodeoredge)   #for each variable in the original node or edge
        #     var = JuMP.Variable()
        for (key,var) in getnodevariables(nodeoredge)
            var2 = nodeoredge2[key]
            vals = JuMP.getvalue(var)
            setvalue(var2,vals)  #the dimensions have to line up
        end
    end
end

setvalue(jarr1::JuMP.JuMPArray,jarr2::JuMP.JuMPArray) = setvalue(jarr1.innerArray,jarr2.innerArray)
setvalue(jarr1::JuMP.JuMPArray,jarr2::Array) = setvalue(jarr1.innerArray,jarr2)
setvalue(jarr1::Array,jarr2::JuMP.JuMPArray) = setvalue(jarr1,jarr2.innerArray)

#need to get solution data back into original models
function solve(m::JuMP.Model,graph::PlasmoGraph)
    println("Creating flattened graph model...")
    m_flat = create_flat_graph_model(m)
    println("Finished model instantiation")
    m_flat.solver = m.solver
    status = JuMP.solve(m_flat)
    if status == :Optimal
        #Now get our solution data back into the original model
        set_solution!(getgraph(m_flat),graph)
    end
end

# function Plasmo.solve(m::JuMP.Model)
#     if is_graphmodel(m)
#         println("Creating flattened graph model...")
#         m_flat = create_flat_graph_model(m)
#         println("Finished model instantiation")
#         m_flat.solver = m.solver
#         status = JuMP.solve(m_flat)
#         if status == :Optimal
#             #Now get our solution data back into the original model
#             _copy_solution!(getgraph(m_flat),graph)
#         end
#     else
#         JuMP.solve(m)
#     end
# end
