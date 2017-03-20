using JuMP
import JuMP.getvariable
import MathProgBase

type JuMPData
    objective             #Individual objective function for a node.
    vardict::Dict         #Dictionary of symbols to JuMP variables
    constraints::Vector   #Vector of model constraints on this node
    model::JuMP.Model     #The actual JuMP model the user supplies for a node
    solndict::Dict        #Current solution held in memory
end

JuMPData() = JuMPData(0,Dict{Any,Any}(),[],nothing,Dict{Any,Any}())
JuMPData(model::JuMP.Model) = JuMPData(0,Dict{Any,Any}(),[],model,Dict{Any,Any}())

#setmodel!(ne::NodeOrEdge,m::JuMP.Model) = addattribute!(ne,:JuMPData,JuMPData(m))
function setmodel!(ne::NodeOrEdge,m::JuMP.Model)
    if hasattribute(ne,:JuMPData)
        getattribute(ne,:JuMPData).model = m
    else
        addattribute!(ne,:JuMPData,JuMPData(m))
    end
end

getmodel(ne::NodeOrEdge) = getattribute(ne,:JuMPData).model

function hasmodel(ne::NodeOrEdge)
    bool = false
    if :JuMPData in getattributenames(ne)
        bool = getmodel(ne) != nothing
    else
        bool = false
    end
    return bool
end

vardict(ne::NodeOrEdge) = getattribute(ne,:JuMPData).vardict
constraints(ne::NodeOrEdge) = getattribute(ne,:JuMPData).constraints
objective(ne::NodeOrEdge) = getattribute(ne,:JuMPData).objective
getvariable(ne::NodeOrEdge,s::Symbol) = vardict(ne)[s]

function addnode!(graph::Graph,m::JuMP.Model)
    node = Node(graph)
    setmodel!(node,m)
    return node
end

function addnode!(node::Node,m::JuMP.Model)
    subnode = addnode!(node)
    setmodel!(subnode,m)
    return subnode
end

# function copymodel(m::JuMP.Model)
#     num_vars = m.numCols
# end

function buildnodemodel!(ne::NodeOrEdge,m::JuMP.Model)
    if hasmodel(ne)   #if the node or edge has a JuMP model attribute
        node_model = getattribute(ne,:JuMPData).model
        num_vars = node_model.numCols
        v_map = Dict()              #this dict will map linear index of the node model variables to the new model JuMP variables {index => JuMP.Variable}
        node_map = Dict()           #nodemap. {varkey => [var1,var2,...]}

        #add the node model variables to the new model
        for i = 1:num_vars
            @variable(m,node_model.colLower[i] <= x <= node_model.colUpper[i])   #just call the new variable x
            var_name = string(Variable(node_model,i))
            setname(x,"$(ne.name)$(ne.id)_"*var_name)       #rename the variable to the node model variable name plus the node or edge name
            setcategory(x,node_model.colCat[i])                                  #set the variable to the same category
            setvalue(x,node_model.colVal[i])                                     #set the variable to the same value
            v_map[i] = x                                                         #map the linear index of the node model variable to the new variable
            m.varDict[Symbol("$(ne.name)$(ne.id)_"*var_name)] = x
        end
        delete!(m.varDict,:x)

        #setup the node_map dictionary.  This maps the node model's variable keys to variables in the newly constructed model.
        #It's basically a copy of the varDict.  It probably will be the same as varDict eventually.
        #this also depends on the user not messing up node_model.varDict
        for key in keys(node_model.varDict)
            #node_map[key] = []
            if isa(node_model.varDict[key],Union{JuMP.JuMPArray,Array})     #if the JuMP variable is an array or a JuMPArray
                vars = node_model.varDict[key]
                if isa(vars,JuMP.JuMPArray)
                    vars = vars.innerArray
                end
                dims = JuMP.size(vars)
                node_map[key] = Array{JuMP.Variable}(dims)
                for j = 1:length(vars)
                    var = vars[j]
                    node_map[key][j] = v_map[var.col]
                end
                #for var in node_model.varDict[key][:]                       #need to keep array structure somehow
                #    push!(node_map[key],v_map[var.col])
                #end
            else
                node_map[key] =  [v_map[node_model.varDict[key].col]]       #else it's a single variable
                #push!(node_map[key],v_map[node_model.varDict[key].col])
            end
        end
        getattribute(ne,:JuMPData).vardict = node_map  #set the node vardict attribute

        #copy the linear constraints to the new model
        for i = 1:length(node_model.linconstr)
            con = node_model.linconstr[i]
            t = collect(linearterms(con.terms))
            @constraint(m,reference,con.lb <= sum{t[i][1]*v_map[linearindex(t[i][2])],i = 1:length(t)} + con.terms.constant <= con.ub)
            push!(constraints(ne),reference)
        end

        #Copy the non-linear constraints to the new model
        d = JuMP.NLPEvaluator(node_model)         #Get the NLP evaluator object.  Initialize the expression graph
        MathProgBase.initialize(d,[:ExprGraph])
        num_cons = MathProgBase.numconstr(node_model)
        for i = 1:num_cons
            if !(MathProgBase.isconstrlinear(d,i))    #if it's not a linear constraint
                expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
                _splicevars!(expr,v_map)              #splice the variables from v_map into the expression
                JuMP.addNLconstraint(m,expr)          #raw expression input for non-linear constraint
            end
        end
        #If the objective is linear, store it as a node object
        if MathProgBase.isobjlinear(d)                #get the node model objective.  set it to a node attribute
            t = collect(linearterms(node_model.obj.aff))
            #Make the objective a minimization
            getobjectivesense(node_model) == :Min? sense = 1: sense = -1
            obj = @objective(m,Min, sense*sum{t[i][1]*v_map[linearindex(t[i][2])],i = 1:length(t)})
            getattribute(ne,:JuMPData).objective = obj
        end
        #I don't have the non-linear objective yet, but it shouldn't be any different than the constraints.
    end
end

#splice variables into a constraint expression
function _splicevars!(expr::Expr,v_map::Dict)
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splicevars!(expr.args[i],v_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]   #this is the actual index in x[1], x[2], etc...
                new_var = :($(v_map[var_index]))   #get the JuMP variable from v_map using the index
                expr.args[i] = new_var             #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end

#generate an aggregate model that equates common variable names across models
function aggregate_model!(nodes::Array{Node,1})  #This should be Array{Node}, but it's easy to end up as Any
    anode = aggregate!(nodes)
    graph = anode.subgraph
    m = generatemodel!(graph)
    setmodel!(anode,m)
    return anode
end

#run coupling functions for nodes and edges
function _runcouplings(m::JuMP.Model,graph::Graph)
    for node in getnodes(graph)
        if node.subgraph != nothing
            _runcouplings(m,node.subgraph)
        end
    end
    #coupling arguments are custom arguments to pass to coupling functions
    for (node, func) in graph.node_coupling_functions
        coupling_args = getcouplingargs(graph,node)
        if coupling_args != nothing
            func(m,node,coupling_args...)
        else
            func(m,node)
        end
    end
    for (edge,func) in graph.edge_coupling_functions
        coupling_args = getcouplingargs(graph,edge)
        if coupling_args != nothing
            func(m,edge,coupling_args...)
        else
            func(m,edge)
        end
    end
end
#set individual node model value to the current graph solution
function storecurrentsolution!(graph)
    nodelist = getnodelist(graph)
    edgelist = getedgelist(graph)
    for ne in [nodelist;edgelist]
        if hasattribute(ne,:JuMPData)
            vdict = vardict(ne)
            sdict = Dict()
            for (key,vars) in vdict
                soln = [getvalue(var) for var in vars]
                sdict[key] = soln
            end
            getattribute(ne,:JuMPData).solndict = sdict
        end
    end
end

function setnodevalues!(graph::Graph)
    for ne in [getnodelist(graph);getedgelist(graph)]
        if hasmodel(ne)
            node_model = getmodel(ne)
            solndict = getattribute(ne,:JuMPData).solndict
            for (var_name,values) in solndict
                model_var = getvariable(node_model,var_name)
                if isa(model_var,Array) || isa(model_var,JuMP.JuMPArray)
                    for i = 1:length(model_var)
                        setvalue(model_var[:][i],values[i])
                    end
                else
                    setvalue(model_var,values[1])
                end
            end
        end
    end
end

function generatemodel!(graph::Graph)
    m = Model()
    nodelist = getnodelist(graph)  #all nodes
    edgelist = getedgelist(graph)  #all edges
    #nodelist = getbottomnodes(graph) #just bottom nodes
    #edgelist = getbottomedges(graph) #just bottom edges
    #nodelist = getnodes(graph)
    #edgelist = getedges(graph)
    #build each model
    objs = []
    for edge in edgelist
        buildnodemodel!(edge,m)
        if hasmodel(edge)
            push!(objs,objective(edge))
        end
    end
    for node in nodelist
        buildnodemodel!(node,m)
        if hasmodel(node)
            push!(objs,objective(node))
        end
    end
    _runcouplings(m,graph)
    #println(objs)
    @objective(m,Min,sum(objs))
    addattribute!(graph,:current_model,m)
    return m
end

function generatetoplevelmodel!(graph::Graph)
    #Generate a model based on the top level node models
    m = Model()
    nodelist = getnodes(graph)
    edgelist = getedges(graph)
    #build each model
    objs = []
    for edge in edgelist
        if hasmodel(edge)
            buildnodemodel!(edge,m)
            push!(objs,objective(edge))
        end
    end
    for node in nodelist
        if hasmodel(node)
            buildnodemodel!(node,m)
            push!(objs,objective(node))
            #redirect node variables to the aggregated model
            if node.subgraph != nothing
                for subnode in node.subgraph.nodes
                    for key in keys(vardict(subnode))
                        new_vars = []
                        vars = vardict(subnode)[key]
                        for var in vars
                            new_var = vardict(node)[symbol(var)][1]
                            push!(new_vars,new_var)
                        end
                        vardict(subnode)[key] = new_vars
                    end
                end
                for subedge in node.subgraph.edges
                    if hasmodel(subedge)
                        for key in keys(vardict(subedge))
                            new_vars = []
                            vars = vardict(subedge)[key]
                            for var in vars
                                new_var = vardict(node)[symbol(var)][1]
                                push!(new_vars,new_var)
                            end
                            vardict(subedge)[key] = new_vars
                        end
                    end
                end

            end
        end
    end
    _runcouplings(m,graph)
    @objective(m,Min,sum(objs))
    addattribute!(graph,:current_model,m)
    return m
end


function getmodel(graph::Graph)
    if hasattribute(graph,:current_model)
        return getattribute(graph,:current_model)
    else
        return generatemodel!(graph)
    end
end
