# Plasmo Code
# Yankai Cao, Victor Zavala
# UW-Madison, 2016
# objective cannot be nonlinear	(can be quadratic)

module PlasmoOld
using JuMP
export NetModel, @addNode, getNode, @Linkingconstraint
export getparent, getchildren, getchildrenDict 
export Ipopt_solve
export PipsNlp_solve
export getsumobjectivevalue, RandomStochasticModel, StochasticModel, copyStoModel, copyModel, copyNLModel, extensiveSimplifiedModel, addNLconstraint2, getData, _splicevars! 
using Base.Meta
using MathProgBase      
using Distributions

type NetData
    children::Vector{JuMP.Model}
    parent
    childrenDict::Dict{String, JuMP.Model}
end
NetData() = NetData(JuMP.Model[], nothing, Dict{String, JuMP.Model}())
#NetData() = NetData(nothing, Dict{String, JuMP.Model}())

function NetModel(buildType="serial")
    m = JuMP.Model()
    m.ext[:Net] = NetData()
    m.ext[:BuildType] = buildType
    m.ext[:linkingId] = []
    return m
end

function getNet(m::JuMP.Model)
    if haskey(m.ext, :Net)
        return m.ext[:Net]
    else
        error("This functionality is only available to NetModel")
    end
end


function getData(m::JuMP.Model)
    if haskey(m.ext, :Data)
        return m.ext[:Data]
    else
        error("This functionality is only available")
    end
end

getparent(m::JuMP.Model)       = getNet(m).parent
getchildrenDict(m::JuMP.Model)     = getNet(m).childrenDict

getname(c::Symbol) = c
getname(c::Void) = ()
getname(c::Expr) = c.args[1]

function getchildren(m::JuMP.Model)
    if haskey(m.ext, :Net)
        return getNet(m).children
    else
        return []
    end
end

function getNode(m::JuMP.Model, modelname::String)
    if !haskey(getNet(m).childrenDict, modelname)
        error("No model with name $modelname")
    elseif getNet(m).childrenDict[modelname] === nothing
        error("Multiple models with name $modelname")
    else
        return getNet(m).childrenDict[modelname]
    end
end


function registermodel(m::JuMP.Model, modelname::String, value::JuMP.Model)
    if haskey(getNet(m).childrenDict, modelname)
        getNet(m).childrenDict[modelname] = nothing # indicate duplicate
        error("Multiple models with name $modelname")
    else
        getNet(m).childrenDict[modelname] = value
    end
end

macro Linkingconstraint(m, args...)
      expr = quote
      	  id_start = length($(m).linconstr) + 1
          @constraint($(m),$(args...))
	  id_end = length($(m).linconstr)
	  $(m).ext[:linkingId] = [$(m).ext[:linkingId]; id_start:id_end]   
      end
      return esc(expr)
end

macro addNode(m, node)
    if isa(node, Symbol)
       return quote       	    
    	   if haskey($(esc(node)).ext, :Net)
	        getNet($(esc(node))).parent = $(esc(m))
    	   else
		$(esc(node)).ext[:Net] = NetData(JuMP.Model[], $(esc(m)),Dict{Symbol, JuMP.Model}())
           end
	   push!(getNet($(esc(m))).children, $(esc(node)))
    	   registermodel($(esc(m)), string($(quot(node))), $(esc(node)))
       end
    else
	error("not supported")
    end
end

macro addNode(m, node, nodename)
       return quote
           if haskey($(esc(node)).ext, :Net)
                getNet($(esc(node))).parent = $(esc(m))
           else
                $(esc(node)).ext[:Net] = NetData(JuMP.Model[], $(esc(m)),Dict{Symbol, JuMP.Model}())
           end
           push!(getNet($(esc(m))).children, $(esc(node)))
           registermodel($(esc(m)), $(esc(nodename)), $(esc(node)))
       end
end


# functions from this line are all only for stochastic problems
function getsumobjectivevalue(m::JuMP.Model) 
	 children = getchildren(m)
         leafModelList = children
         modelList = [m; leafModelList]
	 objVal = 0
         for (idx,node) in enumerate(modelList)
	    objVal += node.objVal     
         end
	 return objVal
end


function StochasticModel(createModel, data)
    nscen = length(data)
    master=NetModel()
    node1 = createModel(data[1])
    nfirst = length(node1.ext[:firstVarsId])
    @variable(master, first[1:nfirst])
    for i in 1:nscen
    	node = createModel(data[i])
    	@addNode(master, node, "s$i")
	firstVarsId = node.ext[:firstVarsId]
	for j = 1:nfirst
	    if firstVarsId[j] > 0
	        @constraint(master, first[j] == Variable(node, firstVarsId[j]))
	    end			
	end	    
    end
end

#0.5  2.0
function RandomStochasticModel(createModel, nscen=100, nfirst=5, nparam=5, rdl=0.0, rdu=2.0, adl=-10, adu=10 )
    srand(1234)	 
    master=NetModel()
    @variable(master, first[1:nfirst])
    firstVarsId = collect(1:nfirst)

    for i in 1:nscen
        node = createModel()
        @addNode(master, node, "s$i")
	firstVarsId = 1:nfirst
	node.ext[:firstVarsId] = firstVarsId
        for j = 1:nfirst
            if firstVarsId[j] > 0
                @constraint(master, first[j] == Variable(node, firstVarsId[j]))
            end
        end
	nmodified = 0
	if i == 1
	    continue
	end
	
        if (typeof(nparam) == Int64) && nmodified < floor(nparam/2)
            for c = 1:length(node.quadconstr)
                con = node.quadconstr[c]
                terms = con.terms
                aff = terms.aff
                vars = [terms.qvars1;terms.qvars2; aff.vars]
                varsId = zeros(Int, length(vars))
                for k in 1:length(vars)
                    varsId[k] = vars[k].col
                end
                varsId = sort(union(varsId))
                if length(findin(varsId, firstVarsId))  == length(varsId)
                    continue
                end
                if typeof(nparam) == Int64
                    if nmodified >= nparam
                        break
                    end
                end

            	if  con.sense == :(==)
                    connew = copy(con, node)
		    connew.terms.aff.constant = addnoise(connew.terms.aff.constant, 0, adu, rdl, rdu)
		    connew.sense = :(>=)			    
                    con.terms.aff.constant = addnoise(con.terms.aff.constant, adl, 0, -rdu, -rdl)
                    con.sense = :(<=)
                    push!(node.quadconstr, connew)
                elseif con.sense == :(>=)
                    aff.constant = addnoise(aff.constant, 0, adu, rdl, rdu)
            	elseif con.sense == :(<=)
                    aff.constant = addnoise(aff.constant, adl, 0, -rdu, -rdl)
            	end
		nmodified += 1
            end
        end
	

	for c in 1:length(node.linconstr)
	    con = node.linconstr[c]    	    
	    terms = con.terms
	    vars = terms.vars
	    varsId = zeros(Int, length(vars))
	    for k in 1:length(vars)
	    	varsId[k] = vars[k].col
	    end
	    varsId = sort(union(varsId))
	    if length(findin(varsId, firstVarsId)) == length(varsId)
	        continue
	    end

	    if typeof(nparam) == Int64
                if nmodified >= nparam
                    break
            	end
	    elseif typeof(nparam) == UnitRange{Int64}
	        if !in(c, nparam)
                    continue
                end	
	    end

	    if  con.lb == con.ub
	    	connew = copy(con, node)	
	    	connew.lb = addnoise(connew.lb, adl, 0, -rdu, -rdl)
		connew.ub = Inf
                con.ub = addnoise(con.ub, 0, adu, rdl, rdu)
		con.lb = -Inf
		push!(node.linconstr, connew)
	    elseif con.lb == -Inf
	        con.ub = addnoise(con.ub, 0, adu, rdl, rdu)
	    elseif con.ub  == -Inf
	        con.lb = addnoise(con.lb, adl, 0, -rdu, -rdl)
            end	    
            nmodified += 1
	end

	#=
	if (typeof(nparam) == Int64) && nmodified < nparam
	    for c = 1:length(node.quadconstr)
                con = node.quadconstr[c]
                terms = con.terms
		aff = terms.aff
                vars = [terms.qvars1;terms.qvars2; aff.vars]
                varsId = zeros(Int, length(vars))
                for k in 1:length(vars)
                    varsId[k] = vars[k].col
                end
		varsId = sort(union(varsId))
                if length(findin(varsId, firstVarsId))  == length(varsId)
                    continue
                end
		if typeof(nparam) == Int64
                    if nmodified >= nparam
                        break
                    end
		end
                if  con.sense == :(==)
                    connew = copy(con, node)
                    connew.terms.aff.constant = addnoise(connew.terms.aff.constant, 0, adu, rdl, rdu)
                    connew.sense = :(>=)
                    con.terms.aff.constant = addnoise(con.terms.aff.constant, adl, 0, -rdu, -rdl)
                    con.sense = :(<=)
                    push!(node.quadconstr, connew)
                elseif con.sense == :(>=)
                    aff.constant = addnoise(aff.constant, 0, adu, rdl, rdu)
                elseif con.sense == :(<=)
                    aff.constant = addnoise(aff.constant, adl, 0, -rdu, -rdl)
                end
		#aff.constant = addnoise(aff.constant, adl, adu, rdl, rdu)
		nmodified += 1
	    end
	end
	=#

        for v in 1:node.numCols
            if typeof(nparam) == Int64
                if nmodified >= nparam
                    break
                end
            end
	    if length(findin(v, firstVarsId)) != 0
	        continue
	    end
            node.colLower[v] = addnoise(node.colLower[v], adl, 0, -rdu, -rdl)
            node.colUpper[v] = addnoise(node.colUpper[v], 0, adu, rdl, rdu)
            nmodified += 1
        end

	if (i == 1) && (typeof(nparam) == Int64) && (nmodified < nparam)
	    println("warning: the number of linear/quadratic second stage constraint  ", nmodified, " is less than nparam ",nparam) 	     
	end
    end
    return master
end


function addnoise(a, adl, adu, rdl, rdu)
    if a == 0
        d = Uniform(adl, adu)
        a = a + rand(d)      	 
    else
	d = Uniform(rdl, rdu)
	a = a + abs(a) * rand(d)		 
    end
    return a
end


function Base.copy(v::Array{Variable}, new_model::Model, indent::Int)
    ret = similar(v, Variable, size(v))
    for I in eachindex(v)
        ret[I] = Variable(new_model, v[I].col+indent)
    end
    ret
end

function Base.copy(v::Array{Variable}, new_model::Model, v_map::Array{Int, 1})
    ret = similar(v, Variable, size(v))
    for I in eachindex(v)
        ret[I] = Variable(new_model, v_map[v[I].col])
    end
    ret
end

#=
function extensiveModel(P::JuMP.Model)
        m = Model()
	m.ext[:v_map] = []
        children = getchildren(P)
        nscen = length(children)
        NLobj = false
        if P.nlpdata != nothing
            if P.nlpdata.nlobj != nothing
                NLobj = true
            end
        else
            for scen = 1:nscen
                node = children[scen]
                if node.nlpdata != nothing
                    if node.nlpdata.nlobj != nothing
                        NLobj = true
                        break
                    end
                end
            end
        end
        num_vars = P.numCols
        #add the node model variables to the new model
        for i = 1:num_vars
            x = JuMP.@variable(m)            #create an anonymous variable
            setlowerbound(x, P.colLower[i])
            setupperbound(x, P.colUpper[i])
            var_name = string(Variable(P,i))
            setname(x,var_name)                            #rename the variable to the node model variable name plus the node or edge name
            setcategory(x, P.colCat[i])                            #set the variable to the same category
            setvalue(x, P.colVal[i])                               #set the variable to the same value
        end

        m.obj = copy(P.obj, m)
        m.objSense = P.objSense
        for scen = 1:nscen
            modelname = "s$(scen)"
            node = children[scen]
            nodeobj = addnodemodel!(m, node, modelname)
            m.obj.qvars1 = [m.obj.qvars1; nodeobj.qvars1]
            m.obj.qvars2 = [m.obj.qvars2; nodeobj.qvars2]
            m.obj.qcoeffs = [m.obj.qcoeffs; nodeobj.qcoeffs]
            m.obj.aff.vars = [m.obj.aff.vars; nodeobj.aff.vars]
            m.obj.aff.coeffs = [m.obj.aff.coeffs; nodeobj.aff.coeffs]
            m.obj.aff.constant += nodeobj.aff.constant
        end

        m.linconstr  = [map(c->copy(c, m), P.linconstr); m.linconstr]
        children = getchildren(P)
        for i = 1:length(P.linconstr)
            Pcon = P.linconstr[i]
            Pvars = Pcon.terms.vars
            mcon = m.linconstr[i]
            mvars = mcon.terms.vars
            for (j, Pvar) in enumerate(Pvars)
                if (Pvar.m != P)
		    scen = -1
                    for k in 1:length(children)
                        if Pvar.m == children[k]
                            scen = k
                        end
                    end
                    new_id =  m.ext[:v_map][scen][Pvar.col]
                    mvars[j] = Variable(m, new_id)
                end
            end
        end	
        return m
end
=#

function extensiveSimplifiedModel(P::JuMP.Model)
    ncols_first = P.numCols
    scenarios = PlasmoOld.getchildren(P)
    for (idx,scenario) in enumerate(scenarios)
        scenario.ext[:firstVarsId] = zeros(Int, ncols_first)
        scenario.ext[:firstVarsId][1:end]= -1
    end
    for c in 1:length(P.linconstr)
        coeffs = P.linconstr[c].terms.coeffs
        vars   = P.linconstr[c].terms.vars
        firstVarId = 0
        for (it,ind) in enumerate(coeffs)
            if (vars[it].m) == P
                firstVarId = vars[it].col
                break
            end
        end
        for (it,ind) in enumerate(coeffs)
            if (vars[it].m) != P
               scenario = vars[it].m
               scenario.ext[:firstVarsId][firstVarId] = vars[it].col
            end
        end
    end
    for (idx,scenario) in enumerate(scenarios)
        firstVarsId = scenario.ext[:firstVarsId]
	for i in 1:ncols_first
            if firstVarsId[i] > 0
	        if scenario.colCat[firstVarsId[i]] == :Bin ||  scenario.colCat[firstVarsId[i]] == :Int
		    P.colCat[i] = scenario.colCat[firstVarsId[i]]
		end    
                if scenario.colLower[firstVarsId[i]] >  P.colLower[i]
                    P.colLower[i] = scenario.colLower[firstVarsId[i]]
                end
                if scenario.colUpper[firstVarsId[i]] <  P.colUpper[i]
                    P.colUpper[i] = scenario.colUpper[firstVarsId[i]]
                end
            end
        end
    end

        m = Model()
	m.ext[:v_map] = []
        NLobj = false
        children = getchildren(P)
        nscen = length(children)
        NLobj = false
        if P.nlpdata != nothing
            if P.nlpdata.nlobj != nothing
                NLobj = true
            end
        else
            for scen = 1:nscen
                node = children[scen]
                if node.nlpdata != nothing
                    if node.nlpdata.nlobj != nothing
                        NLobj = true
                        break
                    end
                end
            end
        end
	if NLobj
	    error("sorry, currently not support nonlinear objective function, please formulate it as a constraint!")
	end

        num_vars = P.numCols
        #add the node model variables to the new model
        for i = 1:num_vars
            x = JuMP.@variable(m)            #create an anonymous variable
            setlowerbound(x, P.colLower[i])
            setupperbound(x, P.colUpper[i])
            var_name = string(Variable(P,i))
            setname(x,var_name)				       #rename the variable to the node model variable name plus the node or edge name
            setcategory(x, P.colCat[i])                        #set the variable to the same category
            setvalue(x, P.colVal[i])                           #set the variable to the same value
         end
 
	m.obj = copy(P.obj, m)
        m.objSense = P.objSense

	children = getchildren(P)
	nscen = length(children)
	ncols = Array{Int}(nscen+1)
	nlinconstrs = Array{Int}(nscen+1)
	ncols[1] = m.numCols
	nlinconstrs[1] = length(m.linconstr)
	for scen in 1:nscen
	    node = children[scen]
	    modelname = "s$(scen)"
            nodeobj = addnodeSimplifiedmodel!(m, node, modelname, scen)
            m.obj.qvars1 = [m.obj.qvars1; nodeobj.qvars1]
            m.obj.qvars2 = [m.obj.qvars2; nodeobj.qvars2]
            m.obj.qcoeffs = [m.obj.qcoeffs; nodeobj.qcoeffs]
            m.obj.aff.vars = [m.obj.aff.vars; nodeobj.aff.vars]
            m.obj.aff.coeffs = [m.obj.aff.coeffs; nodeobj.aff.coeffs]
            m.obj.aff.constant += nodeobj.aff.constant
	    ncols[scen + 1] = m.numCols
	    nlinconstrs[scen + 1] = length(m.linconstr)
	end
	m.ext[:ncols] = ncols
	m.ext[:nlinconstrs] = nlinconstrs
	return m    
end


function addnodeSimplifiedmodel!(m::JuMP.Model,node::JuMP.Model, nodename, scen)
        old_numCols = m.numCols
        num_vars = node.numCols
	firstVarsId = node.ext[:firstVarsId] 

        v_map = Array{Int}(num_vars)		#this dict will map linear index of the node model variables to the new model JuMP variables {index => JuMP.Variable}
        #add the node model variables to the new model
        for i = 1:num_vars
	    first = findin(firstVarsId, i);
	    if length(first) == 0  
                x = JuMP.@variable(m)            #create an anonymous variable
                setlowerbound(x, node.colLower[i])
                setupperbound(x, node.colUpper[i])  	    
            	var_name = string(Variable(node,i)) 
            	setname(x,nodename*var_name)                              #rename the variable to the node model variable name plus the node or edge name
		m.colCat[x.col] = node.colCat[i]
            	#setcategory(x, node.colCat[i])                            #set the variable to the same category
            	setvalue(x,node.colVal[i])                                #set the variable to the same value
            	v_map[i] = x.col                                              #map the linear index of the node model variable to the new variable
	    else
		v_map[i] = first[1]
	    end
        end
	push!(m.ext[:v_map], v_map)

        for i = 1:length(node.linconstr)
	    con = copy(node.linconstr[i], node)
	    varsId = extractVarsId(con.terms.vars)
	    # if it is a first stage constraints (a constraints has only only first stage variables)
	    if length(findin(varsId, firstVarsId))  == length(varsId) && scen != 1
	        continue
	    end
	    con.terms.vars = copy(con.terms.vars, m, v_map)
	    m.linconstr = [m.linconstr; con]
        end

        for i = 1:length(node.quadconstr)
            con = copy(node.quadconstr[i], node)
            terms = con.terms
            varsId = extractVarsId([terms.qvars1;terms.qvars2;terms.aff.vars])
            if length(findin(varsId, firstVarsId))  == length(varsId) && scen	!= 1
                continue
            end
            terms.qvars1 = copy(terms.qvars1, m, v_map)
            terms.qvars2 = copy(terms.qvars2, m, v_map)
            terms.aff.vars = copy(terms.aff.vars, m, v_map)
            m.quadconstr = [m.quadconstr; con]
        end

        #Copy the non-linear constraints to the new model
        d = JuMP.NLPEvaluator(node)         #Get the NLP evaluator object.  Initialize the expression graph
        MathProgBase.initialize(d,[:ExprGraph])
        num_cons = MathProgBase.numconstr(node)
        for i = (1+length(node.linconstr)+length(node.quadconstr)):num_cons
            if !(MathProgBase.isconstrlinear(d,i))    #if it's not a linear constraint
                expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
            	varsId = extractVarsId(expr)
            	if length(findin(varsId, firstVarsId))  == length(varsId) && scen != 1
                    continue
                end
                _splicevars!(expr, m, v_map)              #splice the variables from v_map into the expression
                addNLconstraint2(m,expr)          	  #raw expression input for non-linear constraint
            end
        end
        nodeobj = copy(node.obj)
        nodeobj.qvars1 = copy(nodeobj.qvars1, m, v_map)
        nodeobj.qvars2 = copy(nodeobj.qvars2, m, v_map)
        nodeobj.aff.vars = copy(nodeobj.aff.vars, m, v_map)
        return nodeobj
end



#=
function addnodemodel!(m::JuMP.Model,node::JuMP.Model, nodename)
        old_numCols = m.numCols
        num_vars = node.numCols
        v_map = Array{Int}(num_vars)               #this dict will map linear index of the node model variables to the new model JuMP variables {index => JuMP.Variable}
        #add the node model variables to the new model
        for i = 1:num_vars
            x = JuMP.@variable(m)            #create an anonymous variable
            setlowerbound(x, node.colLower[i])
            setupperbound(x, node.colUpper[i])
            var_name = string(Variable(node,i))
            setname(x,nodename*var_name)                             #rename the variable to the node model variable name plus the node or edge name
            setcategory(x, node.colCat[i])                            #set the variable to the same category
            setvalue(x,node.colVal[i])                               #set the variable to the same value
            v_map[i] = x.col                                         #map the linear index of the node model variable to the new variable
        end
	push!(m.ext[:v_map], v_map)

        for i = 1:length(node.linconstr)
	    con	= copy(node.linconstr[i], node)
            con.terms.vars = copy(con.terms.vars, m, v_map)
            m.linconstr	   = [m.linconstr; con]
        end

	for i = 1:length(node.quadconstr)
            con = copy(node.quadconstr[i], node)
	    terms = con.terms
            terms.qvars1 = copy(terms.qvars1, m, old_numCols)
            terms.qvars2 = copy(terms.qvars2, m, old_numCols)
            terms.aff.vars = copy(terms.aff.vars, m, old_numCols)
	    m.quadconstr = [m.quadconstr; con]
	end

        #Copy the non-linear constraints to the new model
        d = JuMP.NLPEvaluator(node)         #Get the NLP evaluator object.  Initialize the expression graph
        MathProgBase.initialize(d,[:ExprGraph])
        num_cons = MathProgBase.numconstr(node)
        for i = (1+length(node.linconstr)+length(node.quadconstr)):num_cons
            if !(MathProgBase.isconstrlinear(d,i))    #if it's not a linear constraint
                expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
                _splicevars!(expr,m, v_map)              #splice the variables from v_map into the expression
                addNLconstraint2(m,expr)          #raw expression input for non-linear constraint
            end
        end
        nodeobj = copy(node.obj)
        nodeobj.qvars1 = copy(nodeobj.qvars1, m, old_numCols)
        nodeobj.qvars2 = copy(nodeobj.qvars2, m, old_numCols)
        nodeobj.aff.vars = copy(nodeobj.aff.vars, m, old_numCols)
        return nodeobj
end
=#


function copyModel(P::JuMP.Model)
    m = Model()	
    # Variables
    m.numCols = P.numCols
    m.colNames = P.colNames[:]
    m.colNamesIJulia = P.colNamesIJulia[:]
    m.colLower = P.colLower[:]
    m.colUpper = P.colUpper[:]
    m.colCat = P.colCat[:]
    m.colVal = P.colVal[:]
    m.linconstrDuals = P.linconstrDuals[:]
    m.redCosts = P.redCosts[:]
    # Constraints
    m.linconstr  = map(c->copy(c, m), P.linconstr)
    m.quadconstr = map(c->copy(c, m), P.quadconstr)

    
    m.objDict = Dict{Symbol,Any}()
    m.varData = ObjectIdDict()
    for (symb,o) in P.objDict
        newo = copy(o, m)
        m.objDict[symb] = newo
        if haskey(P.varData, o)
            m.varData[newo] = P.varData[o]
            #m.varData[newvar] = copy(P.varData[v]) # should we copy this too ? We need to define copy(::JuMPContainerData) too then
        end
    end
    # Objective
    m.obj = copy(P.obj, m)
    m.objSense = P.objSense
    #Copy the non-linear constraints to the new model
    d = JuMP.NLPEvaluator(P)     #Get the NLP evaluator object.  Initialize the expression graph
    MathProgBase.initialize(d,[:ExprGraph])
    num_cons = MathProgBase.numconstr(P)
    for i = (1+length(P.linconstr)+length(P.quadconstr)):num_cons
            expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
	    _splicevars!(expr, m)
            addNLconstraint2(m,expr)          #raw expression input for non-linear constraint
    end
    if  !MathProgBase.isobjlinear(d) && !MathProgBase.isobjquadratic(d)   #P.nlpdata.nlobj != nothing 
        expr = MathProgBase.obj_expr(d)
        _splicevars!(expr, m)
        JuMP.setNLobjective(m, P.objSense, expr)
    end
    
    if !isempty(P.ext)
        m.ext = similar(P.ext)
        for (key, val) in P.ext
            m.ext[key] = try
                copy(P.ext[key])
            catch
                continue;  #error("Error copying extension dictionary. Is `copy` defined for all your user types?")
            end
        end
    end
    return m
end


function copyNLModel(P::JuMP.Model)  # copy model and convert quadratic consriants and obj to nonlinear
    m = Model()
    # Variables
    m.numCols = P.numCols
    m.colNames = P.colNames[:]
    m.colNamesIJulia = P.colNamesIJulia[:]
    m.colLower = P.colLower[:]
    m.colUpper = P.colUpper[:]
    m.colCat = P.colCat[:]
    m.colVal = P.colVal[:]
    m.linconstrDuals = P.linconstrDuals[:]
    m.redCosts = P.redCosts[:]
    # Constraints
    m.linconstr  = map(c->copy(c, m), P.linconstr)
    #m.quadconstr = map(c->copy(c, m), P.quadconstr)

    m.objDict = Dict{Symbol,Any}()
    m.varData = ObjectIdDict()
    for (symb,o) in P.objDict
        newo = copy(o, m)
        m.objDict[symb] = newo
        if haskey(P.varData, o)
            m.varData[newo] = P.varData[o]
        end
    end

    # Objective
    m.obj = copy(P.obj, m)
    m.objSense = P.objSense
    #Copy the non-linear constraints to the new model
    d = JuMP.NLPEvaluator(P)     #Get the NLP evaluator object.  Initialize the expression graph
    MathProgBase.initialize(d,[:ExprGraph])
    num_cons = MathProgBase.numconstr(P)
    for i = (1+length(P.linconstr)):num_cons
            expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
	    _splicevars!(expr, m) 
            addNLconstraint2(m,expr)          #raw expression input for non-linear constraint
    end
    if !MathProgBase.isobjlinear(d)  #isa(P.nlpdata.nlobj, NonlinearExprData) || MathProgBase.isobjquadratic(d)  #P.nlpdata.nlobj != nothing 
        expr = MathProgBase.obj_expr(d)
	_splicevars!(expr, m)
        JuMP.setNLobjective(m, P.objSense, expr)
    end
    if !isempty(P.ext)
        m.ext = similar(P.ext)
        for (key, val) in P.ext
            m.ext[key] = try
                copy(P.ext[key])
            catch
                continue;  #error("Error copying extension dictionary. Is `copy` defined for all your user types?")
            end
        end
    end
    return m
end


function copyStoModel(P::JuMP.Model)
    if haskey(P.ext, :Net)	 
        m = NetModel() 
	children = getchildren(P)
	nscen = length(children)
    	for scen = 1:nscen 
	    modelname = "s$(scen)" 
    	    nodecopy = copyModel(children[scen])    	
	    @addNode(m, nodecopy, modelname)   
        end
    	m.numCols = P.numCols
    	m.colNames = P.colNames[:]
    	m.colNamesIJulia = P.colNamesIJulia[:]
    	m.colLower = P.colLower[:]
    	m.colUpper = P.colUpper[:]
    	m.colCat = P.colCat[:]
    	m.colVal = P.colVal[:]
	m.linconstr  = map(c->copy(c, m), P.linconstr)
	for i = 1:length(P.linconstr)
	    Pcon = P.linconstr[i]	
	    Pvars = Pcon.terms.vars
	    mcon = m.linconstr[i]
	    mvars = mcon.terms.vars 
            for (j, Pvar) in enumerate(Pvars)
                if (Pvar.m != P)
                    scen = -1
                    for k in 1:length(children)
                        if Pvar.m == children[k]
                            scen = k
                        end
                    end
		    mvars[j] = Variable(getchildren(m)[scen], Pvar.col)
                end
            end
        end
        m.obj = copy(P.obj, m)
    	m.objSense = P.objSense
    	#Copy the non-linear objective
	if P.nlpdata != nothing
    	if P.nlpdata.nlobj != nothing  # !isobjquadratic(d)&& ! isobjlinear(d)
	    removed = removeConnection(P)
            d = JuMP.NLPEvaluator(P)     #Get the NLP evaluator object.  Initialize the expression graph
            MathProgBase.initialize(d,[:ExprGraph])
            expr = MathProgBase.obj_expr(d)
            JuMP.setNLobjective(m, P.objSense, expr)
	    P.linconstr = [P.linconstr; removed]
    	end
	end
    else
	m = copy(P)	
    end
    return m
end


#splice variables into a constraint expression
function _splicevars!(expr::Expr, m::JuMP.Model, v_map=nothing)
if v_map !=nothing
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splicevars!(expr.args[i], m, v_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]   #this is the actual index in x[1], x[2], etc...
                new_var = :($(Variable(m, v_map[var_index])))   #get the JuMP variable from v_map using the index
                expr.args[i] = new_var             #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
else
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splicevars!(expr.args[i], m, v_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]   #this is the actual index in x[1], x[2], etc...
                new_var = :($(Variable(m, var_index)))   #get the JuMP variable from v_map using the index
                expr.args[i] = new_var             #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end
end


function extractVarsId(expr::Expr)
    varsId = Int[]	 
    if expr.head == :ref 
       var_index = expr.args[2]   #this is the actual index in x[1], x[2], etc...
       return [var_index]
    end
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            varsId = [varsId; extractVarsId(expr.args[i])]
        end
    end
    return union(varsId) 
end

function extractVarsId(vars::Array{Variable,1})
    varsId = Int[]
    for i = 1:length(vars)
    	push!(varsId, vars[i].col)
    end
    return union(varsId)
end




function removeConnection(master::JuMP.Model)
    deleterow = []
    for row in 1:length(master.linconstr)
        coeffs = master.linconstr[row].terms.coeffs
        vars   = master.linconstr[row].terms.vars
        for (it,ind) in enumerate(coeffs)
            if (vars[it].m) != master
                push!(deleterow, row)
                break
            end
        end
    end
    deleted = master.linconstr[deleterow]
    deleteat!(master.linconstr,deleterow)
    return deleted
end


# args[3] is sure 0 #a number (prabably 0.0)
# Ex: addNLconstraint(m, :($x + $y^2 <= 1))
function addNLconstraint2(m::Model, ex::Expr)
     @assert ex.args[3] == 0.0 	 
     @assert JuMP.isexpr(ex, :call) 
     #=
     if (ex.args[3] != 0.0) || !JuMP.isexpr(ex, :call)
     	 return addNLconstraint(m, ex)
     end
     =#
    JuMP.initNLP(m)
    m.internalModelLoaded = false
        # Simple comparison - move everything to the LHS
        op = ex.args[1]
        if op == :(==)
            lb = 0.0
            ub = 0.0
        elseif op == :(<=) || op == :(≤)
            lb = -Inf
            ub = 0.0
        elseif op == :(>=) || op == :(≥)
            lb = 0.0
            ub = Inf
        else
            error("in addNLconstraint ($ex): expected comparison operator (<=, >=, or ==).")
        end
        lhs = :($(ex.args[2]))
        c = JuMP.NonlinearConstraint(JuMP.NonlinearExprData(m, lhs), lb, ub)
        push!(m.nlpdata.nlconstr, c)
        return ConstraintRef{Model,NonlinearConstraint}(m, length(m.nlpdata.nlconstr))
end


end
#include("NetParPipsNlp.jl")
include("NetIpopt.jl")
#end



