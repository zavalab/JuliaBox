module Plasmo
import JuMP
export NetModel, @addNode, getNode
export getparent, getchildren 
#getchildrenDict,  organizeData
export Ipopt_solve
export PipsNlp_solve
using Base.Meta

type NetData
    children::Vector{JuMP.Model}
    parent
    childrenDict::Dict{ASCIIString, JuMP.Model}
end
NetData() = NetData(JuMP.Model[], nothing, Dict{ASCIIString, JuMP.Model}())


function NetModel()
    m = JuMP.Model()
    m.ext[:Net] = NetData()
    return m
end

function getNet(m::JuMP.Model)
    if haskey(m.ext, :Net)
        return m.ext[:Net]
    else
        error("This functionality is only available")
    end
end

getparent(m::JuMP.Model)       = getNet(m).parent
#getchildren(m::JuMP.Model)     = getNet(m).children
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

function getNode(m::JuMP.Model, modelname::ASCIIString)
    if !haskey(getNet(m).childrenDict, modelname)
        error("No model with name $modelname")
    elseif getNet(m).childrenDict[modelname] === nothing
        error("Multiple models with name $modelname")
    else
        return getNet(m).childrenDict[modelname]
    end
end

function registermodel(m::JuMP.Model, modelname::ASCIIString, value::JuMP.Model)
    if haskey(getNet(m).childrenDict, modelname)
        getNet(m).childrenDict[modelname] = nothing # indicate duplicate
        error("Multiple models with name $modelname")
    else
        getNet(m).childrenDict[modelname] = value
    end
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

function organizeData(raw_data)
        relation=Array(ASCIIString, length(raw_data), 2)
        for i in 1:length(raw_data)
            relation[i,1:2] = raw_data[i][1:2]
        end

        master_index = find(relation[1:end,2].=="")
        assert(length(master_index)==1)
        master_index=master_index[1]
        function getData(parent_index)
                 data = []
                 local_data = raw_data[parent_index][3:end]
                 children_data = []
                 parent_name = relation[parent_index,1]
                 children_index = find(relation[1:end,2].==parent_name)
                 for (idx,child_index) in enumerate(children_index)
                     child_data = getData(child_index)
                     push!(children_data, child_data)
                 end

                 if (length(children_index) !=0)
                    push!(data, local_data)
                    push!(data, children_data)
                    push!(data, raw_data[parent_index][1])
                    return data
                 else
                    push!(local_data, raw_data[parent_index][1])
                    return local_data
                 end
        end
        return getData(master_index)
end

end

include("NetParPipsNlp.jl")
include("NetClusterIPM.jl")
include("NetIpopt.jl")