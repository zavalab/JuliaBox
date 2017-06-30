# Plasmo Code
# Yankai Cao, Victor Zavala
# UW-Madison, 2016

module Plasmo
import JuMP
export NetModel, GraphModel, @addNode, getNode, @Linkingconstraint
export getparent, getchildren
export Ipopt_solve
export PipsNlp_solve
using Base.Meta

type NetData
    children::Vector{JuMP.Model}
    parent
    childrenDict::Dict{String, JuMP.Model}
end
NetData() = NetData(JuMP.Model[], nothing, Dict{String, JuMP.Model}())


function NetModel(buildType="serial")
    m = JuMP.Model()
    m.ext[:Net] = NetData()
    m.ext[:BuildType] = buildType
    m.ext[:linkingId] = []
    return m
end

GraphModel() = NetModel()

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

end
include("NetParPipsNlp.jl")
include("NetIpopt.jl")
