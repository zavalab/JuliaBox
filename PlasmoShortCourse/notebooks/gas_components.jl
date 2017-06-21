#Gas components file with functions for:
# 1. Building JuMP models for specific gas network components
# 2. Types that facilitate modeling gas systems
# 3. Coupling Functions for junctions and pipes at multiple levels
using JuMP
# physical property data
eps= 0.025		             # pipe rugosity - [mm]
z= 0.80        			     # gas compressibility  - []
rhon=0.72         		     # density of air at normal conditions - [kg/m3]
R=8314.0       			     # universal gas constant [J/kgmol-K]
M=18.0    			         # gas molar mass [kg/kgmol]
pi=3.14         		     # pi
Tgas = 293.15      		     # reference temperature [K]
Cp = 2.34        		     # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85        		     # heat capacity @ constant volume [kJ/kg-K]
U = 1.0*0.1     		     # pipe heat transfer coefficient [J/m2-s-K]
Tamb = 20+273.15   		     # soil temperature [K]
Tsup = 30+273.15   		     # supply temperature [K]

#scaling factors
ffac=(1e+6*rhon)/(24*3600)                     # from scmx10-6/day to kg/s
ffac2=(3600)/(1e+4*rhon)                       # from kg/s to scmx10-4/hr
pfac=1e+5                                      # from bar to Pa
pfac2=1e-5                                     # from Pa to bar
dfac=1e-3                                      # from mm to m
lfac=1e+3                                      # from km to m

#calculated constants
gam = Cp/Cv       		     	# expansion coefficient [-]
nu2 = gam*z*R*Tgas/M  			# gas speed of sound
om = (gam-1.0)/gam 		     	# aux constant [-]
c4 = (1/ffac2)*(Cp*Tgas)

################################
# Gas Node Component
################################
# type JunctionData
#     time_grid
#     pressure_lower
#     pressure_upper
#     pressure_start
# end
#
# function gasjunction(data::JunctionData)
#     m = Model()
#     time_grid = data.time_grid
#     @variable(m,data.pressure_lower <= pressure[time_grid] <= data.pressure_upper, start = data.pressure_start)
#     @variable(m,supply[time_grid] >= 0)  #supply flow to gasnode
#     @variable(m,deliver[time_grid] >= 0) #delivered flow from gas node
#     return m
# end

type NodeData
    time_grid
    n_supplies
    supply_upper
    n_demands
    pressure_lower
    pressure_upper
    demand_cost
    gen_cost
end

function gasnode(data::NodeData)
    m = Model()
    time_grid = data.time_grid
    supplies = collect(1:data.n_supplies)
    demands = collect(1:data.n_demands)
    @variable(m,data.pressure_lower <= pressure[time_grid] <= data.pressure_upper,start = 50)
    @variable(m,0 <= fgen[supplies,time_grid] <= data.supply_upper, start = 10)
    @variable(m,demand[demands,time_grid] >= 0)
    @variable(m,fdeliver[demands,time_grid] >= 0,start = 100)
    @constraint(m, flowLimit[d = demands,t = time_grid], fdeliver[d,t] <= demand[d,t])

    @variable(m,total_delivered[time_grid] >= 0)
    @variable(m,total_supplied[time_grid] >= 0)
    @constraint(m,sum_deliver[t = time_grid],total_delivered[t] == sum(fdeliver[d,t] for d in demands))
    @constraint(m,sum_supply[t = time_grid],total_supplied[t] == sum(fgen[s,t] for s in supplies))

    @variable(m,delivercost[demands]) #each demand has a deliver cost (revenue)
    @variable(m,total_delivercost)
    @variable(m,supplycost[supplies])
    @variable(m,total_supplycost)

    @constraint(m,[d = demands],delivercost[d] == sum(data.demand_cost*fdeliver[d,t] for t = time_grid))
    @constraint(m, integratedDemandCost, total_delivercost == sum(delivercost[d] for d = demands))
    @constraint(m,[s = supplies],supplycost[s] == sum(data.gen_cost*fgen[s,t] for t = time_grid))
    @constraint(m, integratedSupplyCost, total_supplycost == sum(supplycost[s] for s = supplies))

    @objective(m, Min, total_delivercost + total_supplycost)
    return m
end
# ####################################
# # Supply Component
# ####################################
# type SupplyData
#     time_grid
#     cost
#     fgen_lower
#     fgen_upper
# end
#
# function gassupply(data::SupplyData)
#     cost = data.cost
#     time_grid = data.time_grid
#     m = Model()
#     @variable(m, data.fgen_lower <= fgen[time_grid] <= data.fgen_upper, start = 10)
#     @variable(m, gencost[time_grid])
#     @constraint(m,costconstraint[t = time_grid], gencost[t] == cost*fgen[t])
#     return m
# end
# ####################################
# # Demand Component
# ####################################
# type DemandData
#     time_grid
#     cost
# end
#
# function gasdemand(data::DemandData)
#     cost = data.cost
#     time_grid = data.time_grid
#     m = Model()
#     @variable(m,fdeliver[time_grid] >= 0, start = 100)
#     @variable(m, demandcost)
#     @variable(m, fdemand[time_grid] >= 0)
#     @constraint(m, flowLimit[t = time_grid], fdeliver[t] <= fdemand[t])
#     @constraint(m, integratedGasCost, demandcost == sum(cost*fdeliver[t] for t = time_grid))
#     @objective(m, Min, demandcost)
#     return m
# end
#
# ############################################
# Gas pipeline components
############################################
# The macros are useful here for defining base components and then extending them with other macros
macro basepipeeqns(m)
    expr = quote
        @variable(m,pin[time_grid] >= 0, start = 60)
        @variable(m,pout[time_grid] >= 0, start = 60)
        @variable(m,fin[time_grid] >=0, start = 100)
        @variable(m,fout[time_grid] >= 0, start = 100)
        @variable(m, min_pressure <= px[time_grid,x_grid] <= max_pressure, start = 60)#start = 0.5*($data.min_pressure + $data.max_pressure)) #60   # link pressure profile - [bar]
        @variable(m, min_flow <= fx[time_grid,x_grid] <= max_flow, start = 10)#100 #start = 0.5*($data.min_flow + $data.max_flow))
        @constraint(m, flow_in[t = time_grid],  fx[t,1] == fin[t])
        @constraint(m, flow_out[t = time_grid], fx[t,x_grid[end]] == fout[t])
    end
    return esc(expr)
end
#
macro steadystatestart(m)
    expr = quote
        @constraint(m, flow_ss[t = 1, x = x_grid[1:end-1]], (fx[t,x+1] - fx[t,x]) == 0)
        @constraint(m, pres_ss[t = 1, x = x_grid[1:end-1]], -c2*(px[t,x+1] - px[t,x])/dx - slack[t,x] == 0)
    end
    return esc(expr)
end
#
# pressure boundary without a compressor at the inlet
macro passivepressureboundary(m)
    expr = quote
        @constraint(m,press_in[t = time_grid], px[t,1] == pin[t])
        @constraint(m,press_out[t = time_grid], px[t,x_grid[end]] == pout[t])
    end
    return esc(expr)
end
#
# pressure boundary with a compressor at the inlet (there is a dp term)
macro activepressureboundary(m)
    expr = quote
        @constraint(m,press_in[t = time_grid],  px[t,1] == pin[t] + dp[t])
        @constraint(m,press_out[t = time_grid], px[t,x_grid[end]] == pout[t])
    end
    return esc(expr)
end
#
# dynamic transport equations for flow and pressure
macro isothermaleuler(m)
    expr = quote
        @variable(m, slack2[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable
        @variable(m, slack3[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable
        @variable(m, slack[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable for friction loss term
        @NLconstraint(m, slackeq[t = time_grid, x = x_grid],  slack[t,x]*px[t,x] - c3*fx[t,x]*fx[t,x] == 0)
        @NLconstraint(m, slackeq2[t = time_grid, x = x_grid],  slack2[t,x]*px[t,x] - 2*c1*fx[t,x] == 0)
        @NLconstraint(m, slackeq3[t = time_grid, x = x_grid],  slack3[t,x]*px[t,x]*px[t,x] - c1*fx[t,x]*fx[t,x] == 0)
        @constraint(m, press[t = time_grid[1:end-1], x = x_grid[1:end-1]], (px[t+1,x]-px[t,x])/dt + c1*(fx[t+1,x+1]-fx[t+1,x])/dx == 0 )
        @constraint(m, flow[t = time_grid[1:end-1], x = x_grid[1:end-1]], (fx[t+1,x]-fx[t,x])/dt == -slack2[t+1,x]*(fx[t+1,x+1]-fx[t+1,x])/dx +
                                        slack3[t+1,x]*(px[t+1,x+1]-px[t+1,x])/dx -c2*(px[t+1,x+1]-px[t+1,x])/dx - slack[t+1,x])
    end
    return esc(expr)
end

macro weymouthapprox(m)
    expr = quote
        @variable(m, slack[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable for friction loss term
        @NLconstraint(m, slackeq[t = time_grid, x = x_grid],  slack[t,x]*px[t,x] - c3*fx[t,x]*fx[t,x] == 0)
        @constraint(m, press[t = time_grid[1:end-1], x = x_grid[1:end-1]], (px[t+1,x]-px[t,x])/dt + c1*(fx[t+1,x+1]-fx[t+1,x])/dx == 0 )
        @constraint(m, flow[t = time_grid[1:end-1], x = x_grid[1:end-1]], (fx[t+1,x]-fx[t,x])/dt == -c2*(px[t+1,x+1]-px[t+1,x])/dx - slack[t+1,x])
    end
    return esc(expr)
end

#Working on this.  I think I need to go back to the original pipe equations and get rid of the x_grid index.  This approximation
#integrates out the entire time profile.
macro simoneapprox(m)
    expr = quote
        L = len
        xL = x_grid[end]
        @variable(m, slack2[time_grid] >= 0, start = 1)  #auxiliary variable
        @variable(m, slack3[time_grid] >= 0, start = 1)  #auxiliary variable
        @constraint(m,mass[t = time_grid[1:end-1]], (L/2)*(px[t+1,1] - px[t,1])/dt + (px[t+1,xL] - px[t,xL])/dt + c1*(fx[t+1,xL]-fx[t+1,1]) == 0)
        @constraint(m, momentum[t = time_grid[1:end-1]],(L/2)*((fx[t+1,1] + fx[t,1])/dt + (fx[t+1,xL] + fx[t,xL])/dt) ==
        -c2*(px[t+1,xL]-px[t+1,1])-c3*(L/2)*(slack2[t+1] + slack3[t+1]))

        @constraint(m,slack2cons[t = time_grid],slack2[t]*px[t,1] == fx[t,1]*fx[t,1])
        @constraint(m,slack3cons[t = time_grid],slack3[t]*px[t,xL] == fx[t,xL]*fx[t,xL])
        #Constant velocity in pipe
        @constraint(m,constVel[t = time_grid, x = x_grid[1:end-1]], px[t,x+1]*(fx[t,x+1]-fx[t,x])/dx == fx[t,x+1]*(px[t,x+1]-px[t,x])/dx)
    end
    return esc(expr)
end
#Working on this
macro simonesteadystatestart(m)
    expr = quote
        @constraint(m, mass_ss, c1*(fx[1,xL]-fx[1,1]) == 0)
        @constraint(m, momentum_ss,-c2*(px[1,xL]-px[1,1])-c3*(L/2)*(slack2[1] + slack3[1]) == 0)
    end
    return esc(expr)
end
#
macro adiabaticapprox(m)
    expr = quote
        @variable(m, slack[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable for friction loss term
        @NLconstraint(m, slackeq[t = time_grid, x = x_grid],  slack[t,x]*px[t,x] - c3*fx[t,x]*fx[t,x] == 0)
        @constraint(m, mass[t = time_grid[1:end-1], x = x_grid[1:end-1]], (px[t+1,x]-px[t,x])/dt + c1*(fx[t+1,x+1]-fx[t+1,x])/dx == 0 )
        @constraint(m, momentum[t = time_grid[1:end-1], x = x_grid[1:end-1]], 0 == -c2*(px[t+1,x+1]-px[t+1,x])/dx - slack[t+1,x])
    end
    return esc(expr)
end
# steady state transport equations for flow and pressure
macro steadystateflow(m)
    expr = quote
        @variable(m, slack[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable for friction loss term
        @NLconstraint(m, slackeq[t = time_grid, x = x_grid],  slack[t,x]*px[t,x] - c3*fx[t,x]*fx[t,x] == 0)
        @constraint(m, press[t = time_grid,x = x_grid[1:end-1]], c1*(fx[t,x+1]-fx[t,x])/dx == 0 )
        @constraint(m, flow[t = time_grid,x = x_grid[1:end-1]],-c2*(px[t,x+1]-px[t,x])/dx - slack[t,x] == 0)
    end
    return esc(expr)
end
# linepack corresponds to the mass of gas in the piping
macro linepackeqns(m)
    expr = quote
        @variable(m,linepack[time_grid])
        @constraint(m,linepack_def[t = time_grid],linepack[t] == sum(fx[t,x] for x in x_grid)*dx)
        #Periodic terminal constraint
        @constraint(m,linepack_cons, linepack[time_grid[end]] >= linepack[time_grid[1]])
    end
    return esc(expr)
end
# comperssor power and dp
macro compressoreqns(m)
    expr = quote
        @variable(m, dp_min <= dp[time_grid] <= dp_max, start = 10)
        @variable(m, min_power <= pow[time_grid] <= max_power, start = 500)
        @variable(m, powercost)
        @NLconstraint(m, powereqn[t = time_grid], pow[t] == c4*fin[t]*(((pin[t]+dp[t])/pin[t])^om-1))
        @constraint(m,boostcosteqn, powercost == sum(cost*pow[t] for t = time_grid)*dt/3600)
        @objective(m, Min, powercost)
    end
    return esc(expr)
end
#Pipe data struct
type PipeData
    len
    diameter
    time_grid
    x_grid
    min_pressure
    max_pressure
    min_flow
    max_flow
end

#Compressor data struct
type CompData
    dp_min
    dp_max
    min_power
    max_power
    cost
end

macro compressordata(data)
    expr = quote
        dp_min = $data.dp_min
        dp_max = $data.dp_max
        min_power = $data.min_power
        max_power = $data.max_power
        cost = $data.cost
    end
    return esc(expr)
end

# I didn't feel like throwing all of this in the pipe function, so here's a macro that unpacks all of the data
macro unpackpipedata(data)
    expr = quote
        diameter = $data.diameter
        len = $data.len
        x_grid = $data.x_grid
        time_grid = $data.time_grid
        min_pressure = $data.min_pressure
        max_pressure = $data.max_pressure
        min_flow = $data.min_flow
        max_flow = $data.max_flow
        area = (1/4)*pi*diameter*diameter
        lam = (2*log10(3.7*diameter/(eps*dfac)))^(-2)   #pipe friction coefficient
        c1 = (pfac2/ffac2)*(nu2/area)
        c2 = area*(ffac2/pfac2)
        c3 = area*(pfac2/ffac2)*(8*lam*nu2)/(pi*pi*diameter^5)
    end
    return esc(expr)
end

#Functions to build pipeline models using the defined macros
#Build a steady state pipe without no compressor boundary
function sspassivelink(data::PipeData)
    @unpackpipedata(data)
    dx = len / (length(x_grid) - 1)
    m = Model()
    @basepipeeqns(m)
    @steadystateflow(m)
    @passivepressureboundary(m)
    return m
end

#Build a steady state pipe model with a compressor boundary
function ssactivelink(pdata::PipeData,cdata::CompData)
    @unpackpipedata(pdata)
    @compressordata(cdata)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @compressoreqns(m)
    @activepressureboundary(m)
    @steadystateflow(m)
    return m
end

#Build a dynamic pipe model without a compressor boundary
function weymouthpassivelink(data::PipeData)
    @unpackpipedata(data)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @passivepressureboundary(m) #boundary conditions
    @weymouthapprox(m)  #isothermal euler PDE
    #@steadystatestart(m)  #Initial condition
    @linepackeqns(m)  #terminal constraint
    return m
end

function steadystatestart(m::Model,data::PipeData)
    diameter = data.diameter
    len = data.len
    area = (1/4)*pi*diameter*diameter
    x_grid = data.x_grid
    dx = len / (length(x_grid) - 1)
    c2 = area*(ffac2/pfac2)
    fx = m[:fx]
    px = m[:px]
    slack = m[:slack]
    @constraint(m, flow_ss[t = 1, x = x_grid[1:end-1]], (fx[t,x+1] - fx[t,x]) == 0)
    @constraint(m, pres_ss[t = 1, x = x_grid[1:end-1]], -c2*(px[t,x+1] - px[t,x])/dx - slack[t,x] == 0)
    return m
end

#Build a dynamic pipe model with a compressor boundary
function weymouthactivelink(pdata::PipeData,cdata::CompData)
    @unpackpipedata(pdata)
    @compressordata(cdata)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @compressoreqns(m)
    @activepressureboundary(m)
    @weymouthapprox(m)
    #@steadystatestart(m)
    @linepackeqns(m)
    return m
end

function adiabaticpassivelink(data::PipeData)
    @unpackpipedata(data)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @passivepressureboundary(m) #boundary conditions
    @adiabaticapprox(m)  #isothermal euler PDE
    #@steadystatestart(m)  #Initial condition
    @linepackeqns(m)  #terminal constraint
    return m
end

function adiabaticactivelink(pdata::PipeData,cdata::CompData)
    @unpackpipedata(pdata)
    @compressordata(cdata)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @compressoreqns(m)
    @activepressureboundary(m)
    @adiabaticapprox(m)
    #@steadystatestart(m)
    @linepackeqns(m)
    return m
end


function euleractivelink(pdata::PipeData,cdata::CompData)
    @unpackpipedata(pdata)
    @compressordata(cdata)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @compressoreqns(m)
    @activepressureboundary(m)
    @isothermaleuler(m)
    #@steadystatestart(m)
    @linepackeqns(m)
    return m
end

#Build a dynamic pipe model without a compressor boundary
function eulerpassivelink(data::PipeData)
    @unpackpipedata(data)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @passivepressureboundary(m) #boundary conditions
    @isothermaleuler(m)  #isothermal euler PDE
    #@steadystatestart(m)  #Initial condition
    @linepackeqns(m)  #terminal constraint
    return m
end

function simoneactivelink(pdata::PipeData,cdata::CompData)
    @unpackpipedata(pdata)
    @compressordata(cdata)
    dx = len / (length(x_grid) - 1)
    dt = horizon / length(time_grid)
    m = Model()
    @basepipeeqns(m)
    @compressoreqns(m)
    @activepressureboundary(m)
    @simoneapprox(m)
    @simonesteadystatestart(m)
    @linepackeqns(m)
    return m
end
###################################################
# Coupling functions for gas network
###################################################
#coupling functions - gas node
#coupling assumes equal time grid points for now
#couple a gas junction to its supplies and demands
# function couplegasjunction!(m::JuMP.Model,node::Node)
#     #graph = node.graph
#     graph = getgraph(node) #couple in the node's graph context
#     supply = getnodevariable(node,:supply)
#     deliver = getnodevariable(node,:deliver)
#
#     #couple supply and delivery to junction variables for supply and delivery
#     @constraint(m,fsup[t = time_grid], supply[t]  ==  sum(getnodevariable(neighbors_in(graph,node)[i],:fgen)[t] for i = 1:n_neighbors_in(graph,node)))
#     #@constraint(m,fsup[t = time_grid], supply[t]  ==  sum(getnodevardict(neighbors_in(graph,node)[i])[:fgen][t] for i = 1:n_neighbors_in(graph,node)))
#     @constraint(m,fdel[t = time_grid], deliver[t] ==  sum(getnodevariable(neighbors_out(graph,node)[i],:fdeliver)[t] for i = 1:n_neighbors_out(graph,node)))
# end
#
# #couple a gas junction to other junctions
# function couplegasnode!(m::JuMP.Model,node::Node)
#     topgraph = node.graph.parent.graph  #couple in the higher graph context (I will have syntax to do this better soon)
#     #graph = getgraph(node)
#     #topgraph = getparentgraph(graph)
#     supply1 = getnodevariable(node,:supply)
#     deliver1 = getnodevariable(node,:deliver)
#     links_in = edges_in(topgraph,node)  #the links coming in at the higher level graph
#     links_out = edges_out(topgraph,node)
#
#     flow = @constraint(m,[t = time_grid], 0 == sum(getnodevariable(links_in[i],:fout)[t] for i = 1:length(links_in)) -
#     sum(getnodevariable(links_out[i],:fin)[t] for i = 1:length(links_out)) + supply1[t] - deliver1[t])
# end
#
# #couple a link to two junctions
# function couplelink!(m::JuMP.Model,edge::Edge)
#     node_from = getconnectedfrom(edge)
#     node_to = getconnectedto(edge)
#     pin = getnodevariable(edge,:pin)
#     pout = getnodevariable(edge,:pout)
#     @constraint(m,pressure_in[t = time_grid],pin[t] == getnodevariable(node_from,:pressure)[t])
#     @constraint(m,pressure_out[t = time_grid],pout[t] == getnodevariable(node_to,:pressure)[t])
# end

#####################################################################################################################################
# These are convenience types for creating gas network systems.  They use the PLASMO Graph types and coupling functions.
# This shows how the graph can be extended to specific modeling applications
#####################################################################################################################################
# type GasSystem
#     graph::Graph
#     junction::Node
#     supplies::Vector{Node}
#     demands::Vector{Node}
# end
#
# function GasSystem(data::JunctionData)
#     graph = Graph()
#     junction = addnode(graph,gasjunction(data))  #add a new node and set the model to the gasjunction model
#     addattribute!(junction,:jdata,data)
#     setcouplingfunction(graph,junction,couplegasjunction!) #set the coupling function
#     GasSystem(graph,junction,Node[],Node[])
# end
#
# function addsupply!(gassystem::GasSystem,data::SupplyData)
#     graph = gassystem.graph
#     supply = addnode(graph,gassupply(data))
#     addattribute!(supply,:sdata,data)
#     push!(gassystem.supplies,supply)
#     addedge(supply,gassystem.junction)
#     return supply
# end
#
# function adddemand!(gassystem::GasSystem,data::DemandData)
#     graph = gassystem.graph
#     demand = addnode(graph,gasdemand(data))
#     addattribute!(demand,:ddata,data)
#     push!(gassystem.demands,demand)
#     addedge(gassystem.junction,demand)
#     return demand
# end
#
# getsystemgraph(gassystem::GasSystem) = gassystem.graph
