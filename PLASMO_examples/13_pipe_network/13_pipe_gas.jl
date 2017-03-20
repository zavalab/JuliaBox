push!(LOAD_PATH, pwd())
include("../src/PlaSMO_debug.jl")
include("../gas_components.jl")

#fix supply pressures
function fix_supply_pressure!(m::JuMP.Model,supply_nodes::Any,spressure::Number)
    for supply in supply_nodes
        junction = neighbors_out(supply)[1]
        pressure = getvariable(junction,:pressure)
        @constraint(m,fix_pressure[t = time_grid],pressure[t] == spressure)
    end
end

function fix_demands!(m::JuMP.Model,demand_nodes,demand_matrix)
    i = 1
    for node in demand_nodes
        stoch_data = demand_matrix[i,:]
        i += 1
        demand = getvariable(node,:fdemand)
        @constraint(m, fix_dem[t = time_grid],demand[t] == stoch_data[t])
    end
end

function create_13_pipe_network(node_data,link_data,supply_data,demand_data,time_grid,x_grid,horizon)
    #Create the top level gas network
    network = Graph()

    #create a dictionary to keep track of nodes {node_name:node_object}
    nodeDict = Dict()

    #Create the individual gas systems
    for row in node_data
        pmin = row["pmin"]
        pmax = row["pmax"]
        pstart = (pmin + pmax)/2
        jdata = JunctionData(time_grid,pmin,pmax,pstart)  #time grid, p_low,p_high,p_start
        node = addnode!(network)
        #addattribute!(node,:jdata,jdata)
        gassystem = GasSystem(jdata)
        nodeDict[row["NODE"]] = gassystem
        addgraph!(node,getgraph(gassystem))  #add a gas junction system to the node
        setcouplingfunction(network,gassystem.junction,couplegasnode!)
    end

    #Create gas supplies and connect them to their corresponding gas nodes
    #supply_locations = []
    supply_nodes = []
    for row in supply_data
        supply_location = row["sloc"]                #get the supply node location
        s_min =  row["smin"]
        s_max = row["smax"]*ffac*ffac2
        gassystem = nodeDict[supply_location]
        sdata = SupplyData(time_grid,0,s_min,s_max)
        snode = addsupply!(gassystem,sdata)
        addattribute!(snode,:sdata,sdata)
        push!(supply_nodes,snode)
    end

    #create the gas demands and connect them
    demand_nodes = []
    for row in demand_data
        demand_location = row["dloc"]
        gassystem = nodeDict[demand_location]
        ddata = DemandData(time_grid,-1000)
        dnode = adddemand!(gassystem,ddata)
        addattribute!(dnode,:ddata,ddata)
        push!(demand_nodes,dnode)
    end

    #create all the links
    for row in link_data
        start_junction = nodeDict[row["lstartloc"]]  #get the starting gas node
        end_junction = nodeDict[row["lendloc"]]      #get the ending gas node
        diameter = row["ldiam"]*dfac
        llength = row["llength"]*lfac
        min_pressure = 0
        max_pressure = 100
        min_flow = 0
        max_flow = 100
        dpmin = 0
        dpmax = 100
        min_power = 0
        max_power = 3000
        comp_cost = 0.1

        pdata = PipeData(llength,diameter,time_grid,x_grid,min_pressure,max_pressure,min_flow,max_flow)
        cdata = CompData(dpmin,dpmax,min_power,max_power,comp_cost)
        link = Edge(network,start_junction.junction,end_junction.junction)  #connect junctions with pipelines
        #add attributes to the link to reference them easily later
        addattribute!(link,:pdata,pdata)
        if row["ltype"] == "a"
            setmodel!(link,ssactivelink(pdata,cdata))
            addattribute!(link,:cdata,cdata)
        elseif row["ltype"] == "p"
            setmodel!(link,sspassivelink(pdata))
        else
            error("link type should be either a (active) or p (passive)")
        end
        setcouplingfunction(network,link,couplelink!)
    end
    return network,supply_nodes,demand_nodes
end
