using JuMP
using Gurobi
using DelimitedFiles
using CSV
using DataFrames
using Statistics

function National_Upcycling_MIP(path, budget, time_limit)

######################################################################
# PATH DEFINITION
######################################################################

# Define a file path to folder containing all relevant csv files
cd(path)

######################################################################
# FUNCTION DEFINITIONS
######################################################################

# Define a function to initialize dicts
function DictInit(OrderedMeyList,InitValue)
    """
    Initializes a dictionary with the keys in OrderedMeyList and assigns each
    key the value in InitValue
    Inputs:
        OrderedMeyList - a list of lists (i.e., list of lists of keys)
        InitValue - probably either 0 or false
    Outputs:
        Out - a dictionary
    Note: replaces the initialization loop for dictionaries that require full
    population; i.e., in the case of set intersections
    """
    Out = Dict()
    if length(OrderedMeyList) == 1
        for key in OrderedMeyList[1]
            Out[key] = InitValue
        end
    else
        SplatMeys = collect(Iterators.product(OrderedMeyList...))
        for key in SplatMeys
            Out[key] = InitValue
        end
    end
    return Out
end

# Define a function to calculate arc length
function GreatCircle(ref_lon, ref_lat, dest_lon, dest_lat)
    """
    Calculates the great circle distance between a reference location (ref)
    and a destination (dest) and returns the great circle distances in km
    """
    # Haversine formula
    dist = 2*R*asin(sqrt(((sind((dest_lat - ref_lat)/2))^2) +
        cosd(dest_lat)*cosd(ref_lat)*((sind((dest_lon - ref_lon)/2))^2)))
    return(dist)
end

# Define an output formatting function
function FilePrint(Variable,OrderedIndexList,filename;Header="*"^50,DataName="",VarName="")
    """
    Generates a list of strings and prints them to file with nice formatting;
    reduces required script code clutter.
    > Variable: the result of a JuMP getvalue() call; a Dict().
    > OrderedIndexList: a list of indices in the same order as the indices of
      the data in Variale; each element is a list of index elements. As an
      example: [A,B,C] where A = [a1,a2,...,aN], B = ....
    > filename: the file name for printing
    > Header: a string to be used as a header above the printed data
    > DataName: a header to appear above the printed data
    > VarName: the desired output name of the variable on screen; a string
    """
    # Header and DataName:
    print(filename,"\n"*Header*"\n"*DataName)

    # Collect indices via splatting to create all permutations; print each permuted index and value to file
    if length(OrderedIndexList) == 1
        SplatIndex = OrderedIndexList[1]
        for index in SplatIndex
            print(filename,"\n"*VarName*"(\""*string(index)*"\") = "*string(Variable[index]))
        end
    else
        SplatIndex = collect(Iterators.product(OrderedIndexList...))
        for index in SplatIndex
            print(filename,"\n"*VarName*string(index)*" = "*string(Variable[[i for i in index]...]))
        end
    end
end

######################################################################
# DATA READING AND PREPARATION
######################################################################

# Import, define and pretreat data
arc_matrix = CSV.read(pwd()*"/arc_matrix.csv", DataFrame)
node_matrix = CSV.read(pwd()*"/node_matrix.csv", DataFrame)

prod_matrix = CSV.read(pwd()*"/product_matrix.csv", DataFrame)
demand_matrix = CSV.read(pwd()*"/demand_matrix.csv", DataFrame)
supply_matrix = CSV.read(pwd()*"/supply_matrix.csv", DataFrame)
technology_matrix = CSV.read(pwd()*"/technology_matrix.csv", DataFrame)

alpha_matrix = readdlm(pwd()*"/alpha_matrix.csv",',');
site_matrix = CSV.read(pwd()*"/site_matrix.csv", DataFrame)

# Set a constant used for distance calculation
R = 6335.439

# Nodes, products, customers, suppliers and technologies
NODES = node_matrix[:,1]; # all nodes
PRODS = prod_matrix[:,1] # all products
DEMS  = demand_matrix[:,1] # all demands
SUPS  = supply_matrix[:,1] # all supply
TECHS = technology_matrix[:,1] # all technologies
ARCS  = arc_matrix[:,1] # all arcs
SITES = site_matrix[:,1] # all technology providers

# Node properties
node_alia = Dict(zip(NODES, node_matrix[:,2])); # node alias
node_lat = Dict(zip(NODES, node_matrix[:,3])); # node longitude
node_long  = Dict(zip(NODES, node_matrix[:,4])); # node latitude

# Product properties
prod_name = Dict(zip(PRODS, prod_matrix[:,2])); # product names
prod_trans_vc_truck= Dict(zip(PRODS, prod_matrix[:,3])); # product transportation bids
prod_trans_fc_truck= Dict(zip(PRODS, prod_matrix[:,4]));

# Demand properties
dem_node  = Dict(zip(DEMS, demand_matrix[:,2])); # demand locations
dem_prod  = Dict(zip(DEMS, demand_matrix[:,3])); # demand products
dem_bid   = Dict(zip(DEMS, demand_matrix[:,4])); # demand bids
dem_cap   = Dict(zip(DEMS, demand_matrix[:,5])); # demand capacities

# Supplier properties
sup_node  = Dict(zip(SUPS, supply_matrix[:,2])); # supply locations
sup_prod  = Dict(zip(SUPS, supply_matrix[:,3])); # supply products
sup_bid   = Dict(zip(SUPS, supply_matrix[:,4])); # supply bids
sup_cap   = Dict(zip(SUPS, supply_matrix[:,5])); # supply capacities

# technology properties
tech_cap  = Dict(zip(TECHS, technology_matrix[:,2]));
tech_inv  = Dict(zip(TECHS, technology_matrix[:,4])); # technology investment costs (not useful here)
tech_bid  = Dict(zip(TECHS, technology_matrix[:,5])); # technology operational costs (bids) per unit reference prod
tech_refprod = Dict(zip(TECHS, technology_matrix[:,3])); # technology reference products

# Technology properties
#   - MODIFICATION: tech_size_3 and tech_cost_3 are eliminated, tech_cost_1 and tech_cost_2 now reference different columns
tech_cap  = Dict(zip(TECHS, technology_matrix[:,2])); # technology capacities with regard to the reference product

tech_size_1  = Dict(zip(TECHS, technology_matrix[:,7]));
tech_size_2  = Dict(zip(TECHS, technology_matrix[:,8]));
# tech_size_3  = Dict(zip(TECHS, technology_matrix[:,9]));

tech_cost_1  = Dict(zip(TECHS, technology_matrix[:,9]));
tech_cost_2  = Dict(zip(TECHS, technology_matrix[:,10]));
# tech_cost_3  = Dict(zip(TECHS, technology_matrix[:,12]));

tp_site = Dict(zip(SITES, site_matrix[:,2])); # node location of the technology provider
tp_tech = Dict(zip(SITES, site_matrix[:,3])); # technology type that the provider can provide
tp_indicator = Dict(zip(SITES, ones(length(SITES)))); # technology indicator



# Define allowed combinations of nodes and tech types (dictated by sites)
site_tuples = Vector{Tuple{Int,String}}(undef, length(SITES))
for s in SITES
    site_tuples[s] = (site_matrix[s,2], site_matrix[s,3])
end

# Define vector of tuple combinations (nodes & tech types) that do not correspond to allowed sites
NONSITES = Vector{Tuple{Int,String}}()
for n in NODES
    for t in TECHS
        if !((n,t) in site_tuples)
            push!(NONSITES, (n,t))
        end
    end
end

# Define a two-key dictionary for transformation factors
transfer = Dict((TECHS[1],PRODS[1]) => 0.5);
for i in 1:length(TECHS)
    for k in 1: length(PRODS)
        transfer[(TECHS[i], PRODS[k])] = alpha_matrix[i,k];
    end
end



######################################################################

# Define a subset of arcs in forward and reverse directions
Ain = Dict() # given node s(n), provides all arcs a in A directed towards node s
Aout = Dict() # given node s(n), provides all arcs a in A directed out of node s
[Ain[n] = Vector{Int}(undef,0) for n in NODES]
[Aout[n] = Vector{Int}(undef,0) for n in NODES]
for a = 1:length(ARCS)
    push!(Ain[arc_matrix[a,3]], ARCS[a])
    push!(Aout[arc_matrix[a,2]], ARCS[a])
end

# Define distances traveled by arcs between nodes (using the Haversine formula)
distance = Dict();
for a in ARCS
    lat_in = node_matrix[arc_matrix[a,3], 4]
    lat_out = node_matrix[arc_matrix[a,2], 4]
    lon_in = node_matrix[arc_matrix[a,3], 3]
    lon_out = node_matrix[arc_matrix[a,2], 3]
    distance[a] = GreatCircle(lon_in, lat_in, lon_out, lat_out)
end

######################################################################
# VARIABLE AND CONSTRAINT DEFINITIONS
######################################################################

# Define the "Parameter for Big_M formulation"
M=100000000

# Define the annualized factor (converts total capital cost to annual capital cost)
af=0.117

# Print a message as solving preparation begins
println("Solving...")

# m = Model(with_optimizer(Gurobi.Optimizer))
m = JuMP.Model(Gurobi.Optimizer)

@variable(m, f[a in ARCS, pr in PRODS]>= 0);

# Define demand and supply variables
@variable(m, dem[DEMS] >= 0);
@variable(m, d[NODES,PRODS] >= 0);
@variable(m, sup[SUPS] >= 0);
@variable(m, s[NODES,PRODS] >= 0);
@variable(m, t[NODES] >= 0);

# Define generation & consumption variables by technologies
@variable(m, x[NODES,PRODS,TECHS]);
@variable(m, p[NODES, PRODS]);

@constraint(m, techflow[i in NODES, t in TECHS], x[i,tech_refprod[t],t] <= 0);

# Define cost variables 
@variable(m,transcost);
@variable(m,opcost);
@variable(m,capcost);
@variable(m,demrevn);
@variable(m,supcost);

# Define demand and supply contraints
@constraint(m, demeq[n in NODES, pr in PRODS], d[n,pr] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr
                && dem_node[dd]==n));
@constraint(m, supeq[n in NODES, pr in PRODS], s[n,pr] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr
                && sup_node[ss]==n));


# Define balance constraint - rewritten to accomodate a subset of arcs
@constraint(m, balance[i in NODES, pr in PRODS], s[i,pr] + p[i,pr] + sum(f[a,pr] for a in Ain[i]) == d[i,pr] + sum(f[a,pr] for a in Aout[i]))

# Define conversion constraints
@constraint(m, process[i in NODES, pr in PRODS], p[i,pr] == sum(x[i,pr,t] for t in TECHS));
@constraint(m, transfer_pr[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] ==
                                    transfer[t,pr]/transfer[t,tech_refprod[t]]*x[i,tech_refprod[t],t]);
              

# Define demand capacity constraints
@constraint(m, demand_capacity[i in DEMS], dem[i] <= dem_cap[i]);

# Define supply capacity constraints
@constraint(m, supply_capacity[i in SUPS], sup[i] <= sup_cap[i]);

# Define z as a binary variable
#   - MODIFICATION: 1:3 -> 1:2
@variable(m, z[SITES, 1:2], Bin, start = 0)

# Define constraint: consumption at allowed techsites must be less than the sum of capacities for technologies placed there.
@constraint(m, bin_tech_cap[s in SITES], - x[tp_site[s],tech_refprod[tp_tech[s]],tp_tech[s]] <= z[s,1]*tech_size_1[tp_tech[s]]
                                                                                              + z[s,2]*tech_size_2[tp_tech[s]]);

# Define added constraint: consumption at nodes with no allowed sites must be zero.
@constraint(m, [ns in NONSITES], x[ns[1], tech_refprod[convert(String3, ns[2])], convert(String3, ns[2])] == 0)



# Define capital cost constraint
@constraint(m, capcost == sum(z[s,1]*tech_cost_1[tp_tech[s]] for s in SITES)
                        + sum(z[s,2]*tech_cost_2[tp_tech[s]] for s in SITES));

# MODIFICATION: Define budget
@constraint(m, capcost <= budget);

# Define expressions used in the objective
opcost = @expression(m, -sum(x[i,tech_refprod[t],t]*tech_bid[t] for i in NODES for t in TECHS));
transcost = @expression(m, sum(prod_trans_vc_truck[pr]*distance[a]*f[a,pr] + f[a,pr]*prod_trans_fc_truck[pr] for a in ARCS for pr in PRODS));
demrevn = @expression(m, sum(dem[i]*dem_bid[i] for i in DEMS));
supcost = @expression(m, sum(sup[i]*sup_bid[i] for i in SUPS));

# Set time limit
set_time_limit_sec(m, time_limit)

# Define social welfare objective
@objective(m, Max, demrevn - supcost - opcost - transcost - af*capcost)
optimize!(m)

######################################################################
# SOLUTION OUTPUT
######################################################################

# Output JuMP variables
z_out = JuMP.value.(z)
f_out = JuMP.value.(f)
p_out = JuMP.value.(p)
drevn_out = JuMP.value.(demrevn)
scost_out = JuMP.value.(supcost)
opcost_out = JuMP.value.(opcost)
trcost_out = JuMP.value.(transcost)
afcapcost_out = af * JuMP.value.(capcost)
swf_out = drevn_out - scost_out - opcost_out - trcost_out - afcapcost_out
d_out = JuMP.value.(d)
dem_out = JuMP.value.(dem)
s_out = JuMP.value.(s)
sup_out = JuMP.value.(sup)

# Format output variables for printing to a text file
solution = open(pwd()*"/VMain_Solution_File.txt", "w")
    FilePrint(swf_out, [], solution; Header="*"^50, DataName="Social Welfare", VarName="swf")
    FilePrint(drevn_out, [], solution; Header="*"^50, DataName="Demand Revenue", VarName="drevn")
    FilePrint(scost_out, [], solution; Header="*"^50, DataName="Supply Cost", VarName="scost") 
    FilePrint(opcost_out, [], solution; Header="*"^50, DataName="Operating Cost", VarName="opcost")
    FilePrint(trcost_out, [], solution; Header="*"^50, DataName="Transportation Cost", VarName="trcost") 
    FilePrint(afcapcost_out, [], solution; Header="*"^50, DataName="Annualized Capital Cost", VarName="afcapcost") 
    FilePrint(z_out, [SITES,1:2], solution; Header="*"^50, DataName="Technology Placement: (Tech Site, Scale)", VarName="z")
    FilePrint(f_out, [ARCS,PRODS], solution; Header="*"^50, DataName="Product Flows: (Arc, Product)", VarName="f")
    FilePrint(p_out, [NODES,PRODS], solution; Header="*"^50, DataName="Technology Flows: (Node, Product)", VarName="p")
    FilePrint(d_out, [NODES,PRODS], solution; Header="*"^50, DataName="Demand: (Node, Product)", VarName="d")
    FilePrint(dem_out, [DEMS], solution; Header="*"^50, DataName="Demand: (ID)", VarName="dem")
    FilePrint(s_out, [NODES,PRODS], solution; Header="*"^50, DataName="Supply: (Node, Product)", VarName="s")
    FilePrint(sup_out, [SUPS], solution; Header="*"^50, DataName="Supply: (ID)", VarName="sup")

close(solution)





##################################################################
# Interpretation

# Read in relevant matrix.csv data
##################################################################

# Define a working directory path
solution_file = "/VMain_Solution_File.txt"

# Read in and organize site data
site_matrix = CSV.read(pwd()*"/site_matrix.csv", DataFrame)
site_ids = site_matrix[:,1]
site_nodes = site_matrix[:,2]
site_techs = site_matrix[:,3]
total_sites = length(site_ids)

# Read in and organize node data
node_matrix = CSV.read(pwd()*"/node_matrix.csv", DataFrame)
node_ids = node_matrix[:,1]
node_lons = node_matrix[:,3]
node_lats = node_matrix[:,4]

# Read in and organize arc data
arc_matrix = CSV.read(pwd()*"/arc_matrix.csv", DataFrame)
arc_ids = arc_matrix[:,1]
arc_1sts = arc_matrix[:,2]
arc_2nds = arc_matrix[:,3]
arc_dists = arc_matrix[:,5]

# Read in product data
product_matrix = CSV.read(pwd()*"/product_matrix.csv", DataFrame)
product_ids = product_matrix[:,1]
product_names = product_matrix[:,2]
product_transvc = product_matrix[:,3]

# Read in tech data
technology_matrix = CSV.read(pwd()*"/technology_matrix.csv", DataFrame)
tech_ids = technology_matrix[:,1]
tech_names = technology_matrix[:,6]
tech_capacities = Matrix(technology_matrix[:,7:9])

# Read in supply data
supply_matrix = CSV.read(pwd()*"/supply_matrix.csv", DataFrame)
supply_ids = supply_matrix[:,1]
supply_nodes = supply_matrix[:,2]
supply_products = supply_matrix[:,3]
supply_capacities = supply_matrix[:,5]
total_supply_capacity = sum(supply_capacities)

# Read in demand data
demand_matrix = CSV.read(pwd()*"/demand_matrix.csv", DataFrame)
demand_ids = demand_matrix[:,1]
demand_nodes = demand_matrix[:,2]
demand_products = demand_matrix[:,3]
demand_capacities = demand_matrix[:,5]



##################################################################
# Read in solution placement data
##################################################################

# Read solution text file into initial matrix (each line forming one row)
solution_data = readdlm(pwd()*solution_file, '\n')

# Transfer technology placements to dedicated vector
placement_data = Vector{String}(undef, 0)
for i = 1:size(solution_data, 1)
    # Include line only if it corresponds to technology placement data
    if solution_data[i][1:2] == "z("
        push!(placement_data, solution_data[i])
    end
end

# Split the placement integer from each line of placement_data off into a new vector
#   - MODIFICATION: 3 -> 2
placement_vector = getindex.(split.(placement_data, "= "), 2)

# Convert each number string to an int and reorganize into a matrix
placement_matrix = zeros(trunc(Int, length(placement_data)/2), 2)
for i = 1:2
    for j = 1:convert(Int, length(placement_data)/2)
        # Define an index to check in placement_vector depending on site and scale of technology
        ind = convert(Int, j + (length(placement_data)/2) * (i - 1))
        placement_matrix[j,i] = convert(Int, round(parse(Float64, placement_vector[ind])))
    end
end
placement_matrix = convert(Matrix{Int}, placement_matrix)

##################################################################
# Summarize tech & scale placement
##################################################################

# Initialize a matrix to summarize total tech placement of each type and each scale
summary_matrix = zeros(length(tech_ids), 2)

# Add up placed technologies of each type and scale combination
for i = 1:2
    for j = 1:convert(Int, length(placement_vector)/2)
        ind = convert(Int, j + (length(placement_data)/2) * (i - 1))
        for k = 1:length(tech_ids)
            if site_techs[j] == tech_ids[k]
                summary_matrix[k,i] += round(parse(Float64, placement_vector[ind]))
            end
        end
    end
end
summary_matrix = convert(Matrix{Int}, summary_matrix)

# Initialize a DataFrame for saving the summary_matrix to CSV
#   MODIFICATION: Removed "Scale 3" column
summary_headers = ["Tech Type", "Scale 1", "Scale 2"]
summary_df = DataFrame(A=String[], B=Int[], C=Int[])

# Fill the DataFrame with the rows of summary_matrix
for i = 1:size(summary_matrix, 1)
    push!(summary_df, (tech_ids[i], summary_matrix[i,1], summary_matrix[i,2]))
end

# Write the complete file
CSV.write(pwd()*"/summary_matrix.csv", summary_df, header=summary_headers)

##################################################################
# Write placements to Techmap File
##################################################################

# Prepare and initialize a techmap DataFrame
techmap_headers = ["1. Tech location reference ID", "2. Node ID", "3. Time ID", "4. Tech ID"]
techmap_df = DataFrame(A=Int[], B=Int[], C=Int[], D=String[])

# Add a tech ID for each nonzero value in placement_matrix
global techmap_counter = 0
for i = 1:size(placement_matrix, 1)
    for j = 1:size(placement_matrix, 2)
        if placement_matrix[i,j] != 0
            for k = 1:placement_matrix[i,j]
                global techmap_counter += 1
                node_id = site_nodes[i]
                tech_id = string(site_techs[i], ".", j)     # ex. id = "tA1.2"
                push!(techmap_df, (techmap_counter, node_id, 1, tech_id))
            end
        end
    end
end

# Write a file. This will be converted to non-integer form (to string IDs) for use with CSC.jl separately.
CSV.write(pwd()*"/techmap_matrix.csv", techmap_df, header=techmap_headers)

end