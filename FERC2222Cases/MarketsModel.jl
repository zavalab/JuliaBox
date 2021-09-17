################################################################################
### NULIA PACMAGE IMPORTS ###
using JuMP
using Clp
using DelimitedFiles

#******************************************************************************#
#********** FILE OPTIONS, FOLDER PATHWAYS, AND HARD-CODED PARAMETERS **********#
#******************************************************************************#
################################################################################
### NOTES ###
#=
=#
################################################################################
### FILE OPTIONS ###
# Print the model in the REPL? (Note: model is always output as a text file.)
PrintModel = false

# Print output to REPL? (Note: model output is always printed to text file.)
PrintOutput = true

# Specify separator text for printing
PrintSpacer = "*"^50

################################################################################
### SET DIRECTORY ###
# Where data files are stored; presumably a folder in ModelDirectory
DataFolder = "Post 171-69"

# Output model file name
ModelOutputFileName = "_Model.txt"

# Output solution data file name
SolutionOutputFileName = "_SolutionData.txt"

#******************************************************************************#
#******************************************************************************#
#******************************************************************************#

################################################################################
### USER FUNCTIONS ###
function DataCheck(MaybeEmptyFile,Delimitor,Type,Comments)
    """
    Given:
    >MaybeEmptyFile: a .csv file that may or may not be empty
    >Delimitor: the delimitor for the .csv file
    >Type: output readdlm() data type; i.e., Any
    >Comments: Boolean; are there comment lines in the file?
    Returns: data contained in the file (if any) otherwise an empty array: []
    """
    data = try
        readdlm(MaybeEmptyFile,Delimitor,Type,comments=Comments)
    catch
        []
    end
    return(data)
end

function TextListFromCSV(IDs,DataCol)
    """
    This function recovers comma-separated lists of text labels stored as
    single entries in a csv file (separated by a non-comma separator) and
    returns them as a dictionary indexed by the given IDs
    i.e., |P01,P02,P03| -> ["P01","P02","P03"]
    Inputs:
        IDs - a list of labels to be used as dictionary keys
        DataCol - The column of data with rows corresponding to the labels in
            IDS, with each containing one or more text labels separated by commas
    Outputs:
        Out - a dictionary mapping the keys in IDs to the values in DataCol

    Note: This function replaces code of the form:
    for s = 1:length(asset_id) # AssetData[s,4] is a comma-separated list
        asset_inputs[asset_id[s]] = [String(i) for i in split(AssetData[s,4],",")]
    end
    """
    Out = Dict()
    for i = 1:length(IDs)
        Out[IDs[i]] = [String(i) for i in split(DataCol[i],",")]
    end
    return Out
end

function NumericListFromCSV(IDs,ID2s,DataCol)
    """
    For data that may be stored as a list of numeric values within CSV data
    i.e., |0.1,0.2,0.3,0.4| in a "|"-separated data file. This function parses
    these data as Float64 values and assigns them to a dictionary using given
    keys; the dictionary is returned.
    Inputs:
        IDs - a list of IDs corresponding to the rows of the data in DataCol; also used as dict keys
        ID2s - a list of secondary IDs for use as keys in the dict; a list of lists
        DataCol - a column of data from a .csv file

    Note: This function replaces code like this:
    for s = 1:length(asset_id) # AssetData[s,7] is a comma-separated list
        check = typeof(AssetData[s, 7])
        if check != Float64 && check != Int64 # then it's a string of numbers
            temp = split(asset_data[s, 7],",")
            values = [parse(Float64,temp[i]) for i = 1:length(temp)] # Float64 array
        else # in the case it was a Float64 or an Int64
            values = asset_data[s, 7]
        end
        key1 = asset_id[s]
        key2s = asset_inputs[key1]
        for i = 1:length(key2s)
            asset_input_stoich[key1,key2s[i]] = values[i]
        end
    end
    """
    Out = Dict()
    L = length(IDs)
    for l = 1:L
        check = typeof(DataCol[l])
        if check != Float64 && check != Int64 # then it's a string of numbers
            temp = split(DataCol[l],",") # separate by commas into temporary list
            values = [parse(Float64,temp[i]) for i = 1:length(temp)] # parse as Float64 array
        else # in the case it was already a Float64 or an Int64
            values = DataCol[l]
        end
        key1 = IDs[l]
        key2s = ID2s[key1]
        for i = 1:length(key2s)
            Out[key1,key2s[i]] = values[i]
        end
    end
    return Out
end

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

function PrettyPrint(Data, OrderedIndexList, TruthTable=nothing; Header="*"^50, DataName="", VarName="")
    """
    Prints Data values indexed by OrderedIndexList based on whether or not
    the corresponding index value of TruthTable is true; Data and TruthTable
    share the same index pattern.
    Inputs:
        Data - an indexed data structure; a dictionary or JuMP variable
        OrderedIndexList - a list of the indices corresponding to Data and TruthTable
        TruthTable - a dictionary indexed on OrderedindexList which outputs Boolean values (true/false)
        Header - A string to print above any data
        DataName - A string to be printed as a header for the data
        VarName - A string to be printed before each line
    """
    # check truthtable; if not, default to true
    if TruthTable == nothing
        TruthTable = DictInit(OrderedIndexList, true)
    end
    # Start printing headers
    println(Header)
    println(DataName)
    # Check index index list length for proper index handling
    if length(OrderedIndexList) == 1
        SplatIndex = OrderedIndexList[1]
        for index in SplatIndex
            if TruthTable[index]
                println(VarName*"(\""*string(index)*"\"): "*string(Data[index]))
            end
        end
    else
        SplatIndex = collect(Iterators.product(OrderedIndexList...))
        for index in SplatIndex
            if TruthTable[index...]
                println(VarName*string(index)*": "*string(Data[index...]))
            end
        end
    end
end

function Nonzeros(Data, OrderedIndexList; Threshold=1E-9)
    """
    Given a dictionary of data indexed by the labels in OrderedIndexList
    returns a Boolean dictionary pointing to nonzero indices in Data, where
    nonzero is subject to a threshold value Threshold, defaulting to 1E-9.
    """
    OutDict = Dict()
    if length(OrderedIndexList) == 1
        SplatIndex = OrderedIndexList[1]
        for index in SplatIndex
            if abs(Data[index]) > Threshold
                OutDict[index] = true
            else
                OutDict[index] = false
            end
        end
    else
        SplatIndex = collect(Iterators.product(OrderedIndexList...))
        for index in SplatIndex
            if abs(Data[[i for i in index]...]) > Threshold
                OutDict[[i for i in index]...] = true
            else
                OutDict[[i for i in index]...] = false
            end
        end
    end
    return(OutDict)
end

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

function RawDataPrint(data,filename;Header="*"^50,DataName="")
    """
    Prints raw data to file for record-keeping purposes;
    reduces required script code clutter.
    > data: an array of data read from a .csv file.
    > filename: the file name for printing
    > Header: a string to be used as a header above the printed data
    > DataName: a header to appear above the printed data
    """
    # Header and DataName:
    print(filename,"\n"*Header*"\n"*DataName)

    # number of rows
    n = size(data)[1]

    # print rows of data array to file
    for i = 1:n
        print(filename,"\n")
        print(filename,data[i,:])
    end
end

################################################################################
### DATA IMPORT ###
node_data = DataCheck(DataFolder*"/node_data.csv",'|',Any,true)
arc_data = readdlm(DataFolder*"/arc_data.csv",'|',Any,comments=true)
product_data = DataCheck(DataFolder*"/product_data.csv",'|',Any,true)
demand_data = DataCheck(DataFolder*"/demand_data.csv",'|',Any,true)
supply_data = DataCheck(DataFolder*"/supply_data.csv",'|',Any,true)
tech_data = DataCheck(DataFolder*"/technology_data.csv",'|',Any,true)
techsite_data = DataCheck(DataFolder*"/techsite_data.csv",'|',Any,true)

################################################################################
### DATA SETUP ###
# Node data parsing
node_id = node_data[:, 1]
[node_id[i] = string(node_id[i]) for i = 1:length(node_id)]
node_alias = Dict(zip(node_id, node_data[:, 2])) # node name
#node_lon = Dict(zip(node_id, node_data[:, 3])) # longitude
#node_lat  = Dict(zip(node_id, node_data[:, 4])) # latitude

# Arc data parsing
arc_id = arc_data[:, 1]
[arc_id[a] = string(arc_id[a]) for a = 1:length(arc_id)] # convert to string from substring
arc_n = Dict(zip(arc_id, arc_data[:, 2])) # first node "n"
arc_m  = Dict(zip(arc_id, arc_data[:, 3])) # second node "m"
arc_cap  = Dict(zip(arc_id, arc_data[:, 4])) # arc capacity
arc_len = Dict(zip(arc_id, arc_data[:, 5])) # arc length

# Product data parsing
product_id = product_data[:, 1]
[product_id[i] = String(product_id[i]) for i = 1:length(product_id)]
product_alias = Dict(zip(product_id, product_data[:, 2])) # product name
product_trspt = Dict(zip(product_id, product_data[:, 3])) # product transport cost

# Demand data parsing
dem_id = demand_data[:, 1]
[dem_id[i] = string(dem_id[i]) for i = 1:length(dem_id)]
dem_node  = Dict(zip(dem_id, demand_data[:, 2])); # demand node
dem_prod  = Dict(zip(dem_id, demand_data[:, 3])); # demand product
dem_bid   = Dict(zip(dem_id, demand_data[:, 4])); # demand bid
dem_cap   = Dict(zip(dem_id, demand_data[:, 5])); # demand capacity

# Supply data parsing
sup_id = supply_data[:, 1]
[sup_id[i] = string(sup_id[i]) for i = 1:length(sup_id)]
sup_node  = Dict(zip(sup_id, supply_data[:, 2])); # demand node
sup_prod  = Dict(zip(sup_id, supply_data[:, 3])); # demand product
sup_bid   = Dict(zip(sup_id, supply_data[:, 4])); # demand bid
sup_cap   = Dict(zip(sup_id, supply_data[:, 5])); # demand capacity

# Technology transformation parsing
tech_id = tech_data[:, 1]
[tech_id[i] = string(tech_id[i]) for i = 1:length(tech_id)] # de-substring-ification
tech_output = TextListFromCSV(tech_id, tech_data[:,2]) # products made
tech_input = TextListFromCSV(tech_id, tech_data[:, 3]) # product consumed
tech_output_stoich = NumericListFromCSV(tech_id, tech_output, tech_data[:, 4]) # product yield factor
tech_input_stoich = NumericListFromCSV(tech_id, tech_input, tech_data[:, 5]) # reference yield factor
tech_bid_ref = Dict(zip(tech_id, tech_data[:, 6])) # reference product for bids
[tech_bid_ref[i] = string(tech_bid_ref[i]) for i in keys(tech_bid_ref)] # de-substring-ification
tech_bid = Dict(zip(tech_id, tech_data[:, 7])) # bid value; interpreted as operating cost
tech_cap = NumericListFromCSV(tech_id, tech_input, tech_data[:, 8]) # capacity value
tech_defn = Dict(zip(tech_id, tech_data[:, 9])) # technology description

# Technology site data
techsite_id = techsite_data[:, 1]
[techsite_id[i] = String(techsite_id[i]) for i = 1:length(techsite_id)]
techsite_node = Dict(zip(techsite_id, techsite_data[:, 2])) # technology node
techsite_tech = Dict(zip(techsite_id, techsite_data[:, 3])) # technology existing at specified site

################################################################################
### SET INDICES AND ALIASES ###
A = arc_id
N = node_id
M = node_id
P = product_id
Q = product_id
D = dem_id
S = sup_id
T = tech_id
L = techsite_id

# Truth table for available arcs; nodes N and M are connected by A
ANM = DictInit([N,M],false) # Truth table for available arcs; nodes N and M are connected by A
NMA = Dict() # Reverse-dictionary for arc index from nodes
NMP = [] # List of node pairs in arcs, forward and reverse, with products
for a in A
    ANM[arc_n[a],arc_m[a]] = true # forward
    ANM[arc_m[a],arc_n[a]] = true # reverse
    NMA[arc_n[a],arc_m[a]] = a
    NMA[arc_m[a],arc_n[a]] = a
    for p in P
        push!(NMP,(arc_n[a],arc_m[a],p))
        push!(NMP,(arc_m[a],arc_n[a],p))
    end
end
NMPTruthTable = DictInit([N,M,P],false)
for key in NMP
    NMPTruthTable[key] = true
end

LTN = []
for l in L
    push!(LTN,(techsite_tech[l],techsite_node[l]))
end

# Mapping set: product P can be used to produce product Q
TPQ = []
for t in tech_id
    for p in tech_output[t]
        for q in tech_input[t]
            push!(TPQ,(t,p,q)) # (prod,ref) pairs
        end
    end
end
TPQ = unique(TPQ) # remove duplicates

# Mapping set: node n produces product p from technology t (truth table and index list)
NPT = DictInit([N,P,T],false)
NPTlist = []
# Mapping set: node n consumes product q for technology t (truth table and index list)
NQT = DictInit([N,Q,T],false)
NQTlist = []
for ts in techsite_id
    n = techsite_node[ts] # element
    t = techsite_tech[ts] # element
    peas = tech_output[t] # list of elements
    ques = tech_input[t] # list of elements
    for p in peas
        NPT[n,p,t] = true
        push!(NPTlist,(n,p,t))
    end
    for q in ques
        NQT[n,q,t] = true
        push!(NQTlist,(n,q,t))
    end
end

# given n,p, provide the associated t's
NPt = Dict()
for n in node_id
    for p in product_id
        teas = [] # empty list of t indices
        for t in tech_id
            if NPT[n,p,t] # is true
                push!(teas,t)
            end
        end
        NPt[n,p] = teas
    end
end
# given n,q, provide the associated t's
NQt = Dict()
for n in node_id
    for q in product_id
        teas = [] # empty list of t indices
        for t in tech_id
            if NQT[n,q,t] # is true
                push!(teas,t)
            end
        end
        NQt[n,q] = teas
    end
end

################################################################################
### CALCULATED PARAMETERS ###
# Inter-node distances
node_dist = Dict()
for a in A
    n = arc_n[a]
    m = arc_m[a]
    node_dist[n,m] = arc_len[a]
    node_dist[m,n] = arc_len[a]
end

# Maximum demand by node and product
dMAX = DictInit([N,P],0)
for n in node_id
    for p in product_id
        for j in dem_id
            if dem_node[j] == n && dem_prod[j] == p
                dMAX[n,p] += dem_cap[j]
            end
        end
    end
end

# Maximum supply by node and product
sMAX = DictInit([N,P],0)
for n in node_id
    for p in product_id
        for i in sup_id
            if sup_node[i] == n && sup_prod[i] == p
                sMAX[n,p] += sup_cap[i]
            end
        end
    end
end

α = Dict() # yield of p from q in technology t
for t in tech_id
    for p in tech_output[t] # list of products made by t
        for q in tech_input[t] # list of products consumed by t
            α[t,p,q] = tech_output_stoich[t,p]/tech_input_stoich[t,q]
        end
    end
end

# Maximum production levels
gMAX = DictInit([T,N,P],0)
for ts in techsite_id # techsite produces (t,n) pairs
    for p in tech_output[techsite_tech[ts]] # list of output products from t
        for q in tech_input[techsite_tech[ts]] # list of input products for t
            # get technology capacity output yields cooresponding to all inputs q to t
            # here we are adding entries to gMAX with an extra index, and afer this loop, we will find the minimum stoichiometric yield
            gMAX[techsite_tech[ts],techsite_node[ts],p,q] = α[techsite_tech[ts],p,q]*tech_cap[techsite_tech[ts],q]
        end
        # get technology capacity and assign to gMAX by site and product index
        gMAX[techsite_tech[ts],techsite_node[ts],p] = minimum(gMAX[techsite_tech[ts],techsite_node[ts],p,q] for q in tech_input[techsite_tech[ts]])
    end
end

xMAX = DictInit([T,N,Q],0)
for ts in techsite_id # produces (t,n) pairs
    for q in tech_input[techsite_tech[ts]] # a list of product id's
        xMAX[techsite_tech[ts],techsite_node[ts],q] = tech_cap[techsite_tech[ts],q] # capacity for technology at ts
    end
end

# Require specialized set of big M constraints for flow variable f
# All n=m scenarios should be zero
fMAX = DictInit([N,M,P],0)
for a in arc_id
    n = arc_n[a]
    m = arc_m[a]
    if ANM[n,m]
        for p in product_id
            #fMAX[n,m,p] = minimum([sMAX[n,p] + sum(gMAX[t,n,p] for t in tech_id), arc_cap[a]])
            fMAX[n,m,p] = arc_cap[a]
            fMAX[m,n,p] = fMAX[n,m,p]
        end
    end
end

################################################################################
### MODEL STATEMENT ###
MOD = JuMP.Model(optimizer_with_attributes(Clp.Optimizer))

################################################################################
### VARIABLES ###
# Model variables
@variable(MOD, 0 <= s[n=N,p=P] <= sMAX[n,p]) # total supply by product and node
@variable(MOD, 0 <= d[n=N,p=P] <= dMAX[n,p]) # total demand by product and node
@variable(MOD, 0 <= sl[i=S] <= sup_cap[i]) # individual supplies by supply list
@variable(MOD, 0 <= dl[j=D] <= dem_cap[j]) # individual demands by demand list
@variable(MOD, 0 <= f[n=N,m=M,p=P] <= fMAX[n,m,p]) # transport from n to m
@variable(MOD, 0 <= x[t=T,n=N,q=Q] <= xMAX[t,n,q]) # consumption, standard: P committed for production of something else
@variable(MOD, 0 <= g[t=T,n=N,p=P] <= gMAX[t,n,p]) # generation, standard:  P produced from something else
# Objective function parts
@variable(MOD, 0 <= transport_cost) # objective transport cost term
@variable(MOD, 0 <= operating_cost) # objective technology operating cost term

################################################################################
### EQUATIONS ###
# supply and demand total equal total of individual supplies and individual demands
@constraint(MOD, SupplyBalance[n=N,p=P], s[n,p] == sum(sl[i] for i in S if sup_node[i] == n && sup_prod[i] == p))
@constraint(MOD, DemandBalance[n=N,p=P], d[n,p] == sum(dl[j] for j in D if dem_node[j] == n && dem_prod[j] == p))

# System mass balance
@constraint(MOD, Balance[n=N,p=P], s[n,p] + sum(f[m,n,p] for m in M if ANM[m,n]) + sum(g[t,n,p] for t in NPt[n,p])
    == d[n,p] + sum(f[n,m,p] for m in M if ANM[n,m]) + sum(x[t,n,p] for t in NQt[n,p]))

# Conversion relationships (yield-based)
@constraint(MOD, Conversion[n=N,(t,p,q) in TPQ], g[t,n,p] == α[t,p,q]*x[t,n,q])

################################################################################
### OBJECTIVE ###
demand_revenue = @expression(MOD, sum(dem_bid[dd]*dl[dd] for dd in D))
supply_revenue = @expression(MOD, sum(sup_bid[ss]*sl[ss] for ss in S))
#transport_revenue = @expression(MOD, sum(product_trspt[p]*node_dist[n,m]*f[n,m,p] for n in N,m in M,p in P))
transport_revenue = @expression(MOD, sum(product_trspt[p]*node_dist[n,m]*f[n,m,p] for (n,m,p) in NMP))
operating_revenue = @expression(MOD, sum(x[t,n,q]*tech_bid[t] for t in T, n in N, q in Q if (NQT[n,q,t] && tech_bid_ref[t] == q)))

# Full objective
@objective(MOD, Max, demand_revenue - supply_revenue - transport_revenue - operating_revenue)

################################################################################
### DISPLAY MODEL FORMULATION ###
filename = open(DataFolder*"/"*ModelOutputFileName,"w")
print(filename, MOD)
close(filename)
if PrintModel == true
    print(MOD)
end

################################################################################
### SOLVE AND DATA RETRIEVAL ###
# Display statistics
println(PrintSpacer*"\nModel statistics:")
println("Variables: "*string(length(all_variables(MOD))))
println("Total inequality constraints: "*string(num_constraints(MOD,AffExpr, MOI.LessThan{Float64})+num_constraints(MOD,AffExpr, MOI.GreaterThan{Float64})+num_constraints(MOD,VariableRef, MOI.LessThan{Float64})+num_constraints(MOD,VariableRef, MOI.GreaterThan{Float64})))
println("Total equality constraints: "*string(num_constraints(MOD,VariableRef, MOI.EqualTo{Float64})+num_constraints(MOD,AffExpr, MOI.EqualTo{Float64})))
println("Variable bounds: "*string(num_constraints(MOD,VariableRef, MOI.LessThan{Float64})+num_constraints(MOD,VariableRef, MOI.GreaterThan{Float64})))
println("Model inequality constraints: "*string(num_constraints(MOD,AffExpr, MOI.LessThan{Float64})+num_constraints(MOD,AffExpr, MOI.GreaterThan{Float64})))
println("Model equality constraints: "*string(num_constraints(MOD,AffExpr, MOI.EqualTo{Float64})))
println(PrintSpacer)

# Solve
println("Solving original problem...")
JuMP.optimize!(MOD)
println("Primal status: "*string(primal_status(MOD)))
println("Dual status: "*string(dual_status(MOD)))
z_out = JuMP.objective_value(MOD)
s_out = JuMP.value.(s)
d_out = JuMP.value.(d)
sl_out = JuMP.value.(sl)
dl_out = JuMP.value.(dl)
f_out = JuMP.value.(f)
x_out = JuMP.value.(x)
g_out = JuMP.value.(g)
cp_out = JuMP.dual.(Balance)

################################################################################
### PROFIT CALCULATIONS ###

PhiSupply = Dict()
for i in S
    PhiSupply[i] = (cp_out[sup_node[i],sup_prod[i]] - sup_bid[i])*sl_out[i]
end

PhiDemand = Dict()
for j in D
    PhiDemand[j] = (dem_bid[j] - cp_out[dem_node[j],dem_prod[j]])*dl_out[j]
end

PhiTransport = Dict()
cp_transport = Dict()
for n in N
    for m in M
        for p in P
            # transport prices (notation: n -> m for p)
            cp_transport[n,m,p] = cp_out[m,p] - cp_out[n,p]
            # transport profits
            PhiTransport[n,m,p] = (cp_transport[n,m,p] - product_trspt[p])*f_out[n,m,p]
        end
    end
end

PhiTech = DictInit([T,N],0)
cp_tech = Dict()
for t in T
    for n in N
        # technology prices
        cp_tech[t,n] = sum(cp_out[n,p]*tech_output_stoich[t,p] for p in tech_output[t]) - sum(cp_out[n,q]*tech_input_stoich[t,q] for q in tech_input[t])
    end
end
for ts in techsite_id
    n = techsite_node[ts]
    t = techsite_tech[ts]
    p = tech_bid_ref[t]
    PhiTech[t,n] = (cp_tech[t,n] - tech_input_stoich[t,p]*tech_bid[t])*x_out[t,n,tech_bid_ref[t]]
end

################################################################################
### WRITE OUTPUT TO FILE AND DISPLAY TO REPL ###
filename = open(DataFolder*"/"*SolutionOutputFileName,"w")
print(filename, PrintSpacer*"\nModel statistics:")
print(filename, "\nVariables: "*string(length(all_variables(MOD))))
print(filename, "\nTotal inequality constraints: "*string(num_constraints(MOD,AffExpr, MOI.LessThan{Float64})+num_constraints(MOD,AffExpr, MOI.GreaterThan{Float64})+num_constraints(MOD,VariableRef, MOI.LessThan{Float64})+num_constraints(MOD,VariableRef, MOI.GreaterThan{Float64})))
print(filename, "\nTotal equality constraints: "*string(num_constraints(MOD,VariableRef, MOI.EqualTo{Float64})+num_constraints(MOD,AffExpr, MOI.EqualTo{Float64})))
print(filename, "\nVariable bounds: "*string(num_constraints(MOD,VariableRef, MOI.LessThan{Float64})+num_constraints(MOD,VariableRef, MOI.GreaterThan{Float64})))
print(filename, "\nModel inequality constraints: "*string(num_constraints(MOD,AffExpr, MOI.LessThan{Float64})+num_constraints(MOD,AffExpr, MOI.GreaterThan{Float64})))
print(filename, "\nModel equality constraints: "*string(num_constraints(MOD,AffExpr, MOI.EqualTo{Float64})))
print(filename, "\n"*PrintSpacer*"\nObjective value: "*string(z_out))
FilePrint(s_out,[N,P],filename,DataName="Supply values:",VarName="s")
FilePrint(d_out,[N,P],filename,DataName="Demand values:",VarName="d")
FilePrint(f_out,[N,M,P],filename,DataName="Transport values:",VarName="f")
FilePrint(x_out,[T,N,P],filename,DataName="Consumption values:",VarName="x")
FilePrint(g_out,[T,N,P],filename,DataName="Generation values:",VarName="g")
FilePrint(cp_out,[N,P],filename,DataName="Nodal clearing prices:",VarName="π")
FilePrint(cp_transport,[N,M,P],filename,DataName="Transport clearing prices:",VarName="π_f")
FilePrint(cp_tech,[T,N],filename,DataName="Technology clearing prices:",VarName="π_t")
FilePrint(PhiDemand,[D],filename,DataName="Demand profits:",VarName="Φ_D")
FilePrint(PhiSupply,[S],filename,DataName="Supply profits:",VarName="Φ_S")
FilePrint(PhiTransport,[N,M,P],filename,DataName="Transport profits:",VarName="Φ_f")
FilePrint(PhiTech,[T,N],filename,DataName="Technology profits:",VarName="Φ_t")
# To add source data to files:
print(filename, "\n"*PrintSpacer*"\nRaw data:")
RawDataPrint(node_data,filename,DataName="Node data:")
RawDataPrint(arc_data,filename,DataName="Arc data:")
RawDataPrint(product_data,filename,DataName="Product data:")
RawDataPrint(demand_data,filename,DataName="Demand data:")
RawDataPrint(supply_data,filename,DataName="Supply data:")
RawDataPrint(tech_data,filename,DataName="Technology data:")
RawDataPrint(techsite_data,filename,DataName="Technology site data:")
close(filename)

# And print to REPL
if PrintOutput
    println(PrintSpacer)
    println("Objective value: ", z_out)
    # NOTE to display all values, replace Nonzeros() with DictInit([A,B,...Z],true)
    PrettyPrint(s_out, [N,P], Nonzeros(s_out, [N,P]), DataName="Supply values:", VarName="s")
    PrettyPrint(d_out, [N,P], Nonzeros(d_out, [N,P]), DataName = "Demand values:", VarName="d")
    PrettyPrint(f_out, [N,M,P], Nonzeros(f_out, [N,M,P]), DataName = "Transport values:", VarName="f")
    #PrettyPrint(f_out, [N,M,P], NMPTruthTable, DataName = "Transport values:", VarName="f")
    PrettyPrint(x_out, [T,N,P], Nonzeros(x_out, [T,N,P]), DataName = "Consumption values:", VarName="x")
    PrettyPrint(g_out, [T,N,P], Nonzeros(g_out, [T,N,P]), DataName = "Generation values:", VarName="g")
    PrettyPrint(cp_out, [N,P], Nonzeros(cp_out, [N,P]), DataName = "Nodal clearing prices:", VarName="π")
    PrettyPrint(cp_transport, [N,M,P], Nonzeros(cp_transport, [N,M,P]), DataName = "Transport clearing prices:", VarName="π_f")
    #PrettyPrint(cp_transport, [N,M,P], NMPTruthTable, DataName = "Transport clearing prices:", VarName="π_f")
    PrettyPrint(cp_tech, [T,N], Nonzeros(cp_tech, [T,N]), DataName = "Technology clearing prices:", VarName="π_t")
    PrettyPrint(PhiDemand, [D], Nonzeros(PhiDemand, [D]), DataName = "Demand profits:", VarName="Φ_D")
    PrettyPrint(PhiSupply, [S], Nonzeros(PhiSupply, [S]), DataName = "Supply profits:", VarName="Φ_S")
    PrettyPrint(PhiTransport, [N,M,P], Nonzeros(PhiTransport, [N,M,P]), DataName = "Transport profits:", VarName="Φ_f")
    PrettyPrint(PhiTech, [T,N], Nonzeros(PhiTech, [T,N]), DataName = "Technology profits:", VarName="Φ_t")
end

################################################################################
### END OF CODE ###
println(PrintSpacer,"\nAll done\n",PrintSpacer)
