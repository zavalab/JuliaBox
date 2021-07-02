################################################################################
### JULIA PACMAGE IMPORTS ###
using JuMP
using Gurobi
using DelimitedFiles

#******************************************************************************#
#********** FILE OPTIONS, FOLDER PATHWAYS, AND HARD-CODED PARAMETERS **********#
#******************************************************************************#
################################################################################
### NOTES ###
#=

=#
################################################################################
### CONTROL FLOW ###
# Print the model in the REPL? (Note: model is always output as a text file.)
PrintModel = false

# Print output to REPL? (Note: model output is always printed to text file.)
PrintOutput = false

# Specify separator text for printing
PrintSpacer = "*"^50

# Use custom arc lengths or determine from longitude/latitude coordinates?
CustomLengths = true

################################################################################
### SET DIRECTORY ###
# Where data files are stored
DataFolder = "TestingFiles/WisconsinV3" # Base case
#DataFolder = "TestingFiles/WisconsinV3NoStorage" # To illustrate steady-state comparison
#DataFolder = "TestingFiles/WisconsinV3TripleManure" # Motivating additional study case
#DataFolder = "TestingFiles/WisconsinV3UnlimitedStorage" # Second attempt at illustrating price volatility

# Output model file name
ModelOutputFileName = "_Model.txt"

# Output solution data file name
SolutionOutputFileName = "_SolutionData.txt"

################################################################################
### PARAMETERIZED VALUES ###
R = 6335.439 # Earth radius used for distance calculation

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

function GreatCircle(ref_lon, ref_lat, dest_lon, dest_lat)
    """
    Calculates the great circle distance between a reference location (ref)
    and a destination (dest) and returns the great circle distances in km
    """
    # haversine formula
    dist = 2*R*asin(sqrt(((sind((dest_lat - ref_lat)/2))^2) +
        cosd(dest_lat)*cosd(ref_lat)*((sind((dest_lon - ref_lon)/2))^2)))
    return(dist)
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
storage_data = DataCheck(DataFolder*"/storage_data.csv",'|',Any,true)
tech_data = DataCheck(DataFolder*"/technology_data.csv",'|',Any,true)
techsite_data = DataCheck(DataFolder*"/techsite_data.csv",'|',Any,true)
time_data = DataCheck(DataFolder*"/time_data.csv",'|',Any,true)

################################################################################
### DATA SETUP ###
# Node data parsing
node_id = node_data[:, 1]
[node_id[i] = string(node_id[i]) for i = 1:length(node_id)]
node_alias = Dict(zip(node_id, node_data[:, 2])) # node name
node_lon = Dict(zip(node_id, node_data[:, 3])) # longitude
node_lat  = Dict(zip(node_id, node_data[:, 4])) # latitude

# Arc data parsing
arc_id = arc_data[:, 1]
[arc_id[a] = string(arc_id[a]) for a = 1:length(arc_id)] # convert to string from substring
arc_n = Dict(zip(arc_id, arc_data[:, 2])) # first node "n"
arc_m  = Dict(zip(arc_id, arc_data[:, 3])) # second node "m"
arc_cap  = Dict(zip(arc_id, arc_data[:, 4])) # arc capacity
if CustomLengths # Only need to define this data if working with customized arc lengths
    arc_len  = Dict(zip(arc_id, arc_data[:, 5])) # arc length; optional for distances
end
cardA = length(arc_id)

# Product data parsing
product_id = product_data[:, 1]
[product_id[i] = String(product_id[i]) for i = 1:length(product_id)]
product_alias = Dict(zip(product_id, product_data[:, 2])) # product name
product_trspt = Dict(zip(product_id, product_data[:, 3])) # product transport cost

# Demand data parsing
dem_id = demand_data[:, 1]
[dem_id[i] = string(dem_id[i]) for i = 1:length(dem_id)]
dem_node  = Dict(zip(dem_id, demand_data[:, 2])) # demand node
dem_prod  = Dict(zip(dem_id, demand_data[:, 3])) # demand product
dem_bid   = Dict(zip(dem_id, demand_data[:, 4])) # demand bid
dem_cap   = Dict(zip(dem_id, demand_data[:, 5])) # demand capacity
dem_time  = Dict(zip(dem_id, demand_data[:, 6])) # demand period

# Supply data parsing
sup_id = supply_data[:, 1]
[sup_id[i] = string(sup_id[i]) for i = 1:length(sup_id)]
sup_node  = Dict(zip(sup_id, supply_data[:, 2])) # supply node
sup_prod  = Dict(zip(sup_id, supply_data[:, 3])) # supply product
sup_bid   = Dict(zip(sup_id, supply_data[:, 4])) # supply bid
sup_cap   = Dict(zip(sup_id, supply_data[:, 5])) # supply capacity
sup_time  = Dict(zip(sup_id, supply_data[:, 6])) # supply period
sup_sup = Dict(zip(sup_id, supply_data[:, 7])) # supply supplier entity (i.e., for plotting totals)
#sup_con = Dict(zip(sup_id, supply_data[:, 8])) # supplier contract

# Storage data parsing
sto_id = storage_data[:, 1]
[sto_id[i] = string(sto_id[i]) for i = 1:length(sto_id)]
sto_node  = Dict(zip(sto_id, storage_data[:, 2])) # storage node
sto_prod  = Dict(zip(sto_id, storage_data[:, 3])) # storage product
sto_bid   = Dict(zip(sto_id, storage_data[:, 4])) # storage bid
sto_cap   = Dict(zip(sto_id, storage_data[:, 5])) # storage capacity
#sto_init  = Dict(zip(sto_id, storage_data[:, 6])) # storage opening inventory

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
#tech_cap = NumericListFromCSV(tech_id, tech_input, tech_data[:, 8]) # capacity value
tech_cap = Dict(zip(tech_id, tech_data[:, 8])) # capacity value
tech_defn = Dict(zip(tech_id, tech_data[:, 9])) # technology description

# Technology site data
techsite_id = techsite_data[:, 1]
[techsite_id[i] = String(techsite_id[i]) for i = 1:length(techsite_id)]
techsite_node = Dict(zip(techsite_id, techsite_data[:, 2])) # technology node
techsite_tech = Dict(zip(techsite_id, techsite_data[:, 3])) # technology existing at specified site

# Time data parsing
time_id = time_data[:, 1]
[time_id[i] = String(time_id[i]) for i = 1:length(time_id)]

# Update
println(PrintSpacer*"\nData Load Complete\n"*PrintSpacer)

################################################################################
### SET INDICES AND ALIASES ###
A = arc_id
K = sto_id
N = node_id
M = node_id
P = product_id
Q = product_id
D = dem_id
S = sup_id
T = tech_id
L = techsite_id
Z = time_id

# Create an internal set of arcs to ensure that both forward and backward arcs exist
for a = 1:cardA
    ai = a + cardA # shift index for added arc
    push!(A, "AI"*lpad(a,ndigits(cardA),"0")) # add a label; repeat original numbers, but prefix differently
    arc_n[A[ai]] = arc_m[A[a]] # flip arc at new position
    arc_m[A[ai]] = arc_n[A[a]] # flip arc at new position
    arc_cap[A[ai]] = arc_cap[A[a]] # same capacity (consistent with past code)
    if CustomLengths
        arc_len[A[ai]] = arc_len[A[a]] # same length (consistent with past code)
    end
end
cardA = length(A) # update cardA value

# Truth table for available arcs; nodes N and M are connected by A
AnM = Dict() # given node n, provides all nodes m in M connected by an arc a in A
Ain = Dict() # given node n, provides all arcs a in A directed towards node n
Aout = Dict() # given node n, provides all arcs a in A directed out of node n
[AnM[n] = Vector{String}(undef,0) for n in N]
[Ain[n] = Vector{String}(undef,0) for n in N]
[Aout[n] = Vector{String}(undef,0) for n in N]
for a in A
    push!(AnM[arc_n[a]], arc_m[a])
    push!(Ain[arc_m[a]],a)
    push!(Aout[arc_n[a]],a)
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

# Z subsets for balance constraint definition
Z1 = [Z[1]] # the first time index needs its own balance; create as a 1D array

Zz = Z[1:end-1] # for dual prices, we need a z+1 referencing system

ZZ = Z[2:end] # the "interior" time points share a common structure

Zprior = Dict()
for z = 2:length(Z)
    Zprior[Z[z]] = Z[z-1]
end

Zpost = Dict()
for z = 1:(length(Z)-1)
    Zpost[Z[z]] = Z[z+1]
end

################################################################################
### CALCULATED PARAMETERS ###
# Inter-node distances
node_dist = Dict()
if CustomLengths # Use the distances specified with arcs
    for a in A
        n = arc_n[a]
        m = arc_m[a]
        node_dist[n,m] = arc_len[a]
        #node_dist[m,n] = arc_len[a] # not needed in new formulation
    end
    println("### WARNING - NODE DISTANCES DETERMINED BY USER, SET \"CustomLengths\" TO \"false\" TO USE GREAT CIRCLE DISTANCES! ###")
else # calculate distances from node longitude and latitude coordinates
    for a in A
        node_dist[arc_n[a],arc_m[a]] = GreatCircle(node_lon[arc_n[a]], node_lat[arc_n[a]], node_lon[arc_m[a]], node_lat[arc_m[a]])
    end
end

# Maximum demand by node and product
dMAX = DictInit([N,P,Z],0)
for n in node_id
    for p in product_id
        for z in Z
            for j in dem_id
                if dem_node[j] == n && dem_prod[j] == p && dem_time[j] == z
                    dMAX[n,p,z] += dem_cap[j]
                end
            end
        end
    end
end

# Maximum supply by node and product
sMAX = DictInit([N,P,Z],0)
for n in node_id
    for p in product_id
        for z in Z
            for i in sup_id
                if sup_node[i] == n && sup_prod[i] == p && sup_time[i] == z
                    sMAX[n,p,z] += sup_cap[i]
                end
            end
        end
    end
end

# Maximum storage by node and product
vMAX = DictInit([N,P,Z],0)
for n in node_id
    for p in product_id
        for z in time_id
            for k in sto_id
                if sto_node[k] == n && sto_prod[k] == p #&& z != Z[end]
                    vMAX[n,p,z] += sto_cap[k]
                end
            end
        end
    end
end

# yield of p from q in technology t
α = Dict()
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
        gMAX[techsite_tech[ts],techsite_node[ts],p] = tech_output_stoich[techsite_tech[ts],p]*tech_cap[techsite_tech[ts]]
    end
end

xMAX = DictInit([T,N,Q],0)
for ts in techsite_id # produces (t,n) pairs
    for q in tech_input[techsite_tech[ts]] # a list of product id's
        xMAX[techsite_tech[ts],techsite_node[ts],q] = tech_input_stoich[techsite_tech[ts],q]*tech_cap[techsite_tech[ts]] # capacity for technology at ts
    end
end

# Require specialized set of big M constraints for flow variable f
# All n=m scenarios should be zero
fMAX = DictInit([A,P],0)
for a in A
    for p in P
        fMAX[a,p] = arc_cap[a]
    end
end

# Update
println(PrintSpacer*"\nData Processing Complete\n"*PrintSpacer)

################################################################################
### MODEL STATEMENT ###
MOD = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag"=>0))

################################################################################
### VARIABLES ###
@variable(MOD, 0 <= s[n=N,p=P,z=Z] <= sMAX[n,p,z]) # total supply by product and node
@variable(MOD, 0 <= d[n=N,p=P,z=Z] <= dMAX[n,p,z]) # total demand by product and node
@variable(MOD, 0 <= v[n=N,p=P,z=Z] <= vMAX[n,p,z]) # total demand by product and node
@variable(MOD, 0 <= sl[i=S] <= sup_cap[i]) # individual supplies by supply list
@variable(MOD, 0 <= dl[j=D] <= dem_cap[j]) # individual demands by demand list
@variable(MOD, 0 <= vl[k=K,z=Z] <= sto_cap[k]) # individual storage offering by list
@variable(MOD, 0 <= f[a=A,p=P,z=Z] <= fMAX[a,p]) # transport from n to m
@variable(MOD, 0 <= x[t=T,n=N,q=Q,z=Z] <= xMAX[t,n,q]) # consumption, standard: P committed for production of something else
@variable(MOD, 0 <= g[t=T,n=N,p=P,z=Z] <= gMAX[t,n,p]) # generation, standard:  P produced from something else

################################################################################
### EQUATIONS ###
# supply and demand total equal total of individual supplies and individual demands
@constraint(MOD, SupplyBalance[n=N,p=P,z=Z], s[n,p,z] == sum(sl[i] for i in S if sup_node[i] == n && sup_prod[i] == p && sup_time[i] == z))
@constraint(MOD, DemandBalance[n=N,p=P,z=Z], d[n,p,z] == sum(dl[j] for j in D if dem_node[j] == n && dem_prod[j] == p && dem_time[j] == z))
@constraint(MOD, StoragBalance[n=N,p=P,z=Z], v[n,p,z] == sum(vl[k,z] for k in K if sto_node[k] == n && sto_prod[k] == p))

# System mass balance
@constraint(MOD, Balance1[n=N,p=P,z=Z1], s[n,p,z]                    + sum(f[a,p,z] for a in Ain[n]) + sum(g[t,n,p,z] for t in NPt[n,p])
    == d[n,p,z] + v[n,p,z] + sum(f[a,p,z] for a in Aout[n]) + sum(x[t,n,p,z] for t in NQt[n,p]))
@constraint(MOD, BalanceZ[n=N,p=P,z=ZZ], s[n,p,z] + v[n,p,Zprior[z]] + sum(f[a,p,z] for a in Ain[n]) + sum(g[t,n,p,z] for t in NPt[n,p])
    == d[n,p,z] + v[n,p,z] + sum(f[a,p,z] for a in Aout[n]) + sum(x[t,n,p,z] for t in NQt[n,p]))

#=
OLD FORMULATION - EASIER TO SEE TRANSPORT VARIABLE RELATIONSHIPS
@constraint(MOD, Balance1[n=N,p=P,z=Z1], s[n,p,z] +                  + sum(f[m,n,p,z] for m in AnM[n]) + sum(g[t,n,p,z] for t in NPt[n,p])
    == d[n,p,z] + v[n,p,z] + sum(f[n,m,p,z] for m in AnM[n]) + sum(x[t,n,p,z] for t in NQt[n,p]))
@constraint(MOD, BalanceZ[n=N,p=P,z=ZZ], s[n,p,z] + v[n,p,Zprior[z]] + sum(f[m,n,p,z] for m in AnM[n]) + sum(g[t,n,p,z] for t in NPt[n,p])
    == d[n,p,z] + v[n,p,z] + sum(f[n,m,p,z] for m in AnM[n]) + sum(x[t,n,p,z] for t in NQt[n,p]))
=#

# Conversion relationships (yield-based)
@constraint(MOD, Conversion[n=N,(t,p,q) in TPQ,z=Z], g[t,n,p,z] == α[t,p,q]*x[t,n,q,z]) # generation simultaneous with consumption

################################################################################
### OBJECTIVE ###
# Cost expressions
demand_revenue = @expression(MOD, sum(dem_bid[j]*dl[j] for j in D))
supply_revenue = @expression(MOD, sum(sup_bid[i]*sl[i] for i in S))
transport_cost = @expression(MOD, sum(product_trspt[p]*f[a,p,z] for a in A,p in P,z in Z))
operating_cost = @expression(MOD, sum(x[t,n,q,z]*tech_bid[t] for t in T, n in N, q in Q, z in Z if (NQT[n,q,t] && tech_bid_ref[t] == q)))
storage_cost = @expression(MOD, sum(sto_bid[k]*vl[k,z] for k in K,z in Z))

# Full objective
@objective(MOD, Max, demand_revenue - supply_revenue - transport_cost - operating_cost - storage_cost)

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
v_out = JuMP.value.(v)
sl_out = JuMP.value.(sl)
dl_out = JuMP.value.(dl)
vl_out = JuMP.value.(vl)
f_out = JuMP.value.(f)
x_out = JuMP.value.(x)
g_out = JuMP.value.(g)
cp_out1 = JuMP.dual.(Balance1)
cp_outZ = JuMP.dual.(BalanceZ)
#cp_out = JuMP.dual.(Balance)
# merge the duals together along the time axis
cp_out = DictInit([N,P,Z],0)
for n in N
    for p in P
        for z in Z
            if z == Z[1] # or: z in Z1
                cp_out[n,p,z] = cp_out1[n,p,z]
            else # z in ZZ
                cp_out[n,p,z] = cp_outZ[n,p,z]
            end
        end
    end
end

################################################################################
### PROFIT CALCULATIONS ###

PhiSupply = Dict()
for i in S
    PhiSupply[i] = (cp_out[sup_node[i],sup_prod[i],sup_time[i]] - sup_bid[i])*sl_out[i]
end

PhiDemand = Dict()
for j in D
    PhiDemand[j] = (dem_bid[j] - cp_out[dem_node[j],dem_prod[j],dem_time[j]])*dl_out[j]
end

PhiStorage = Dict()
cp_storage = Dict()
for n in N
    for p in P
        for z in Zz
            cp_storage[n,p,z] = cp_out[n,p,Zpost[z]] - cp_out[n,p,z] # indexes into future; undefined for z = Z[end]
        end
    end
end
for k in K
    n = sto_node[k]
    p = sto_prod[k]
    PhiStorage[k] = sum((cp_storage[n,p,z] - sto_bid[k])*vl_out[k,z] for z in Zz)
end

PhiTransport = Dict()
cp_transport = Dict()
for a in A
    for p in P
        for z in Z
            # transport prices (notation: n -> m for p)
            cp_transport[a,p,z] = cp_out[arc_m[a],p,z] - cp_out[arc_n[a],p,z]
        end
        # transport profits
        PhiTransport[a,p] = sum((cp_transport[a,p,z] - product_trspt[p])*f_out[a,p,z] for z in Z)
    end
end

PhiTech = DictInit([T,N],0)
cp_tech = DictInit([T,N,Z],0)
for t in T
    for n in N
        for z in Z
            # technology prices
            cp_tech[t,n,z] = sum(cp_out[n,p,z]*tech_output_stoich[t,p] for p in tech_output[t]) - sum(cp_out[n,q,z]*tech_input_stoich[t,q] for q in tech_input[t])
        end
    end
end
for ts in techsite_id
    n = techsite_node[ts]
    t = techsite_tech[ts]
    p = tech_bid_ref[t]
    PhiTech[t,n] = sum((cp_tech[t,n,z] - tech_input_stoich[t,p]*tech_bid[t])*x_out[t,n,p,z] for z in Z)
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
FilePrint(s_out,[N,P,Z],filename,DataName="Supply values:",VarName="s")
FilePrint(d_out,[N,P,Z],filename,DataName="Demand values:",VarName="d")
FilePrint(v_out,[N,P,Z],filename,DataName="Storage values:",VarName="v")
FilePrint(f_out,[A,P,Z],filename,DataName="Transport values:",VarName="f")
FilePrint(x_out,[T,N,P,Z],filename,DataName="Consumption values:",VarName="x")
FilePrint(g_out,[T,N,P,Z],filename,DataName="Generation values:",VarName="g")
FilePrint(cp_out,[N,P,Z],filename,DataName="Nodal clearing prices:",VarName="π")
FilePrint(cp_storage,[N,P,Zz],filename,DataName="Storage clearing prices:",VarName="π_k")
FilePrint(cp_transport,[A,P,Z],filename,DataName="Transport clearing prices:",VarName="π_f")
FilePrint(cp_tech,[T,N,Zz],filename,DataName="Technology clearing prices:",VarName="π_t")
FilePrint(PhiDemand,[D],filename,DataName="Demand profits:",VarName="Φ_d")
FilePrint(PhiSupply,[S],filename,DataName="Supply profits:",VarName="Φ_s")
FilePrint(PhiStorage,[K],filename,DataName="Storage profits:",VarName="Φ_v")
FilePrint(PhiTransport,[A,P],filename,DataName="Transport profits:",VarName="Φ_f")
FilePrint(PhiTech,[T,N],filename,DataName="Technology profits:",VarName="Φ_t")
# To add source data to files:
print(filename, "\n"*PrintSpacer*"\nRaw data:")
RawDataPrint(node_data,filename,DataName="Node data:")
RawDataPrint(arc_data,filename,DataName="Arc data:")
RawDataPrint(product_data,filename,DataName="Product data:")
RawDataPrint(demand_data,filename,DataName="Demand data:")
RawDataPrint(supply_data,filename,DataName="Supply data:")
RawDataPrint(storage_data,filename,DataName="Storage data:")
RawDataPrint(tech_data,filename,DataName="Technology data:")
RawDataPrint(techsite_data,filename,DataName="Technology site data:")
RawDataPrint(time_data,filename,DataName="Time index data:")
close(filename)

# And print to REPL
if PrintOutput
    println(PrintSpacer)
    println("Objective value: ", z_out)
    # NOTE to display all values, replace Nonzeros() with DictInit([A,B,...Z],true), or (for specific values) with a TruthTable
    PrettyPrint(s_out, [N,P,Z], Nonzeros(s_out, [N,P,Z]), DataName="Supply values:", VarName="s")
    PrettyPrint(d_out, [N,P,Z], Nonzeros(d_out, [N,P,Z]), DataName = "Demand values:", VarName="d")
    PrettyPrint(v_out, [N,P,Z], Nonzeros(v_out, [N,P,Z]), DataName = "Storage values:", VarName="v")
    PrettyPrint(f_out, [A,P,Z], Nonzeros(f_out, [A,P,Z]), DataName = "Transport values:", VarName="f")
    PrettyPrint(x_out, [T,N,P,Z], Nonzeros(x_out, [T,N,P,Z]), DataName = "Consumption values:", VarName="x")
    PrettyPrint(g_out, [T,N,P,Z], Nonzeros(g_out, [T,N,P,Z]), DataName = "Generation values:", VarName="g")
    PrettyPrint(cp_out, [N,P,Z], Nonzeros(cp_out, [N,P,Z]), DataName = "Nodal clearing prices:", VarName="π")
    PrettyPrint(cp_storage, [N,P,Zz], Nonzeros(cp_storage, [N,P,Zz]), DataName = "Storage clearing prices:", VarName="π")
    PrettyPrint(cp_transport, [A,P,Z], Nonzeros(cp_transport, [A,P,Z]), DataName = "Transport clearing prices:", VarName="π_f")
    PrettyPrint(cp_tech, [T,N,Zz], Nonzeros(cp_tech, [T,N,Zz]), DataName = "Technology clearing prices:", VarName="π_t")
    PrettyPrint(PhiDemand, [D], Nonzeros(PhiDemand, [D]), DataName = "Demand profits:", VarName="Φ_D")
    PrettyPrint(PhiSupply, [S], Nonzeros(PhiSupply, [S]), DataName = "Supply profits:", VarName="Φ_S")
    PrettyPrint(PhiStorage, [K], Nonzeros(PhiStorage, [K]), DataName = "Storage profits:", VarName="Φ_K")
    PrettyPrint(PhiTransport, [A,P], Nonzeros(PhiTransport, [A,P]), DataName = "Transport profits:", VarName="Φ_f")
    PrettyPrint(PhiTech, [T,N], Nonzeros(PhiTech, [T,N]), DataName = "Technology profits:", VarName="Φ_t")
end

# Save records
include("RecordMaker.jl")

################################################################################
### END OF CODE ###
println(PrintSpacer,"\nAll done\n",PrintSpacer)