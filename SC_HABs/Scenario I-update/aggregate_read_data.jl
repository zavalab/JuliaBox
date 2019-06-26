function Getindex(a,b);
    k = 0;
    for i in 1:length(a)
        if a[i] == b;
            k = i;
            break
        end
    end
    return k
end

## Reading data
node_matrix       = readdlm("node_matrix.csv",',');
demand_matrix     = readdlm("demand_matrix.csv",',');
product_matrix    = readdlm("product_matrix.csv",',');
Plimit_matrix     = readdlm("Plimit_matrix.csv",',');
Nlimit_matrix     = readdlm("Nlimit_matrix.csv",',');

arg1 = 1348;
## Define sets
NODES = node_matrix[:,1];
DEMS  = demand_matrix[:,1];
PRODS = product_matrix[:,1];
TIME  = 1:15 # From April to October


## Define dictionaries
node_lat   = Dict(zip(NODES, node_matrix[:,2]));           # latitude of each node
node_long  = Dict(zip(NODES, node_matrix[:,3]));           # longitude of each node


prod_trans = Dict(zip(PRODS, product_matrix[:,3]));        # transportation cost of each product $/km/tonne
prod_P     = Dict(zip(PRODS, product_matrix[:,4]));        # phsphorus release coefficient of each product kg/kg
prod_N     = Dict(zip(PRODS, product_matrix[:,5]));        # nitrogen release coefficient of each product kg/kg

dem_node   = Dict(zip(DEMS, demand_matrix[:,2]));          # node for the demand sink
dem_prod   = Dict(zip(DEMS, demand_matrix[:,3]));          # product type of demanded product
dem_time   = Dict(zip(DEMS, demand_matrix[:,4]));          # demand time
dem_cap    = Dict(zip(DEMS, demand_matrix[:,5]));          # maximum capacity of the required product tonne
dem_price  = Dict(zip(DEMS, demand_matrix[:,6]));          # price of the demanded product $/tonne

R = 6335.439

# distance between nodes
D = Dict(("n1", "n2") => 1.1);
for i in NODES
  for j in NODES
    D[(i, j)] = 2*R*asin(sqrt(sin((node_lat[j] - node_lat[i])*pi/2/180)^2 + cos(node_lat[j]*pi/180)*cos(node_lat[i]*pi/180)*sin((node_long[j] - node_long[i])*pi/2/180)^2)); ## Using the Haversine formula
  end
end

Plimit = Dict(("n1",1) => 0.1);
Nlimit = Dict(("n1",1) => 0.1);
for i in 1:length(NODES)
    for τ in TIME
        Plimit[(NODES[i],τ)] = Plimit_matrix[i,τ+1]/1000;
        Nlimit[(NODES[i],τ)] = Nlimit_matrix[i,τ+1]/1000;
    end
end

FARM1 = node_matrix[1:181,1];
LAND = node_matrix[182:1348,1];
Other = node_matrix[1349:end,1];



nutrient_info     = readdlm("nutrients_analysis.csv",',');
croptype = Dict(zip(LAND, nutrient_info[:,2]));
Pneed    = Dict(zip(LAND, nutrient_info[:,5]));
Nneed    = Dict(zip(LAND, nutrient_info[:,8]));
landarea = Dict(zip(LAND, nutrient_info[:,9]));
