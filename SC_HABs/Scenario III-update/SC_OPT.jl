## Supply Chain Optimization Model for transportation routing and technology placement
## Add time dimension
## Coded by Yicheng Hu 2018-07

## Reading data
node_matrix       = readdlm("ag_node_matrix.csv",',');
supply_matrix     = readdlm("supply_matrix.csv",',');
demand_matrix     = readdlm("ag_demand_matrix.csv",',');
inventory_matrix  = readdlm("inventory_matrix.csv",',');
technology_matrix = readdlm("technology_matrix.csv",',');
product_matrix    = readdlm("product_matrix.csv",',');
alpha_matrix      = readdlm("alpha_matrix.csv",',');
Plimit_matrix     = readdlm("ag_Plimit_matrix.csv",',');
Nlimit_matrix     = readdlm("ag_Nlimit_matrix.csv",',');


## Define sets
NODES = node_matrix[:,1];
arg1 = length(NODES);

SUPS  = supply_matrix[:,1];
DEMS  = demand_matrix[:,1];
TECHS = technology_matrix[:,1];
PRODS = product_matrix[:,1];
INVS   = inventory_matrix[:,1];
TIME  = 1:15 # From April to October
NODES_N = NODES[1:arg1]     # nodes that cannot be candidates for technologies
NODES_C = NODES[arg1+1:end] # nodes that can be candidates for technologies


## Define dictionaries
node_lat   = Dict(zip(NODES, node_matrix[:,2]));           # latitude of each node
node_long  = Dict(zip(NODES, node_matrix[:,3]));           # longitude of each node


prod_trans = Dict(zip(PRODS, product_matrix[:,3]));        # transportation cost of each product $/km/tonne
prod_P     = Dict(zip(PRODS, product_matrix[:,4]));        # phsphorus release coefficient of each product kg/kg
prod_N     = Dict(zip(PRODS, product_matrix[:,5]));        # nitrogen release coefficient of each product kg/kg

tech_ref   = Dict(zip(TECHS, technology_matrix[:,2]));     # reference product of each technology (with alpha = -1)
tech_cap   = Dict(zip(TECHS, 2*technology_matrix[:,3]));     # technlogy capacity (tonne/week)
tech_op    = Dict(zip(TECHS, technology_matrix[:,4]));     # operational cost of each technology $/tonne ref prod
tech_inv   = Dict(zip(TECHS, technology_matrix[:,5]));     # investment cost of each technology $

inv_node   = Dict(zip(INVS, inventory_matrix[:,2]));        # node for the storage system
inv_prod   = Dict(zip(INVS, inventory_matrix[:,3]));        # product type for the storage system
inv_cap    = Dict(zip(INVS, 1.5*inventory_matrix[:,4]/1000));   # inventory capacity tonne

sup_node   = Dict(zip(SUPS, supply_matrix[:,2]));          # node for the supply source
sup_prod   = Dict(zip(SUPS, supply_matrix[:,3]));          # product type of the supply source
sup_time   = Dict(zip(SUPS, supply_matrix[:,5]));          # supply time
sup_value  = Dict(zip(SUPS, 2*1.5*supply_matrix[:,4]/1000));          # amount of the supplied product tonne, 1.5 consider water used for washing
sup_price  = Dict(zip(SUPS, supply_matrix[:,6]));          # price of the supplied product $/tonne

dem_node   = Dict(zip(DEMS, demand_matrix[:,2]));          # node for the demand sink
dem_prod   = Dict(zip(DEMS, demand_matrix[:,3]));          # product type of demanded product
dem_time   = Dict(zip(DEMS, demand_matrix[:,5]));          # demand time
dem_cap    = Dict(zip(DEMS, 2*demand_matrix[:,4]));          # maximum capacity of the required product tonne
dem_price  = Dict(zip(DEMS, demand_matrix[:,6]));          # price of the demanded product $/tonne

R = 6335.439

# stoichiometric coefficient
α = Dict(("tA1","p1") => 0.5);         #Just used as an initiator to set up the dictionary with two keys
for t in 1:length(TECHS)
    for pr in 1: length(PRODS)
        α[(TECHS[t], PRODS[pr])] = alpha_matrix[t,pr];
    end
end

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
        Plimit[(NODES[i],τ)] = Plimit_matrix[i,τ+1];
        Nlimit[(NODES[i],τ)] = Nlimit_matrix[i,τ+1];
    end
end

using JuMP
using Gurobi
#using Cbc
m = Model(solver=GurobiSolver(Threads = 1,MIPGap = 1e-2, NodefileStart=0.25, TimeLimit = 43200));
#m = Model(solver = CbcSolver());

@variable(m, f[NODES,NODES,PRODS,TIME] >= 0);  # product flow
@variable(m, s[NODES,PRODS,TIME] >= 0);        # product supply
@variable(m, d[NODES,PRODS,TIME] >= 0);        # product demand
@variable(m, I[NODES,PRODS,[0;TIME]] >= 0);    # product inventory
#@variable(m, g[NODES_C,PRODS,TECHS,TIME]);     # generation/consumption amount
#@variable(m, y[NODES_C,TECHS], Bin)            # binary variable

# convert data form and set values/bounds
@constraint(m, [n in NODES, pr in PRODS, τ in TIME], s[n,pr,τ] == sum(sup_value[sup]*(sup_prod[sup] == pr)*(sup_node[sup] == n) for sup in SUPS));
@constraint(m, [n in NODES, pr in PRODS, τ in TIME], d[n,pr,τ] <= sum(dem_cap[dem]*(dem_prod[dem] == pr)*(dem_node[dem] == n) for dem in DEMS));
@constraint(m, [n in NODES, pr in PRODS, τ in TIME], I[n,pr,τ] <= sum(inv_cap[inv]*(inv_prod[inv] == pr)*(inv_node[inv] == n) for inv in INVS));
@constraint(m, [pr in PRODS, n in NODES], I[n,pr,0] == 0.63*sum(inv_cap[inv]*(inv_prod[inv] == pr)*(inv_node[inv] == n) for inv in INVS));  # set initial values of inventory levels
@constraint(m, [pr in PRODS, n in NODES], I[n,pr,15] <= 0.05*sum(inv_cap[inv]*(inv_prod[inv] == pr)*(inv_node[inv] == n) for inv in INVS)); # set ending constraints

#@variable(m, ss[SUPS] >= 0);    # supply source flow
#@variable(m, dd[DEMS] >= 0);    # demand sink flow
#@variable(m, ii[INVS] >= 0);    # inventory flow
#@constraint(m, [n in NODES, pr in PRODS, τ in TIME], s[n,pr,τ] == sum(ss[sup] for sup in SUPS if sup_prod[sup] == pr && sup_node[sup] == n && sup_time[sup] == τ));
#@constraint(m, [n in NODES, pr in PRODS, τ in TIME], d[n,pr,τ] == sum(dd[dem] for dem in DEMS if dem_prod[dem] == pr && dem_node[dem] == n && dem_time[dem] == τ));
#@constraint(m, [n in NODES, pr in PRODS, τ in TIME], I[n,pr,τ] <= sum(ii[inv] for inv in INVS if inv_prod[inv] == pr && inv_node[inv] == n));
#@constraint(m, [sup in SUPS], ss[sup] == sup_value[sup]);   # fix supply values
#@constraint(m, [dem in DEMS], dd[dem] <= dem_cap[dem]);     # set demand upper bounds
#@constraint(m, [inv in INVS], ii[inv] <= inv_cap[inv]);     # set inventory upper bounds

# balance constraints
@constraint(m, balance1[n in NODES_N, pr in PRODS, τ in TIME], I[n,pr,τ-1] + s[n,pr,τ] + sum(f[j,n,pr,τ] for j in NODES) == I[n,pr,τ] + d[n,pr,τ] + sum(f[n,j,pr,τ] for j in NODES));
#@constraint(m, balance2[n in NODES_C, pr in PRODS, τ in TIME], I[n,pr,τ-1] + s[n,pr,τ] + sum(f[j,n,pr,τ] for j in NODES) + sum(g[n,pr,t,τ] for t in TECHS) == I[n,pr,τ] + d[n,pr,τ] + sum(f[n,j,pr,τ] for j in NODES));

# transformation constraints
#@constraint(m, transfer[n in NODES_C, pr in PRODS, τ in TIME, t in TECHS], g[n,pr,t,τ] == α[t,pr]/α[t,tech_ref[t]]*g[n,tech_ref[t],t,τ]);
#@constraint(m, trasnfer_bounds_1[n in NODES_C, τ in TIME, t in TECHS], g[n,tech_ref[t],t,τ] <= 0);
#@constraint(m, transfer_bounds_2[n in NODES_C, τ in TIME, t in TECHS], g[n,tech_ref[t],t,τ] >= -y[n,t]*tech_cap[t]);

# economic metic calculation
@variable(m, Cinv >= 0);    # total investment cost
@variable(m, Cop >= 0);     # total oprational cost
@variable(m, Ctrans >= 0);  # total transportation cost
@variable(m, Revenue);      # total revenue
#@constraint(m, Cinv == 1/40*sum(y[n,t]*tech_inv[t] for n in NODES_C for t in TECHS));
#@constraint(m, Cop == sum(g[n,tech_ref[t],t,τ]*tech_op[t] for n in NODES_C for t in TECHS for τ in TIME));
@constraint(m, Ctrans == sum(f[i,j,pr,τ]*D[i,j]*prod_trans[pr] for i in NODES for j in NODES for pr in PRODS for τ in TIME));
@constraint(m, Revenue == sum(d[n,pr,τ]*sum(dem_price[dem]*(dem_node[dem] == n)*(dem_prod[dem] == pr) for dem in DEMS) for n in NODES for pr in PRODS for τ in TIME) - sum(s[n,pr,τ]*sum(sup_price[sup]*(sup_node[sup] == n)*(sup_prod[sup] == pr) for sup in SUPS) for n in NODES for pr in PRODS for τ in TIME));

# estimate nutrient release
@variable(m, TP[NODES, TIME] >= 0);  # net P release at each node at each time
@variable(m, TN[NODES, TIME] >= 0);  # net N release at each node at each time
@variable(m, ferP[NODES,TIME] >= 0); # P supplemented by fertilizers
@variable(m, ferN[NODES,TIME] >= 0); # N supplemented by fertilizers
@variable(m, sumferP >= 0);
@variable(m, sumferN >= 0);
@variable(m, sumTP);            # total net P release in the area
@variable(m, sumTN);            # total net N release in the area
@constraint(m, P_estimator[n in NODES, τ in TIME], TP[n,τ] >= ferP[n,τ] + sum(prod_P[pr]*d[n,pr,τ] for pr in PRODS) - Plimit[n,τ]);
@constraint(m, N_estimator[n in NODES, τ in TIME], TN[n,τ] >= ferN[n,τ] + sum(prod_N[pr]*d[n,pr,τ] for pr in PRODS) - Nlimit[n,τ]);
@constraint(m, [n in NODES, t in TIME], ferP[n,t] + sum(prod_P[pr]*d[n,pr,t] for pr in PRODS) - Plimit[n,t] >= 0);
@constraint(m, [n in NODES, t in TIME], ferN[n,t] + sum(prod_N[pr]*d[n,pr,t] for pr in PRODS) - Nlimit[n,t] >= 0);
@constraint(m, sumTP == sum(TP[n,τ] for n in NODES for τ in TIME));
@constraint(m, sumTN == sum(TN[n,τ] for n in NODES for τ in TIME));
@constraint(m, sumferP == sum(ferP[n,τ] for n in NODES for τ in TIME));
@constraint(m, sumferN == sum(ferN[n,τ] for n in NODES for τ in TIME));

# objective function (weights can be flexible)
λ1 = 1;
λ2 = 33.07*1000; # apoorva $/tonne
λ3 = 4.54*1000;  # traci equivalent
@objective(m, Min, λ1*(Cinv + Cop + Ctrans - Revenue + 1.34*1000*sumferP + 0.86*1000*sumferN) + λ2*sumTP + λ3*sumTN);
