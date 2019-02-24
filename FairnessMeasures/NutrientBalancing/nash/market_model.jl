## Market coordination case study
## Market analysis of the rock river area

using JuMP
using Ipopt
using Gurobi

# import, define and pretreat data

node_matrix = readdlm("node_matrix.csv",',');
prod_matrix = readdlm("product_matrix.csv",',');
demand_matrix = readdlm("demand_matrix.csv",',');
supply_matrix = readdlm("supply_matrix.csv",',');
technology_matrix = readdlm("technology_matrix.csv",',');
Plimit_matrix = readdlm("Plimit_matrix.csv",',');

alpha_matrix = readdlm("alpha_matrix.csv",',');
site_matrix = readdlm("site_matrix.csv",',');

# setting constants
R = 6335.439 # used for distance calculation
λ = 0.001      # penalty coefficient ($/tonne P)

# nodes, products, customers, suppliers and technologies
NODES = node_matrix[:,1]; # all nodes
PNODES = Plimit_matrix[:, 1] # nodes accepting P
PRODS = prod_matrix[:,1] # all products
DEMS  = demand_matrix[:,1] # all demands
SUPS  = supply_matrix[:,1] # all supply matrix
TECHS = technology_matrix[:,1] # all technologoes
TECH_PRVD = site_matrix[:,1] # all technology providers

# node properties
node_alia = Dict(zip(NODES, node_matrix[:,2])); # node alias
node_lat = Dict(zip(NODES, node_matrix[:,3])); # node longitude
node_long  = Dict(zip(NODES, node_matrix[:,4])); # node latitude

# product properties
prod_name = Dict(zip(PRODS, prod_matrix[:,2])); # product names
prod_trans= Dict(zip(PRODS, prod_matrix[:,3])); # product transportation bids
prod_Pemit = Dict(zip(PRODS, prod_matrix[:,4])); # emission coefficient of phosphorous of each product

# custormer properties
dem_node  = Dict(zip(DEMS, demand_matrix[:,2])); # demand locations
dem_prod  = Dict(zip(DEMS, demand_matrix[:,3])); # demand products
dem_bid   = Dict(zip(DEMS, demand_matrix[:,4])); # demand bids
dem_cap   = Dict(zip(DEMS, demand_matrix[:,5])); # demand capacities

# supplier properties
sup_node  = Dict(zip(SUPS, supply_matrix[:,2])); # supply locations
sup_prod  = Dict(zip(SUPS, supply_matrix[:,3])); # supply products
sup_bid   = Dict(zip(SUPS, supply_matrix[:,4])); # supply bids
sup_cap   = Dict(zip(SUPS, supply_matrix[:,5])); # supply capacities

# technology properties
tech_cap  = Dict(zip(TECHS, technology_matrix[:,2])); # technology capacities with regard to the reference product
tech_inv  = Dict(zip(TECHS, technology_matrix[:,4])); # technology investment costs (not useful here)
tech_bid  = Dict(zip(TECHS, technology_matrix[:,5])); # technology operational costs (bids) per unit reference prod
tech_refprod = Dict(zip(TECHS, technology_matrix[:,3])); # technology reference products

# phosphorous limit in each nodes
P_limit = Dict(zip(PNODES, Plimit_matrix[:,2]));

# distance between nodes (using the Haversine formula)
distance = Dict((NODES[1], NODES[2]) => 0.5);
for i in NODES
    for j in NODES
    distance[(i,j)] = 2*R*asin(sqrt(sin((node_lat[j] - node_lat[i])*pi/2/180)^2
            + cos(node_lat[j]*pi/180)*cos(node_lat[i]*pi/180)*sin((node_long[j]
                    - node_long[i])*pi/2/180)^2));
    end
end

# modeling
#m = Model(solver=GurobiSolver(Threads=3))
m = Model(solver=IpoptSolver(mu_strategy="adaptive", linear_solver="ma57"))

# flows
@variable(m, f[NODES,NODES,PRODS]>= 0, start = 100);

# demand and supply
@variable(m, dem[DEMS] >= 0, start = 1);
@variable(m, d[NODES,PRODS] >= 0, start = 10);
@variable(m, sup[SUPS] >= 0, start = 500);
@variable(m, s[NODES,PRODS] >= 0, start = 500);
@variable(m, t[PNODES] >= 1e-4, start = 0.5);
@variable(m, Pn[NODES] >= 0, start = 0.5);

# social welfare
@variable(m, swf, start = 10)
@variable(m, nash_obj, start = 10)
#@variable(m, transcost >= 0, start = 1e5)

@constraint(m, [n in PNODES], Pn[n] == sum(d[n,pr]*prod_Pemit[pr] for pr in PRODS));
#@constraint(m, [n in NODES], t[n] >= Pn[n] - P_limit[n]);


@constraint(m, [n in PNODES], t[n] <=  (P_limit[n] - Pn[n])) # using the epigraph formulation

# demand and supply
@constraint(m, demeq[n in NODES, pr in PRODS], d[n,pr] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr
                && dem_node[dd]==n));
@constraint(m, supeq[n in NODES, pr in PRODS], s[n,pr] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr
                && sup_node[ss]==n));

# Eliminating self flows
for n in NODES
	for p in PRODS
		setlowerbound(f[n,n,p], 0.0)
		setupperbound(f[n,n,p], 0.0)
	end
end

# balance and conversion constraints
@constraint(m, balance[i in NODES, pr in PRODS],s[i,pr]+sum(f[j,i,pr] for j in NODES) == sum(f[i,j,pr]  for j in NODES)+d[i,pr]);

# demand capacity constraints
@constraint(m, demand_capacity[i in DEMS], dem[i] <= dem_cap[i]);

# supply capacity constraints
@constraint(m, supply_capacity[i in SUPS], sup[i] == sup_cap[i]);
#@constraint(m, transcost == sum(prod_trans[pr]*distance[i,j]*f[i,j,pr] for i in NODES for j in NODES for pr in PRODS))
#@constraint(m, transcost <= 1e7)

@NLconstraint(m, nash_obj_cons, nash_obj <= sum(log(t[n]) for n in PNODES) )

@objective(m, Max, nash_obj)
