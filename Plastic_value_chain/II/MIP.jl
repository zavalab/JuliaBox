using JuMP
using Gurobi
using DelimitedFiles
using CSV
using DataFrames
# import, define and pretreat data

node_matrix = CSV.read("consumer//node_matrix.csv", DataFrame)

prod_matrix = CSV.read("consumer//product_matrix.csv", DataFrame)
demand_matrix = CSV.read("consumer//demand_matrix.csv", DataFrame)
supply_matrix = CSV.read("consumer//supply_matrix.csv", DataFrame)
technology_matrix = CSV.read("consumer//technology_matrix.csv", DataFrame)


alpha_matrix = readdlm("consumer//alpha_matrix.csv",',');
site_matrix = CSV.read("consumer//site_matrix.csv", DataFrame)


# setting constants
R = 6335.439 # used for distance calculation

# nodes, products, customers, suppliers and technologies
NODES = node_matrix[:,1]; # all nodes
PRODS = prod_matrix[:,1] # all products
DEMS  = demand_matrix[:,1] # all demands
SUPS  = supply_matrix[:,1] # all supply
TECHS = technology_matrix[:,1] # all technologoes
#ARCS  = arc_matrix[:,1] # all arcs
TECH_PRVD = site_matrix[:,1] # all technology providers

# node properties
node_alia = Dict(zip(NODES, node_matrix[:,2])); # node alias
node_lat = Dict(zip(NODES, node_matrix[:,3])); # node longitude
node_long  = Dict(zip(NODES, node_matrix[:,4])); # node latitude

# product properties
prod_name = Dict(zip(PRODS, prod_matrix[:,2])); # product names
prod_trans_vc_truck= Dict(zip(PRODS, prod_matrix[:,3])); # product transportation bids
prod_trans_fc_truck= Dict(zip(PRODS, prod_matrix[:,4]));



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
tech_cap  = Dict(zip(TECHS, technology_matrix[:,2]));
tech_inv  = Dict(zip(TECHS, technology_matrix[:,4])); # technology investment costs (not useful here)
tech_bid  = Dict(zip(TECHS, technology_matrix[:,5])); # technology operational costs (bids) per unit reference prod
tech_refprod = Dict(zip(TECHS, technology_matrix[:,3])); # technology reference products



# technology properties
tech_cap  = Dict(zip(TECHS, technology_matrix[:,2])); # technology capacities with regard to the reference product

tech_size_1  = Dict(zip(TECHS, technology_matrix[:,7]));
tech_size_2  = Dict(zip(TECHS, technology_matrix[:,8]));
tech_size_3  = Dict(zip(TECHS, technology_matrix[:,9]));


tech_cost_1  = Dict(zip(TECHS, technology_matrix[:,10]));
tech_cost_2  = Dict(zip(TECHS, technology_matrix[:,11]));
tech_cost_3  = Dict(zip(TECHS, technology_matrix[:,12]));


tp_site = Dict(zip(TECH_PRVD, site_matrix[:,2])); # node location of the technology provider
tp_tech = Dict(zip(TECH_PRVD, site_matrix[:,5])); # technology type that the provider can provide
tp_indicator = Dict(zip(TECH_PRVD, ones(length(TECH_PRVD)))); # technology indicator

# define two-key dictionaries
# transformation factors
transfer = Dict((TECHS[1],PRODS[1]) => 0.5);
for i in 1:length(TECHS)
    for k in 1: length(PRODS)
        transfer[(TECHS[i], PRODS[k])] = alpha_matrix[i,k];
    end
end

# distance between nodes (using the Haversine formula)
distance = Dict((NODES[1], NODES[2]) => 0.5);
for i in NODES
    for j in NODES
    distance[(i,j)] = 2*R*asin(sqrt(sin((node_lat[j] - node_lat[i])*pi/2/180)^2
            + cos(node_lat[j]*pi/180)*cos(node_lat[i]*pi/180)*sin((node_long[j]
                    - node_long[i])*pi/2/180)^2));
    end
end
######################################################################
C=zeros(length(NODES), length(NODES),length(PRODS))

for i in NODES
    for j in NODES
    if supply_matrix[j,5]>=10000
    C[i,j,2]=1
    end
end
end

for i in NODES
    for j in NODES
        for k in 3:4
    if supply_matrix[j,5]>=10000 && supply_matrix[i,5]>=10000
    C[i,j,k]=1
    end
end
end
end

M=100000000

af=0.117

m = Model(with_optimizer(Gurobi.Optimizer))
@variable(m, f[i in NODES,j in NODES,pr in PRODS;C[i,j,pr]==1]>= 0);

# demand and supply
@variable(m, dem[DEMS] >= 0);
@variable(m, d[NODES,PRODS] >= 0);
@variable(m, sup[SUPS] >= 0);
@variable(m, s[NODES,PRODS] >= 0);
@variable(m, t[NODES] >= 0);

@variable(m, y[TECHS, NODES]);

# generated/consumed amount by technologies
@variable(m, x[NODES,PRODS,TECHS]);
@variable(m, p[NODES, PRODS]);

@constraint(m, techflow[i in NODES, t in TECHS], x[i,tech_refprod[t],t] <= 0);

# transportation cost and operational cost


@variable(m,transcost);
@variable(m,opcost);
@variable(m,capcost);
@variable(m,demrevn);
@variable(m,supcost);

@variable(m, swf)

# technology site
@constraint(m, siteeq[i in NODES, t in TECHS], y[t,i] == sum(tp_indicator[tp] for tp in TECH_PRVD if tp_site[tp] == i &&
                tp_tech[tp] == t));

# demand and supply
@constraint(m, demeq[n in NODES, pr in PRODS], d[n,pr] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr
                && dem_node[dd]==n));
@constraint(m, supeq[n in NODES, pr in PRODS], s[n,pr] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr
                && sup_node[ss]==n));


# balance and conversion constraints
@constraint(m, balance[i in NODES, pr in PRODS],s[i,pr]+p[i,pr]+sum(f[j,i,pr] for j in NODES if C[j,i,pr]==1 ) ==
                                    sum(f[i,j,pr]  for j in NODES if C[i,j,pr]==1 )+d[i,pr]);
@constraint(m, process[i in NODES, pr in PRODS], p[i,pr] == sum(x[i,pr,t] for t in TECHS));
@constraint(m, transfer_pr[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] ==
                                    transfer[t,pr]/transfer[t,tech_refprod[t]]*x[i,tech_refprod[t],t]);




# demand capacity constraints
@constraint(m, demand_capacity[i in DEMS], dem[i] <= dem_cap[i]);

# supply capacity constraints
#*******************************
@constraint(m, supply_capacity[i in SUPS],  sup[i] == sup_cap[i]);
#*******************************

@variable(m, z[NODES,TECHS,1:3], Int);
@constraint(m, [i in NODES, t in TECHS, k in 1:3],  z[i,t,k] <=5);
@constraint(m, [i in NODES, t in TECHS, k in 1:3],  z[i,t,k] >=0);
#@constraint(m, [i in NODES, t in TECHS],  sum(z[i,t,k] for k in 1:3)<=5*y[t,i]);
#@constraint(m, [i in NODES, t in TECHS],  sum(z[i,t,k] for k in 1:3)<=5*y[t,i]);
@constraint(m, [i in NODES, t in TECHS], - x[i,tech_refprod[t],t]  <=   z[i,t,1]*tech_size_1[t]
                                                                         +z[i,t,2]*tech_size_2[t]
                                                                         +z[i,t,3]*tech_size_3[t]);
@constraint(m, capcost ==   sum( z[i,t,1]*tech_cost_1[t] for i in NODES for t in TECHS )
                            + sum( z[i,t,2]*tech_cost_2[t] for i in NODES for t in TECHS )
                            + sum( z[i,t,3]*tech_cost_3[t] for i in NODES for t in TECHS ));
##Objective
@constraint(m, opcost == -sum(x[i,tech_refprod[t],t]*tech_bid[t] for i in NODES for t in TECHS));
@constraint(m, transcost == sum(prod_trans_vc_truck[pr]*distance[i,j]*f[i,j,pr] + f[i,j,pr]*prod_trans_fc_truck[pr] for i in NODES for j in NODES for pr in PRODS if C[i,j,pr]==1));


@constraint(m, demrevn == sum(dem[i]*dem_bid[i] for i in DEMS));
@constraint(m, supcost == sum(sup[i]*sup_bid[i] for i in SUPS));
@constraint(m, swf == demrevn - supcost - opcost - transcost - af*capcost);

@objective(m, Max, swf)
optimize!(m)

