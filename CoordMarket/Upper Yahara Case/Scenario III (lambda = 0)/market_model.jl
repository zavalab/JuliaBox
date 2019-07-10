## Market coordination case study
## Market analysis of the rock river area

using JuMP
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
#flow_matrix = readdlm("flow_matrix.csv",',');
#arc_matrix = readdlm("arc_matrix.csv",',');

# setting constants
R = 6335.439 # used for distance calculation
M = 1e20     # big M
λ = 1e-8      # penalty coefficient ($/tonne P)

# nodes, products, customers, suppliers and technologies
NODES = node_matrix[:,1]; # all nodes
PRODS = prod_matrix[:,1] # all products
DEMS  = demand_matrix[:,1] # all demands
SUPS  = supply_matrix[:,1] # all supply matrix
TECHS = technology_matrix[:,1] # all technologoes
#ARCS  = arc_matrix[:,1] # all arcs
TECH_PRVD = site_matrix[:,1] # all technology providers

#arc_matrix = Matrix(length(NODES)^2*length(PRODS),6);
#n = 1;
#for pr in PRODS
#    for i in NODES
#        for j in NODES
#            arc_matrix[n,1] = "l$n";
#            arc_matrix[n,2] = i;
#            arc_matrix[n,3] = j;
#            arc_matrix[n,4] = pr;
#            arc_matrix[n,5] = 1e20;
#            arc_matrix[n,6] = 0;
#            n = n + 1;
#        end
#    end
#end
#ARCS  = arc_matrix[:,1]; # all arcs

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
P_limit = Dict(zip(NODES, Plimit_matrix[:,2]));

# arc properties
#arc_send = Dict(zip(ARCS, arc_matrix[:,2])); # sending node of each arc
#arc_recv = Dict(zip(ARCS, arc_matrix[:,3])); # recieving node of each arc
#arc_prod = Dict(zip(ARCS, arc_matrix[:,4])); # product of each arc
#arc_cap  = Dict(zip(ARCS, arc_matrix[:,5])); # capacity of each arc
#arc_dist = Dict(zip(ARCS, arc_matrix[:,6])); # either offer distance, or calculate distance
#for l in ARCS
#    arc_dist[l] = 2*R*asin(sqrt(sin((node_lat[arc_recv[l]] - node_lat[arc_send[l]])*pi/2/180)^2
#            + cos(node_lat[arc_recv[l]]*pi/180)*cos(node_lat[arc_send[l]]*pi/180)*sin((node_long[arc_recv[l]]
#                    - node_long[arc_send[l]])*pi/2/180)^2));
#end

# technology provider properties
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

# technology sites
#y = Dict((TECHS[1],NODES[1]) => 0.5);
#for t in TECHS
#    for i in NODES
#        y[(t,i)] = 0;
#    end
#end

#total_site,~ = size(site_matrix);
#for i in 1:total_site
#    y[(site_matrix[i,4],site_matrix[i,1])] = 1;
#end

# flow setting
#ifflow = Dict((NODES[1],NODES[2]) => 0.5)
#for i in NODES
    #for j in NODES
        #ifflow[(i,j)] = flow_matrix[i,j];
    #end
#end

# modeling
m = Model(solver = GurobiSolver(Threads = 1));

# flows
#@variable(m, edge[ARCS] >= 0);
@variable(m, f[NODES,NODES,PRODS]>= 0);

# demand and supply
@variable(m, dem[DEMS] >= 0);
@variable(m, d[NODES,PRODS] >= 0);
@variable(m, sup[SUPS] >= 0);
@variable(m, s[NODES,PRODS] >= 0);
@variable(m, t[NODES] >= 0);
@variable(m, Pn[NODES] >= 0);

@constraint(m, [n in NODES], Pn[n] == sum(d[n,pr]*prod_Pemit[pr] for pr in PRODS));
@constraint(m, [n in NODES], t[n] >= Pn[n] - P_limit[n]);

# technology site
@variable(m, y[TECHS, NODES]);

# generated/consumed amount by technologies
@variable(m, x[NODES,PRODS,TECHS]);
@variable(m, p[NODES, PRODS]);

@constraint(m, techflow[i in NODES, t in TECHS], x[i,tech_refprod[t],t] <= 0);

# transportation cost and operational cost
@variable(m,transcost);
@variable(m,opcost);
@variable(m,demrevn);
@variable(m,supcost);

# social welfare
@variable(m, swf)

# assign edge to flows
#@constraint(m, edgeeq[i in NODES, j in NODES, pr in PRODS], f[i,j,pr] == sum(edge[l] for l in ARCS if arc_send[l]==i &&
#                arc_recv[l] == j && arc_prod[l] == pr));

# technology site
@constraint(m, siteeq[i in NODES, t in TECHS], y[t,i] == sum(tp_indicator[tp] for tp in TECH_PRVD if tp_site[tp] == i &&
                tp_tech[tp] == t));

# demand and supply
@constraint(m, demeq[n in NODES, pr in PRODS], d[n,pr] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr
                && dem_node[dd]==n));
@constraint(m, supeq[n in NODES, pr in PRODS], s[n,pr] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr
                && sup_node[ss]==n));
# bound flows
#@constraint(m, flowbd[i in NODES, j in NODES, pr in PRODS], f[i,j,pr] <= ifflow[(i,j)]*M) #M for unbounded flows

# balance and conversion constraints
@constraint(m, balance[i in NODES, pr in PRODS],s[i,pr]+p[i,pr]+sum(f[j,i,pr] for j in NODES) ==
                                    sum(f[i,j,pr]  for j in NODES)+d[i,pr]);
@constraint(m, process[i in NODES, pr in PRODS], p[i,pr] == sum(x[i,pr,t] for t in TECHS));
@constraint(m, transfer_pr[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] ==
                                    transfer[t,pr]/transfer[t,tech_refprod[t]]*x[i,tech_refprod[t],t]);

# technology capacity constriants
@constraint(m, techonofflb[i in NODES, t in TECHS], x[i,tech_refprod[t],t] >= -y[t,i]*tech_cap[t]);
@constraint(m, techonoffub[i in NODES, t in TECHS], x[i,tech_refprod[t],t] <= +y[t,i]*tech_cap[t]);

# demand capacity constraints
@constraint(m, demand_capacity[i in DEMS], dem[i] <= dem_cap[i]);

# supply capacity constraints
@constraint(m, supply_capacity[i in SUPS], sup[i] <= sup_cap[i]);

# arc capacityt constraints
#@constraint(m, edge_capacity[l in ARCS], edge[l] <= arc_cap[l]);

##Objective
@constraint(m, opcost == - sum(x[i,tech_refprod[t],t]*tech_bid[t] for i in NODES for t in TECHS));
@constraint(m, transcost == sum(prod_trans[pr]*distance[i,j]*f[i,j,pr] for i in NODES for j in NODES for pr in PRODS));
#@constraint(m, transcost == sum(prod_trans[arc_prod[l]]*arc_dist[l]*edge[l] for l in ARCS));
@constraint(m, demrevn == sum(dem[i]*dem_bid[i] for i in DEMS));
@constraint(m, supcost == sum(sup[i]*sup_bid[i] for i in SUPS));
@constraint(m, swf == demrevn - supcost - opcost - transcost - λ*(sum(t)-t["n1371"]-t["n1372"]));

@objective(m, Max, swf)
