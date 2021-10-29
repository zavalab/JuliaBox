
using JuMP
using Gurobi
using DelimitedFiles
using CSV
using Distributions
using DataFrames
# import, define and pretreat data



node_matrix = CSV.read("time_test//node_matrix.csv", DataFrame)
prod_matrix = CSV.read("time_test//product_matrix.csv", DataFrame)
demand_matrix = CSV.read("time_test//demand_matrix.csv", DataFrame)
supply_matrix = CSV.read("time_test//supply_matrix.csv", DataFrame)
technology_matrix = CSV.read("time_test//technology_matrix.csv", DataFrame)
alpha_matrix = readdlm("time_test//alpha_matrix.csv",',');
site_matrix = CSV.read("time_test//site_matrix.csv", DataFrame)

######################################################################

RecordLB=zeros(3,10)
N = [1000,10000,100000]

for q in 1:3

for w in 1:10
@time begin
# setting constants
R = 6335.439 # used for distance calculation
#M = 1e20     # big M


# nodes, products, customers, suppliers and technologies
NODES = node_matrix[:,1]; # all nodes
PRODS = prod_matrix[:,1] # all products
DEMS  = demand_matrix[:,1] # all demands
SUPS  = supply_matrix[:,1] # all supply matrix
TECHS = technology_matrix[:,1] # all technologoes
#ARCS  = arc_matrix[:,1] # all arcs
TECH_PRVD = site_matrix[:,1] # all technology providers

# node properties
node_alia = Dict(zip(NODES, node_matrix[:,2])); # node alias
node_lat = Dict(zip(NODES, node_matrix[:,3])); # node longitude
node_long  = Dict(zip(NODES, node_matrix[:,4])); # node latitude

# product properties
prod_name = Dict(zip(PRODS, prod_matrix[:,2])); # product names
prod_trans= Dict(zip(PRODS, prod_matrix[:,3])); # product transportation bids

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

tech_size_1  = Dict(zip(TECHS, technology_matrix[:,7]));
tech_size_2  = Dict(zip(TECHS, technology_matrix[:,8]));
tech_size_3  = Dict(zip(TECHS, technology_matrix[:,9]));

tech_cost_1  = Dict(zip(TECHS, technology_matrix[:,10]));
tech_cost_2  = Dict(zip(TECHS, technology_matrix[:,11]));
tech_cost_3  = Dict(zip(TECHS, technology_matrix[:,12]));

tech_inv  = Dict(zip(TECHS, technology_matrix[:,4])); # technology investment costs (not useful here)
tech_bid  = Dict(zip(TECHS, technology_matrix[:,5])); # technology operational costs (bids) per unit reference prod
tech_refprod = Dict(zip(TECHS, technology_matrix[:,3])); # technology reference products


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
######################################################################
#############################lowerbound###############################
#############################sample edges#############################


lb_size_of_sample=N[q]

A_LB = Vector(1:length(NODES)*length(NODES))
S_sample=sample(A_LB,lb_size_of_sample,replace=false)

C=zeros(length(NODES), length(NODES))

for k in S_sample
    i=(k-1)÷length(NODES)+1
    j=(k-1)%length(NODES)+1
    C[i,j]=1
end

CC=zeros(1372, 1372)
for i in 1:length(NODES)
    for j in 1:length(NODES)
        CC[NODES[i],NODES[j]]=C[i,j]
    end
end


M=100000000


m_l = Model(with_optimizer(Gurobi.Optimizer,MIPGap=0.003))
@variable(m_l, f[i in NODES,j in NODES,PRODS;CC[i,j]==1]>= 0);
# demand and supply
@variable(m_l, dem[DEMS] >= 0);
@variable(m_l, d[NODES,PRODS] >= 0);
@variable(m_l, sup[SUPS] >= 0);
@variable(m_l, s[NODES,PRODS] >= 0);


# generated/consumed amount by technologies
@variable(m_l, x[NODES,PRODS,TECHS]);
@variable(m_l, p[NODES, PRODS]);

@constraint(m_l, techflow[i in NODES, t in TECHS], x[i,tech_refprod[t],t] <= 0);

# transportation cost and operational cost
@variable(m_l,transcost);
@variable(m_l,opcost);
@variable(m_l,capcost);
@variable(m_l,demrevn);
@variable(m_l,supcost);

# social welfare
@variable(m_l, swf)

# demand and supply
@constraint(m_l, demeq[n in NODES, pr in PRODS], d[n,pr] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr
                && dem_node[dd]==n));
@constraint(m_l, supeq[n in NODES, pr in PRODS], s[n,pr] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr
                && sup_node[ss]==n));

# balance and conversion constraints
@constraint(m_l, balance[i in NODES, pr in PRODS],s[i,pr]+p[i,pr]+sum(f[j,i,pr] for j in NODES if CC[j,i]==1) ==
                                    sum(f[i,j,pr]  for j in NODES if CC[i,j]==1)+d[i,pr]);
@constraint(m_l, process[i in NODES, pr in PRODS], p[i,pr] == sum(x[i,pr,t] for t in TECHS));
@constraint(m_l, transfer_pr[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] ==
                                    transfer[t,pr]/transfer[t,tech_refprod[t]]*x[i,tech_refprod[t],t]);


# demand capacity constraints
@constraint(m_l, demand_capacity[i in DEMS], dem[i] <= dem_cap[i]);

# supply capacity constraints
@constraint(m_l, supply_capacity[i in SUPS], sup[i] <= sup_cap[i]);


##Objective
@constraint(m_l, opcost == - sum(x[i,tech_refprod[t],t]*tech_bid[t] for i in NODES for t in TECHS));

#modular constraint on the size of equipment

@variable(m_l, z[NODES,TECHS,1:3], Bin);
@constraint(m_l, [i in NODES, t in TECHS], - x[i,tech_refprod[t],t]  <=   z[i,t,1]*tech_size_1[t]
                                                                         +z[i,t,2]*tech_size_2[t]
                                                                         +z[i,t,3]*tech_size_3[t]);
@constraint(m_l, capcost ==   sum( z[i,t,1]*tech_cost_1[t] for i in NODES for t in TECHS )
                            + sum( z[i,t,2]*tech_cost_2[t] for i in NODES for t in TECHS )
                            + sum( z[i,t,3]*tech_cost_3[t] for i in NODES for t in TECHS ));

@constraint(m_l, transcost == sum(prod_trans[pr]*distance[i,j]*f[i,j,pr] for i in NODES for j in NODES for pr in PRODS if CC[i,j]==1));
@constraint(m_l, demrevn == sum(dem[i]*dem_bid[i] for i in DEMS if dem_node[i] in NODES));
@constraint(m_l, supcost == sum(sup[i]*sup_bid[i] for i in SUPS if sup_node[i] in NODES));
@constraint(m_l, swf == demrevn - supcost - opcost - 0.05*capcost - transcost);

@objective(m_l, Max, swf)
optimize!(m_l)


RecordLB[q,w]=JuMP.value.(swf)
end
end
end
println(RecordLB)
