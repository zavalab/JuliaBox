# P Recovery Model
# Victor M Zavala, 2016, UW-Madison
# Apoorva M Sampat, 2016, UW-Madison

using JuMP;
using Gurobi;

m = Model(solver=GurobiSolver(Threads = 2, MIPGap = 0.02))
#, MIPGap = 0.0002

#Importing Data
technology_matrix 	= 	readdlm("./InputData/technology_matrix.csv"	,',')
node_matrix 		= 	readdlm("./InputData/node_matrix.csv"	,',')
product_matrix 		= 	readdlm("./InputData/product_matrix.csv"	,',')
supply_matrix 		= 	readdlm("./InputData/supply_matrix.csv"	,',')
demand_matrix		=	readdlm("./InputData/demand_matrix.csv"	,',')
alpha_matrix 		= 	readdlm("./InputData/alpha_matrix.csv"	,',')
supply_values		= 	readdlm("./InputData/supply_values.csv"	,',')
node_distance 		=	readdlm("./InputData/node_distance.csv", ',')
node_dist_correc	=	readdlm("./InputData/node_dist.csv", ',')

# Sets
# Imporatant that these are defined in ascending order
TECHS = technology_matrix[:,1]		# set of technologies
NODES = node_matrix[:,1]		# set of nodes
PRODS = product_matrix[:,1]		# set of products
SUPS  = supply_matrix[:,1]		# set of supplies
DEMS  = demand_matrix[:,1]		# set of demands
SUPS_VALUES = supply_values[:,1]

NODES_dist = node_dist_correc[:,1]

# Parameters 	(no need to define separately for julia model)
tech_cap	= 	Dict(zip(TECHS, technology_matrix[:, 3]))             	 # technology capacity
tech_alias 	=	Dict(zip(TECHS, technology_matrix[:, 2])) 		 # technology name alias
tech_invcost 	= 	Dict(zip(TECHS, technology_matrix[:, 4]))	         # technology investment cost
tech_refprod 	= 	Dict(zip(TECHS, technology_matrix[:, 5]))
node_lat	= 	Dict(zip(NODES, node_matrix[:, 3]))		 # node latitude
node_long	= 	Dict(zip(NODES, node_matrix[:, 4]))   	 	 # node longitude
node_alias 	= 	Dict(zip(NODES, node_matrix[:, 2]))     	 # node alias name
prod_alias 	= 	Dict(zip(PRODS, product_matrix[:, 2]))     	 # product alias name
prod_transcost 	= 	Dict(zip(PRODS, product_matrix[:, 3]))           # product transportation cost
sup_node 	= 	Dict(zip(SUPS, supply_matrix[:, 2]))        	 # supply node
sup_prod 	= 	Dict(zip(SUPS, supply_matrix[:, 3]))        	 # supply product
sup_cap 	= 	Dict(zip(SUPS, supply_matrix[:, 4])) 		 # supply flow capacity
sup_price 	= 	Dict(zip(SUPS, supply_matrix[:, 5]))        	 # supply price
dem_node 	= 	Dict(zip(DEMS, demand_matrix[:, 2])) 		 # demand node
dem_prod 	= 	Dict(zip(DEMS, demand_matrix[:, 3]))        	 # demand product
dem_cap 	= 	Dict(zip(DEMS, demand_matrix[:, 4]))          	 # demand flow capacity
dem_price 	= 	Dict(zip(DEMS, demand_matrix[:, 5]))    	 # demand price
sup_value	=	Dict(zip(SUPS_VALUES, supply_values[:, 2]))
tech_opcost	=	Dict(zip(TECHS, technology_matrix[:, 6]))	# technology operating cost
p_score 	= 	Dict(zip(NODES, node_matrix[:, 5]))
prod_p		=	Dict(zip(PRODS, product_matrix[:, 4]))		# P content of each product
prod_num 	= 	Dict(zip(product_matrix[:, 2], PRODS)) 		# Helps in calling results by the product name, instead of prod no.
prod_units	=	Dict(zip(PRODS, product_matrix[:, 6]))
node_dist	=	node_distance[:, 8]

sup_tot		= 	sum(supply_values[:, 2])
# Emissions Metric
co2_per_km	=	0.1 						 # 0.1g of CO2 emitted per km per kg of freight

# Unit conversions
kW_to_MW = 1e-3
rec_value = 500 # Usd per MWh

# TradeOff Analysis Parameters
cost_min	=	0
cost_max	=	5000.000000000002			# Obtained without introducing epsilon constraint
budget		=	190e6				# Max value corressponding to siting the most expensive tech at all nodes
epsilon		=	"combined_5000AU"



# CVaR Parameters
#alpha = 0.8

# Haversine Formula Parameters
R = 6335.439

## Defining Two Variable Dictionaries ##
transfer = Dict(("t1","p1") => 0.5) 				 #Just used as an initiator to set up the dictionary

for i in 1: length(TECHS)
	for j in 1: length(PRODS)
		transfer[(TECHS[i], PRODS[j])] = alpha_matrix[i, j]
		j = j +1
	end
i = i +1
end

distance = Dict(("n1", "n2") => 1.1)

# Creating the distance dictionary, all distances are in km
counter = 1
for i in NODES_dist
       for j in NODES_dist
	       distance[(i,j)] = node_dist[counter]
	       counter = counter + 1
       end
end

M  = 1e15;                        			 	# bigM
Mflow = 10e6;

# variables
@variable(m, flow[NODES,NODES,PRODS] >= 0) 	# product flow (snd,rec,prod)
@variable(m, sup[SUPS] >= 0) 			# supply flow
@variable(m, dem[DEMS] >=0)
@variable(m, y[TECHS,NODES])			# Used for solving the LP
#@variable(m, y[TECHS,NODES], Bin)
@variable(m, flowin[NODES,PRODS] >=0)
@variable(m, flowout[NODES,PRODS] >=0)
@variable(m, suptot[NODES,PRODS] >=0)
@variable(m, demtot[NODES,PRODS] >=0)
@variable(m, gentot[NODES,PRODS,TECHS])

# Basic Variables
@variable(m, swf)
@variable(m, invcost >= 0)
@variable(m, transcost >= 0)
@variable(m, obj)
@variable(m, opcost >= 0)

@variable(m, co2_total >= 0)

# CVaR Variables
#@variable(m, dissatisfaction[NODES])
#@variable(m, phi[NODES] >= 0 )
#@variable(m, cvar)
#@variable(m, risk)
#@variable(m, eta>=0)

@variable(m, processed[NODES, PRODS] >= 0)
@variable(m, unprocessed[NODES, PRODS] >= 0)

@variable(m, umanure >= 0)
#@variable(m, p_tot >= 0)
@variable(m, prod_revenue )
@variable(m, demtot_price[NODES, PRODS])

@variable(m, rec_revenue)
@variable(m, p_credit)
@variable(m, rin_revenue)

#@variable(m, roi)

# Capacity Constraints for the variables
#@constraint(m, supply_capacity[i = SUPS], sup[i] <= sup_cap[i])
#@constraint(m, demand_capacity[i = DEMS], dem[i] <= dem_cap[i])

# variable in objective function
@constraint(m, swf == sum(sup_price[s]*sup[s] for s in SUPS) - sum(dem_price[d]*dem[d] for d in DEMS))
@constraint(m, invcost == sum(tech_invcost[t]*y[t,n] for t in TECHS for n in NODES));
@constraint(m, transcost == sum(prod_transcost[p]*flow[n,nn,p]*distance[n,nn] for n in NODES for nn in NODES for p in PRODS));

# Objective Function
#@constraint(m, obj ==  (invcost/20/365) + transcost)
#@constraint(m, obj ==  cost_total) 		# Obj2
#@constraint(m, obj ==  -struvite_total)
#@constraint(m, obj == cvar)
#@constraint(m,obj==eta)
#@constraint(m, obj == -sum{gentot[n,"p2",t], t in TECHS, n in NODES})
#@constraint(m, obj == sum{dissatisfaction[n], n in NODES})
# @constraint(m, objcons, obj == umanure)		# Obj1
#@constraint(m, obj_p_cons, obj == -prod_revenue + (transcost + opcost/365))
@constraint(m, obj_p_cons, obj == -prod_revenue + (transcost + opcost/365) + invcost/20/365 -p_credit)
@objective(m, Min, obj)


# Product balances
@constraint(m, flowineq[n in NODES, p in PRODS], flowin[n,p] == sum(flow[nn,n,p] for nn in NODES))
@constraint(m, flowouteq[n in NODES, p in PRODS], flowout[n,p] == sum(flow[n,nn,p] for nn in NODES));

@constraint(m, supeq[n in NODES, p in PRODS], suptot[n,p] == sum(sup[s] for s in SUPS if sup_prod[s] == p && sup_node[s] == n));
@constraint(m, demeq[n in NODES, p in PRODS], demtot[n,p] == sum(dem[d] for d in DEMS if dem_prod[d]==p && dem_node[d]==n)) ;

@constraint(m, geneqlb[n in NODES, p in PRODS, t in TECHS], gentot[n,p,t] - (processed[n,tech_refprod[t]] )*transfer[t,p] >= -(1 - y[t,n])*Mflow);
@constraint(m, genequb[n in NODES, p in PRODS, t in TECHS], gentot[n,p,t] - (processed[n,tech_refprod[t]] )*transfer[t,p] <= +(1-y[t,n])*Mflow);

@constraint(m, floweq[n in NODES, p in PRODS],     +  flowin[n,p]
					           - flowout[n,p]
           	        		           +  sum(gentot[n,p,t] for t in TECHS)
			                           +    suptot[n,p]
			                           -    demtot[n,p] == 0) ;

# Adding Split Equations #
@constraint(m, spliteq[n in NODES, p in PRODS], processed[n, p] + unprocessed[n, p] == flowin[n, p] + suptot[n, p]);

# techonology capacities
#techcapeq{n in NODES, t in TECHS}:
@constraint(m, techcapeq[n in NODES, t in TECHS], processed[n,tech_refprod[t]]  <= y[t,n]*tech_cap[t] + (1-y[t,n])*Mflow)

# logic constraint (at most one technology per node)
@constraint(m, onetecheq2[n in NODES], sum(y[t,n] for t in TECHS) <= 1)

# logic constraint (if no technology installed)

@constraint(m, techonofflb[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] >= -y[t,n]*tech_cap[t])
@constraint(m, techonoffub[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] <= +y[t,n]*tech_cap[t])

# Assigning 0 value to processed variable (if no technology installed)

@constraint(m, processedlb[n in NODES, p in PRODS], processed[n,p] >= -sum(y[t,n]*tech_cap[t] for t in TECHS))
@constraint(m, processedub[n in NODES, p in PRODS], processed[n,p] <= +sum(y[t,n]*tech_cap[t] for t in TECHS))

#Eliminating Self-Flows
@constraint(m, fix[n in NODES,p in PRODS], flow[n,n,p] == 0)

## Fixing Supply Values
@constraint(m, fix_supply_values[i in SUPS_VALUES], sup[i] == sup_value[i])

@constraint(m, co2_total == co2_per_km*sum(flow[n,nn,p]*distance[n,nn] for n in NODES for nn in NODES for p in PRODS));
#@constraint(m, struvite_total == sum{demtot[n,"p2"], n in NODES})

#@constraint(m, transcost ==  epsilon*cost_max)
#@constraint(m, cost_total <= epsilon*cost_max)
@expression(m, daily_cost,  (invcost/20/365) + transcost + opcost )
#@constraint(m, invcost <= epsilon*budget)

#Defining unprocessed manure varaible
@constraint(m, unpro_manure_cons, umanure == sum(demtot[n, "p1"] for n in NODES))

# Adding the constraint for operating cost
@constraint(m, opcost_cons, opcost == sum(tech_opcost[t]*y[t, n] for t in TECHS for n in NODES))		# Annual operating cost


@expression(m, dist_flow[p in PRODS], sum(distance[n,nn]*flow[n,nn,p] for n in NODES for nn in NODES))

@expression(m, tot_flow[p in PRODS], sum(flow[n,nn,p] for n in NODES for nn in NODES))
@expression(m, prod_recovered[p in PRODS], sum(demtot[n,p] for n in NODES))

@expression(m, p_recovered[p in PRODS], prod_p[p]*prod_recovered[p])

### Introducing beta_scores/geographical P scores ranging from 1 to 10
beta_score = Dict(("d1") => 1.0)

for i in 1: length(DEMS)
       if dem_price[DEMS[i]] >= 0
       		beta_score[DEMS[i]] = 1
       else
       		beta_score[DEMS[i]] = p_score[dem_node[DEMS[i]]]/10
       end
       i = i +1
end

# Total product revenue
@constraint(m, prod_revenue == sum(dem[d]*dem_price[d]*beta_score[d] for d in DEMS) ) # including the geogrpahical p concentration

# Defining expressions for reporting the solution
#@expression(m, roi, (prod_revenue - (invcost/365/20) - (opcost) - transcost )*365/(invcost)*100)

## Incentives for renewable energy


# Renewable energy credits (REC)
# Renewable energy credits (REC)
# Assuming value of USD 1/MWh
#@constraint(m, rec_cons, rec_revenue == prod_recovered[prod_num["Electricity"]]*rec_value*kW_to_MW)
# @expression(m, rec_revenue[n in NODES], rec_value*rec[n])

# Nutrient Tax Credits (S.988 Scenario)
# Assuming Virginia's Value of USD 10.10 / lb P = USD 22.37 /kg P
p_credit_value = 22.04
@constraint(m, p_credit_cons, p_credit == p_credit_value*sum(p_recovered[p] for p in [prod_num["Cake1"], prod_num["Cake2"], prod_num["Cake3"], prod_num["Struvite"], prod_num["Struvite+Solids"]]) )
epsilon		=	"p_credit_$(p_credit_value)"

# Electricity to the grid
