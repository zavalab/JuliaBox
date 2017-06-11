# P Recovery Model
# Victor M Zavala, 2016, UW-Madison
# Apoorva M Sampat, 2016, UW-Madison

using JuMP;
using Gurobi;

m = Model(solver=GurobiSolver(Threads = 2))
#, MIPGap = 0.0002

#Importing Data
technology_matrix 	= 	readdlm("technology_matrix.csv"	,',')
node_matrix 		= 	readdlm("node_matrix.csv"	,',')
product_matrix 		= 	readdlm("product_matrix.csv"	,',')
supply_matrix 		= 	readdlm("supply_matrix.csv"	,',')
demand_matrix		=	readdlm("demand_matrix.csv"	,',')
alpha_matrix 		= 	readdlm("alpha_matrix.csv"	,',')          			 # product transfer matrix
supply_values		= 	readdlm("supply_values.csv"	,',')
#ideal_values		= 	readdlm("ideal_values.csv"	,',')
#weight_matrix		=	readdlm("weights.csv"		,',')

# Sets
# Imporatant that these are defined in ascending order
TECHS = technology_matrix[:,1]		# set of technologies
NODES = node_matrix[:,1]		# set of nodes
PRODS = product_matrix[:,1]		# set of products
SUPS  = supply_matrix[:,1]		# set of supplies
DEMS  = demand_matrix[:,1]		# set of demands
SUPS_VALUES = supply_values[:,1]
#IDEAL = ideal_values[:,1]

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
#ide_value	=	Dict(zip(IDEAL, ideal_values[:, 2]))
#weight		=	Dict(zip(NODES, weights[:,2]))
tech_opcost	=	Dict(zip(TECHS, technology_matrix[:, 6]))	# technology operating cost
p_score		=	Dict(zip(NODES, node_matrix[:, 5]))

# Emissions Metric
co2_per_km	=	0.1 						 # 0.1g of CO2 emitted per km per kg of freight

# TradeOff Analysis Parameters
cost_min	=	0
cost_max	=	5000.000000000002			# Obtained without introducing epsilon constraint
budget		=	5e5
#cost_max	=	152189.715716414			# Obtained with epsilon = 1 for the above value of cost_max
#struvite_max	=	828488.5492
#struvite_min	=	0
epsilon		=	0.009


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

#=
weight = Dict(("n1","n1") => 0.5)

for i in 1: length(NODES)
	for j in 1: length(NODES)
		weight[(NODES[i], NODES[j])] = weight_matrix[i, j+1]
		j = j +1
	end
i = i +1
end
=#
#=
ide_value = Dict(("n1","n1") => 1.0)
for i in 1: length(NODES)
	for j in 1: length(NODES)
		ide_value[(NODES[i], NODES[j])] = ideal_values[i, j+1]
		j = j +1
	end
i = i +1
end
=#

distance = Dict(("n1", "n2") => 1.1)

# Using the Haversine formula
for i in NODES
	for j in NODES
		#distance[(NODES[i], NODES[j])] = sqrt((node_matrix[i,3] - node_matrix[j,3])^2 + (node_matrix[i,4] - node_matrix[j,4])^2)
		distance[(i, j)] = 2*R*asin(sqrt(sin((node_lat[j] - node_lat[i])*pi/2/180)^2 + cos(node_lat[j]*pi/180)*cos(node_lat[i]*pi/180)*sin((node_long[j] - node_long[i])*pi/2/180)^2))
	end
end

M  = 1e15;                        			 	# bigM
Mflow = 10e6;
						
# variables
@variable(m, flow[NODES,NODES,PRODS] >= 0) 	# product flow (snd,rec,prod)
@variable(m, sup[SUPS] >= 0) 			# supply flow
@variable(m, dem[DEMS] >=0) 
#@defVar(m, y[TECHS,NODES])			# Used for solving the LP
@variable(m, y[TECHS,NODES], Bin)
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
@variable(m, struvite_total >= 0)
#@variable(m, cost_total >= 0)

# CVaR Variables
#@variable(m, dissatisfaction[NODES])
#@variable(m, phi[NODES] >= 0 )
#@variable(m, cvar)
#@variable(m, risk)
#@variable(m, eta>=0)

@variable(m, processed[NODES, PRODS] >= 0)
@variable(m, unprocessed[NODES, PRODS] >= 0)

@variable(m, umanure >= 0)

# Capacity Constraints for the variables
#@constraint(m, supply_capacity[i = SUPS], sup[i] <= sup_cap[i]) 
#@constraint(m, demand_capacity[i = DEMS], dem[i] <= dem_cap[i])

# variable in objective function
@constraint(m, swf == sum{sup_price[s]*sup[s], s in SUPS} - sum{dem_price[d]*dem[d], d in DEMS})
@constraint(m, invcost == sum{tech_invcost[t]*y[t,n], t in TECHS, n in NODES});
@constraint(m, transcost == sum{prod_transcost[p]*flow[n,nn,p]*distance[n,nn], n in NODES,nn in NODES,p in PRODS});

# Objective Function
#@constraint(m, obj ==  (invcost/20/365) + transcost) 
#@constraint(m, obj ==  cost_total) 		# Obj2
#@constraint(m, obj ==  -struvite_total) 
#@constraint(m, obj == cvar)
#@constraint(m,obj==eta)
#@constraint(m, obj == -sum{gentot[n,"p2",t], t in TECHS, n in NODES})
#@constraint(m, obj == sum{dissatisfaction[n], n in NODES})
@constraint(m, objcons, obj == umanure)		# Obj1
@objective(m, Min, obj)  
    
# Product balances
@constraint(m, flowineq[n in NODES, p in PRODS], flowin[n,p] == sum{flow[nn,n,p], nn in NODES})
@constraint(m, flowouteq[n in NODES, p in PRODS], flowout[n,p] == sum{flow[n,nn,p], nn in NODES});

@constraint(m, supeq[n in NODES, p in PRODS], suptot[n,p] == sum{sup[s], s in SUPS; sup_prod[s] == p && sup_node[s] == n});
@constraint(m, demeq[n in NODES, p in PRODS], demtot[n,p] == sum{dem[d], d in DEMS; dem_prod[d]==p && dem_node[d]==n}) ;

@constraint(m, geneqlb[n in NODES, p in PRODS, t in TECHS], gentot[n,p,t] - (processed[n,tech_refprod[t]] )*transfer[t,p] >= -(1 - y[t,n])*Mflow);
@constraint(m, genequb[n in NODES, p in PRODS, t in TECHS], gentot[n,p,t] - (processed[n,tech_refprod[t]] )*transfer[t,p] <= +(1-y[t,n])*Mflow);

@constraint(m, floweq[n in NODES, p in PRODS],     +  flowin[n,p] 
					           - flowout[n,p]
           	        		           +  sum{gentot[n,p,t], t in TECHS}
			                           +    suptot[n,p]
			                           -    demtot[n,p] == 0) ;

# Adding Split Equations #
@constraint(m, spliteq[n in NODES, p in PRODS], processed[n, p] + unprocessed[n, p] == flowin[n, p] + suptot[n, p]);


													      
# techonology capacities
#techcapeq{n in NODES, t in TECHS}:
@constraint(m, techcapeq[n in NODES, t in TECHS], processed[n,tech_refprod[t]]  <= y[t,n]*tech_cap[t] + (1-y[t,n])*Mflow)

# logic constraint (at most one technology per node)
@constraint(m, onetecheq2[n in NODES], sum{y[t,n], t in TECHS} <= 1)

# logic constraint (if no technology installed)

@constraint(m, techonofflb[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] >= -y[t,n]*tech_cap[t])
@constraint(m, techonoffub[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] <= +y[t,n]*tech_cap[t])

# Assigning 0 value to processed variable (if no technology installed)

@constraint(m, processedlb[n in NODES, p in PRODS], processed[n,p] >= -sum{y[t,n]*tech_cap[t], t in TECHS})
@constraint(m, processedub[n in NODES, p in PRODS], processed[n,p] <= +sum{y[t,n]*tech_cap[t], t in TECHS})


#@addConstraint(m, techonofflb[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] >= -y[t,n]*Mflow)
#@addConstraint(m, techonoffub[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] <= +y[t,n]*Mflow)

# logic constraint (no flow to technology if not installed)
#  noflowineq{n in NODES, t in TECHS, p in PRODS}:  flowin[n,p,t] <= y[t,n]*M;

#Eliminating Self-Flows
@constraint(m, fix[n in NODES,p in PRODS], flow[n,n,p] == 0)

## Fixing Supply Values
@constraint(m, fix_supply_values[i in SUPS_VALUES], sup[i] == sup_value[i])

@constraint(m, co2_total == co2_per_km*sum{flow[n,nn,p]*distance[n,nn], n in NODES,nn in NODES,p in PRODS});
@constraint(m, struvite_total == sum{demtot[n,"p2"], n in NODES})

#@constraint(m, transcost ==  epsilon*cost_max)
#@constraint(m, cost_total <= epsilon*cost_max)
@expression(m, daily_cost,  (invcost/20/365) + transcost)
@constraint(m, daily_cost <= epsilon*budget)		

#@constraint(m, Dissatisfaction[n in NODES], dissatisfaction[n] == - sum{gentot[n, "p2", t], t in TECHS}) 

#@constraint(m, worstcons[n in NODES], dissatisfaction[n]<= eta)

#@constraint(m, Phi_constraint[n in NODES], phi[n] >= dissatisfaction[n] - risk) 
#@constraint(m, cvar == (1/100)*sum{risk + (phi[n]/(1 - alpha)), n in NODES})

#Defining unprocessed manure varaible
@constraint(m, unpro_manure_cons, umanure == sum{p_score[n]*demtot[n, "p1"], n in NODES})

# Adding the constraint for operating cost
@constraint(m, opcost_cons, opcost == sum{tech_opcost[t]*y[t, n], t in TECHS, n in NODES})		# Annual operating cost


@expression(m, dist_flow[p in PRODS], sum{distance[n,nn]*flow[n,nn,p], n in NODES, nn in NODES})
#@expression(m, dist_flow_p2, sum{distance[n,nn]*flow[n,nn,"p2"], n in NODES, nn in NODES}) 
#@expression(m, dist_flow_p3, sum{distance[n,nn]*flow[n,nn,"p3"], n in NODES, nn in NODES})

@expression(m, tot_flow[p in PRODS], sum{flow[n,nn,p], n in NODES, nn in NODES})


