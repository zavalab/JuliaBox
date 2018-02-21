# P Recovery Model
# Victor M Zavala, 2016, UW-Madison
# Apoorva M Sampat, 2016, UW-Madison

using JuMP;
using Gurobi;

m = Model(solver=GurobiSolver(Threads = 4, MIPGap = 0.05))
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

const co2_per_km	=	0.1e-3 						# 0.1g of CO2 emitted per km per kg of freight
const kg_to_tons = 1e-3
const tons_to_kg = 1e3
const dollar_to_thousand = 1e-3
const thoudand_to_dollar = 1e3

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

sup_tot		= 	sum(supply_values[:, 2]) #kg/day
# Emissions Metric


# Unit conversions
kW_to_MW = 1e-3
rec_value = 500 # Usd per MWh

# TradeOff Analysis Parameters
epsilon		=	"combined_5000AU"

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
# Units: [tons/day]
@variable(m, flow[NODES,NODES,PRODS] >= 0) 	# product flow (snd,rec,prod)
# Units: [tons/day]
@variable(m, sup[SUPS] >= 0) 			# supply flow
# Units: [tons/day]
@variable(m, dem[DEMS] >=0)
#@defVar(m, y[TECHS,NODES])			# Used for solving the LP

# Binary
@variable(m, y[TECHS,NODES], Bin)
# Units: [tons/day]
@variable(m, flowin[NODES,PRODS] >=0)
# Units: [tons/day]
@variable(m, flowout[NODES,PRODS] >=0)
# Units: [tons/day]
@variable(m, suptot[NODES,PRODS] >=0)
# Units: [tons/day]
@variable(m, demtot[NODES,PRODS] >=0)
# Units: [tons/day]
@variable(m, gentot[NODES,PRODS,TECHS])

# Basic Variables
# Units: [$/day]
@variable(m, swf)
# Units: [Thousand $]
@variable(m, invcost >= 0)
# Units: [Thousand $/day]
@variable(m, transcost >= 0)
# Units: [Thousand $]
@variable(m, obj)
# Units: [Thousand $/day]
@variable(m, opcost >= 0)
# Units: [Thousand $/day]
@variable(m, daily_cost >= 0)
# Units: [kg/day]
@variable(m, co2_total >= 0)

# CVaR Variables
#@variable(m, dissatisfaction[NODES])
#@variable(m, phi[NODES] >= 0 )
#@variable(m, cvar)
#@variable(m, risk)
#@variable(m, eta>=0)
# Units: [tons/day]
@variable(m, processed[NODES, PRODS] >= 0)
# Units: [tons/day]
@variable(m, unprocessed[NODES, PRODS] >= 0)
# Units: [tons/day]
@variable(m, umanure >= 0)
#@variable(m, p_tot >= 0)

# Units: [Thousand $/day]
@variable(m, prod_revenue )

# Units: [Thousand $/day]
@variable(m, demtot_price[NODES, PRODS])

# Units: [Thousand $/day]
@variable(m, rec_revenue)
# Units: [Thousand $/day]
@variable(m, p_credit)
# Units: [Thousand $/day]
@variable(m, rin_revenue)

#@variable(m, roi)

# Capacity Constraints for the variables
#@constraint(m, supply_capacity[i = SUPS], sup[i] <= sup_cap[i])
#@constraint(m, demand_capacity[i = DEMS], dem[i] <= dem_cap[i])

# variable in objective function
#@constraint(m, swf == sum(sup_price[s]*sup[s] for s in SUPS) - sum(dem_price[d]*dem[d] for d in DEMS))

# LHS: [Thousand $]
# RHS: [$]*[Binary]*dollar_to_thousand = [Thousand $]
@constraint(m, invcost == sum(tech_invcost[t]*y[t,n] for t in TECHS for n in NODES)*dollar_to_thousand);

# LHS: [Thousand $/day]*thoudand_to_dollar = [$/day]
# RHS: [$/kg/km]*[ton/day]*[km]*tons_to_kg = [$/day]
@constraint(m, transcost*thoudand_to_dollar == sum(prod_transcost[p]*flow[n,nn,p]*distance[n,nn] for n in NODES for nn in NODES for p in PRODS)*tons_to_kg);

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

# LHS: [Thousand $/day]
# RHS: [Thousand $/day]
@constraint(m, obj_p_cons, 20*365*obj == -20*365*prod_revenue + (20*365*transcost + 20*opcost) + invcost - 20*365*p_credit)
@objective(m, Min, obj)


# Product balances
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, flowineq[n in NODES, p in PRODS], flowin[n,p] == sum(flow[nn,n,p] for nn in NODES))
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, flowouteq[n in NODES, p in PRODS], flowout[n,p] == sum(flow[n,nn,p] for nn in NODES));
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, supeq[n in NODES, p in PRODS], suptot[n,p] == sum(sup[s] for s in SUPS if sup_prod[s] == p && sup_node[s] == n));
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, demeq[n in NODES, p in PRODS], demtot[n,p] == sum(dem[d] for d in DEMS if dem_prod[d]==p && dem_node[d]==n)) ;
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, geneqlb[n in NODES, p in PRODS, t in TECHS], gentot[n,p,t] - (processed[n,tech_refprod[t]] )*transfer[t,p] >= -(1 - y[t,n])*Mflow*kg_to_tons);
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, genequb[n in NODES, p in PRODS, t in TECHS], gentot[n,p,t] - (processed[n,tech_refprod[t]] )*transfer[t,p] <= +(1-y[t,n])*Mflow*kg_to_tons);
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, floweq[n in NODES, p in PRODS],     +  flowin[n,p]
					           - flowout[n,p]
           	        		           +  sum(gentot[n,p,t] for t in TECHS)
			                           +    suptot[n,p]
			                           -    demtot[n,p] == 0) ;

# Adding Split Equations #
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, spliteq[n in NODES, p in PRODS], processed[n, p] + unprocessed[n, p] == flowin[n, p] + suptot[n, p]);

# techonology capacities
#techcapeq{n in NODES, t in TECHS}:
# LHS: [tons/day]
# RHS: [Binary]*[tons/day] = [tons/day]
@constraint(m, techcapeq[n in NODES, t in TECHS], processed[n,tech_refprod[t]]  <= (y[t,n]*tech_cap[t] + (1-y[t,n])*Mflow)*kg_to_tons)

# logic constraint (at most one technology per node)
# LHS: [Binary]
# RHS: [unitless]
@constraint(m, onetecheq2[n in NODES], sum(y[t,n] for t in TECHS) <= 1)

# logic constraint (if no technology installed)
# LHS: [tons/day]
# RHS: [Binary]*[tons/day]
@constraint(m, techonofflb[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] >= -y[t,n]*tech_cap[t]*kg_to_tons)
# LHS: [tons/day]
# RHS: [Binary]*[tons/day]
@constraint(m, techonoffub[n in NODES, t in TECHS, p in PRODS], gentot[n,p,t] <= +y[t,n]*tech_cap[t]*kg_to_tons)

# Assigning 0 value to processed variable (if no technology installed)
# LHS: [tons/day]
# RHS: [Binary]*[tons/day]
@constraint(m, processedlb[n in NODES, p in PRODS], processed[n,p] >= -sum(y[t,n]*tech_cap[t] for t in TECHS)*kg_to_tons)
# LHS: [tons/day]
# RHS: [Binary]*[tons/day]
@constraint(m, processedub[n in NODES, p in PRODS], processed[n,p] <= +sum(y[t,n]*tech_cap[t] for t in TECHS)*kg_to_tons)

#Eliminating Self-Flows
#@constraint(m, fix[n in NODES,p in PRODS], flow[n,n,p] == 0)

for n in NODES
	for p in PRODS
		setlowerbound(flow[n,n,p], 0.0)
		setupperbound(flow[n,n,p], 0.0)
	end
end

## Fixing Supply Values
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, fix_supply_values[i in SUPS_VALUES], sup[i] == sup_value[i]*kg_to_tons)

# LHS: [kg/day]
# RHS: [kg/km/kg]*[tons/day]*[km]*tons_to_kg = [kg/day]
@constraint(m, co2_total == co2_per_km*sum(flow[n,nn,p]*distance[n,nn] for n in NODES for nn in NODES for p in PRODS)*tons_to_kg);
#@constraint(m, struvite_total == sum{demtot[n,"p2"], n in NODES})

#@constraint(m, transcost ==  epsilon*cost_max)
#@constraint(m, cost_total <= epsilon*cost_max)

# LHS: [Thoudand $/day]
# RHS: [Thousand $/day]
@constraint(m, daily_cost_cons, 20*365*daily_cost ==  (invcost) + 20*365*transcost + 20*365*opcost )
#@constraint(m, invcost <= epsilon*budget)

#Defining unprocessed manure varaible
# LHS: [tons/day]
# RHS: [tons/day]
@constraint(m, unpro_manure_cons, umanure == sum(demtot[n, "p1"] for n in NODES))

# Adding the constraint for operating cost
# LHS: [Thousand $/day]
# RHS: [$/day]*[Binary]*dollar_to_thousand = [Thousand $/day]
@constraint(m, opcost_cons, opcost == sum(tech_opcost[t]*y[t, n] for t in TECHS for n in NODES)*dollar_to_thousand)		# Annual operating cost


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

# LHS: [Thoudand $/day]
# RHS: [(tons/day)*($/kg)*(unitless) =  Thousand $/day]
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
p_credit_value = 22.04*1.5
# LHS: [Thoudand $/day]
# RHS: [($/kg)*(tons/day) =  Thousand $/day]
@constraint(m, p_credit_cons, p_credit == p_credit_value*sum(p_recovered[p] for p in [prod_num["Cake1"], prod_num["Cake2"], prod_num["Cake3"], prod_num["Struvite"], prod_num["Struvite+Solids"]]) )
epsilon		=	"p_credit_$(p_credit_value)"

# Electricity to the grid
