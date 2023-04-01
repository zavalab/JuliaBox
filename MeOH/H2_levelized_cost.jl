using JuMP
using Gurobi
using Distributions
using Plots
using CSV
using DataFrames
using DelimitedFiles
########production scale of 1000ton/day MeOH,########################
######## unconstant production of steam reforming####################
####################################################s#################
##############define electricity market##############################
rtm_matrix = CSV.read("RTM_pan_price.csv", DataFrame)
Interval=35136
T=1:Interval
Elec_price=rtm_matrix[:,5]
#####################################################################
##############define parameters######################################
#utility
Water_price=1.35 #1.35/ton
O2_price=40  #40$/ton
#capital investment:
Af=0.117 #annualized factor
#Physical parameters
#Electrolyzer:
X_e_H2O=9 # making one tone of H2 needs 9 ton of water
X_e_O2=8 # making one tone of H2, generating 8 tons of oxygen

Coef_E_H2=54 # making one ton of H2  requires 54Mwh of electricity
Coef_BOP_cost=258000*54*4 # Balance of Plant (BOP) investment is 258000$/MW, making one ton h2 uses 54mwh of power, one hour has four 15mins intervals
#Coef_Stack_cost=90000*54 # stack investment is 90000/MW,  making one ton h2 uses 54mwh of power, one hour contains 12 5mins
Coef_Stack_cost=81 # stack cost is 81$/t of hydrogen
Coef_Tank_cost_fix=116552 #116552+84*H2 in kg
Coef_Tank_cost_var=84000 #11677+84*H2 in kg

D_H2=2.5 #Demand of H2: 10t/hr
UB_grid=4.62# when maximal power is 1000mw, the max H2 prodcution of electrolzer is 4.62ton/15mins
Min_load = 0.03 #minimum part-load
#####################################################################

m = Model(Gurobi.Optimizer)
#mass balance
@variable(m, fe[T]>=0)  #H2 produced by electrolyzer at time t
@variable(m, fr[T]>=0) #H2 released by tank at time t
@variable(m, fs[T]>=0) #H2 stored at tank at time t
@variable(m, fc[T]>=0) #H2 consumed at tank at time t
@variable(m, s[T]>=0)  #storage level of tank at time t
@variable(m, fw[T]>=0) #Water consumed by electrolyzer at time t
@variable(m, fo[T]>=0) #O2 produced by electrolyzer at time t
@variable(m, fe_max >= 0)  #the full capacity of electrolyzer

@constraint(m, [i in T], fe[i]==fs[i]+fc[i]  ); #PEM balance
@constraint(m, [i in T], fw[i]==fe[i]*X_e_H2O ); #PEM balance
@constraint(m, [i in T], fo[i]==fe[i]*X_e_O2  ); #PEM balance
@constraint(m, [i in T], fc[i]+fr[i] == D_H2 );  # constant demand on Hydrogen
@constraint(m, [i in 1:Interval], fr[i]<=s[i]  );   # gas supplemented by tank is less than the current storage level
@constraint(m, [i in T], fe[i]>=fs[i]  );  # gas stored in each time interval is less than gas produced

@constraint(m, [i in T], fe[i] <= fe_max );  # capacity limit
@constraint(m, [i in T], fe[i] >= Min_load * fe_max  );  # the minimum partial load of the electrolyzer
@constraint(m,  fe_max <= UB_grid  );  # the upper limit of power grid

@constraint(m, s[Interval]==0  );    # all the gas in the tank has to be used up by the end of the year
@constraint(m, s[1]==0  ); # inital storage level is zero
@constraint(m, [i in 1:Interval-1], s[i+1]==s[i]-fr[i]+fs[i] );  #tank balance
#Operation cost
#electrolyzer
@variable(m,Cost_Water>=0) #water cost
@constraint(m, Cost_Water == sum( fe[i]*X_e_H2O*Water_price  for i in T ));

@variable(m,Cost_O2>=0) #oxygen revenue
@constraint(m, Cost_O2 == sum( fe[i]*X_e_O2*O2_price  for i in T ));

@variable(m,Cost_E_op) #electrolyzer power cost
@constraint(m, Cost_E_op == sum( fe[i]*Coef_E_H2*Elec_price[i]  for i in T ));

@variable(m,Cap_Stack>=0) #investment of Stack (Frequent stack replacements are necessary, and there is no need to factor in annualization)
@constraint(m, Cap_Stack== Coef_Stack_cost*sum(fe[i] for i in T) );

@variable(m,OPEX)
@constraint(m, OPEX==Cost_Water+Cost_E_op+Cap_Stack);

@variable(m,Cap_BOP>=0) #investment of BOP
@constraint(m, [i in T],Cap_BOP >= Af*Coef_BOP_cost*fe[i] );

@variable(m,Cap_Tank>=0) #investment of tank
@constraint(m, [i in T],Cap_Tank >= Af*(Coef_Tank_cost_var*s[i]+Coef_Tank_cost_fix) );

@variable(m,CAPEX>=0)
@constraint(m, CAPEX ==Cap_Tank+Cap_BOP+Cap_Stack);

@variable(m,TAC)
@constraint(m, OPEX+CAPEX==TAC);

@objective(m, Min,TAC)
optimize!(m)
#


First_Investment=JuMP.value.(Cap_Tank/Af+Cap_BOP/Af)*1.12*1.45    #installation factor:1.12, contigency, engineering design etc cost factor: 1.45
Annual_OPEX=JuMP.value.(Cost_Water+Cost_E_op+Cap_Stack)+First_Investment*0.05+2000000  #maintenance cost = First_Investment*0.05, 200,0000 salary


H2_produced=sum(JuMP.value.(fe))
D=First_Investment/1.1+sum(Annual_OPEX/((1.1)^t) for t in 1:20) # denominator
N=sum(H2_produced/((1.1)^t) for t in 1:20) # numerator
L=D/N # levelized cost
