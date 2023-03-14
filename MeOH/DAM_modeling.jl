using JuMP
using Gurobi
using Distributions
using Plots
using CSV
using DataFrames
using DelimitedFiles
########production scale of 1000ton/day MeOH,########################
######## unconstant production of steam reforming####################
#####################################################################
##############define electricity market##############################
DAM_matrix = CSV.read("DAM_pan_price.csv", DataFrame)
Interval=8784;
T=1:Interval;
Elec_price=DAM_matrix[:,3];
#####################################################################
##############define parameters######################################
#utility
Water_price=1.35 #1.35/ton
O2_price=40  #40$/ton
#capital investment:
Af=0.117 #annualized factor
#Physical parameters
#Electrolyzer:
X_e_H2O=9 # making one ton of H2 consumes 9 ton of water
X_e_O2=8 # making one ton of H2 generates 8 tons of oxygen

Coef_E_H2=54 # making one ton of H2 requires 54Mwh of electricity
Coef_BOP_cost=258000*54 # Balance of Plant (BOP) investment is 258000$/MW, making one ton h2 uses 54mwh of power
Coef_Stack_cost=81 # stack cost is 81$/t of hydrogen
Coef_Tank_cost_fix=116552 #116552+84*H2 in kg
Coef_Tank_cost_var=84000 #11677+84*H2 in kg

D_H2=11 #Demand of H2: 11t/hr
UB_grid=18.48 # when maximal power from gird is 1000mw, the max H2 prodcution of electrolzer is 18.48t/hr

#####################################################################
##############solve the probelms#####################################
m = Model(Gurobi.Optimizer)
#mass balance
@variable(m, fe[T] >= 0)  #H2 produced by electrolyzer at time t
@variable(m, fr[T] >= 0) #H2 released by tank at time t
@variable(m, fs[T] >= 0) #H2 stored at tank at time t
@variable(m, fc[T] >= 0) #H2 consumed at tank at time t
@variable(m, s[T] >= 0)  #storage level of tank at time t
@variable(m, fw[T] >= 0)  #Water consumed by electrolyzer at time t
@variable(m, fo[T] >= 0)  #O2 produced by electrolyzer at time t
@constraint(m, [i in T], fe[i] == fs[i]+fc[i]  ); #PEM balance
@constraint(m, [i in T], fw[i] == fe[i]*X_e_H2O ); #PEM balance
@constraint(m, [i in T], fo[i] == fe[i]*X_e_O2  ); #PEM balance
@constraint(m, [i in T], fc[i]+fr[i] == D_H2 );  # constant demand on Hydrogen
@constraint(m, [i in 1:Interval], fr[i] <= s[i]  );   # gas supplemented by tank is less than the current storage level
@constraint(m, [i in T], fe[i] >= fs[i]  );  # gas stored in each time interval is less than gas produced

@constraint(m, [i in T], fe[i] <= UB_grid  );  # the upper limit of power grid

@constraint(m, s[Interval] == 0  );   # all the gas in the tank has to be used up by the end of the year
@constraint(m, s[1] == 0  ); # inital storage level is zero
@constraint(m, [i in 1:Interval-1], s[i+1] == s[i]-fr[i]+fs[i] );  #tank balance

#Operation cost
#electrolyzer
@variable(m,Cost_Water>=0) #water cost
@constraint(m, Cost_Water == sum( fe[i]*X_e_H2O*Water_price  for i in T ));
@variable(m,Cost_O2>=0) #oxygen revenue
@constraint(m, Cost_O2 == sum( fe[i]*X_e_O2*O2_price  for i in T ));
@variable(m,Cost_E_op) #electrolyzer power cost
@constraint(m, Cost_E_op == sum( fe[i]*Coef_E_H2*Elec_price[i]  for i in T ));

@variable(m,OPEX)
@constraint(m, OPEX == Cost_Water+Cost_E_op-Cost_O2);

@variable(m,Cap_BOP >= 0) #investment of BOP
@constraint(m, [i in T],Cap_BOP >= Af*Coef_BOP_cost*fe[i] );

@variable(m,Cap_Stack >= 0) #investment of Stack ((Frequent stack replacements are necessary, and there is no need to factor in annualization))
@constraint(m, Cap_Stack == Coef_Stack_cost*sum(fe[i] for i in T) );

@variable(m,Cap_Tank >= 0) #investment of tank
@constraint(m, [i in T],Cap_Tank >= Af*(Coef_Tank_cost_var*s[i]+Coef_Tank_cost_fix) );

@variable(m,CAPEX >= 0)
@constraint(m, CAPEX == Cap_Tank+Cap_BOP+Cap_Stack);

@variable(m,TAC)
@constraint(m, OPEX+CAPEX == TAC);
@objective(m, Min, TAC)

optimize!(m)
#
Record_OPEX=JuMP.value.(OPEX)
Record_CAPEX=JuMP.value.(CAPEX)
Record_C_Tank=JuMP.value.(Cap_Tank)
Record_C_BOP=JuMP.value.(Cap_BOP)
Record_C_Stack=JuMP.value.(Cap_Stack)
Record_E_op=JuMP.value.(Cost_E_op)
Record_O2=JuMP.value.(Cost_O2)
Record_Water=JuMP.value.(Cost_Water)

Record_s=zeros(Interval)
Record_fe=zeros(Interval)
Record_fr=zeros(Interval)
Record_fs=zeros(Interval)
Record_fc=zeros(Interval)
Record_fw=zeros(Interval)
Record_fo=zeros(Interval)

Record_s[:]=JuMP.value.(s)
Record_fe[:]=JuMP.value.(fe)
Record_fr[:]=JuMP.value.(fr)
Record_fs[:]=JuMP.value.(fs)
Record_fc[:]=JuMP.value.(fc)
Record_fw[:]=JuMP.value.(fw)
Record_fo[:]=JuMP.value.(fo)
