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
#for laptop format
rtm_matrix = CSV.read("/MeOH//RTM_pan_price.csv")
Interval=35136
T=1:Interval
Elec_price=rtm_matrix[:,5]
#####################################################################
##############define parameters######################################
#utility
Water_price=1.35 #1.35/ton
O2_price=40  #40$/ton
#capital investment:
Af=0.117
#Physical parameters
#Electrolyzer:
X_e_H2O=9 # making one tone of H2 needs 9 ton of water
X_e_O2=8 # making one tone of H2, generating 8 tons of oxygen

Coef_E_H2=54 # making one ton of H2  requires 54Mwh of electricity
Coef_BOP_cost=258000*54*4 # BOP investment is 258000/MW, making one ton h2 uses 54mwh of power, one hour contains 12 5mins
#Coef_Stack_cost=90000*54 # stack investment is 90000/MW,  making one ton h2 uses 54mwh of power, one hour contains 12 5mins
Coef_Stack_cost=81 # stack cost is 81$/t of hydrogen
Coef_Tank_cost=84000 #11677+84*H2 in kg

life_e=600000  # stack lifetime: 600000 intervel, each intervel is 5mins
D_H2=11/4 #Demand of H2: 11/4t / 15mins
UB_grid=4.62# when maximal power load is 1000mw, the max H2 prodcution of electrolzer is 4.62ton/15mins

#####################################################################
##############solve the probelms#####################################
m = Model(with_optimizer(Gurobi.Optimizer))
#mass balance
@variable(m, fe[T]>=0)  #H2 produced by electrolyzer at time t
@variable(m, fr[T]>=0) #H2 released by tank at time t
@variable(m, fs[T]>=0) #H2 stored at tank at time t
@variable(m, fc[T]>=0) #H2 consumed at tank at time t
@variable(m, s[T]>=0)  #storage level of tank at time t
@variable(m, fw[T]>=0)  #Water consumed by electrolyzer at time t
@variable(m, fo[T]>=0)  #O2 produced by electrolyzer at time t
@constraint(m, [i in T], fe[i]==fs[i]+fc[i]  ); #PEM balance
@constraint(m, [i in T], fw[i]==fe[i]*X_e_H2O ); #PEM balance
@constraint(m, [i in T], fo[i]==fe[i]*X_e_O2  ); #PEM balance
@constraint(m, [i in T], fc[i]+fr[i] == D_H2 );  # satisfies the demand on Hydrogen
@constraint(m, [i in 1:Interval], fr[i]<=s[i]  );   # gas tank release is less than storage
@constraint(m, [i in T], fe[i]>=fs[i]  );  # gas stored is less than gas produced

@constraint(m, [i in T], fe[i]<=UB_grid  );  # the upper limit of power grid

@constraint(m, s[Interval]==0  );   # al the gas will be released at the end of year
@constraint(m, s[1]==0  );
@constraint(m, [i in 1:Interval-1], s[i+1]==s[i]-fr[i]+fs[i] );  #gas stored in the tank

#Operation cost
#electrolyzer
@variable(m,Cost_Water>=0) #water cost
@constraint(m, Cost_Water == sum( fe[i]*X_e_H2O*Water_price  for i in T ));
@variable(m,Cost_O2>=0) #oxygen revenue
@constraint(m, Cost_O2 == sum( fe[i]*X_e_O2*O2_price  for i in T ));
@variable(m,Cost_E_op) #electrolyzer power cost
@constraint(m, Cost_E_op == sum( fe[i]*Coef_E_H2*Elec_price[i]  for i in T ));

@variable(m,OPEX)
@constraint(m, OPEX==Cost_Water+Cost_E_op-Cost_O2);

@variable(m,Cap_BOP>=0) #investment of BOP
@constraint(m, [i in T],Cap_BOP >= Af*Coef_BOP_cost*fe[i] );

@variable(m,Cap_Stack>=0) #investment of Stack
@constraint(m, Cap_Stack== Coef_Stack_cost*sum(fe[i] for i in T) );

@variable(m,Cap_Tank>=0) #investment of tank
@constraint(m, [i in T],Cap_Tank >= Af*Coef_Tank_cost*s[i]+116552*Af );

@variable(m,CAPEX>=0)
@constraint(m, CAPEX ==Cap_Tank+Cap_BOP+Cap_Stack);

@variable(m,TAC)
@constraint(m, OPEX+CAPEX==TAC);
@objective(m, Min, TAC+0.0001*sum(fr[:]))

optimize!(m)



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
