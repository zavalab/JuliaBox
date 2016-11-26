# CHP model for computing Utopian points
# Yankai Cao, Siyu Chen, Luis Fuentes
# UW-Madison, 2016

function get_scenario_model(S,n,NS)

m=Model()

@variable(m, QCHP[1:T,S]>=0,start=99.75)		#heat from CHP
@variable(m, GCHP[1:T,S]>=0,start=85.4)   		#hot water from CHP
@variable(m, 70<=TST[1:T,S]<=105,start=100)  		#temperature of storage tank
@variable(m, GST[i=1:T,j in S]>=0,start=GD[i,j]/2)  	#hot water from storage tank   
@variable(m, GC[i=1:T,j in S]>=0,start=GD[i,j]/2)   	#cold water from utilities
@variable(m, GDU[1:T,S]>=0,start=0)  			#waste water
@variable(m, QL[1:T,S]>=0,start=0)   			#heat lost
@variable(m, NG[1:T,S]>=0,start=210) 			#nature gas /kWh
@variable(m, 0<=WCHP[1:T,S]<=100,start=50)  		#electricity from CHP
@variable(m, WB[i=1:T,j in S]>=0,start=WD[i,j]-50)  	#electricity purchase
@variable(m, 0<=WS[1:T,S]<=100,start=0)  		#electricity sell 
@variable(m, 0<=WC[1:T,S]<=100,start=50) 		#electricity for building
@variable(m, VST[1:T,S]>=0.1,start=12)   		#storage level of tank
@variable(m, VST0>=0.1,start=12)  			#initial storage level     
@variable(m, 0.35<=PL[1:T,S]<=1,start=0.95)     	#partical load
@variable(m, 0<=efPL[1:T,S]<=0.3725,start=0.36) 	#electrical efficiency
@variable(m, Ar[1:T,S]>=0,start=9)       		#heat transfer area
@variable(m, 0<=Fun[1:T,S]<=1,start=0.9) 		#function about PL
@variable(m, EI[1:T,S],start=22.15)      		#costs
@variable(m, GHGE[1:T,S]>=0,start=0)     		#emissions
@variable(m, SW[1:T,S]>=0,start=14.15)   		#water supply
@variable(m, TCost>=0,start=2)            		#total cost
@variable(m, GHGT>=0,start=2)             		#total emissions
@variable(m, SWT>=0,start=2)              		#total water supply
@variable(m, 0<=WMAX<=200, start=110)  			#total capacity
@variable(m, CostCHP>=0,start=19)      			#fixed cost of CHP system
@variable(m, 0<=ST<=50, start=20)      			#storage tank level
@variable(m, CostST>=0,start=1)        			#fixed cost of storage tank
@variable(m, CostFix>=0,start=30)       		#fixed cost of whole sytem
@variable(m, CostOpe>=0,start=500)       		#operation cost of whole system

@NLconstraint(m, constr1[i=1:T,j in S], WCHP[i,j]==WC[i,j]+WS[i,j])     
@NLconstraint(m, constr2[i=1:T,j in S], WD[i,j]==WC[i,j]+WB[i,j])
@NLconstraint(m, constr3[i=1:T,j in S], QCHP[i,j]==efQCHP*(NG[i,j]))
@NLconstraint(m, constr4[i=1:T,j in S], WCHP[i,j]==efPL[i,j]*(NG[i,j]))
@NLconstraint(m, constr5[i=1:T,j in S], PL[i,j]*WMAX==WCHP[i,j])
@NLconstraint(m, constr6[i=1:T,j in S], QCHP[i,j]==GCHP[i,j]*CpW*(TCHP-TAMB[i,j]))
@NLconstraint(m, constr7[j in S], rho*(VST[1,j]-VST0)==GCHP[1,j]-GST[1,j])
@NLconstraint(m, constr8[i=1:(T-1),j in S], rho*(VST[i+1,j]-VST[i,j])==GCHP[i+1,j]-GST[i+1,j])
@NLconstraint(m, constr9[i=1:T,j in S], GD[i,j]+GDU[i,j]==GST[i,j]+GC[i,j])
@NLconstraint(m, constrAr[i=1:T,j in S], ((Ar[i,j]^3)/6)==((VST[i,j])^(2)))
@NLconstraint(m, constr9a[i=1:T,j in S], QL[i,j]==U*Ar[i,j]*(TST[i,j]-TAMB[i,j]))
@NLconstraint(m, constrmix[i=1:T,j in S], ((GD[i,j]+GDU[i,j])*CpW*TS)==(GST[i,j]*CpW*TST[i,j])+(GC[i,j]*CpW*TAMB[i,j]))
@NLconstraint(m, constrtank0[j in S], (rho*CpW*((VST[1,j]*TST[1,j])-(VST0*TST0)))==(GCHP[1,j]*CpW*TCHP)-(GST[1,j]*CpW*TST[1,j])-QL[1,j])
@NLconstraint(m, constrtank[i=1:(T-1),j in S], (rho*CpW*((VST[i+1,j]*TST[i+1,j])-(VST[i,j]*TST[i,j])))==(GCHP[i+1,j]*CpW*TCHP)-(GST[i+1,j]*CpW*TST[i+1,j])-QL[i+1,j])
@NLconstraint(m, constrFun[i=1:T,j in S],(0.3249*(PL[i,j]^2))+(0.445*PL[i,j])+0.2353==Fun[i,j])
@NLconstraint(m, constrfun[i=1:T,j in S],Fun[i,j]*efWCHP==efPL[i,j])
@NLconstraint(m, constr[i=1:T,j in S], GHGE[i,j]==(fNG*NG[i,j]))
@NLconstraint(m, constr[i=1:T,j in S], SW[i,j]==GD[i,j]+GDU[i,j])
@NLconstraint(m, CostCHP==WMAX*CkWe*PT*kF)
@NLconstraint(m, CostST==ST*CST*PT*kF)
@NLconstraint(m, constrSC[i=1:T,j in S], WMAX>=WCHP[i,j])
@NLconstraint(m, constrST[i=1:T,j in S], ST>=VST[i,j])
@NLconstraint(m, constrST[i=1:T,j in S], ST>=VST0)
@NLconstraint(m, constr[i=1:T,j in S], EI[i,j]==(-PSE*(WS[i,j]+WC[i,j]))-(PHW*(QCHP[i,j]-QL[i,j]-(GDU[i,j]*CpW*(TS-TAMB[i,j]))))+(CostNG*NG[i,j])+(CostE[i]*WB[i,j])+(CMCHP*WCHP[i,j]))
@NLconstraint(m, constr[j in S],VST[T,j]==VST0)
@NLconstraint(m, VST0>=0.75*ST)

@NLconstraint(m, CostOpe==sum{EI[i,j],i=1:T,j in S})
SPB=length(S)
@NLconstraint(m, CostFix==(CostST+CostCHP)*SPB)
@NLconstraint(m, TCost==CostOpe+CostFix)
@NLconstraint(m, GHGT==sum{GHGE[i,j],i=1:T,j in S})
@NLconstraint(m, SWT==sum{SW[i,j],i=1:T,j in S})


#=
if n==1
   @NLobjective(m,Min,TCost)
elseif n==2
   @NLobjective(m,Min,GHGT)
elseif n==3
   @NLobjective(m,Min,SWT)
elseif n==4
   @NLobjective(m,Min,(TCost*NS/SPB-CostUP)/(CostNP-CostUP)+(GHGT*NS/SPB-GHGEUP)/(GHGENP-GHGEUP)+(SWT*NS/SPB-SWUP)/(SWNP-SWUP))
end
=#

return m
end