# solution of stochMOO problem (formulation D)
# Yankai Cao, Siyu Chen, Luis Fuentes
# UW-Madison, 2016

push!(LOAD_PATH,pwd())
using Ipopt
using Plasmo
using JuMP

MPI.Init()          #initialize MPI
include("data.jl")  #load chp data
NS=188              # number of scenarios
NB =188             # number of block
SPB=NS/NB           # number of scenarios per block
ao=0                # alpha in objective space
as=0                # alpha in scenario space
CostUP = 339.198    # utopia point for cost
GHGEUP = 1.489      # utopia point for emissions
SWUP = 48785.538    # utopia point for water
No = 3              # number of objectives

include("chp_model.jl")  #scenario model building function

# create two-stage graph model
CHP=NetModel()

# define design variables
@variable(CHP, 0<=ST0<=50,start=20)    #storage tank level
@variable(CHP,0<=WMAX0<=200,start=110) #total capacity
@variable(CHP, 0<=z<=10, start=0)      #auxiliary variable

# create array of scenario models
CHPap=Array(JuMP.Model, NB)
for j in 1:NB  
           # get scenario model    
           CHPap[j] = get_scenario_model((j-1)*SPB+1:j*SPB,1,NS)
           # add first-stage variables and constraints    
           @addNode(CHP, CHPap[j], "s$j")
           @constraint(CHP, getvariable(CHPap[j],:ST)==ST0)
           @constraint(CHP, getvariable(CHPap[j],:WMAX)==WMAX0)

           # add second-stage variables and constraints
           @variable(CHPap[j], N>=0, start=2)	   
           @variable(CHPap[j], 0<= t <=10, start=0)
           @variable(CHPap[j], u[1:No] >= 0, start=2)
           @variable(CHPap[j], F, start=2)
	   @constraint(CHP, N >= F - z)
           @constraint(CHPap[j], Fdef, (1-ao)*No*(F-t) == sum{u[i], i in 1:No})
           @constraint(CHPap[j], u[1] >= (getvariable(CHPap[j],:TCost)) - t )
           @constraint(CHPap[j], u[2] >= (getvariable(CHPap[j],:GHGT))- t)
           @constraint(CHPap[j], u[3] >= (getvariable(CHPap[j],:SWT))-t)
end
@variable(CHP, ob, start=2)
@constraint(CHP, NS*(1-as)* (ob - z) == sum{getvariable(CHPap[j], :N),  j in 1:NS})
@objective(CHP, Min, 100*ob)

# solve two-stage program with PIPSNLP
ParPipsNlp_solve(CHP)

# solve with Ipopt
#Ipopt_solve(CHP)

# process and analyze results          
Costmean = 0
CostWorst = 0
GHGEmean = 0
GHGEWorst = 0
SWmean = 0
SWWorst = 0
Cost = Array(Float64, NB)
GHGE = Array(Float64, NB)
SW = Array(Float64, NB)
Worst_per_scenario = []
mean_per_scenario = []

for j in 1:NB
   Cost0 = getvalue(getvariable(CHPap[j], :TCost))
   Costmean += Cost0
   
   GHGE0 = getvalue(getvariable(CHPap[j], :GHGT))
   GHGEmean += GHGE0

   SW0 = getvalue(getvariable(CHPap[j], :SWT))
   SWmean += SW0
   
   Cost[j] = Cost0
   GHGE[j] = GHGE0
   SW[j] = SW0  

   push!(Worst_per_scenario, max(Cost0, GHGE0, SW0))
   push!(mean_per_scenario, (Cost0 + GHGE0 + SW0)/3)

   if CostWorst < Cost0
      CostWorst = Cost0
   end

    if GHGEWorst < GHGE0
      GHGEWorst = GHGE0
   end

    if SWWorst < SW0
      SWWorst = SW0
   end
end
Costmean = Costmean/NB
GHGEmean = GHGEmean/NB
SWmean = SWmean/NB

# report different objective statistics
obj1 = (Costmean + GHGEmean + SWmean)/3
obj2 = (CostWorst + GHGEWorst + SWWorst)/3
obj3 = mean(Worst_per_scenario)
obj4 = maximum(mean_per_scenario)
obj5 = max(Costmean, GHGEmean, SWmean)
obj6 = maximum(Worst_per_scenario)

comm = MPI.COMM_WORLD
if(MPI.Comm_rank(comm) == 0)
    println("z: ", getvalue(z))
    println("Objective Value ", getobjectivevalue(CHP))
    println("tank level", getvalue(ST0))
    println("capacity", getvalue(WMAX0))
    println("Costmean",Costmean)
    println("GHGEmean",GHGEmean)
    println("SWmean",SWmean)
    println("GHGEWorst",GHGEWorst)
    println("SWWorst",SWWorst)
    println("CostWorst",CostWorst)

    println("obj1", obj1)
    println("obj2", obj2)
    println("obj3", obj3)
    println("obj4", obj4)	
    println("obj5", obj5)
    println("obj6", obj6)

    filename = string("Cost_ao", ao, "_as",as, ".txt")
    writedlm(filename, Cost)
    filename = string("GHGE_ao", ao, "_as",as, ".txt")
    writedlm(filename, GHGE)
    filename = string("SW_ao", ao, "_as",as, ".txt")
    writedlm(filename, SW)
end
MPI.Finalize()




