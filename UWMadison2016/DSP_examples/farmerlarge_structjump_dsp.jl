push!(LOAD_PATH, pwd())
using MPI, JuMP, StructJuMP, DSPsolver, Distributions
ENV["LD_LIBRARY_PATH"] = "DSP/lib"

NS = 1000;                   # number of scenarios
NP = 10;                    # number of products
S = collect(1:NS)           # scenario set
P = collect(1:NP)           # set of crops (1=wheat,2=corn,3=beets)

prcost = zeros(NP)
d = Uniform(100,200)
prcost = rand(d,NP)

pcost = zeros(NP)
d = Uniform(100,200)
pcost = rand(d,NP)

scost = zeros(NP)
scost = pcost - 50

demand = zeros(NP)
d = Uniform(100,300)
demand = rand(d,NP)/NP

# assign random data
yield = zeros(length(S),NP)
d = Uniform(5,20)
for j in 1:(NP-1)
    yield[S,j] = rand(d,1)[1]
end
d = Uniform(10,30)
yield[S,NP] = rand(d,NS)

sellub = zeros(NP)
d = Uniform(2000,8000)
sellub[P] = rand(d,NP)

#Model Expression
m = StructuredModel(num_scenarios=NS);

# first-stage 
@variable(m, x[i=P] >= 0)
@variable(m, s2 >= 0)
@objective(m, Min, sum{prcost[i] * x[i], i=P})
@constraint(m, capacity, sum{x[i], i=P} + s2 == 200)

# second-stage
for i in 1:NS
    sb = StructuredModel(parent=m, id=i, prob=1/NS);
    @variable(sb, y[P] >= 0)
    @variable(sb, 0<= w[j in P] <= sellub[j in P])
    @variable(sb, s[P] >= 0)
    @objective(sb, Min, sum{pcost[j]*y[j]-scost[j]*w[j], j=P})
    @constraint(sb, bal[j=P], yield[i,j]*x[j] + y[j] - w[j] - s[j] == demand[j])
end

#print(m)
#print(getchildren(m)[1])
#print(getchildren(m)[2])

option = "BD"

# load problem to model object
DSPsolver.loadProblem(m);

# set parameters
DSPsolver.setIntParam("LOG_LEVEL",1);
DSPsolver.setIntParam("ITER_LIM",100);
DSPsolver.setDblParam("SCIP/GAP_TOL",0.0);

# solve problem
if option == "DE"
   DSPsolver.solve(DSP_SOLVER_DE);
elseif option == "BD"
       DSPsolver.setIntParam("BD/NUM_CORES",2);
       DSPsolver.solve(DSP_SOLVER_BD);
elseif option == "DD"
       DSPsolver.setIntParam("DD/FEAS_CUTS",1);
       DSPsolver.setIntParam("DD/OPT_CUTS",1);
       DSPsolver.setIntParam("DD/EVAL_UB",1);
       DSPsolver.setIntParam("DD/MASTER_ALGO",IPM_FEAS);
       DSPsolver.solve(DSP_SOLVER_DD);
end

if MPI.Initialized() == false || MPI.Comm_rank(MPI.COMM_WORLD) == 0
   # print some results
   println("Solution Status: ", DSPsolver.getSolutionStatus());
   println("Primal Bound   : ", DSPsolver.getPrimalBound());
   println("Dual Bound     : ", DSPsolver.getDualBound());
end

if option == "DD"
   MPI.Finalize();
end


DSPsolver.getSolution(m);
#println(m.colNames)	
#println(m.colVal)
#println(getvalue(getvariable(m,:x)))
#println("w scenario 1 ")
#println(getvalue(getvariable(getchildren(m)[1], :w)))
#println("w scenario 2 ")
#println(getvalue(getvariable(getchildren(m)[2], :w)))