# Farmer example from Birge and Louveaux book.

# Model Data
# STOCHASTIC MODELING FRAMEWORK

NS = 3;                        # number of scenarios
probability = [1/3, 1/3, 1/3]; # probability


# FIRST-STAGE MODEL
CROPS = 1:3; # set of crops (wheat, corn and sugar beets, resp.)
Cost = [150 230 260]; # cost of planting crops
Budget = 500; # budget capacity


# SECOND-STAGE MODELS
PURCH = 1:2; # set of crops to purchase (wheat and corn, resp.)
SELL  = 1:4; # set of crops to sell (wheat, corn, sugar beets under 6K and those over 6K)
Purchase = [238 210;
            238 210;
            238 210];   # purchase price
Sell = [170 150 36 10;
        170 150 36 10;
        170 150 36 10]; # selling price
Yield = [3.0 3.6 24.0;
         2.5 3.0 20.0;
         2.0 2.4 16.0];
Minreq = [200 240 0;
          200 240 0;
          200 240 0]; # minimum crop requirement
println("")


#Model Expression
push!(LOAD_PATH, pwd())
using MPI, JuMP, StructJuMP, DSPsolver
# CREATE STOCHASTIC MODEL
ENV["LD_LIBRARY_PATH"] = "DSP/lib"
m = StructuredModel(num_scenarios=NS);


# FIRST-STAGE MODEL

# first-stage variables
@variable(m, x[i=CROPS] >= 0, Int)

# first-stage objective
@objective(m, Min, sum{Cost[i] * x[i], i=CROPS})

# first-stage constraint
@constraint(m, const_budget,
               sum{x[i], i=CROPS} <= Budget)


# SECOND-STAGE MODELS

for s in 1:NS
    # stochastic block
    sb = StructuredModel(parent=m, id = s, prob = probability[s]);

    # second-stage variables
    @variable(sb, y[j=PURCH] >= 0)
    @variable(sb, w[k=SELL] >= 0)

    # objective
    @objective(sb, Min,
                  sum{Purchase[s,j] * y[j], j=PURCH}
                  - sum{Sell[s,k] * w[k], k=SELL})
    # constraints
    @constraint(sb, const_minreq[j=PURCH],
                   Yield[s,j] * x[j] + y[j] - w[j] >= Minreq[s,j])
    @constraint(sb, const_minreq_beets,
                   Yield[s,3] * x[3] - w[3] - w[4] >= Minreq[s,3])
    @constraint(sb, const_aux, w[3] <= 6000)
end

print(m)

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
println(m.colNames)	
println(m.colVal)
println(getchildren(m)[1].colNames)
println(getchildren(m)[1].colVal)
println("y scenario 2 ")
println(getvalue(getvariable(getchildren(m)[2], :y)))