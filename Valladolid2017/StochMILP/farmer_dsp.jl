using JuMP, Dsp, MPI;
MPI.Init()

# Model parameters

NS = 5;                  # number of scenarios
S = collect(1:NS);       # scenario set
P = collect(1:3);        # set of crops (1=wheat,2=corn,3=beets)

# Data
prcost = zeros(3)          # production (planting) cost
prcost[1] = 150;
prcost[2] = 230;
prcost[3] = 260;

pcost = zeros(3)           # purchase cost
pcost[1] = 238;
pcost[2] = 210;
pcost[3] = 0;

scost = zeros(3)           # sales cost
scost[1] = 170;
scost[2] = 150;
scost[3] = 36;

demand = zeros(3)          # demand
demand[1] = 200;
demand[2] = 240;
demand[3] = 0;

# assign random data
yield = zeros(length(S),3) # yields
yield[S,1] = 2.5;
yield[S,2] = 3.0;
yield[1,3] = 10;
yield[2,3] = 15;
yield[3,3] = 20;
yield[4,3] = 25;
yield[5,3] = 30;

# Create DSP model
m = Model(NS)

@variable(m, x[P] >= 0)    # acres devoted to crops
@constraint(m, cap, sum(x[j] for j in P) <= 500)

for s = 1:NS
    q = Model(m,s,1/NS)
    @variable(q, y[P] >= 0)    # crops purchase
    @variable(q, w[P] >= 0)    # crops sold
    @variable(q, cost)         # per scenario cost
    @constraint(q, varcost, cost == sum(prcost[j]*x[j] + pcost[j]*y[j] - scost[j]*w[j] for  j in P))
    @constraint(q, bal[j in P], yield[s,j]*x[j]+y[j]-w[j] >= demand[j])
    @constraint(q, sellb, w[3] <= 6000)
    @constraint(q, buyb, y[3] <= 0)
    @objective(q, Min, cost)
end

# Dsp solve types
solve_types = [:Dual, :Benders]

status = solve(m, solve_type = solve_types[2])

# Results
println(getvalue(x))
println("")
# get the children models
ch= m.ext[:DspBlocks].children

for s in 1:length(ch)
    println(getvalue(ch[s][:cost]))
end

println("")
println("obj ", getobjectivevalue(m))

MPI.Finalize();
