# Example on using linking constraints in PLASMO
# Yankai Cao and Victor M. Zavala
# University of Wisconsin-Madison, 2016

push!(LOAD_PATH, pwd())
using JuMP
using Distributions
using Ipopt
using Plasmo
MPI.Init()  # Initialize MPI

srand(123)
NS = 100;                   # number of scenarios
NP = 1000;
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

# create plasmo model
m = GraphModel()
@variable(m, x[P] >= 0)    # acres devoted to crops
@variable(m, s2 >= 0)
@constraint(m, cap, (sum(x[j] for j in P) + s2) == 200)
@objective(m, Min, sum(prcost[j]*x[j] for j in P))
for i in 1:NS
    bl = Model()
    @variable(bl, y[P] >= 0)    # crops purchase
    @variable(bl, 0<=w[j in P] <= sellub[j in P])    # crops sold
    @variable(bl, s[P] >= 0)
    @constraint(m, bal[j in P], yield[i,j]*x[j]+y[j]-w[j] - s[j] == demand[j])
    @variable(bl, cost)
    @constraint(bl, cost ==sum(pcost[j]*y[j] - scost[j]*w[j] for j in P))
    @objective(bl, Min, 1.0/NS*cost)
    @addNode(m, bl, "s$i")
end

# impose expected value constraint on cost
@constraint(m, sum(getNode(m,"s$i")[:cost] for i in 1:NS) >= 50000)

# solve with Ipopt
Ipopt_solve(m)

# solve with PIPS-NLP
#ParPipsNlp_solve(m)

# access solution
#println(getvalue(getvariable(m, :x)))
#println(getvalue(getvariable(getNode(m,"s1"), :w)))
#println(getvalue(getvariable(getNode(m,"s1"), :cost)))
#println(getvalue(getvariable(getNode(m,"s2"), :cost)))

MPI.Finalize()
