# Large farm problem to test scalability of PIPS-NLP
# Yankai Cao and Victor M. Zavala
# University of Wisconsin-Madison, 2016

push!(LOAD_PATH, pwd())

#import MPI


using JuMP
using Distributions
using Ipopt
using Plasmo
 # Initialize MPI
MPI.Init()

srand(123)
NS = 1000;                  # number of scenarios
NP = 100;                   # number of products
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

# create PLASMO model and solve with PIPS-NLP
m = GraphModel()
@variable(m, x[P] >= 0)
@variable(m, s2 >= 0)
@constraint(m, cap, (sum(x[j] for j in P) + s2) == 200)
@objective(m, Min, sum(prcost[j]*x[j] for j in P))
for i in 1:NS
    bl = Model()
    @variable(bl, y[P] >= 0)    # crops purchase
    @variable(bl, 0<=w[j in P] <= sellub[j in P])    # crops sold
    @variable(bl, s[P] >= 0)
    @constraint(m, [j in P], yield[i,j]*x[j]+y[j]-w[j] - s[j] == demand[j])
    @objective(bl, Min, 1.0/NS*sum(pcost[j]*y[j] - scost[j]*w[j] for j in P))
    @addNode(m, bl, "s$i")
end
#ParPipsNlp_solve(m)

#println(getvalue(getvariable(m, :x)))
#println(getvalue(getvariable(getNode(m,"s1"), :w)))

MPI.Finalize()
