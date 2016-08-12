
push!(LOAD_PATH, pwd())
using Ipopt
using NetJuMP
include("NetParPipsNlp.jl")
using JuMP
import MPI
using Distributions

NS = 10;                   # number of scenarios
S = collect(1:NS)          # scenario set

data = zeros(length(S),1)
srand(123)
μ = 0; σ = 10
d = Uniform(0,10)
d = Normal(μ,σ)
data = rand(d,NS)

# build the model and solve with IPOPT
m = Model(solver=IpoptSolver())
   @variable(m,   x >= 0)
   @variable(m,y[S] >= 0)
 @constraint(m, firststage, x <= 1)
 @constraint(m, secondstage[s in S], x + y[s] >= data[s])
@NLobjective(m, Min, (x-1)^2 + (1/NS)*sum{(y[s]-2)^2,s in S})
solve(m)

# Results
println(getvalue(x))
println("")
println(getvalue(y))
println("")
println("obj ",getObjectiveValue(m))

# build with PLASMO and solve with PIPS-NLP
m = NetModel()
@variable(m, xg>=0)
for i in 1:NS
    bl = Model()
    @variable(bl, x >= 0)         
    @variable(bl, y)    
    @constraint(bl, firststage, x  <= 1)
    @NLobjective(bl, Min, (1/NS)*(x-1)^2 + (1/NS)*(y-2)^2)
    @constraint(bl, secondstage, x+ y >= data[i])
    @addNode(m, bl, "s$i")
    @constraint(m,  xg ==  getvariable(bl, :x))
end
ParPipsNlp_solve(m)

println(getvalue(getvariable(m, :xg)))
println(getvalue(getvariable(getNode(m,"s2"), :y)))
data


