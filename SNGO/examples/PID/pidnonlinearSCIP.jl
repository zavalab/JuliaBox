# Yankai Cao, Victor Zavala
# UW-Madison, 2018

push!(LOAD_PATH, ENV["SINGODIR"])
using JuMP, SCIP, Ipopt, PlasmoOld 
using Distributions

srand(0)
N=20
Tf=15
NS=10
h=Tf/N
T=collect(1:N)
Tm=collect(1:(N-1))
mT=collect(2:(N))
S=collect(1:NS)
time = zeros(N)
for t = T
    time[t]=h*(t-1)
end


x0=0
K = zeros(NS)
Kd = zeros(NS)
tau = zeros(NS)
xsp = zeros(NS)
d = Array(Float64, NS , N)


K[1]=5.0
K[2]=1.0
K[3]=2.0
K[4]=1.0
K[5]=1.0
K[6]=2.0
K[7]=3.0
K[8]=4.0
K[9]=5.0
K[10]=1.5

Kd[1]=0.5
Kd[2]=0.4
Kd[3]=0.3
Kd[4]=0.5
Kd[5]=0.8
Kd[6]=0.5
Kd[7]=0.4
Kd[8]=0.3
Kd[9]=0.7
Kd[10]=0.25

tau[1]=0.5
tau[2]=0.5
tau[3]=0.3
tau[4]=0.4
tau[5]=0.7
tau[6]=0.5
tau[7]=0.5
tau[8]=0.3
tau[9]=0.4
tau[10]=0.6

xsp[1]=-1.0
xsp[2]= 1.0
xsp[3]= -2.0
xsp[4]= 1.5
xsp[5]= 2.0
xsp[6]=-1.5
xsp[7]= 1.5
xsp[8]= 0
xsp[9]= 1.5
xsp[10]= 2.0

d[1,:]= -1
d[2,:]= +1
d[3,:]=  2
d[4,:]= -2
d[5,:]= -1
d[6,:]= -0.5
d[7,:]= +0.5
d[8,:]=  2.5
d[9,:]= -1.5
d[10,:]=  -1

for j in 11:NS
    dis = Uniform(1.0,5.0)
    K[j] = rand(dis,1)[1]
    dis = Uniform(0.2,0.8)
    Kd[j] = rand(dis,1)[1]
    dis = Uniform(0.2,0.8)
    tau[j] = rand(dis,1)[1]
    dis = Uniform(-2.3,2.3)
    xsp[j] = rand(dis,1)[1]
    dis = Uniform(-2.5,2.5)
    d[j,:] = rand(dis,1)[1]
end

function scenario(S,NS)
    m=Model()
    @variable(m, -10<= Kc <=10, start=1)
    @variable(m, -100<= tauI <=100, start=1)
    @variable(m, -100<= tauD <=1000, start=1)

    @variable(m, -2.5<=x[j = S, t=T]<=2.5)
    @variable(m, -5.0<=u[j = S, t=mT]<=5.0)
    @variable(m, int[j = S, t=T])

    @constraint(m, [s in S,t in Tm], (1/tau[s])*(x[s,t+1]-x[s,t])/h + x[s,t+1]*x[s,t+1] == K[s]*u[s,t+1]+Kd[s]*d[s,t])
    @constraint(m, [s in S,t in Tm], u[s,t+1] == Kc*(xsp[s]-x[s,t])+ tauI*int[s,t+1] + tauD*(x[s,t+1]-x[s,t])/h)
    @constraint(m, [s in S,t in Tm], (int[s,t+1]-int[s,t])/h == (xsp[s]-x[s,t+1]))

    @constraint(m, [s in S], x[s,1]==x0)
    @constraint(m, [s in S], int[s,1]==0)

    @variable(m, cost[s = S, t=T])
    @constraint(m, [s in S], cost[s,1] == (10*(x[s,1]-xsp[s])^2))
    @constraint(m, [s in S,t in mT], cost[s,t] == (10*(x[s,t]-xsp[s])^2 + 0.01*u[s,t]^2))

    @variable(m, costS[s = S])
    @constraint(m, [s in S], costS[s] == 100/length(T)/NS * sum{cost[s,t], t=T})
    @objective(m, Min, sum{costS[s], s = S})
    return m
end


master=NetModel()
@variable(master, -10<= Kc <=10, start=1)
@variable(master, -100<= tauI <=100, start=1)
@variable(master, -100<= tauD <=1000, start=1)
nodes = Array(JuMP.Model, NS)

for j in 1:NS
    nodes[j] = scenario(j:j, NS)
    @addNode(master, nodes[j], "s$j")
    @constraint(master, getvariable(nodes[j],:Kc)==Kc)
    @constraint(master, getvariable(nodes[j],:tauI)==tauI)
    @constraint(master, getvariable(nodes[j],:tauD)==tauD)
end


mingap = 0.01 
m= extensiveSimplifiedModel(master)
m.solver = SCIPSolver("limits/gap", mingap, "limits/absgap", mingap, "limits/time", 43200.0)
solve(m)
