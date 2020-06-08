using Plasmo

T = 100          #number of time points
d = sin.(1:T)    #disturbance vector

graph = OptiGraph()
@optinode(graph,state[1:T])
@optinode(graph,control[1:T-1])

for node in state
    @variable(node,x)
    @constraint(node, x >= 0)
    @objective(node,Min,x^2)
end
for node in control
    @variable(node,u)
    @constraint(node, u >= -1000)
    @objective(node,Min,u^2)
end

@linkconstraint(graph,[i = 1:T-1],state[i+1][:x] == state[i][:x] + control[i][:u] + d[i])
n1 = state[1]
@constraint(n1,n1[:x] == 0)
