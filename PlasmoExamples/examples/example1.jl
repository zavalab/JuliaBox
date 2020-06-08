using Plasmo
using Plots;pyplot()
using GLPK

graph1 = OptiGraph()

@optinode(graph1,n1)
@variable(n1, y >= 2)
@variable(n1,x >= 0)
@constraint(n1,x + y >= 3)
@objective(n1, Min, y)

@optinode(graph1,n2)
@variable(n2, y)
@variable(n2,x >= 0)
@constraint(n2,x + y >= 3)
@objective(n2, Min, y)

@optinode(graph1,n3)
@variable(n3, y )
@variable(n3,x >= 0)
@constraint(n3,x + y >= 3)
@objective(n3, Min, y)

#Create a link constraint linking the 3 models
@linkconstraint(graph1, n1[:x] + n2[:x] + n3[:x] == 5)
optimize!(graph1,GLPK.Optimizer)

#Query Solution
println("n1[:x] = ",value(n1,n1[:x]))
println("n2[:x] = ",value(n2,n2[:x]))
println("n3[:x] = ",value(n3,n3[:x]))
println("objective value = ",objective_value(graph1))


plt_graph1 = Plots.plot(graph1,node_labels = true,markersize = 60,labelsize = 30,
linewidth = 4,layout_options = Dict(:tol => 0.01,:iterations => 2));
plt_matrix1 = Plots.spy(graph1,node_labels = true,markersize = 30);
