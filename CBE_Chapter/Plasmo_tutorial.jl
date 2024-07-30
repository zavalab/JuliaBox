using Plasmo, Ipopt

# Define OptiGraph
graph = OptiGraph()

# Define OptiNodes
@optinode(graph, n1)
@optinode(graph, n2)
@optinode(graph, n3)

# Add variables, constraint, and objective to OptiNode 1
@variable(n1, x >= 0)
@variable(n1, y >= 2)
@constraint(n1, x + y >= 3)
@objective(n1, Min, 2 * x + y)

# Add variables, constraint, and objective to OptiNode 2
@variable(n2, x >= 0)
@variable(n2, y >= 1)
@constraint(n2, x + 2 * y >= 4)
@objective(n2, Min, 2 * x + y)

# Add variables, constraint, and objective to OptiNode 3
@variable(n3, x >= 1)
@variable(n3, y >= 0)
@constraint(n3, x + 3 * y >= 2)
@objective(n3, Min, 2 * x + y)

# Define linking constraints
@linkconstraint(graph, n1[:x] + n2[:x] == 1)
@linkconstraint(graph, n1[:x] + n2[:x] + n3[:x] == 4)

# Optimize
set_optimizer(graph, Ipopt.Optimizer)
set_to_node_objectives(graph)
optimize!(graph)
