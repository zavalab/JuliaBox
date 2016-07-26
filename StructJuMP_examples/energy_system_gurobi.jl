# Created by Alex Dowling at the University of Wisconsin-Madison in 2016

using Gurobi

function solve_with_Gurobi(pd::ProblemData)

	pd.flat = true

	# Create flattened model
	m = create_struct_model(pd)

	setsolver(m, GurobiSolver())

	#JuMP.solve(m, ignore_solve_hook=true)
	solve(m)
	
	return m	

end

##### Script

# Assemble problem data
# Note: To use flattened model, number of blocks must equal 1
p = ProblemData(Array(1:10), 1, [0.9 0.1; 0.3 0.7])

# Solve with DSP
m = solve_with_Gurobi(p)

# Assemble solution object
s = Solution(m, p)

# Plot solution
# plot_solution(s, p)

println("Done.")

