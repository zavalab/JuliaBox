# Created by Alex Dowling at the University of Wisconsin-Madison in 2016

function solve_with_PIPS(pd::ProblemData,my_solver::ASCIIString="PipsNlp")

	m = create_struct_model(pd, true)
	
	solve(m, solver=my_solver)
	
	return m	

end

##### Script

# Assemble problem data
p = ProblemData(Array(1:10), 2, [0.9 0.1; 0.3 0.7])

# Solve with DSP
m = solve_with_PIPS(p,"PipsNlp")

# Assemble solution object
s = Solution(m, p)

# Plot solution
# plot_solution(s, p)

println("Done.")

