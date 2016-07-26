# Created by Alex Dowling at the University of Wisconsin-Madison in 2016

using MPI, DSPsolver

function solve_with_DSP(solver::ASCIIString, pd::ProblemData)

	# Create structured model
	m = create_struct_model(pd)


	DSPsolver.loadProblem(m);

	if(solver == "DE")
		DSPsolver.solve(DSP_SOLVER_DE);
	end

	if(solver == "BD")
	#	DSPsolver.setInitParam("BD/NUM_CORES",1);
		DSPsolver.solve(DSP_SOLVER_BD);
	end

	if(solver == "DD")
		DSPsolver.solve(DSP_SOLVER_DD);
	end

	println(" ")
	if MPI.Initialized() == false || MPI.Comm_rank(MPI.COMM_WORLD) == 0
	   # print some results
	   println("Solution Status: ", DSPsolver.getSolutionStatus());
	   println("Primal Bound   : ", DSPsolver.getPrimalBound());
	   println("Dual Bound     : ", DSPsolver.getDualBound());
	end
	println(" ")

	# Finalize MPI

	# Report results
	DSPsolver.getSolution(m)
	
	# Warning: Returned solution may not be meaningful when DD is used
	
	return m	

end

##### Script

# Assemble problem data
p = ProblemData(Array(1:10), 2, [0.9 0.1; 0.3 0.7])

# Solve with DSP
m = solve_with_DSP("DD", p)

# Assemble solution object
s = Solution(m, p)

# Plot solution
# plot_solution(s, p)

println("Done.")