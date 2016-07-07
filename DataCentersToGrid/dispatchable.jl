# Julia script
# - study impact of adding data centers to the system
# Kibaek Kim - ANL/MCS 10012015

using DSPsolver, StochJuMP, MPI;

# Initialize MPI
MPI.Init();

SEASON           = ARGS[1];
penetration      = float(ARGS[2]);
numDispatchables = parse(Int,ARGS[3]);
DLoadShedPenalty = parse(Int,ARGS[4]);
OUT_DIR          = ARGS[5];

penetrationPerCent = round(Int,penetration * 100);
DSPsolver.readSmps("./smps-NoPenalty/dispatchable-$SEASON-$penetrationPerCent-$numDispatchables");
#DSPsolver.readSmps("./smps/dispatchable-$SEASON-$penetrationPerCent-$numDispatchables-$DLoadShedPenalty");
#DSPsolver.readSmps("./smps/test-$SEASON-$penetrationPerCent-$numDispatchables-$DLoadShedPenalty");

DSPsolver.setIntParam("LOG_LEVEL",2);
DSPsolver.setIntParam("BD/NUM_CORES",16);
DSPsolver.setDblParam("WALL_LIM", 3300);
DSPsolver.setDblParam("SCIP/GAP_TOL", 0.0001);

DSPsolver.solveBdMpi(10);

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
	if DSPsolver.getSolutionStatus() == :Optimal || DSPsolver.getSolutionStatus() == :StoppedTime || DSPsolver.getSolutionStatus() == :IterOrTimeLimit || DSPsolver.getSolutionStatus() == :StoppedGap
		# File name
		filename = "$OUT_DIR/$SEASON-$penetration-$numDispatchables-$DLoadShedPenalty";
		
		# Save Objective Value
		fp = open("$filename-Objval.csv", "w");
		print(fp, DSPsolver.getPrimalBound());
		close(fp);
		
		# Save solution
		solution = DSPsolver.getSolution();
		writecsv("$filename-Solution.csv",solution);
	end
end

MPI.Finalize();
