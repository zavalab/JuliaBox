# Julia script (Deterministic equivalent version)
# - A economic dispatch model with allocation of dispatchable loads
# - formulated as a two-stage stochastic programming problem with integer variables in the first stage.
# - The model allows curtailments from imports, renewables and wind productions.
# - The problem should be alwasy feasible.
# Kibaek Kim - ANL/MCS 12/11/2015

using JuMP;


SEASON           = ARGS[1];
penetration      = float(ARGS[2]);
numDispatchables = parse(Int,ARGS[3]);

nScenarios       = 1000;
maxDispatchable  = 200; # 200 MW for each dispatchable load
DLoadShedPenalty = 0;

# scale penetration level to account for the additional data center load
original_penetration = penetration;
penetration = penetration * (655+96) / 655;

println("reading data ...");
include("ZC_data.jl");

penetration = original_penetration;

tic();
println("creating model ...");
m = Model();

@defVar(m, InstallDispatchable[n=BUSES] >= 0, Int);    # number of dispatchable loads at bus n
@defVar(m, DLoad[n=BUSES,t=PERIODS] >= 0); #
@defVar(m, 0 <= Gen0[i=GENERATORS] <= max_gen[i]);                # Initial power generation
@defVar(m, 0 <= Gen[i=GENERATORS,t=PERIODS] <= max_gen[i]);       # Power generation
@defVar(m, -flowmax[l] <= Flow[l=LINES,t=PERIODS] <= flowmax[l]); # transmission line capacity
@defVar(m, -360 <= Angle[n=BUSES,t=PERIODS] <= 360);              # phase angle
@defVar(m, 0 <= LoadShed[j=LOADS,t=PERIODS] <= dictDemand[j][t]); # load shedding
@defVar(m, 0 <= ImportCurtailment[i=IMPORTS,t=PERIODS] <= dictImportProduction[i][t]);     # import curtailment
@defVar(m, 0 <= RenewableCurtailment[r=REGENERATORS,t=PERIODS] <= dictREProduction[r][t]); # renewable curtailment
@defVar(m, WindCurtailment[w=WINDS,t=PERIODS] >= 0);        # wind curtailment

@setObjective(m, Min,
        DLoadShedPenalty * sum{maxDispatchable * InstallDispatchable[n], n=BUSES, t=PERIODS}
        - DLoadShedPenalty * sum{DLoad[n,t], n=BUSES, t=PERIODS}
	+ sum{gen_cost[i] * Gen[i,t], i=GENERATORS, t=PERIODS}
	+ 100 * sum{WindCurtailment[w,t], w=WINDS, t=PERIODS}
	+ 1000 * sum{LoadShed[j,t], j=LOADS, t=PERIODS}
	+ 1000 * sum{ImportCurtailment[i,t], i=IMPORTS, t=PERIODS}
	+ 2000 * sum{RenewableCurtailment[r,t], r=REGENERATORS, t=PERIODS});

@addConstraint(m, NUMDISPLOAD, 
	sum{InstallDispatchable[n], n=BUSES} <= numDispatchables);

# Dispatchable load
@addConstraint(m, DISPLOAD[n=BUSES,t=PERIODS], 
	 maxDispatchable * InstallDispatchable[n] - DLoad[n,t] >= 0);

# Ramping rate in normal operating status
@addConstraint(m, RAMP_DOWN[i=GENERATORS,t=1], Gen0[i] - Gen[i,1] <= ramp_down[i]);
@addConstraint(m, RAMP_DOWN[i=GENERATORS,t=2:length(PERIODS)], Gen[i,t-1] - Gen[i,t] <= ramp_down[i]);
@addConstraint(m, RAMP_UP[i=GENERATORS,t=1], Gen[i,1] - Gen0[i] <= ramp_up[i]);
@addConstraint(m, RAMP_UP[i=GENERATORS,t=2:length(PERIODS)], Gen[i,t] - Gen[i,t-1] <= ramp_up[i]);

# Power flow equation
@addConstraint(m, POWER_FLOW[l=LINES,t=PERIODS],
	Flow[l,t] - susceptance[l] * (Angle[tobus[l],t] - Angle[frombus[l],t]) == 0);

# Power balance constraints for system
@addConstraint(m, POWER_BALANCE[n=BUSES,t=PERIODS],
	sum{Flow[l,t], l=LINES; tobus[l] == n} 
	- sum{Flow[l,t], l=LINES; frombus[l] == n}
	+ sum{Gen[i,t], i=GENERATORS; gen2bus[i] == n} 
	- DLoad[n,t]
	+ sum{LoadShed[j,t], j=LOADS; dBusLoads[j] == n}
	- sum{ImportCurtailment[i,t], i=IMPORTS; dBusImportPoints[i] == n}
	- sum{RenewableCurtailment[r,t], r=REGENERATORS; dBusREGenerators[r] == n}
	- sum{WindCurtailment[w,t], w=WINDS; dBusWindGenerators[w] == n}
	== sum{dictDemand[j][t], j=LOADS; dBusLoads[j] == n}
	- sum{dictImportProduction[i][t], i=IMPORTS; dBusImportPoints[i] == n}
	- sum{dictREProduction[r][t], r=REGENERATORS; dBusREGenerators[r] == n}
	- sum{dictWindProduction[w][t,1], w=WINDS; dBusWindGenerators[w] == n});

@addConstraint(m, WINDCURT_BOUND[w=WINDS,t=PERIODS],
	WindCurtailment[w,t] <= dictWindProduction[w][t,1]);

toc();

# file name
penetrationPerCent = round(Int,penetration * 100);
fname = "smps-NoPenalty/dispatchable-$SEASON-$penetrationPerCent-$numDispatchables";
#fname = "smps/dispatchable-$SEASON-$penetrationPerCent-$numDispatchables-$DLoadShedPenalty";
#fname = "smps/test-$SEASON-$penetrationPerCent-$numDispatchables-$DLoadShedPenalty";

# write mps file
println("writing $fname.cor ...");
writeMPS(m,"$fname.cor");

# write other files

# write tim file
startVar2 = length(BUSES) + 1; # start of second-stage variables
println("writing $fname.tim ...");
fp = open("$fname.tim", "w");
print(fp,"*23456789 123456789 123456789 123456789 123456789 123456789\n");
print(fp,"TIME          JuMPModel\n");
print(fp,"PERIODS       IMPLICIT\n");
print(fp,"    VAR1      CON1                     PER1\n");
print(fp,"    VAR$startVar2      CON2                     PER2\n");
print(fp,"ENDATA");
close(fp);

tic();
# write sto file
ncons = MathProgBase.numlinconstr(m);
nbalance = length(BUSES) * length(PERIODS);
nwindcurt = length(WINDS) * length(PERIODS);
startofBalanceConstraints = ncons - nwindcurt - nbalance;
startOfWindCurtailment = ncons - nwindcurt;

println("writing $fname.sto ...");
fp = open("$fname.sto", "w");
print(fp,"*23456789 123456789 123456789 123456789 123456789 123456789\n");
print(fp,"STOCH         JuMPModel\n");
print(fp,"SCENARIOS\n");
nscenarios = 1000;
prob = 1/nscenarios;
for s = 1:nscenarios
	print(fp," SC SCEN$s     ROOT      $prob       PER2\n");
	for n = 1:length(BUSES)
		if in(BUSES[n], values(dBusWindGenerators)) == false
			continue;
		end
		for t = 1:length(PERIODS)
			i = startofBalanceConstraints + (n-1) * length(PERIODS) + t;
			rhs = 0.0;
			for j = LOADS
				if dBusLoads[j] == BUSES[n]
					rhs += dictDemand[j][PERIODS[t]];
				end
			end
			for j = IMPORTS
				if dBusImportPoints[j] == BUSES[n]
					rhs -= dictImportProduction[j][PERIODS[t]];
				end
			end
			for j = REGENERATORS
				if dBusREGenerators[j] == BUSES[n]
					rhs -= dictREProduction[j][PERIODS[t]];
				end
			end
			for j = WINDS
				if dBusWindGenerators[j] == BUSES[n]
					rhs -= dictWindProduction[j][PERIODS[t],s];
				end
			end
			print(fp,"    rhs    CON$i    $rhs\n");
		end
	end
	for w = 1:length(WINDS)
		for t = 1:length(PERIODS)
			i = startOfWindCurtailment + (w-1) * length(PERIODS) + t;
			rhs = dictWindProduction[WINDS[w]][PERIODS[t],s];
			print(fp,"    rhs    CON$i   $rhs\n");
		end
	end
end
print(fp,"ENDATA");
close(fp);
toc();

println("done.");
