# Julia script 
# - A economic dispatch model with data center allocations
# - The model allows curtailments from imports, renewable and wind productions.
# - The problem is always feasible.
# Kibaek Kim - ANL/MCS 12/04/2015

m = Model();

@defVar(m, 0 <= Gen0[i=GENERATORS] <= max_gen[i]);                 # Initial power generation
@defVar(m, 0 <= Gen[i=GENERATORS, t=PERIODS] <= max_gen[i]);       # Power generation
@defVar(m, -flowmax[l] <= Flow[l=LINES, t=PERIODS] <= flowmax[l]); # transmission line capacity
@defVar(m, -360 <= Angle[n=BUSES, t=PERIODS] <= 360);              # phase angle
@defVar(m, 0 <= LoadShed[j=LOADS,t=PERIODS] <= dictDemand[j][t]); # load shedding
@defVar(m, 0 <= ImportCurtailment[i=IMPORTS,t=PERIODS] <= dictImportProduction[i][t]);     # import curtailment
@defVar(m, 0 <= RenewableCurtailment[r=REGENERATORS,t=PERIODS] <= dictREProduction[r][t]); # renewable curtailment
@defVar(m, 0 <= WindCurtailment[w=WINDS,t=PERIODS] <= dictWindProduction[w][t,scenario]);  # wind curtailment

@setObjective(m, Min, 
	sum{gen_cost[i] * Gen[i,t], i=GENERATORS, t=PERIODS}  
	+ 100 * sum{WindCurtailment[w,t], w=WINDS, t=PERIODS}
	+ 1000 * sum{LoadShed[j,t], j=LOADS, t=PERIODS}
	+ 1000 * sum{ImportCurtailment[i,t], i=IMPORTS, t=PERIODS}
	+ 2000 * sum{RenewableCurtailment[r,t], r=REGENERATORS, t=PERIODS});

# Ramping rate in normal operating status
@addConstraint(m, RAMP_DOWN0[i=GENERATORS], Gen[i,1] >= Gen0[i] - ramp_down[i]);
@addConstraint(m, RAMP_DOWN[i=GENERATORS, t=2:length(PERIODS)], Gen[i,t-1] - Gen[i,t] <= ramp_down[i]);
@addConstraint(m, RAMP_UP0[i=GENERATORS], Gen[i,1] <= ramp_up[i] + Gen0[i]);
@addConstraint(m, RAMP_UP[i=GENERATORS, t=2:length(PERIODS)], Gen[i,t] - Gen[i,t-1] <= ramp_up[i]);

# Power flow equation
@addConstraint(m, POWER_FLOW[l=LINES, t=PERIODS],
	Flow[l,t] == susceptance[l] * (Angle[tobus[l],t] - Angle[frombus[l],t]));

# Power balance constraints for system
@addConstraint(m, POWER_BALANCE[n=BUSES, t=PERIODS],
	sum{Flow[l,t], l=LINES; tobus[l] == n} 
	- sum{Flow[l,t], l=LINES; frombus[l] == n}
	+ sum{Gen[i,t], i=GENERATORS; gen2bus[i] == n} 
	+ sum{LoadShed[j,t], j=LOADS; dBusLoads[j] == n} 
	== sum{dictDemand[j][t], j=LOADS; dBusLoads[j] == n}
	- sum{dictImportProduction[i][t] - ImportCurtailment[i,t], i=IMPORTS; dBusImportPoints[i] == n}
	- sum{dictREProduction[r][t] - RenewableCurtailment[r,t], r=REGENERATORS; dBusREGenerators[r] == n}
	- sum{dictWindProduction[w][t,scenario] - WindCurtailment[w,t], w=WINDS; dBusWindGenerators[w] == n}
	+ dataCenter[n][t]);

