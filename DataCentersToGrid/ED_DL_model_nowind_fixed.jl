# Julia script
# - A economic dispatch model with fixed dispatchable load allocations at no wind
# Kibaek Kim - ANL/MCS 1/11/2015

maxDispatchable = 200; # 200 MW for each dispatchable load

m = Model();

#@defVar(m, InstallDispatchable[n=BUSES] >= 0, Int);    # number of dispatchable loads at bus n
@defVar(m, DLoad[n=BUSES,t=PERIODS] >= 0);            # dispatch from dispachable loads
@defVar(m, 0 <= Gen0[i=GENERATORS] <= max_gen[i]);                # Initial power generation
@defVar(m, 0 <= Gen[i=GENERATORS,t=PERIODS] <= max_gen[i]);       # Power generation
@defVar(m, -flowmax[l] <= Flow[l=LINES,t=PERIODS] <= flowmax[l]); # transmission line capacity
@defVar(m, -360 <= Angle[n=BUSES,t=PERIODS] <= 360);              # phase angle
@defVar(m, 0 <= LoadShed[j=LOADS,t=PERIODS] <= dictDemand[j][t]); # load shedding
@defVar(m, 0 <= ImportCurtailment[i=IMPORTS,t=PERIODS] <= dictImportProduction[i][t]);     # import curtailment
@defVar(m, 0 <= RenewableCurtailment[r=REGENERATORS,t=PERIODS] <= dictREProduction[r][t]); # renewable curtailment
@defVar(m, 0 <= WindCurtailment[w=WINDS,t=PERIODS] <= 0);        # wind curtailment

@setObjective(m, Min,
	sum{gen_cost[i] * Gen[i,t], i=GENERATORS, t=PERIODS}
	+ 100 * sum{WindCurtailment[w,t], w=WINDS, t=PERIODS}
	+ 1000 * sum{LoadShed[j,t], j=LOADS, t=PERIODS}
	+ 1000 * sum{ImportCurtailment[i,t], i=IMPORTS, t=PERIODS}
	+ 2000 * sum{RenewableCurtailment[r,t], r=REGENERATORS, t=PERIODS});

#@addConstraint(m, NUMDISPLOAD, 
#	sum{InstallDispatchable[n], n=BUSES} <= numDispatchables);

# Dispatchable load
@addConstraint(m, DISPLOAD[n=BUSES,t=PERIODS], 
	DLoad[n,t] <= maxDispatchable * InstallDispatchableFixed[n]);

# Ramping rate in normal operating status
@addConstraint(m, RAMP_DOWN[i=GENERATORS,t=1], Gen0[i] - Gen[i,1] <= ramp_down[i]);
@addConstraint(m, RAMP_DOWN[i=GENERATORS,t=2:length(PERIODS)], Gen[i,t-1] - Gen[i,t] <= ramp_down[i]);
@addConstraint(m, RAMP_UP[i=GENERATORS,t=1], Gen[i,1] - Gen0[i] <= ramp_up[i]);
@addConstraint(m, RAMP_UP[i=GENERATORS,t=2:length(PERIODS)], Gen[i,t] - Gen[i,t-1] <= ramp_up[i]);

# Power flow equation
@addConstraint(m, POWER_FLOW[l=LINES,t=PERIODS],
	Flow[l,t] == susceptance[l] * (Angle[tobus[l],t] - Angle[frombus[l],t]));

# Power balance constraints for system
@addConstraint(m, POWER_BALANCE[n=BUSES,t=PERIODS],
	sum{Flow[l,t], l=LINES; tobus[l] == n} 
	- sum{Flow[l,t], l=LINES; frombus[l] == n}
	+ sum{Gen[i,t], i=GENERATORS; gen2bus[i] == n} 
	- DLoad[n,t]
	== sum{dictDemand[j][t] - LoadShed[j,t], j=LOADS; dBusLoads[j] == n}
	- sum{dictImportProduction[i][t] - ImportCurtailment[i,t], i=IMPORTS; dBusImportPoints[i] == n}
	- sum{dictREProduction[r][t] - RenewableCurtailment[r,t], r=REGENERATORS; dBusREGenerators[r] == n});

