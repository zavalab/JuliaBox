# Julia script for ZC Cloud
# Kibaek Kim - 2015 ANL MCS

if VERSION < VersionNumber(0,4)
	typealias AbstractString String
	typealias UInt8 Uint8
end

# define some common functions
include("CAISO_data_common_functions.jl");

@assert(isdefined(:SEASON));

# Set paths
DATA_DIR = "/lcrc/project/NEXTGENOPT/data/caiso";
WINDCASE = "Case33";
WIND_DIR = "$DATA_DIR/$WINDCASE/$SEASON";

# ---------------
# Read data files
# ---------------

# BUSES
BUSES              = readdlm("$DATA_DIR/Buses.txt");              # list of buses
gen2bus            = readDict("$DATA_DIR/BusGenerators.txt");     # bus to generators
dBusImportPoints   = readDict("$DATA_DIR/BusImportPoints.txt");   # bus to import points
dBusLoads          = readDict("$DATA_DIR/BusLoads.txt");          # bus to loads
dBusREGenerators   = readDict("$DATA_DIR/BusREGenerators.txt");   # bus to RE generators
dBusWindGenerators = readDict("$DATA_DIR/BusWindGenerators.txt"); # bus t wind generators

# Generators
GENERATORS = readdlm("$DATA_DIR/Generators.txt");      # list of generators
max_gen    = readDict("$DATA_DIR/MaxRunCapacity.txt"); # max generation capacity
min_gen    = readDict("$DATA_DIR/MinRunCapacity.txt"); # min generation capacity
ramp_down  = readDict("$DATA_DIR/RampDown.txt");       # ramp down limit
ramp_up    = readDict("$DATA_DIR/RampUp.txt");         # ramp up limit
gen_cost   = readDict("$DATA_DIR/FuelPrice.txt"); # generation cost

# Calculated Netdemand load
LOADS      = readdlm("$DATA_DIR/Loads.txt"); # list of loads
dDemand    = readdlm("$DATA_DIR/Demand$SEASON.txt");
dictDemand = create2DDict(dDemand);

# IMPORTS
IMPORTS              = readdlm("$DATA_DIR/ImportPoints.txt");
dImportProduction    = readdlm("$DATA_DIR/ImportProduction$SEASON.txt");
dictImportProduction = create2DDict(dImportProduction);
dImportGroups        = readdlm("$DATA_DIR/ImportGroups.txt");
dImportLimit         = readdlm("$DATA_DIR/ImportLimit.txt");
dPolarity            = readdlm("$DATA_DIR/Polarity.txt"); # polarity of line for each import group

# Non-wind renewable production
REGENERATORS     = readdlm("$DATA_DIR/REGenerators.txt");
dREProduction    = readdlm("$DATA_DIR/REProduction$SEASON.txt");
dictREProduction = create2DDict(dREProduction);

# Network
LINES       = readdlm("$DATA_DIR/Lines.txt"); # list of lines
frombus     = readDict("$DATA_DIR/FromBus.txt");
tobus       = readDict("$DATA_DIR/ToBus.txt");
flowmax     = readDict("$DATA_DIR/TC.txt"); # line capacity?
susceptance = readDict("$DATA_DIR/Susceptance.txt");

# WINDS
WINDS                    = readdlm("$DATA_DIR/WindGenerators.txt"); # list of wind generators
dWindProductionSamples   = readdlm("$WIND_DIR/WindProductionSamples.txt");

# ADDITIONAL PARAMETERS
nBuses     = length(BUSES);
nPeriods   = 24;
if isdefined(:nScenarios) == false
    nScenarios = 1;
end
prob       = ones(nScenarios) / nScenarios; # equal probabilities
voll       = 1000;                          # value of lost load ($/MWh)

# ADDITIONAL SETS
PERIODS    = 1:nPeriods;
SCENARIOS  = 1:nScenarios;

# Wind production scenarios
# TODO: not really generic
reshapedWind = reshape(dWindProductionSamples[2:24001,2:6], nPeriods, 1000, length(WINDS));
dictWindProduction = Dict{AbstractString,Array{Float64,2}}();
for w in 1:length(WINDS)
    wind = zeros(nPeriods, nScenarios);
    for t in PERIODS
        for s in SCENARIOS
            # This calculates the production level scaled by a given penetration.
            wind[t,s] = reshapedWind[t,s,w] * penetration / 0.15;
        end
    end
    setindex!(dictWindProduction, wind, WINDS[w]);
end

