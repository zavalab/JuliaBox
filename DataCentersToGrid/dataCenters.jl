# Julia script
# - study impact of adding data centers to the system
# Kibaek Kim - ANL/MCS 10012015

using JuMP;

SEASON      = ARGS[1];
penetration = float(ARGS[2]);
scenario    = int(ARGS[3]);
ndatacenters = int(ARGS[4]);

nScenarios = 1000;

# scale penetration level to account for the additional data center load
original_penetration = penetration;
penetration = penetration * 0.15 / 0.13;

include("ZC_data.jl");

penetration = original_penetration;

# data center
DC = readdlm("DCLOCS.txt")[1:ndatacenters];
dataCenter = Dict{AbstractString, Array{Any,1}}();
for n in BUSES
    if in(n,DC)
        setindex!(dataCenter, ones(nPeriods) * 200, n);
    else
        setindex!(dataCenter, zeros(nPeriods), n);
    end
end

include("ED_DC_model.jl");

status = solve(m);

if status == :Optimal

    # File name
    filename = "/home/kibaekkim/ZCCloudResults/DataCenters/$SEASON-$penetration-$ndatacenters-$scenario";

    # Save Objective Value
    fp = open("$filename-Objval.csv", "w");
    print(fp, getObjectiveValue(m));
    close(fp);

    # Save Gen0 and Gen
    fp = open("$filename-Gen.csv", "w");
    for i in GENERATORS
        print(fp, getValue(Gen0[i]));
        for t in PERIODS
            print(fp, ",", getValue(Gen[i,t]));
        end
        print(fp, "\n");
    end
    close(fp);
    
    # Save LoadShed
    fp = open("$filename-Loadshed.txt", "w");
    for i in LOADS
        print(fp, getValue(LoadShed[i,1]));
        for t in 2:nPeriods
            print(fp, ",", getValue(LoadShed[i,t]));
        end
        print(fp, "\n");
    end
    close(fp);
   
    # Save ImportCurtailment
    fp = open("$filename-ImportCurtailment.txt", "w");
    for i in IMPORTS
        print(fp, getValue(ImportCurtailment[i,1]));
        for t in 2:nPeriods
            print(fp, ",", getValue(ImportCurtailment[i,t]));
        end
        print(fp, "\n");
    end
    close(fp);

    # Save RenewableCurtailment
    fp = open("$filename-RenewableCurtailment.txt", "w");
    for i in REGENERATORS
        print(fp, getValue(RenewableCurtailment[i,1]));
        for t in 2:nPeriods
            print(fp, ",", getValue(RenewableCurtailment[i,t]));
        end
        print(fp, "\n");
    end
    close(fp);

    # Save WindCurtailment
    fp = open("$filename-WindCurtailment.txt", "w");
    for i in WINDS
        print(fp, getValue(WindCurtailment[i,1]));
        for t in 2:nPeriods
            print(fp, ",", getValue(WindCurtailment[i,t]));
        end
        print(fp, "\n");
    end
    close(fp);

    # Save LMP
    fp = open("$filename-LMP.txt", "w");
    for i in BUSES
        print(fp, getDual(POWER_BALANCE[i,1]));
        for t in 2:nPeriods
            print(fp, ",", getDual(POWER_BALANCE[i,t]));
        end
        print(fp, "\n");
    end
    close(fp);
end
