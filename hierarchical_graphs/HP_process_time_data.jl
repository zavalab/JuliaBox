using Plasmo
using DelimitedFiles, CSV, DataFrames

###############################################################################
# The script below creates a CSV of the generators in the different layers
# It also provides the scripts for converting the load data from regional data
# to bus data based on the bus participation factor
###############################################################################

# Read in data
line_data = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Lines.csv"))
bus_data = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Buses.csv"))
gen_data = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/All_Generators.csv"))

# Separate generator datasets into subsets that will be used in modeling
gen_DA = gen_data[((gen_data[:, "DAUC Generator?"] .== 1) .&& (gen_data[:, "Must Use?"] .== 0)), :]
gen_r = gen_data[gen_data[:, "Must Use?"] .== 1, :]
gen_ST = gen_data[gen_data[:, "DAUC Generator?"] .!= 1, :]

#CSV.write((@__DIR__)*"/datasets/nrel118/Generators_DAUC.csv", gen_DA)
#CSV.write((@__DIR__)*"/datasets/nrel118/Generators_STUC.csv", gen_ST)
#CSV.write((@__DIR__)*"/datasets/nrel118/Generators_R.csv", gen_r)

alpha = .5

# Load in the demand data; this data is hourly for three regions
DA_demand_3R = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/DA_load_3R.csv"))
RT_demand_3R = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/RT_load_3R.csv"))
ST_demand_3R = DA_demand_3R .* (1 - alpha) .+ RT_demand_3R .* alpha

# Define function for converting the load demand data from a regional level to a bus level
function get_demand_df(data_df)
    demand_df = DataFrame()
    demand_df.Bus = bus_data[:, "Bus Name"]
    demand_df.Region = bus_data[:, "Region"]

    for i in 1:size(DA_demand_3R, 1)
        demand_df[!, "$i"] .= 0.0
    end

    for i in 1:size(bus_data, 1)
        region_val = bus_data[i, "Region"]
        load_pf = bus_data[i, 3]
        if region_val == "R1"
            demand_df[i, 3:end] = data_df[:, 2] .* load_pf
        elseif region_val == "R2"
            demand_df[i, 3:end] = data_df[:, 3] .* load_pf
        elseif region_val == "R3"
            demand_df[i, 3:end] = data_df[:, 4] .* load_pf
        else
            println("REGION NOT SPECIFIED")
        end
    end

    return demand_df
end

DA_demand = get_demand_df(DA_demand_3R)
RT_demand = get_demand_df(RT_demand_3R)
ST_demand = get_demand_df(ST_demand_3R)

#CSV.write((@__DIR__)*"/datasets/nrel118/DA_bus_demand.csv", DA_demand)
#CSV.write((@__DIR__)*"/datasets/nrel118/RT_bus_demand.csv", RT_demand)
#CSV.write((@__DIR__)*"/datasets/nrel118/ST_bus_demand.csv", ST_demand)

solar_DA = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Solar/DA.csv"))
solar_RT = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Solar/RT.csv"))
solar_ST = (solar_DA .+ solar_RT) ./ 2

#CSV.write((@__DIR__)*"/datasets/nrel118/Solar/ST.csv", solar_ST)

wind_DA = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Wind/DA.csv"))
wind_RT = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Wind/RT.csv"))
wind_ST = (wind_DA .+ wind_RT) ./ 2

#CSV.write((@__DIR__)*"/datasets/nrel118/Wind/ST.csv", wind_ST)
