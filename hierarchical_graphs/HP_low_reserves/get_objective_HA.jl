using Plots
using Plasmo
using DelimitedFiles, CSV, DataFrames
using Gurobi

line_data = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Lines.csv"))
bus_data = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Buses.csv"))
gen_data = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/All_Generators.csv"))
gen_r = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Generators_R.csv"))
gen_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Generators_DAUC.csv"))
gen_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Generators_STUC.csv"))

gen_data_conv = vcat(gen_DA, gen_ST)

gen_data_DA = vcat(gen_DA, gen_r)

D_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/DA_bus_demand.csv"))
D_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/ST_bus_demand_15min.csv"))
D_RT = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/RT_bus_demand_15min.csv"))

wind_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Wind/DA.csv"))
wind_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Wind/ST_15min.csv"))
wind_RT = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Wind/RT_15min.csv"))

solar_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Solar/DA.csv"))
solar_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Solar/ST_15min.csv"))
solar_RT = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Solar/RT_15min.csv"))

xi_DA = Dict()
xi_ST = Dict()
xi_RT = Dict()

wind_names = names(wind_DA)[2:end]
solar_names = names(solar_DA)[2:end]

for i in 1:length(wind_names)
    xi_DA[wind_names[i]] = wind_DA[:, wind_names[i]]
    xi_ST[wind_names[i]] = wind_ST[:, wind_names[i]]
    xi_RT[wind_names[i]] = wind_RT[:, wind_names[i]]
end

for i in 1:length(solar_names)
    xi_DA[solar_names[i]] = solar_DA[:, solar_names[i]]
    xi_ST[solar_names[i]] = solar_ST[:, solar_names[i]]
    xi_RT[solar_names[i]] = solar_RT[:, solar_names[i]]
end

B = line_data[:, "Susceptance"]

phi_c = 25
phi_o = 25
phi_u = 5000

include((@__DIR__)*"/layer_construction.jl")
include((@__DIR__)*"/link_solutions.jl")

# Build Dummy graph (used for matching variables to saved solutions)
graph = OptiGraph()

graph_DAUC = OptiGraph()
day_num = 1

build_bus_over_time(graph_DAUC, 25, D_DA[:, 3:27], xi_DA, 1, gen_DA, gen_data_DA, gen_DA)

add_subgraph!(graph, graph_DAUC)

for i in 1:8
    local graph_STUC = OptiGraph()

    t_range = (3 + (i - 1) * 12 + (day_num - 1) * 96):(2 + 16 + (i - 1) * 12 + (day_num - 1) * 96)

    offset = (i - 1) * 12 + (day_num - 1) * 96
    build_bus_over_time(graph_STUC, 16, D_ST[:, t_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = offset)
    add_subgraph!(graph, graph_STUC)

end

gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

for i in 1:96
    local graph_HAED = OptiGraph()
    build_bus_over_time(graph_HAED, 5, D_RT[:, (3 + (i - 1)):(2 + (i - 1) + 5)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1), HA = true)
    add_subgraph!(graph, graph_HAED)
end

all_vars = all_variables(graph)

subgraph_set = getsubgraphs(graph)

# Get the objective based on the performance of the HA-ED levels
# Use only the first data point of each HAED subgraph of each day
# File names are strings corresponding to parts of the file name
function get_objective_HA_level(graph, file_name1, file_name2, file_name3, day_start, day_end; dt = .25, offset = 0)

    obj_val = [0.0]
    # Loop through the number of days you want the objective for
    # Each day, the saved solutions are queried and multiplied by
    # the costs of generation, curtailment, shedding, etc.
    for i in day_start:day_end
        HAED_dict = Dict()
        println("running day $i")
        for j in 1:96
            HAED_vars = all_variables(subgraph_set[9 + j])
            HAED_vals = readdlm((@__DIR__)*file_name1 * "$(j + offset)" * file_name2 * "$(i)" * file_name3, ',')
            HA_subgraph = getsubgraphs(subgraph_set[9 + j])[1]
            for k in 1:length(HAED_vals)
                HAED_dict[HAED_vars[k]] = HAED_vals[k]
            end

            for k in 1:118
                bus_name = bus_data[k, "Bus Name"]
                D_shed_val = HAED_dict[HA_subgraph[:bus][bus_name][:D_shed]]
                obj_val[1] += phi_u * D_shed_val * dt
            end

            for k in 1:size(gen_data, 1)
                var_gen_cost = gen_data[k, "Variable Cost (\$/MWh)"]
                bus_name = gen_data[k, "Bus of Connection"]
                gen_name = gen_data[k, "Generator Name"]
                gen_p_val = HAED_dict[HA_subgraph[:bus][bus_name][:Gen_p][gen_name]]
                gen_m_val = HAED_dict[HA_subgraph[:bus][bus_name][:Gen_m][gen_name]]
                obj_val[1] += var_gen_cost * gen_p_val * dt
                obj_val[1] += phi_o * gen_m_val * dt
            end
        end
    end
    return obj_val
end

monolith_obj_HA = get_objective_HA_level(graph, "/results/32d_monolith/mipgap5_results", "_day", ".csv", 1, 32, offset = 9)
serial_obj_HA   = get_objective_HA_level(graph, "/results/32d_serial/HAED_results", "_day", ".csv", 1, 32)
