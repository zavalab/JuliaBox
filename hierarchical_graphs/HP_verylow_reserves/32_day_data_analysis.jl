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

demand_DA = sum(Matrix(D_DA[:, 3:end]), dims = 1)[:]
demand_ST = sum(Matrix(D_ST[:, 3:end]), dims = 1)[:]
demand_RT = sum(Matrix(D_RT[:, 3:end]), dims = 1)[:]

num_days = 32
DA_len = 24 * num_days
DA_gen_num = size(gen_DA, 1)
com_gens_DA = zeros(DA_len)
gen_m_DA = zeros(DA_len)
gen_p_DA = zeros(DA_len)
D_shed_DA = zeros(DA_len)
on_gens_DA = zeros(DA_len)
off_gens_DA = zeros(DA_len)

for i in 1:num_days
    DAUC_dict = Dict()
    DAUC_vars = all_variables(subgraph_set[1])
    DAUC_vals = readdlm((@__DIR__)*"/results/32d_serial/DAUC_results_day$(i).csv", ',')
    for k in 1:length(DAUC_vars)
        DAUC_dict[DAUC_vars[k]] = DAUC_vals[k]
    end
    for j in 1:24
        DA_subgraph = getsubgraphs(subgraph_set[1])[j]
        for k in 1:DA_gen_num
            bus_name = gen_DA[k, "Bus of Connection"]
            gen_name = gen_DA[k, "Generator Name"]

            x_val = DAUC_dict[DA_subgraph[:bus][bus_name][:x][gen_name]]
            s_val = DAUC_dict[DA_subgraph[:bus][bus_name][:s][gen_name]]
            z_val = DAUC_dict[DA_subgraph[:bus][bus_name][:z][gen_name]]

            index = (i - 1) * 24 + j
            com_gens_DA[index] += x_val
            on_gens_DA[index] += s_val
            off_gens_DA[index] += z_val
        end
        for k in 1:size(bus_data, 1)
            bus_name = bus_data[k, "Bus Name"]
            D_shed_val = DAUC_dict[DA_subgraph[:bus][bus_name][:D_shed]]

            index = (i - 1) * 24 + j
            D_shed_DA[index] += D_shed_val
        end
        for k in 1:size(gen_data_DA, 1)
            bus_name = gen_data_DA[k, "Bus of Connection"]
            gen_name = gen_data_DA[k, "Generator Name"]

            gen_m_val = DAUC_dict[DA_subgraph[:bus][bus_name][:Gen_m][gen_name]]
            gen_p_val = DAUC_dict[DA_subgraph[:bus][bus_name][:Gen_p][gen_name]]

            index = (i - 1) * 24 + j
            gen_m_DA[index] += gen_m_val
            gen_p_DA[index] += gen_p_val
        end
    end
end

gen_tot_DA = gen_m_DA + gen_p_DA

ST_len = 96 * num_days
ST_gen_num = size(gen_ST, 1)
com_gens_ST = zeros(ST_len)
gen_m_ST = zeros(ST_len)
gen_p_ST = zeros(ST_len)
D_shed_ST = zeros(ST_len)
on_gens_ST = zeros(ST_len)
off_gens_ST = zeros(ST_len)

for i in 1:num_days
    for j in 1:8
        STUC_dict = Dict()
        STUC_vals = readdlm((@__DIR__)*"/results/32d_serial/STUC_results$(j)_day$(i).csv", ',')
        STUC_vars = all_variables(subgraph_set[1 + j])
        for k in 1:length(STUC_vars)
            STUC_dict[STUC_vars[k]] = STUC_vals[k]
        end

        for k in 1:12
            ST_subgraph = getsubgraphs(subgraph_set[1 + j])[k]
            for l in 1:ST_gen_num
                bus_name = gen_ST[l, "Bus of Connection"]
                gen_name = gen_ST[l, "Generator Name"]

                x_val = STUC_dict[ST_subgraph[:bus][bus_name][:x][gen_name]]
                s_val = STUC_dict[ST_subgraph[:bus][bus_name][:s][gen_name]]
                z_val = STUC_dict[ST_subgraph[:bus][bus_name][:z][gen_name]]

                index = (i - 1) * 96 + (j - 1) * 12 + k
                com_gens_ST[index] += x_val
                on_gens_ST[index] += s_val
                off_gens_ST[index] += z_val

            end
            for l in 1:size(bus_data, 1)
                bus_name = bus_data[l, "Bus Name"]
                D_shed_val = STUC_dict[ST_subgraph[:bus][bus_name][:D_shed]]

                index = (i - 1) * 96 + (j - 1) * 12 + k
                D_shed_ST[index] += D_shed_val
            end
            for l in 1:size(gen_data, 1)
                bus_name = gen_data[l, "Bus of Connection"]
                gen_name = gen_data[l, "Generator Name"]

                gen_m_val = STUC_dict[ST_subgraph[:bus][bus_name][:Gen_m][gen_name]]
                gen_p_val = STUC_dict[ST_subgraph[:bus][bus_name][:Gen_p][gen_name]]

                index = (i - 1) * 96 + (j - 1) * 12 + k
                gen_m_ST[index] += gen_m_val
                gen_p_ST[index] += gen_p_val
            end
        end
    end
end

gen_tot_ST = gen_m_ST + gen_p_ST

HA_len = 96 * num_days
gen_m_HA = zeros(HA_len)
gen_p_HA = zeros(HA_len)
D_shed_HA = zeros(HA_len)

for i in 1:num_days
    for j in 1:96
        HAED_dict = Dict()
        HAED_vals = readdlm((@__DIR__)*"/results/32d_serial/HAED_results$(j)_day$(i).csv", ',')
        HAED_vars = all_variables(subgraph_set[j + 9])
        for k in 1:length(HAED_vals)
            HAED_dict[HAED_vars[k]] = HAED_vals[k]
        end

        HA_subgraph = getsubgraphs(subgraph_set[9 + j])[1]

        for k in 1:size(bus_data, 1)
            bus_name = bus_data[k, "Bus Name"]
            D_shed_val = HAED_dict[HA_subgraph[:bus][bus_name][:D_shed]]

            index = (i - 1) * 96 + j
            D_shed_HA[index] += D_shed_val
        end
        for k in 1:size(gen_data, 1)
            bus_name = gen_data[k, "Bus of Connection"]
            gen_name = gen_data[k, "Generator Name"]

            gen_m_val = HAED_dict[HA_subgraph[:bus][bus_name][:Gen_m][gen_name]]
            gen_p_val = HAED_dict[HA_subgraph[:bus][bus_name][:Gen_p][gen_name]]

            index = (i - 1) * 96 + j
            gen_m_HA[index] += gen_m_val
            gen_p_HA[index] += gen_p_val
        end
    end
end

gen_tot_HA = gen_m_HA + gen_p_HA

t_DA = 1:1:(num_days * 24)
t_ST = .25:.25:(num_days * 24)
t_HA = .25:.25:(num_days * 24)
DA_results = DataFrame(["time" => t_DA, "demand" => demand_DA, "com_gens" => com_gens_DA, "started" => on_gens_DA,
    "off" => off_gens_DA, "gen_m" => gen_m_DA, "gen_p" => gen_p_DA, "gen_tot" => gen_tot_DA, "D_shed" => D_shed_DA])
#CSV.write((@__DIR__)*"/DA_results_serial_res.csv", DA_results)

ST_results = DataFrame(["time" => t_ST, "demand" => demand_ST, "com_gens" => com_gens_ST, "started" => on_gens_ST,
    "off" => off_gens_ST, "gen_m" => gen_m_ST, "gen_p" => gen_p_ST, "gen_tot" => gen_tot_ST, "D_shed" => D_shed_ST])
#CSV.write((@__DIR__)*"/ST_results_serial_res.csv", ST_results)

HA_results = DataFrame(["time" => t_HA, "demand" => demand_RT, "gen_m" => gen_m_HA, "gen_p" => gen_p_HA, "gen_tot" => gen_tot_HA, "D_shed" => D_shed_HA])
#CSV.write((@__DIR__)*"/HA_results_serial_res.csv", HA_results)

################################################################################################################################################################
################################################################################################################################################################

###################### MONOLITHIC PROBLEM #######################################

################################################################################################################################################################
################################################################################################################################################################

num_days = 32
DA_len = 24 * num_days
DA_gen_num = size(gen_DA, 1)
com_gens_DA_m = zeros(DA_len)
gen_m_DA_m = zeros(DA_len)
gen_p_DA_m = zeros(DA_len)
D_shed_DA_m = zeros(DA_len)
on_gens_DA_m = zeros(DA_len)
off_gens_DA_m = zeros(DA_len)

for i in 1:num_days
    DAUC_dict = Dict()
    DAUC_vars = all_variables(subgraph_set[1])
    DAUC_vals = readdlm((@__DIR__)*"/results/32d_monolith/mipgap10_results1_day$(i).csv", ',')
    for k in 1:length(DAUC_vars)
        DAUC_dict[DAUC_vars[k]] = DAUC_vals[k]
    end
    for j in 1:24
        DA_subgraph = getsubgraphs(subgraph_set[1])[j]
        for k in 1:DA_gen_num
            bus_name = gen_DA[k, "Bus of Connection"]
            gen_name = gen_DA[k, "Generator Name"]

            x_val = DAUC_dict[DA_subgraph[:bus][bus_name][:x][gen_name]]
            s_val = DAUC_dict[DA_subgraph[:bus][bus_name][:s][gen_name]]
            z_val = DAUC_dict[DA_subgraph[:bus][bus_name][:z][gen_name]]

            index = (i - 1) * 24 + j
            com_gens_DA_m[index] += x_val
            on_gens_DA_m[index] += s_val
            off_gens_DA_m[index] += z_val
        end
        for k in 1:size(bus_data, 1)
            bus_name = bus_data[k, "Bus Name"]
            D_shed_val = DAUC_dict[DA_subgraph[:bus][bus_name][:D_shed]]

            index = (i - 1) * 24 + j
            D_shed_DA_m[index] += D_shed_val
        end
        for k in 1:size(gen_data_DA, 1)
            bus_name = gen_data_DA[k, "Bus of Connection"]
            gen_name = gen_data_DA[k, "Generator Name"]

            gen_m_val = DAUC_dict[DA_subgraph[:bus][bus_name][:Gen_m][gen_name]]
            gen_p_val = DAUC_dict[DA_subgraph[:bus][bus_name][:Gen_p][gen_name]]

            index = (i - 1) * 24 + j
            gen_m_DA_m[index] += gen_m_val
            gen_p_DA_m[index] += gen_p_val
        end
    end
end

gen_tot_DA_m = gen_m_DA_m + gen_p_DA_m

ST_len = 96 * num_days
ST_gen_num = size(gen_ST, 1)
com_gens_ST_m = zeros(ST_len)
gen_m_ST_m = zeros(ST_len)
gen_p_ST_m = zeros(ST_len)
D_shed_ST_m = zeros(ST_len)
on_gens_ST_m = zeros(ST_len)
off_gens_ST_m = zeros(ST_len)

for i in 1:num_days
    for j in 1:8
        STUC_dict = Dict()
        STUC_vals = readdlm((@__DIR__)*"/results/32d_monolith/mipgap10_results$(j + 1)_day$(i).csv", ',')
        STUC_vars = all_variables(subgraph_set[1 + j])
        for k in 1:length(STUC_vars)
            STUC_dict[STUC_vars[k]] = STUC_vals[k]
        end

        for k in 1:12
            ST_subgraph = getsubgraphs(subgraph_set[1 + j])[k]
            for l in 1:ST_gen_num
                bus_name = gen_ST[l, "Bus of Connection"]
                gen_name = gen_ST[l, "Generator Name"]

                x_val = STUC_dict[ST_subgraph[:bus][bus_name][:x][gen_name]]
                s_val = STUC_dict[ST_subgraph[:bus][bus_name][:s][gen_name]]
                z_val = STUC_dict[ST_subgraph[:bus][bus_name][:z][gen_name]]

                index = (i - 1) * 96 + (j - 1) * 12 + k
                com_gens_ST_m[index] += x_val
                on_gens_ST_m[index] += s_val
                off_gens_ST_m[index] += z_val

            end
            for l in 1:size(bus_data, 1)
                bus_name = bus_data[l, "Bus Name"]
                D_shed_val = STUC_dict[ST_subgraph[:bus][bus_name][:D_shed]]

                index = (i - 1) * 96 + (j - 1) * 12 + k
                D_shed_ST_m[index] += D_shed_val
            end
            for l in 1:size(gen_data, 1)
                bus_name = gen_data[l, "Bus of Connection"]
                gen_name = gen_data[l, "Generator Name"]

                gen_m_val = STUC_dict[ST_subgraph[:bus][bus_name][:Gen_m][gen_name]]
                gen_p_val = STUC_dict[ST_subgraph[:bus][bus_name][:Gen_p][gen_name]]

                index = (i - 1) * 96 + (j - 1) * 12 + k
                gen_m_ST_m[index] += gen_m_val
                gen_p_ST_m[index] += gen_p_val
            end
        end
    end
end

gen_tot_ST_m = gen_m_ST_m + gen_p_ST_m

HA_len = 96 * num_days
gen_m_HA_m = zeros(HA_len)
gen_p_HA_m = zeros(HA_len)
D_shed_HA_m = zeros(HA_len)

for i in 1:num_days
    for j in 1:96
        HAED_dict = Dict()
        HAED_vals = readdlm((@__DIR__)*"/results/32d_monolith/mipgap10_results$(j + 9)_day$(i).csv", ',')
        HAED_vars = all_variables(subgraph_set[j + 9])
        for k in 1:length(HAED_vals)
            HAED_dict[HAED_vars[k]] = HAED_vals[k]
        end

        HA_subgraph = getsubgraphs(subgraph_set[9 + j])[1]

        for k in 1:size(bus_data, 1)
            bus_name = bus_data[k, "Bus Name"]
            D_shed_val = HAED_dict[HA_subgraph[:bus][bus_name][:D_shed]]

            index = (i - 1) * 96 + j
            D_shed_HA_m[index] += D_shed_val
        end
        for k in 1:size(gen_data, 1)
            bus_name = gen_data[k, "Bus of Connection"]
            gen_name = gen_data[k, "Generator Name"]

            gen_m_val = HAED_dict[HA_subgraph[:bus][bus_name][:Gen_m][gen_name]]
            gen_p_val = HAED_dict[HA_subgraph[:bus][bus_name][:Gen_p][gen_name]]

            index = (i - 1) * 96 + j
            gen_m_HA_m[index] += gen_m_val
            gen_p_HA_m[index] += gen_p_val
        end
    end
end

gen_tot_HA_m = gen_m_HA_m + gen_p_HA_m

t_DA = 1:1:(num_days * 24)
t_ST = .25:.25:(num_days * 24)
t_HA = .25:.25:(num_days * 24)


DA_results_m = DataFrame(["time" => t_DA, "demand" => demand_DA[1:length(t_DA)], "com_gens" => com_gens_DA_m, "started" => on_gens_DA_m,
    "off" => off_gens_DA_m, "gen_m" => gen_m_DA_m, "gen_p" => gen_p_DA_m, "gen_tot" => gen_tot_DA_m, "D_shed" => D_shed_DA_m])
#CSV.write((@__DIR__)*"/DA_results_monolith_int.csv", DA_results_m)

ST_results_m = DataFrame(["time" => t_ST, "demand" => demand_ST[1:length(t_ST)], "com_gens" => com_gens_ST_m, "started" => on_gens_ST_m,
    "off" => off_gens_ST_m, "gen_m" => gen_m_ST_m, "gen_p" => gen_p_ST_m, "gen_tot" => gen_tot_ST_m, "D_shed" => D_shed_ST_m])
#CSV.write((@__DIR__)*"/ST_results_monolith_int.csv", ST_results_m)

HA_results_m = DataFrame(["time" => t_HA, "demand" => demand_RT[1:length(t_HA)], "gen_m" => gen_m_HA_m, "gen_p" => gen_p_HA_m, "gen_tot" => gen_tot_HA_m, "D_shed" => D_shed_HA_m])
#CSV.write((@__DIR__)*"/HA_results_monolith_int.csv", HA_results_m)



plot(t_DA, com_gens_DA, label = "Receding Horizon")
plot!(t_DA, com_gens_DA_m, label = "Monolith")
title!("DAUC Layer")
xlabel!("Time (hr)")
ylabel!("Number of Committed Generators")
#savefig((@__DIR__)*"/June2023_plots/ST_serial_DAUC_com_gens.png")

plot(t_ST, com_gens_ST, label = "Receding Horizon")
plot!(t_ST, com_gens_ST_m, label = "Monolith")
title!("STUC Layer")
xlabel!("Time (hr)")
ylabel!("Number of Committed Generators")
