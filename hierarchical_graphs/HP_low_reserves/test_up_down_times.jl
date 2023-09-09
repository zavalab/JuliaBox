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

# Define functions for testing up and down times
function test_up_time(s_vector, z_vector, up_time, gen_num; dt = 1)
    t_points = length(s_vector)
    up_time_points = Int(up_time / dt)

    for i in 1:(t_points - up_time_points)
        if s_vector[i] == 1 && z_vector[i] == 0
            if sum(z_vector[i:(i + up_time_points - 1)]) != 0
                println("$gen_num has a wrong up time starting at time point $(i)!")
            end
        end
    end
end

function test_down_time(s_vector, z_vector, down_time, gen_num; dt = 1)
    t_points = length(s_vector)
    down_time_points = Int(down_time / dt)

    for i in 1:(t_points - down_time_points)
        if z_vector[i] == 1 && s_vector == 0
            if sum(s_vector[i:(i + down_time_points - 1)]) != 0
                println("$gen_num has a wrong down time  starting at time point $(i)!")
            end
        end
    end
end

# Query the s, z, x solutions for DA layer in receding horizon problem
num_days = 32
DA_len = 24 * num_days
DA_gen_num = size(gen_DA, 1)
s_array = zeros(DA_gen_num, DA_len)
z_array = zeros(DA_gen_num, DA_len)
x_array = zeros(DA_gen_num, DA_len)

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
            x_array[k, index] = x_val
            s_array[k, index] = s_val
            z_array[k, index] = z_val
        end
    end
end

# Query the s, z, x solutions for ST layer in receding horizon problem
ST_len = 96 * num_days
ST_gen_num = size(gen_ST, 1)
s_array_ST = zeros(ST_gen_num, ST_len)
z_array_ST = zeros(ST_gen_num, ST_len)
x_array_ST = zeros(ST_gen_num, ST_len)

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
                x_array_ST[l, index] = x_val
                s_array_ST[l, index] = s_val
                z_array_ST[l, index] = z_val
            end
        end
    end
end

# Test up/down times on DA solutions
for i in 1:DA_gen_num
    up_time   = gen_DA[i, "Min Up Time (h)"]
    down_time = gen_DA[i, "Min Down Time (h)"]
    s_vec = s_array[i, :]
    z_vec = z_array[i, :]
    test_up_time(s_vec, z_vec, up_time, i; dt = 1)
    test_down_time(s_vec, z_vec, down_time, i; dt = 1)
end
println("DONE TESTING DAUC SERIAL")

# Test up/down times on ST solutions
for i in 1:ST_gen_num
    up_time = gen_ST[i, "Min Up Time (h)"]
    down_time = gen_ST[i, "Min Down Time (h)"]
    s_vec = s_array_ST[i, :]
    z_vec = z_array_ST[i, :]
    test_up_time(s_vec, z_vec, up_time, i; dt = .25)
    test_down_time(s_vec, z_vec, down_time, i, dt = .25)
end
println("DONE TESTING STUC SERIAL")

# Query the s, z, x solutions for DA layer in monolithic problem
num_days = 32
DA_len = 24 * num_days
DA_gen_num = size(gen_DA, 1)
s_array = zeros(DA_gen_num, DA_len)
z_array = zeros(DA_gen_num, DA_len)
x_array = zeros(DA_gen_num, DA_len)

for i in 1:num_days
    DAUC_dict = Dict()
    DAUC_vars = all_variables(subgraph_set[1])
    DAUC_vals = readdlm((@__DIR__)*"/results/32d_monolith/mipgap5_results1_day$(i).csv", ',')
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
            x_array[k, index] = x_val
            s_array[k, index] = s_val
            z_array[k, index] = z_val
        end
    end
end

# Query the s, z, x solutions for ST layer in monolithic problem
ST_len = 96 * num_days
ST_gen_num = size(gen_ST, 1)
s_array_ST = zeros(ST_gen_num, ST_len)
z_array_ST = zeros(ST_gen_num, ST_len)
x_array_ST = zeros(ST_gen_num, ST_len)

for i in 1:num_days
    for j in 1:8
        STUC_dict = Dict()
        STUC_vars = all_variables(subgraph_set[1 + j])
        STUC_vals = readdlm((@__DIR__)*"/results/32d_monolith/mipgap5_results$(j+1)_day$(i).csv", ',')
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
                x_array_ST[l, index] = x_val
                s_array_ST[l, index] = s_val
                z_array_ST[l, index] = z_val
            end
        end
    end
end

# Test up/down times on DA solutions
for i in 1:DA_gen_num
    up_time   = gen_DA[i, "Min Up Time (h)"]
    down_time = gen_DA[i, "Min Down Time (h)"]
    s_vec = s_array[i, :]
    z_vec = z_array[i, :]
    test_up_time(s_vec, z_vec, up_time, i; dt = 1)
    test_down_time(s_vec, z_vec, down_time, i; dt = 1)
end
println("DONE TESTING DAUC MONOLITH")

# Test up/down times on ST solutions
for i in 1:ST_gen_num
    up_time = gen_ST[i, "Min Up Time (h)"]
    down_time = gen_ST[i, "Min Down Time (h)"]
    s_vec = s_array_ST[i, :]
    z_vec = z_array_ST[i, :]
    test_up_time(s_vec, z_vec, up_time, i; dt = .25)
    test_down_time(s_vec, z_vec, down_time, i, dt = .25)
end
println("DONE TESTING STUC MONOLITH")
