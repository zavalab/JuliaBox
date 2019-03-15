function ParseFile(file, line_max)
    # Load in data
    data = readdlm(file)

    # find the matrix breaks and parse data
    breaks = find(data[:, 1] .== -999)
    BUSES = data[1:breaks[1] - 1, :]
    GENS = data[breaks[1] + 1:breaks[2] - 1, :]
    LINES = data[breaks[2] + 1:end, :]

    # Get the array dimensions
    num_buses = size(BUSES)[1]

    # Setup line information
    num_lines = length(LINES[:, 1])
    line_set  = collect(1:num_lines)
    snd_bus   = LINES[:, 1]
    rec_bus   = LINES[:, 2]

    # Setup generator information
    bus_gen = GENS[:, 1]
    num_gen = length(bus_gen)
    gen_set = collect(1:num_gen)
    gen_max = GENS[:, 9]
    gen_min = GENS[:, 10]

    # Setup load information
    bus_load = find(BUSES[:, 3] .!= 0)
    num_load = length(bus_load)
    load_set = collect(1:num_load)
    θ_ss = Array{Float64}(BUSES[bus_load, 3])

    # Setup equality matrix
    num_snd_bus = length(snd_bus)
    equal_matrix = zeros(num_buses, num_snd_bus + num_gen + num_load)
    for i = 1:num_buses
        for j = 1:length(snd_bus)
            if snd_bus[j] == i
                equal_matrix[i, j] = -1
            elseif rec_bus[j] == i
                equal_matrix[i, j] = 1
            else
                equal_matrix[i, j] = 0
            end
        end
        for j = 1:num_gen
            if bus_gen[j] == i
                equal_matrix[i, j + num_snd_bus] = 1
            else
                equal_matrix[i, j + num_snd_bus] = 0
            end
        end
        for j = 1:num_load
            if bus_load[j] == i
                equal_matrix[i, j + num_snd_bus + num_gen] = -1
            else
                equal_matrix[i, j + num_snd_bus + num_gen] = 0
            end
        end
    end
    hConsts = zeros(num_buses)

    # Setup inequality matrix
    num_vars = num_lines + num_gen
    inequal_matrix = zeros(2 * num_vars, num_vars)
    for i = 1:num_vars
        inequal_matrix[2 * i - 1, i] = 1
        inequal_matrix[2 * i, i] = -1
    end
    inequal_matrix = hcat(inequal_matrix, zeros(2 * num_vars, num_load))

    # Setup inequality constants
    fConsts = zeros(2 * num_vars)
    for i = 1:num_vars
        if i <= num_lines
            fConsts[2 * i - 1] = -line_max
            fConsts[2 * i] = -line_max
        elseif i <= num_gen + num_lines
            fConsts[2 * i - 1] = -gen_max[i - num_lines]
            fConsts[2 * i] = -gen_min[i - num_lines]
        end
    end

    # Dashboard 1
    Cgen    = gen_set[:] + num_lines
    Rgen    = []
    Cgen    = filter!(e->e∉Rgen, Cgen)
    Cload   = load_set[:] + num_lines + num_gen
    Rload   = load_set[:] + num_lines + num_gen
    Cload   = filter!(e->e∉Rload, Cload)
    control = [line_set, Cgen, Cload]
    sto     = [Rgen, Rload]

    # Prepare output
    hRandoms = [equal_matrix[:, sto[1]] equal_matrix[:, sto[2]]]
    hControls = [equal_matrix[:, control[1]] equal_matrix[:, control[2]] equal_matrix[:, control[3]]]
    fRandoms = [inequal_matrix[:, sto[1]] inequal_matrix[:, sto[2]]]
    fControls = [inequal_matrix[:, control[1]] inequal_matrix[:, control[2]] inequal_matrix[:, control[3]]]

    # return results
    return Dict("fConsts" => fConsts, "fControls" => fControls, "fRandoms" => fRandoms, "hConsts" => hConsts, "hControls" => hControls, "hRandoms" => hRandoms)
end
