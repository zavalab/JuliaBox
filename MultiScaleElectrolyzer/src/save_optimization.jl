using DataFrames, JSON, CSV
include("electrolyzer_struct.jl")

function get_val(model, var::Symbol)
    return vec(float.(Array(value.(model[var]))))
end


struct OptimizationResult
    final_value_results::Dict{Symbol, Any}
    mip_gap::Union{Float64, Nothing} # Use Union for optional fields
    objective_value::Union{Float64, Nothing} # Use Union for optional fields
    status::Symbol
    el::Electrolyzer
    time::Float64
    LCOH::Float64
end

function populate_optimization_result(model, solvetime, θ::Electrolyzer)
    obj_val = nothing
    mip_g = nothing
    time = nothing
    LCOH = nothing
    stat = :Unknown # Default status
    final_vals = Dict{Symbol, Any}()

    try
        # Attempt to get objective value
        obj_val = objective_value(model)
        mip_g = relative_gap(model)
        stat = Symbol(termination_status(model)) # Convert to Symbol
        time = solvetime
        vars = [:A, :e_tot, :h, :z_on, :z_sb, :z_off, :z_start, :z_replace]
        for v in vars
            final_vals[v] = get_val(model, v)
        end
        LCOH = value.(model[:LCOH])
        
    catch e
        @error "Error populating OptimizationResult: $e"
        # Set status to reflect error if something goes wrong during data extraction
        stat = :Error
    end

    return OptimizationResult(final_vals, mip_g, obj_val, stat, θ, time, LCOH)
end

function save_optimization_result_to_json(result::OptimizationResult, filename::String)
    # Extract the Electrolyzer's original constructor parameters for robust reconstruction
    electrolyzer_data = Dict(
        :ϕ => result.el.ϕ,
        :ξ_start => result.el.ξ_start,
        :η_overpotential => result.el.η_overpotential,
        :V_i => result.el.V_i,
        :ℓ => result.el.ℓ,
        :η_startup => result.el.η_startup,
        :λ_H => result.el.λ_H,
        :i => result.el.i,
        :λ_CAPEX_Plant => result.el.λ_CAPEX_Plant,
        :λ_CAPEX_Stack => result.el.λ_CAPEX_Stack
        # Only saving input parameters, derived parameters will be re-calculated on load
    )

    # Convert OptimizationResult to a dictionary for JSON serialization
    data = Dict(
        :final_value_results => result.final_value_results,
        :mip_gap => result.mip_gap,
        :objective_value => result.objective_value,
        :status => String(result.status), # Convert Symbol to String for JSON
        :time => result.time,
        :LCOH => result.LCOH,
        :el => electrolyzer_data # Nested electrolyzer data
    )

    open((@__DIR__)*"/"*filename, "w") do f
        JSON.print(f, data, 4) # Use 4-space indentation for readability
    end
    println("OptimizationResult saved to $filename")
end

function read_optimization_result_from_json(filename::String)::OptimizationResult
    local data_dict # Declare data_dict as local to ensure it's accessible outside try block

    try
        json_string = read(filename, String)
        data_dict = JSON.parse(json_string)
    catch e
        error("Error reading or parsing JSON file '$filename': $e")
    end

    # Extract fields and perform necessary type conversions
    final_value_results = Dict{Symbol, Any}()
    if haskey(data_dict, "final_value_results")
        # Ensure keys are Symbols and values are parsed correctly (e.g., numbers, arrays)
        for (k, v) in data_dict["final_value_results"]
            final_value_results[Symbol(k)] = v
        end
    end

    # Handle Nothing for mip_gap and objective_value
    mip_gap = get(data_dict, "mip_gap", nothing)
    objective_value = get(data_dict, "objective_value", nothing)
    solvetime = get(data_dict, "time", nothing)
    LCOH = get(data_dict, "LCOH", nothing)

    # Convert status string back to Symbol, default to :Unknown if not found
    status_str = get(data_dict, "status", "Unknown")
    status_sym = Symbol(status_str)

    # Reconstruct Electrolyzer
    electrolyzer_data_from_json = get(data_dict, "el", nothing)
    if isnothing(electrolyzer_data_from_json)
        error("Missing 'electrolyzer' data in JSON file '$filename'. Cannot reconstruct Electrolyzer.")
    end
    elec = Electrolyzer(
            electrolyzer_data_from_json["ϕ"],
            electrolyzer_data_from_json["ξ_start"],
            electrolyzer_data_from_json["η_overpotential"],
            electrolyzer_data_from_json["V_i"],
            electrolyzer_data_from_json["ℓ"],
            electrolyzer_data_from_json["η_startup"],
            electrolyzer_data_from_json["λ_H"],
            electrolyzer_data_from_json["i"],
            electrolyzer_data_from_json["λ_CAPEX_Plant"],
            electrolyzer_data_from_json["λ_CAPEX_Stack"]
        )

    return OptimizationResult(final_value_results, mip_gap, objective_value, status_sym, elec, solvetime, LCOH)
end

function pop_and_save(model, solvetime, θ::Electrolyzer, filename::String)
    res = populate_optimization_result(model, solvetime, θ)
    save_optimization_result_to_json(res, filename)
end


function get_model_stats(filename)
    """
    Freq. of on operation & Percentage
    Freq. of off operation & Percentage
    Freq. of sb operation & Percentage
    Cold starts
    Warm starts
    """

    local data_dict # Declare data_dict as local to ensure it's accessible outside try block

    try
        json_string = read(filename, String)
        data_dict = JSON.parse(json_string)
    catch e
        error("Error reading or parsing JSON file '$filename': $e")
    end

    # Extract fields and perform necessary type conversions
    final_value_results = Dict{Symbol, Any}()
    if haskey(data_dict, "final_value_results")
        # Ensure keys are Symbols and values are parsed correctly (e.g., numbers, arrays)
        for (k, v) in data_dict["final_value_results"]
            final_value_results[Symbol(k)] = v
        end
    end

    # Get the frequency and percent of operation modes:
    # keys: [:A, :e_tot, :h, :z_on, :z_sb, :z_off, :z_start, :z_replace]
    local count_on = 0
    local count_off = 0
    local count_sb = 0
    local count_warm = 0
    local count_cold = 0
    counts = Dict(zip([:z_on, :z_sb, :z_off, :z_start], [count_on, count_sb, count_off, count_cold]))
    
    
    for (i, count) in counts
        println(typeof(final_value_results))
        println(typeof(final_value_results[i]))
        println(final_value_results[i][1:10])
        for j in Float64.(final_value_results[i])
            # print(typeof(j))
            # break
            if j >0.5
                counts[i] += 1
            end
            # break
        end
    end

    sb_data = final_value_results[:z_sb]
    on_data = final_value_results[:z_on]
    for t in 2:length(on_data)
        if sb_data[t-1] > 0.5 && on_data[t] > 0.5
            count_warm += 1
        end
    end

    count_on = counts[:z_on]
    count_off = counts[:z_off]
    count_sb = counts[:z_sb]
    count_cold = counts[:z_start]
    total_check = length(final_value_results[:z_on])


    println("On: $count_on, $(count_on/total_check*100)%")
    println("Off: $count_off, $(count_off/total_check*100)%")
    println("Sb: $count_sb, $(count_sb/total_check*100)%")
    println("Warm: $count_warm, $(count_warm/total_check*100)%")
    println("Cold: $count_cold, $(count_cold/total_check*100)%")
end