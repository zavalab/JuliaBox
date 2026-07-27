using JuMP, DataFrames, JLD2
include("electrolyzer.jl")

# Keep your mentor's helper function exactly the same
function get_val(model, var::Symbol)
    return vec(float.(Array(value.(model[var]))))
end

# Keep your mentor's custom data structure exactly the same
struct OptimizationResult
    final_value_results::Dict{Symbol, Any}
    mip_gap::Union{Float64, Nothing} 
    objective_value::Union{Float64, Nothing} 
    status::Symbol
    el::Electrolyzer
    time::Float64
    LCOH::Float64
end

# Keep your mentor's population logic nearly identical
function populate_optimization_result(model, solvetime, θ::Electrolyzer)
    # Use NaN instead of nothing so it stays a Float64 if the try block fails
    obj_val = NaN
    mip_g = NaN
    time = solvetime
    LCOH = NaN
    stat = :Unknown 
    final_vals = Dict{Symbol, Any}()

    try
        obj_val = objective_value(model)
        mip_g = relative_gap(model)
        stat = Symbol(termination_status(model)) 
        
        # Check if model has variables before grabbing them
        vars = [:e_tot, :e_wind, :e_grid, :h, :z_on, :z_sb, :z_off, :z_start]
        for v in vars
            # Fixed to singular 'object_dictionary'
            if haskey(JuMP.object_dictionary(model), v)
                final_vals[v] = get_val(model, v)
            else
                final_vals[v] = Float64[] # Fallback empty array if variable doesn't exist
            end
        end
        
        # Fixed to singular 'object_dictionary'
        if haskey(JuMP.object_dictionary(model), :LCOH)
            LCOH = value(model[:LCOH])
        end
        
    catch e
        @error "Error populating OptimizationResult: $e"
        stat = :Error
    end

    # Now this will never fail due to a type conversion error
    return OptimizationResult(final_vals, mip_g, obj_val, stat, θ, time, LCOH)
end

# ==============================================================================
# THE MAGIC OF JLD2: Massively simplified Save & Read
# ==============================================================================

function save_optimization_result_to_jld2(result::OptimizationResult, filename::String)
    # JLD2 seamlessly saves the entire custom struct and all sub-structs directly!
    jldopen(filename, "w") do file
        file["result_object"] = result
    end
    println("OptimizationResult successfully saved to $filename via JLD2")
end

function read_optimization_result_from_jld2(filename::String)::OptimizationResult
    # Read it back completely assembled without any manual dictionary mapping
    res = jldopen(filename, "r") do file
        return file["result_object"]
    end
    return res
end

# Standard pipeline execution function
function pop_and_save(model, solvetime, θ::Electrolyzer, filename::String)
    res = populate_optimization_result(model, solvetime, θ)
    save_optimization_result_to_jld2(res, filename)
end

# ==============================================================================
# STATS ANALYSIS FUNCTION (Updated for native struct reading)
# ==============================================================================
function get_model_stats(filename)
    # Load the native struct back instantly
    result = read_optimization_result_from_jld2(filename)
    final_value_results = result.final_value_results

    # Setup counters
    counts = Dict(:z_on => 0, :z_sb => 0, :z_off => 0, :z_start => 0)
    count_warm = 0

    # Count operations safely 
    for i in [:z_on, :z_sb, :z_off, :z_start]
        for j in final_value_results[i]
            if j > 0.5
                counts[i] += 1
            end
        end
    end

    # Calculate warm starts
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
    total_check = length(on_data)

    println("\n--- Model Statistics ($filename) ---")
    println("On: $count_on, $(round(count_on/total_check*100, digits=2))%")
    println("Off: $count_off, $(round(count_off/total_check*100, digits=2))%")
    println("Sb: $count_sb, $(round(count_sb/total_check*100, digits=2))%")
    println("Warm Starts: $count_warm")
    println("Cold Starts: $count_cold")
end