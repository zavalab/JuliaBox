using DataFrames
dataframe = DataFrame(mode = [], cost = [], carbon = [], Th = [], Tj = [], Tc = [], status = [])
if length(ARGS) < 10
    for i in 1:10
        push!(ARGS, "")
    end
else
    ARGS[1:10] .= ""
end

for mode in ["uniform", "small_var", "large_var"]
    for cost in ["0", "1", "10","100"]
        for carbon in ["0", "1", "10","100"]
            ARGS[1:8] = ["24", "24", "24", "shifting", "real_time", mode, cost, carbon]
            try
                include("main.jl")
                push!(dataframe, (mode, cost, carbon, "24", "24", "24", "Done"))
            catch
                push!(dataframe, (mode, cost, carbon, "24", "24", "24", "Infeasible"))
            end
        end
    end
end

for mode in ["uniform", "small_var", "large_var"]
    for cost in ["0", "0.1", "1"]
        for carbon in ["0", "0.1", "1"]
            for t in ["1", "2", "3", "6", "12"]
                ARGS[1:8] = [t, "6", "6", "shifting", "real_time", mode, cost, carbon]
                try
                    include("main.jl")
                    push!(dataframe, (mode, cost, carbon, t, "6", "6", "Done"))
                catch
                    push!(dataframe, (mode, cost, carbon, t, "6", "6", "Infeasible"))
                end
                ARGS[1:8] = ["6", t, "6", "shifting", "real_time", mode, cost, carbon]
                try
                    include("main.jl")
                    push!(dataframe, (mode, cost, carbon, "6", t, "6", "Done"))
                catch
                    push!(dataframe, (mode, cost, carbon, "6", t, "6", "Infeasible"))
                end
                ARGS[1:8] = ["6", "6", t, "shifting", "real_time", mode, cost, carbon]
                try
                    include("main.jl")
                    push!(dataframe, (mode, cost, carbon, "6", "6", t, "Done"))
                catch
                    push!(dataframe, (mode, cost, carbon, "6", "6", t, "Infeasible"))
                end

            end
        end
    end
end


for mode in ["uniform", "small_var", "large_var"]
    for t in ["3", "6", "12", "18", "24"]
        ARGS[1:8] = [t, t, t, "shifting", "real_time", mode, "100", "10"]
        try
            include("main.jl")
            push!(dataframe, (mode, "100", "10", t, t, t, "Done"))
        catch
            push!(dataframe, (mode, "100", "10", t, t, t, "Infeasible"))
        end
    end
end

for mode in ["uniform", "small_var", "large_var"]
    ARGS[1:8] = ["24", "24", "24", "shifting", "day_ahead", mode, "100", "10"]
    try 
        include("main.jl")
        push!(dataframe, (mode * "_day_ahead", "100", "10", "24", "24", "24", "Done"))
    catch
        push!(dataframe, (mode * "_day_ahead", "100", "10", "24", "24", "24", "Infeasible"))
    end
end

for mode in ["uniform", "small_var", "large_var"]
    # for c in ["0.0001", "0.001"]
    for c in ["0.001"]
        ARGS[1:10] = ["24", "24", "24", "shifting", "real_time", mode, "0", "0", "", c]
        try
            include("main.jl")            
            push!(dataframe, (mode * "edr" * c, "0", "0", "24", "24", "24", "Done"))
        catch
            push!(dataframe, (mode * "edr" * c, "0", "0", "24", "24", "24", "Infeasible"))
        end
    end
end

ARGS[1:10] .= ""

for mode in ["uniform", "small_var", "large_var"]
    for t in ["9", "12", "24"]

        ARGS[1:9] = [t, t, t, "shifting", "real_time", mode, "0", "0", "vary_cap"]
        try
            include("main.jl")
            push!(dataframe, (mode * "_vary_cap", "0", "0", t, t, t, "Done"))
        catch
            push!(dataframe, (mode * "_vary_cap", "0", "0", t, t, t, "Infeasible"))
        end
    end
end

for mode in ["uniform", "small_var", "large_var"]
    for t in ["6", "9", "12", "18", "24"]

        ARGS[1:9] = [t, t, t, "shifting", "real_time", mode, "100", "10", "vary_cap"]
        try
            include("main.jl")
            push!(dataframe, (mode * "_vary_cap", "100", "10", t, t, t, "Done"))
        catch
            push!(dataframe, (mode * "_vary_cap", "100", "10", t, t, t, "Infeasible"))
        end
    end
end

# CSV.write("results/statuses_noramping_edr.csv", dataframe)