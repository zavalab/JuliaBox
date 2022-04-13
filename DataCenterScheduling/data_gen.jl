using DataFrames, CSV

function give_hour(ϵ::Float64, p::Vector{Float64})::Int64
    for i in 1:length(p)
        if ϵ < p[i]
            return i
        end
    end
end

function generate_data(mode=:uniform)

    data_path = "/home/wzhang483/scheduling/azure_trace_stats/vmtable.csv";
    headers=["vmid","subscriptionid","deploymentid","vmcreated", "vmdeleted", "maxcpu", "avgcpu", "p95maxcpu", "vmcategory", "vmcorecount", "vmmemory"];
    trace_df = DataFrame(CSV.File(data_path, header=headers));

    trace_df[!, :lifetime] = (max.(trace_df[!, :vmdeleted] - trace_df[!, :vmcreated], 300)) ./ 3600;
    trace_df[!, :lifetime_hr] = Int.(round.(trace_df[!, :lifetime]));
    trace_df[!, :corehour] = trace_df[!, :lifetime] .* trace_df[!, :vmcorecount];

    trace_df = trace_df[trace_df[:, :lifetime] .< 24, :]

    pairs = Dict{Tuple{Int64, Int64},Int64}()
    fractionals = Dict()
    Nvm = size(trace_df)[1]
    for i in 1:Nvm
        ncores = trace_df[i, :vmcorecount]
        if trace_df[!, :lifetime_hr][i] < 1
            if ncores in keys(fractionals)
                push!(fractionals[ncores], i)
            else
                fractionals[ncores] = [i]
            end
        else
            time = trace_df[!, :lifetime_hr][i]
            if (ncores, time) in keys(pairs)
                pairs[(ncores, time)] += 1
            else
                pairs[(ncores, time)] = 1
            end
        end
    end
    
    for (i, _) in fractionals
        total_hr = Int(round(sum(trace_df[!, :lifetime][fractionals[i]])))
        if (i, 1) in keys(pairs)
            pairs[(i, 1)] += total_hr
        else
            paris[(i, 1)] = total_hr
        end
    end


    if mode == :uniform
        daily_prob = [i * 1/24 for i in 1:24]
    else
        if mode == :small_var
            p_var = [3,4,4,5,4,5,4,5,5,6,7,8,7,8,8,9,8,8,6,5,5,4,3,2];
        elseif mode == :large_var
            p_var = [2,1,1,2,4,6,8,12,16,18,19,14,15,20,23,16,13,10,11,8,6,4,3,2];
        else
            error("The input mode is not recognized.")
        end
        daily_prob = [p_var[1]]
        for i in 2:24
            push!(daily_prob, p_var[i] + daily_prob[i-1])
        end
        daily_prob = daily_prob / daily_prob[end];
    end

    N = Dict(t => Dict() for t in 1:720) # Job submission profile
    for (i, j) in pairs
        days = rand(1:30, j)
        samples = rand(j)
        hours = [give_hour(i, daily_prob) for i in samples]
        times = (days .- 1)*24 .+ hours
        times_dict = Dict(i => count(x->x==i, times) for i in unique(times))
        for (time, njobs) in times_dict
            # if !(time in keys(N))
            #     N[time] = Dict()
            # end
            if i in keys(N[time])
                N[time][i] += Int(round(njobs * 0.8))
            else
                N[time][i] = Int(round(njobs * 0.8))
            end
        end
    end

    return pairs, N
end