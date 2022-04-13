# carbon emission rate data from: https://dl.acm.org/doi/pdf/10.1145/2342356.2342398
#                                 https://www.nrel.gov/docs/fy13osti/56487.pdf

using JuMP, Gurobi
using Random, Distributions
using DataStructures
using Logging
using CSV, DataFrames
using JLD, HDF5

Random.seed!(53)
include("data_gen.jl")

optimizer = optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0, "MIPGap" => 0.0002)

if ARGS[4] == "shifting"
    shift = true
elseif ARGS[4] == "noshifting"
    shift = false
end

if ARGS[5] == "real_time"
    real_time = true
elseif ARGS[5] == "day_ahead"
    real_time = false
end

mode = Symbol(ARGS[6])
λecost = parse(Float64, ARGS[7])
λcarbon = parse(Float64, ARGS[8])
vary_cap = (ARGS[9] == "vary_cap")
edr = (ARGS[10] != "")
if edr
    λedr = parse(Float64, ARGS[10])
    edr_hrs = Set([13, 36, 61, 83, 104, 133])
    edr_set_points = Dict(13 => 12000, 36 => 18000, 61 => 19000, 83 => 13000, 104 => 19000, 133 => 12000)
    edr_T = 6
end

# Globa range
T = collect(1:24*7) # time frame
Tend = 24*6
Kmax = 20000 # max number of units

N = Dict(t => Dict() for t in T) # Job submission profile

# K level based on price
file = "/home/wzhang483/scheduling/ieee/0725_0801_data.csv"
data = CSV.read(file, DataFrame)
rt_prices = data.rt_prices
da_prices = data.da_prices
fuel_types = [:coal, :gas, :nuclear, :hydro, :wind, :solar]
C_emission_data = Dict(:coal => 968, 
                       :gas => 440, 
                       :nuclear => 15, 
                       :hydro => 13.5, 
                       :wind => 22.5, 
                       :solar => 50
                       ); # unit: kg CO2/MWh
data.cer = [sum(data[t, j]*C_emission_data[j] for j in fuel_types)/sum(data[t, j] for j in fuel_types) for t in 1:size(data)[1]]
cer = data.cer
Random.seed!(21)
dist = Normal(1, 0.11)
deviations = rand(dist, length(data.cer))
cer_pred = data.cer .* deviations

if vary_cap
    drange = 0.1
    stepsizes = rand(MersenneTwister(0), 24 * 7) * 0.3 .- 0.15
    disturbance = zeros(T[end])
    for t in T
        if t == 1
            disturbance[t] = stepsizes[1] * drange * 2
        else
            disturbance[t] = disturbance[t-1] + stepsizes[t] * drange * 2
        end
        disturbance[t] = max(-drange, min(drange, disturbance[t]))
    end
    K = Int.(round.(Kmax * (0.75 .+ disturbance)))
    Random.seed!(36)
    dist = Normal(1, 0.07)
    deviations = rand(dist, length(T))
    K_pred = [min(Int(floor(deviations[i] * K[i])), Kmax) for i in T]
    # K_pred = K
else
    K = [Kmax for _ in T]
end
Ppeak = 100
Pidle = 30

# Randomly distribute jobs to be submitted at the beginning of each week
Random.seed!(53)
total_data, N = generate_data(mode)
KL = collect(keys(total_data))

# gel all possible job lengths
L = sort(unique([l for (k,l) in KL]))
maxL = L[end]

l_to_k = Dict()
for (k,l) in KL
    if !(l in keys(l_to_k))
        l_to_k[l] = [k]
    else
        push!(l_to_k[l], k)
    end
end

# keep N in clear format
for t in T
    for kl in KL
        if !(kl in keys(N[t]))
            N[t][kl] = 0
        end
    end
end

# moving horizon parameters
# we make decision every df hours, with a forecast of incoming load in the next fh hours
# we may also add a capacity forecast component
# and at each time we make decision for the next dh hours
if length(ARGS) > 0
    df = 1
    hl = parse(Int, ARGS[1]) - 1
    jf = parse(Int, ARGS[2]) - 1
    pf = parse(Int, ARGS[3]) - 1
else
    df = 1 # decision frequency (how many hours to realize current schedule)
    pf = 0 # time of power market forecast horizon (for both allocation and price)
    jf = 0 # time of incoming job look-ahead (forecast)
    hl = 0 # moving horizon length
end

#N_h = Int(T[end] / df)
N_h = Tend
last_u = Dict(t => 0 for t in T)
last_hn = Dict((k,l,t) => 0 for (k,l) in KL, t in T)
last_tn = Dict(kl => 0 for kl in KL)
last_n = Dict((kl, t) => 0 for kl in KL, t in T)
last_bn = Dict(kl => 0 for kl in KL)
last_v = Dict((k,l,t) => 0 for (k,l) in KL, t in T)
last_nactive = 0
history_u = Dict()
history_hn = Dict()
history_tn = Dict()
history_n = Dict()
history_bn = Dict()
history_v = Dict()

start_time = time()
canceled_units = []
scheduled_units = []
running_units = [] # number of units running during each interval
carbon_emission = []
electricity_cost = []
lost_prev_timeunits = []
tn_total = []
tn_total_k = []
termination_control = 200

pD = 100
σ = 1
# λecost = 0
# λcarbon = 0
nresults = Dict()

for h in 1:N_h
    start_t = (h-1)*df + 1
    ed_t = start_t + hl
    hT = start_t:ed_t
    hT_ext = start_t:(ed_t+maxL-1)
    @info "============ Start decision horizon $(h):[$(start_t):$(start_t+df-1)], with moving horizon [$(start_t), $(ed_t)]  time: $(round(time() - start_time, digits = 2))============"

    # Formulate and solve the scheduling model
    m = Model(optimizer)

    # Control variables capture decision for the timeframe, although only the first df ones are implemented
    n = @variable(m, [KL, hT], Int, lower_bound = 0, base_name = "n")
    v = Dict()
    for (k,l) in KL, tb in max(1,start_t-l+1):start_t-1
        if tb >= start_t - termination_control ## restrict job termination to jobs started in the latest 24 hours only
            v[(k,l,tb)] = @variable(m, integer=true, lower_bound = 0)
        end
    end

    # State variables capture the states AFTER current time is executed
    u = @variable(m, [1:maxL], Int, lower_bound = 0, base_name = "u")
    fix(u[end], 0, force=true)
    # hn = @variable(m, [KL, T], Int, lower_bound = 0, base_name = "hn")
    hn = Dict()
    for (k,l) in KL, tb in max(1,start_t-l+1):start_t
        hn[(k,l,tb)] = @variable(m, integer=true, lower_bound = 0)
    end
    # hn = @variable(m, [ (k,l,tb) for (k,l) in KL for tb in max(1,start_t-l+1):start_t ], Int, lower_bound = 0, base_name = "hn")
    tn = @variable(m, [KL], Int, lower_bound = 0, base_name = "tn")
    bn = @variable(m, [KL], Int, lower_bound = 0, base_name = "bn")

    # Auxiliary variable for current states (lower-bounded by 0)
    # curr_K = @variable(m, [hT], Int, lower_bound = 0, base_name = "currK")
    nactive = @variable(m, [hT_ext], Int, lower_bound = 0, base_name = "nactive")
    nactive_max = @variable(m, integer = true, base_name = "nactive_max")
    for t in hT
        @constraint(m, nactive_max >= nactive[t])
    end

    # Cost terms
    ecost = @variable(m, [hT_ext], base_name = "ecost")
    carbon = @variable(m, [hT_ext], base_name = "carbon")

    for t in hT_ext
        @constraint(m, ecost[t] == rt_prices[min(t, start_t+pf)] * (nactive[t]/Kmax * (Ppeak-Pidle) + Pidle) )
        if real_time || t == start_t
            @constraint(m, carbon[t] == data.cer[min(t, start_t+pf)] * (nactive[t]/Kmax * (Ppeak-Pidle) + Pidle) )
        else
            @constraint(m, carbon[t] == cer_pred[min(t, start_t+pf)] * (nactive[t]/Kmax * (Ppeak-Pidle) + Pidle) )
        end
    end

    ## job allocation constraints
    for t in hT, kl in KL
        if shift
            @constraint(m, sum(n[kl, s] for s in start_t:t) <= sum(N[s][kl] for s in start_t:min(t, start_t+jf)) + last_tn[kl])
        else
            @constraint(m, n[kl, t] <= N[t][kl]) # case when no shifting (delaying) is allowed
        end
    end
    for kl in KL
        @constraint(m, sum(n[kl, s] for s in hT) >= sum(N[s][kl] for s in start_t:(start_t+min(jf÷2,hl÷2)) ) + last_tn[kl])
        # @constraint(m, sum(n[kl, s] for s in hT) >= N[start_t][kl] + last_tn[kl])
    end

    ## number of active computing units
    # for t in hT
    #     @constraint(m, sum( sum( k * n[(k,l), s] for s in max(start_t, t-l+1):t) for (k,l) in KL)
    #                  + sum(last_u[l] for l in (t-start_t+1):maxL )
    #                  - sum(k*v[(k,l,tb)] for (k,l) in KL, tb in max(1,start_t-l+1,start_t-termination_control):start_t-1 if tb + l - t > 0) == nactive[t])
    # end
    # for t in (ed_t+1):(ed_t+maxL-1)
    #     @constraint(m, sum( sum( k * n[(k,l), s] for s in max(start_t, t-l+1):ed_t) for (k,l) in KL) == nactive[t])
    # end
    for t in hT_ext
        @constraint(m, sum( sum( k * n[(k,l), s] for s in max(start_t, t-l+1):min(ed_t,t)) for (k,l) in KL)
                     + sum(last_u[l] for l in (t-start_t+1):maxL )
                     - sum(k*v[(k,l,tb)] for (k,l) in KL, tb in max(1,start_t-l+1,start_t-termination_control):start_t-1 if tb + l - t > 0) == nactive[t])
    end


    ## Capacity bound
    for t in hT
        if vary_cap && t != start_t
            @constraint(m, nactive[t] <= K_pred[t])
        else
            @constraint(m, nactive[t] <= K[t])
        end
    end

    ## Bound on job cancellation
    # for (k,l) in KL, tb in max(1,start_t-l+1):start_t-1
    for (k,l) in KL, tb in max(1,start_t-l+1,start_t-termination_control):start_t-1
        @constraint(m, v[(k,l,tb)] <= last_hn[(k,l,tb)])
    end

    ## State update
    for l in 1:maxL-1
        expr = AffExpr(0)
        if l+1 in keys(l_to_k)
            expr += sum(k * n[(k,l+1), start_t] for k in l_to_k[l+1])
        end
        # for tb in 1:start_t-1
        for tb in max(1,start_t-termination_control):start_t-1
            if l+start_t-tb+1 in keys(l_to_k)
                expr -= sum(k * v[(k,l+start_t-tb+1,tb)] for k in l_to_k[l+start_t-tb+1])
            end
        end
        @constraint(m, u[l] == last_u[l+1] + expr)
    end

    # for (k,l) in KL, tb in max(1,start_t-l+1):start_t
    for (k,l) in KL, tb in max(1,start_t-l+1,start_t-termination_control):start_t
        if tb + l <= start_t + 1 || l == 1
            fix(hn[(k,l,tb)], 0, force = true)
        else

            if tb == start_t
                @constraint(m, hn[(k,l,tb)] == n[(k,l), tb])
            else
                @constraint(m, hn[(k,l,tb)] == last_hn[(k,l,tb)] - v[(k,l,tb)])
            end
        end
    end

    for (k,l) in KL
        @constraint(m, tn[(k,l)] == last_tn[(k,l)] + N[start_t][(k,l)] - n[(k,l), start_t] + sum(v[(k,l,tb)] for tb in max(1,start_t-l+1,start_t-termination_control):start_t-1) )
    end

    

    for (k,l) in KL
        if l == 1
            @constraint(m, bn[(k,l)] == last_bn[(k,l)] + n[(k,l), start_t])
        # elseif start_t - l + 1 > 0
        elseif start_t - l + 1 > 0
            if start_t - l + 1 > start_t - termination_control
                @constraint(m, bn[(k,l)] == last_bn[(k,l)] + last_hn[(k,l,start_t-l+1)] - v[(k,l,start_t-l+1)])
            else
                @constraint(m, bn[(k,l)] == last_bn[(k,l)] + last_hn[(k,l,start_t-l+1)])
            end
        else
            @constraint(m, bn[(k,l)] == last_bn[(k,l)])
        end
    end

    # @objective(m, Max, (sum( ((ed_t + 1) * k * l - t) * n[(k,l), t] for (k,l) in KL, t in hT) - sum( ((ed_t + 1) * k * l - t) * v[(k,l,t)] for (k,l,t) in keys(v) )) - λecost * nactive_max /10 - λcarbon * sum(carbon)/10);
    if edr
        expr = AffExpr(0)
        for i in 0:edr_T
            if start_t + i in edr_hrs
                edr_hr = start_t + i
                expr -= λedr * (nactive[edr_hr] - edr_set_points[edr_hr])^2
                println("Observing EDR signal for time $(edr_hr)")
                break
            end
        end
        @objective(m, Max, (sum( ((hl + 1) * k * l - (t-start_t)) * n[(k,l), t] for (k,l) in KL, t in hT) - sum( ((hl + 1) * k * l - (t-start_t)) * v[(k,l,t)] for (k,l,t) in keys(v) )) + expr);
    else
        @objective(m, Max, (sum( ((hl + 1) * k * l - (t-start_t)) * n[(k,l), t] for (k,l) in KL, t in hT) - sum( ((hl + 1) * k * l - (t-start_t)) * v[(k,l,t)] for (k,l,t) in keys(v) )) - λecost * nactive_max /10 - λcarbon * sum(carbon)/10);
    end

    @time optimize!(m)

    total_units = 0
    for (k,l) in KL, tb in max(1,start_t-l+1,start_t-termination_control):start_t-1
        total_units += k * value(v[(k,l,tb)])
    end
    @info "  Total number of canceled units: $(total_units) "
    @info "  Termination status: $(termination_status(m))"

    # forward the states
    global last_u = Int.(round.(value.(u)))
    global last_hn = Dict( i => Int(round(value(j))) for (i,j) in hn )
    global last_tn = Int.(round.(value.(tn)))
    global last_bn = Int.(round.(value.(bn)))
    global last_nactive = Int(round(value(nactive[start_t])))

    # record some stats
    if h > 1
        push!(canceled_units, sum(k * Int(round(value(v[k,l,tb]))) for (k,l,tb) in keys(v)))
        push!(lost_prev_timeunits, Int( round( sum( k * (start_t - tb) * value(v[k,l,tb]) for (k,l,tb) in keys(v)) ) ) )
    end
    push!(scheduled_units, sum(k * Int(round(value(n[(k,l),start_t]))) for (k,l) in KL) )
    push!(running_units, Int(round(value(nactive[start_t]))))
    push!(carbon_emission, value(carbon[start_t]))
    push!(electricity_cost, value(ecost[start_t]))
    push!(tn_total, sum(k*l*last_tn[(k,l)] for (k,l) in KL))
    push!(tn_total_k, sum(k*last_tn[(k,l)] for (k,l) in KL))
    nresults[start_t] = [Int(round(value(nactive[i]))) for i in hT]
end

goodput_complete = sum(k * l * last_bn[(k,l)] for (k,l) in KL) / sum(K[1:Tend])
goodput_running = (sum(running_units) - sum(lost_prev_timeunits)) / sum(K[1:Tend])
termination_loss = sum(lost_prev_timeunits) / sum(K[1:Tend])
println("Completed goodput: ", goodput_complete)
println("Running goodput: ", goodput_running)
println("Termination loss: ", termination_loss)
head_string = "miso"
if !shift
    head_string *= "noshift"
end
if !real_time
    head_string *= "day_ahead"
end
if vary_cap
    head_string *= "vary_cap"
end
if edr
    head_string *= "edr$(λedr)_$(edr_T)"
end
dir = "./results/"
jldopen(dir * head_string * "_$(hl+1)_$(jf+1)_$(pf+1)_$(λecost)_$(λcarbon)_$(mode).jld", "w") do f
    write(f, "goodput_complete", goodput_complete)
    write(f, "goodput_running", goodput_running)
    write(f, "termination_loss", termination_loss)
    write(f, "canceled_units", canceled_units)
    write(f, "scheduled_units", scheduled_units)
    write(f, "running_units", running_units)
    write(f, "lost_prev_timeunits", lost_prev_timeunits)
    write(f, "last_bn", last_bn)
    write(f, "last_u", last_u)
    write(f, "last_hn", last_hn)
    write(f, "last_tn", last_tn)
    write(f, "carbon_emission", carbon_emission)
    write(f, "electricity_cost", electricity_cost)
    write(f, "nactive_results", nresults)
    write(f, "K", K)
    write(f, "tn_total", tn_total)
    write(f, "tn_total_k", tn_total_k)
    try
        write(f, "K_pred", K_pred)
    catch
    end
end
