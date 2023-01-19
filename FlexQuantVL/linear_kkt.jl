using JuMP, Gurobi, SCIP

function feasibility(T, buses, branches, gens, loads, vls, vls_in, vls_out)
    # number of components
    nbus = length(buses);
    nbranch = length(branches);
    ngen = length(gens);
    nload = length(loads);
    nvls = length(vls);

    optimizer = optimizer_with_attributes(Gurobi.Optimizer)
    m = Model(optimizer)

    @variable(m, s[1:ngen, 1:T])                      # supply
    @variable(m, f[1:nbranch, 1:T])                   # flow forward
    @variable(m, θ[1:nbus, 1:T])
    @variable(m, δ[1:nvls])

    bals = Dict()
    ineeq_constrs = Dict()
    for tt in 1:T
        # Supply and demand capacity constraints
        for k in 1:ngen
            v = gens[k]
            @constraint(m, s[k,tt] - v["cap"] <= 0)
            @constraint(m, -s[k,tt] <= 0)
            # Ramping constraints
            if tt < T
                @constraint(m, s[k,tt+1] - s[k,tt] - v["ru"] <= 0)
                @constraint(m, s[k,tt] - s[k,tt+1] - v["rd"] <= 0)
            end
        end

        # Flow capacity constraints & DC power flow equations
        for k in 1:nbranch
            v = branches[k]
            @constraint(m, f[k,tt] - v["cap"] <= 0)
            @constraint(m, -f[k,tt] - v["cap"] <= 0)
            @constraint(m, f[k,tt] == v["b"] * (θ[v["f_bus"],tt] - θ[v["t_bus"],tt]))
        end

        # Node balance constraint
        for k in 1:nbus
            bals[(k,tt)] = @constraint(m, sum(s[i,tt] for i in buses[k]["gens"]) - sum(loads[i]["cap"][tt] for i in buses[k]["loads"])
                                + sum(f[i,tt] for i in buses[k]["edges_in"]) - sum(f[i,tt] for i in buses[k]["edges_out"])
                                + sum(δ[v] for v in vls_out[(k,tt)]) - sum(δ[v] for v in vls_in[(k,tt)])
                                == 0)
            
            if length(vls_out[(k,tt)]) + length(vls_in[(k,tt)]) > 0
                @constraint(m, sum(δ[v] for v in vls_out[(k,tt)]) - sum(δ[v] for v in vls_in[(k,tt)]) - sum(loads[i]["cap"][tt] for i in buses[k]["loads"]) <= 0)
            end
        end
    end

    for k in 1:nvls
        v = vls[k]
        @constraint(m, δ[k] - v["cap"] <= 0)
        @constraint(m, -δ[k] <= 0)
    end

    optimize!(m)

    return m
end

function ψ(T, buses, branches, gens, loads, vls, vls_in, vls_out)
    # number of components
    nbus = length(buses);
    nbranch = length(branches);
    ngen = length(gens);
    nload = length(loads);
    nvls = length(vls);

    # number of constraints
    n_s_lb = ngen*T
    n_s_ub = ngen*T
    n_f_lb = nbranch*T
    n_f_ub = nbranch*T
    n_δ_lb = nvls
    n_δ_ub = nvls
    n_d_lb = nbus*T

    n_ineq_constrs = n_s_lb + n_s_ub + n_f_lb + n_f_ub + n_δ_lb + n_δ_ub + n_d_lb

    optimizer = optimizer_with_attributes(Gurobi.Optimizer)
    m = Model(optimizer)

    @variable(m, s[1:ngen, 1:T])                      # supply
    @variable(m, f[1:nbranch, 1:T])                   # flow forward
    @variable(m, θ[1:nbus, 1:T])
    @variable(m, δ[1:nvls])
    @variable(m, t)

    bals = Dict()
    ineeq_constrs = Dict()
    for tt in 1:T
        # Supply and demand capacity constraints
        for k in 1:ngen
            v = gens[k]
            @constraint(m, s[k,tt] - v["cap"] <= t)
            @constraint(m, -s[k,tt] <= t)
            # Ramping constraints
            if tt < T
                @constraint(m, s[k,tt+1] - s[k,tt] - v["ru"] <= t)
                @constraint(m, s[k,tt] - s[k,tt+1] - v["rd"] <= t)
            end
        end

        # Flow capacity constraints & DC power flow equations
        for k in 1:nbranch
            v = branches[k]
            @constraint(m, f[k,tt] - v["cap"] <= t)
            @constraint(m, -f[k,tt] - v["cap"] <= t)
            @constraint(m, f[k,tt] == v["b"] * (θ[v["f_bus"],tt] - θ[v["t_bus"],tt]))
        end

        # Node balance constraint
        for k in 1:nbus
            bals[(k,tt)] = @constraint(m, sum(s[i,tt] for i in buses[k]["gens"]) - sum(loads[i]["cap"][tt] for i in buses[k]["loads"])
                                + sum(f[i,tt] for i in buses[k]["edges_in"]) - sum(f[i,tt] for i in buses[k]["edges_out"])
                                + sum(δ[v] for v in vls_out[(k,tt)]) - sum(δ[v] for v in vls_in[(k,tt)])
                                == 0)
            
            if length(vls_out[(k,tt)]) + length(vls_in[(k,tt)]) > 0
                @constraint(m, sum(δ[v] for v in vls_out[(k,tt)]) - sum(δ[v] for v in vls_in[(k,tt)]) - sum(loads[i]["cap"][tt] for i in buses[k]["loads"]) <= t)
            end
        end
    end

    for k in 1:nvls
        v = vls[k]
        @constraint(m, δ[k] - v["cap"] <= t)
        @constraint(m, -δ[k] <= t)
    end

    @objective(m, Min, t)
    optimize!(m)

    return m, bals
end

function kkt(T, buses, branches, gens, loads, vls, vls_in, vls_out)

    # number of components
    nbus = length(buses);
    nbranch = length(branches);
    ngen = length(gens);
    nload = length(loads);
    nvls = length(vls);

    # println("========================= Problem Data =========================")

    U = 10000

    # number of constraints
    n_s_lb = ngen*T
    n_s_ub = ngen*T
    n_s_ru = ngen*(T-1)
    n_s_rd = ngen*(T-1)
    n_f_lb = nbranch*T
    n_f_ub = nbranch*T
    n_δ_lb = nvls
    n_δ_ub = nvls
    n_d_lb = 0

    ineq_coeffs = Dict("s_lb" => Dict(), "s_ub" => Dict(),
                       "s_ru" => Dict(), "s_rd" => Dict(),
                       "f_lb" => Dict(), "f_ub" => Dict(),
                       "δ_lb" => Dict(), "δ_ub" => Dict(),
                       "d_lb" => Dict())

    ineq_consts = Dict("s_lb" => Dict(), "s_ub" => Dict(),
                       "s_ru" => Dict(), "s_rd" => Dict(),
                       "f_lb" => Dict(), "f_ub" => Dict(),
                       "δ_lb" => Dict(), "δ_ub" => Dict(),
                       "d_lb" => Dict())
        
    for t in 1:T
        for i in 1:ngen
            ineq_coeffs["s_lb"][i,t] = Dict((:p, [i, t]) => -1)
            ineq_coeffs["s_ub"][i,t] = Dict((:p, [i, t]) => 1)
            ineq_consts["s_lb"][i,t] = 0
            ineq_consts["s_ub"][i,t] = -gens[i]["cap"]
            if t < T
                ineq_coeffs["s_ru"][i,t] = Dict((:p, [i, t+1]) => 1, (:p, [i, t]) => -1)
                ineq_coeffs["s_rd"][i,t] = Dict((:p, [i, t+1]) => -1, (:p, [i, t]) => 1)
                ineq_consts["s_ru"][i,t] = -gens[i]["ru"]
                ineq_consts["s_rd"][i,t] = -gens[i]["rd"]
            end
        end
        for i in 1:nbranch
            ineq_coeffs["f_lb"][i,t] = Dict((:f, [i, t]) => -1)
            ineq_coeffs["f_ub"][i,t] = Dict((:f, [i, t]) => 1)
            ineq_consts["f_lb"][i,t] = -branches[i]["cap"]
            ineq_consts["f_ub"][i,t] = -branches[i]["cap"]
        end
        for i in 1:nbus
            if length(vls_out[(i,t)]) + length(vls_in[(i,t)]) > 0
                ineq_coeffs["d_lb"][i,t] = Dict((:δ, [k]) => 1 for k in vls_out[(i,t)])
                for k in vls_in[(i,t)]
                    ineq_coeffs["d_lb"][i,t][(:δ, [k])] = -1
                end
                ineq_consts["d_lb"][i,t] = length(buses[i]["loads"]) > 0 ? -sum(loads[k]["cap"][t] for k in buses[i]["loads"]) : 0
                n_d_lb += 1
            end
        end
    end
    for i in 1:nvls
        ineq_coeffs["δ_lb"][i] = Dict((:δ, [i]) => -1)
        ineq_coeffs["δ_ub"][i] = Dict((:δ, [i]) => 1)
        ineq_consts["δ_lb"][i] = 0
        ineq_consts["δ_ub"][i] = -vls[i]["cap"]
    end

    kvl_coeffs = Dict()
    bals_coeffs = Dict()
    bals_consts = Dict()
    for t in 1:T
        for l in 1:nbranch
            v = branches[l]
            kvl_coeffs[l,t] = Dict((:f, [l,t]) => 1, (:θ, [v["f_bus"],t]) => -v["b"], (:θ, [v["t_bus"],t]) => v["b"])
        end

        for k in 1:nbus
            bals_coeffs[k,t] = Dict()
            for i in buses[k]["gens"]
                bals_coeffs[k,t][(:p, [i,t])] = 1
            end
            for i in buses[k]["edges_in"]
                bals_coeffs[k,t][(:f, [i,t])] = 1
            end
            for i in buses[k]["edges_out"]
                bals_coeffs[k,t][(:f, [i,t])] = -1
            end
            for v in vls_out[(k,t)]
                bals_coeffs[k,t][(:δ, [v])] = 1
            end
            for v in vls_in[(k,t)]
                bals_coeffs[k,t][(:δ, [v])] = -1
            end

            bals_consts[k,t] = length(buses[k]["loads"]) > 0 ? -sum(loads[i]["cap"][t] for i in buses[k]["loads"]) : 0
        end
    end

    n_ineq_constrs = n_s_lb + n_s_ub + n_s_ru + n_s_rd + n_f_lb + n_f_ub + n_δ_lb + n_δ_ub + n_d_lb

    optimizer = optimizer_with_attributes(Gurobi.Optimizer)
    m = Model(optimizer)

    @variable(m, p[1:ngen, 1:T])                      # supply
    @variable(m, f[1:nbranch, 1:T])                   # flow forward
    @variable(m, θ[1:nbus, 1:T])
    @variable(m, δ[1:nvls])
    @variable(m, s[1:n_ineq_constrs] >= 0)
    @variable(m, y[1:n_ineq_constrs], Bin)
    @variable(m, λ[1:n_ineq_constrs] >= 0)
    @variable(m, pai[1:nbus, 1:T])
    @variable(m, μ[1:nbranch, 1:T])
    @variable(m, t)

    bals = Dict()
    ineq_constrs = Dict()
    grad_expr = Dict(i => AffExpr(0) for i in union(p, f, θ, δ))

    ct = 0
    for type in ["s_lb", "s_ub", "s_ru", "s_rd", "f_lb", "f_ub", "d_lb"]
        coeffs = ineq_coeffs[type]
        consts = ineq_consts[type]
        for (i,tt) in keys(coeffs)
            ct += 1
            ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i,tt]) + consts[i,tt] + s[ct] - t == 0)
            for ((var, idx),num) in coeffs[i,tt]
                v = m[var][idx...]
                grad_expr[v] += λ[ct] * num
            end
        end
    end

    for type in ["δ_lb", "δ_ub"]
        coeffs = ineq_coeffs[type]
        consts = ineq_consts[type]
        for i in keys(coeffs)
            ct += 1
            ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i]) + consts[i] + s[ct] - t == 0)
            # ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i]) + consts[i] + s[ct] == 0)
            for ((var, idx),num) in coeffs[i]
                v = m[var][idx...]
                grad_expr[v] += λ[ct] * num
            end
        end
    end

    @constraint(m, sum(λ) == 1)
    # @constraint(m, sum(λ[1:n_ineq_constrs - n_δ_lb - n_δ_ub]) == 1)

    for i in 1:n_ineq_constrs
        @constraint(m, s[i] <= U * (1-y[i]))
        @constraint(m, λ[i] <= y[i])
    end

    for tt in 1:T
        for k in 1:nbranch
            @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in kvl_coeffs[k,tt]) == 0)
            for ((var, idx),num) in kvl_coeffs[k,tt]
                v = m[var][idx...]
                grad_expr[v] += μ[k, tt] * num
            end
        end

        # Node balance constraint
        for k in 1:nbus
            bals[k,tt] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in bals_coeffs[k,tt]) + bals_consts[k,tt] == 0)
            for ((var, idx),num) in bals_coeffs[k,tt]
                v = m[var][idx...]
                grad_expr[v] += pai[k, tt] * num
            end
        end
    end

    for (_,expr) in grad_expr
        @constraint(m, expr == 0)
    end

    # @objective(m, Min, t)

    optimize!(m)

    println("Value of t: ", value(t))

    return m
end

function F(T, buses, branches, gens, loads, vls, vls_in, vls_out;
           economic=nothing, base_cost=0., avg_cost=0., ratio=1.2)

    @assert economic in Set(["total", "average", nothing])

    # number of components
    nbus = length(buses);
    nbranch = length(branches);
    ngen = length(gens);
    nload = length(loads);
    nvls = length(vls);

    U = 10000

    # println("========================= Problem Data =========================")

    # number of constraints
    n_s_lb = ngen*T
    n_s_ub = ngen*T
    n_s_ru = ngen*(T-1)
    n_s_rd = ngen*(T-1)
    n_f_lb = nbranch*T
    n_f_ub = nbranch*T
    n_δ_lb = nvls
    n_δ_ub = nvls
    n_d_lb = 0

    ineq_coeffs = Dict("s_lb" => Dict(), "s_ub" => Dict(),
                       "s_ru" => Dict(), "s_rd" => Dict(),
                       "f_lb" => Dict(), "f_ub" => Dict(),
                       "δ_lb" => Dict(), "δ_ub" => Dict(), 
                       "d_lb" => Dict())

    ineq_consts = Dict("s_lb" => Dict(), "s_ub" => Dict(),
                       "s_ru" => Dict(), "s_rd" => Dict(),
                       "f_lb" => Dict(), "f_ub" => Dict(),
                       "δ_lb" => Dict(), "δ_ub" => Dict(),
                       "d_lb" => Dict())
    
    if !isnothing(economic)
        ineq_coeff_econ_p = Dict()
        ineq_coeff_econ_δ = Dict()
    end

    for t in 1:T
        for i in 1:ngen
            ineq_coeffs["s_lb"][i,t] = Dict((:p, [i, t]) => -1)
            ineq_coeffs["s_ub"][i,t] = Dict((:p, [i, t]) => 1)
            ineq_consts["s_lb"][i,t] = 0
            ineq_consts["s_ub"][i,t] = -gens[i]["cap"]
            if !isnothing(economic)
                ineq_coeff_econ_p[i,t] = gens[i]["bid_cost"][t] 
            end
            if t < T
                ineq_coeffs["s_ru"][i,t] = Dict((:p, [i, t+1]) => 1, (:p, [i, t]) => -1)
                ineq_coeffs["s_rd"][i,t] = Dict((:p, [i, t+1]) => -1, (:p, [i, t]) => 1)
                ineq_consts["s_ru"][i,t] = -gens[i]["ru"]
                ineq_consts["s_rd"][i,t] = -gens[i]["rd"]
            end
        end
        for i in 1:nbranch
            ineq_coeffs["f_lb"][i,t] = Dict((:f, [i, t]) => -1)
            ineq_coeffs["f_ub"][i,t] = Dict((:f, [i, t]) => 1)
            ineq_consts["f_lb"][i,t] = -branches[i]["cap"]
            ineq_consts["f_ub"][i,t] = -branches[i]["cap"]
        end
        for i in 1:nbus
            if length(vls_out[(i,t)]) + length(vls_in[(i,t)]) > 0
                ineq_coeffs["d_lb"][i,t] = Dict((:δ, [k]) => 1 for k in vls_out[(i,t)])
                for k in vls_in[(i,t)]
                    ineq_coeffs["d_lb"][i,t][(:δ, [k])] = -1
                end
                ineq_consts["d_lb"][i,t] = length(buses[i]["loads"]) > 0 ? -sum(loads[k]["cap"][t] for k in buses[i]["loads"]) : 0
                n_d_lb += 1
            end
        end
    end
    for i in 1:nvls
        ineq_coeffs["δ_lb"][i] = Dict((:δ, [i]) => -1)
        ineq_coeffs["δ_ub"][i] = Dict((:δ, [i]) => 1)
        ineq_consts["δ_lb"][i] = 0
        ineq_consts["δ_ub"][i] = -vls[i]["cap"]
        if !isnothing(economic)
            ineq_coeff_econ_δ[i] = vls[i]["bid_cost"]
        end
    end

    kvl_coeffs = Dict()
    bals_coeffs = Dict()
    bals_consts = Dict()
    for t in 1:T
        for l in 1:nbranch
            v = branches[l]
            kvl_coeffs[l,t] = Dict((:f, [l,t]) => 1, (:θ, [v["f_bus"],t]) => -v["b"], (:θ, [v["t_bus"],t]) => v["b"])
        end

        for k in 1:nbus
            bals_coeffs[k,t] = Dict()
            for i in buses[k]["gens"]
                bals_coeffs[k,t][(:p, [i,t])] = 1
            end
            for i in buses[k]["edges_in"]
                bals_coeffs[k,t][(:f, [i,t])] = 1
            end
            for i in buses[k]["edges_out"]
                bals_coeffs[k,t][(:f, [i,t])] = -1
            end
            for v in vls_out[(k,t)]
                bals_coeffs[k,t][(:δ, [v])] = 1
            end
            for v in vls_in[(k,t)]
                bals_coeffs[k,t][(:δ, [v])] = -1
            end

            bals_consts[k,t] = length(buses[k]["loads"]) > 0 ? -sum(loads[i]["cap"][t] for i in buses[k]["loads"]) : 0
        end
    end

    n_ineq_constrs = n_s_lb + n_s_ub + n_s_ru + n_s_rd + n_f_lb + n_f_ub + n_δ_lb + n_δ_ub + n_d_lb
    ineq_names = Dict{Int64, String}()

    optimizer = optimizer_with_attributes(Gurobi.Optimizer)
    m = Model(optimizer)

    @variable(m, p[1:ngen, 1:T])                      # supply
    @variable(m, f[1:nbranch, 1:T])                   # flow forward
    @variable(m, θ[1:nbus, 1:T])
    @variable(m, δ[1:nvls])
    @variable(m, s[1:n_ineq_constrs] >= 0)
    @variable(m, y[1:n_ineq_constrs], Bin)
    @variable(m, λ[1:n_ineq_constrs] >= 0)
    @variable(m, pai[1:nbus, 1:T])
    @variable(m, μ[1:nbranch, 1:T])
    @variable(m, d[1:nbus, 1:T])
    @variable(m, t) # This denotes the parameter for the box constraint now
    if !isnothing(economic)
        @variable(m, s_econ >= 0)
        @variable(m, y_econ, Bin)
        @variable(m, λ_econ >= 0)
    end

    bals = Dict()
    ineq_constrs = Dict()
    grad_expr = Dict(i => AffExpr(0) for i in union(p, f, θ, δ))

    ct = 0
    for type in ["s_lb", "s_ub", "s_ru", "s_rd", "f_lb", "f_ub", "d_lb"]
        coeffs = ineq_coeffs[type]
        consts = ineq_consts[type]
        for (i,tt) in keys(coeffs)
            ct += 1
            ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i,tt]) + consts[i,tt] + s[ct] == 0)
            ineq_names[ct] = type * "[$(i),$(tt)]"
            for ((var, idx),num) in coeffs[i,tt]
                v = m[var][idx...]
                grad_expr[v] += λ[ct] * num
            end
        end
    end

    for type in ["δ_lb", "δ_ub"]
        coeffs = ineq_coeffs[type]
        consts = ineq_consts[type]
        for i in keys(coeffs)
            ct += 1
            ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i]) + consts[i] + s[ct] == 0)
            ineq_names[ct] = type * "[$(i)]"
            for ((var, idx),num) in coeffs[i]
                v = m[var][idx...]
                grad_expr[v] += λ[ct] * num
            end
        end
    end

    if !isnothing(economic)
        @constraint(m, sum(λ) + λ_econ == 1)
    else
        @constraint(m, sum(λ) == 1)
    end
    # @constraint(m, sum(λ[1:n_ineq_constrs - n_δ_lb - n_δ_ub]) == 1)

    for i in 1:n_ineq_constrs
        @constraint(m, s[i] <= U * (1-y[i]))
        @constraint(m, λ[i] <= y[i])
    end
    if economic == "total"
        @constraint(m, s_econ <= 10*base_cost * (1-y_econ))
        @constraint(m, λ_econ <= y_econ)
    elseif economic == "average"
        total_d = sum(sum(v["cap"]) for (_,v) in loads)
        @constraint(m, s_econ <= 10*ratio*avg_cost*total_d * (1-y_econ))
        @constraint(m, λ_econ <= y_econ)
    end

    for tt in 1:T
        for k in 1:nbranch
            @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in kvl_coeffs[k,tt]) == 0)
            for ((var, idx),num) in kvl_coeffs[k,tt]
                v = m[var][idx...]
                grad_expr[v] += μ[k, tt] * num
            end
        end

        # Node balance constraint
        for k in 1:nbus
            # bals[k,tt] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in bals_coeffs[k,tt]) + bals_consts[k,tt] == 0)
            bals[k,tt] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in bals_coeffs[k,tt]) + d[k,tt] == 0)
            @constraint(m, d[k,tt] <= bals_consts[k,tt] * (1-0.5*t))
            @constraint(m, d[k,tt] >= bals_consts[k,tt] * (1+0.5*t))
            for ((var, idx),num) in bals_coeffs[k,tt]
                v = m[var][idx...]
                grad_expr[v] += pai[k, tt] * num
            end
        end
    end

    if economic == "total"
        ineq_constrs[ct+1] = @constraint(m, sum(ineq_coeff_econ_p[i,t] * p[i,t] for i in 1:ngen for t in 1:T) + 
                                            sum(ineq_coeff_econ_δ[v] * δ[v] for v in 1:nvls) - ratio*base_cost + s_econ == 0)
    elseif economic == "average"
        ineq_constrs[ct+1] = @constraint(m, sum(ineq_coeff_econ_p[i,t] * p[i,t] for i in 1:ngen for t in 1:T) + 
                                            sum(ineq_coeff_econ_δ[v] * δ[v] for v in 1:nvls) + ratio*avg_cost*sum(d) + s_econ == 0)
    end

    if !isnothing(economic)
        for t in 1:T
            for i in 1:ngen
                v = m[:p][i,t]
                grad_expr[v] += ineq_coeff_econ_p[i,t] * λ_econ
            end
        end
        for i in 1:nvls
            v = m[:δ][i]
            grad_expr[v] += ineq_coeff_econ_δ[i] * λ_econ
        end
    end

    for (_,expr) in grad_expr
        @constraint(m, expr == 0)
    end

    @objective(m, Min, t)

    # optimize!(m)

    return m, ineq_constrs, ineq_names
end

function F_ellip(T, buses, branches, gens, loads, vls, vls_in, vls_out)

    # number of components
    nbus = length(buses);
    nbranch = length(branches);
    ngen = length(gens);
    nload = length(loads);
    nvls = length(vls);

    U = 10000

    # println("========================= Problem Data =========================")

    # number of constraints
    n_s_lb = ngen*T
    n_s_ub = ngen*T
    n_s_ru = ngen*(T-1)
    n_s_rd = ngen*(T-1)
    n_f_lb = nbranch*T
    n_f_ub = nbranch*T
    n_δ_lb = nvls
    n_δ_ub = nvls
    n_d_lb = 0

    ineq_coeffs = Dict("s_lb" => Dict(), "s_ub" => Dict(),
                       "s_ru" => Dict(), "s_rd" => Dict(),
                       "f_lb" => Dict(), "f_ub" => Dict(),
                       "δ_lb" => Dict(), "δ_ub" => Dict(), 
                       "d_lb" => Dict())

    ineq_consts = Dict("s_lb" => Dict(), "s_ub" => Dict(),
                       "s_ru" => Dict(), "s_rd" => Dict(),
                       "f_lb" => Dict(), "f_ub" => Dict(),
                       "δ_lb" => Dict(), "δ_ub" => Dict(),
                       "d_lb" => Dict())
        
    for t in 1:T
        for i in 1:ngen
            ineq_coeffs["s_lb"][i,t] = Dict((:p, [i, t]) => -1)
            ineq_coeffs["s_ub"][i,t] = Dict((:p, [i, t]) => 1)
            ineq_consts["s_lb"][i,t] = 0
            ineq_consts["s_ub"][i,t] = -gens[i]["cap"]
            if t < T
                ineq_coeffs["s_ru"][i,t] = Dict((:p, [i, t+1]) => 1, (:p, [i, t]) => -1)
                ineq_coeffs["s_rd"][i,t] = Dict((:p, [i, t+1]) => -1, (:p, [i, t]) => 1)
                ineq_consts["s_ru"][i,t] = -gens[i]["ru"]
                ineq_consts["s_rd"][i,t] = -gens[i]["rd"]
            end
        end
        for i in 1:nbranch
            ineq_coeffs["f_lb"][i,t] = Dict((:f, [i, t]) => -1)
            ineq_coeffs["f_ub"][i,t] = Dict((:f, [i, t]) => 1)
            ineq_consts["f_lb"][i,t] = -branches[i]["cap"]
            ineq_consts["f_ub"][i,t] = -branches[i]["cap"]
        end
        for i in 1:nbus
            if length(vls_out[(i,t)]) + length(vls_in[(i,t)]) > 0
                ineq_coeffs["d_lb"][i,t] = Dict((:δ, [k]) => 1 for k in vls_out[(i,t)])
                for k in vls_in[(i,t)]
                    ineq_coeffs["d_lb"][i,t][(:δ, [k])] = -1
                end
                ineq_consts["d_lb"][i,t] = length(buses[i]["loads"]) > 0 ? -sum(loads[k]["cap"][t] for k in buses[i]["loads"]) : 0
                n_d_lb += 1
            end
        end
    end
    for i in 1:nvls
        ineq_coeffs["δ_lb"][i] = Dict((:δ, [i]) => -1)
        ineq_coeffs["δ_ub"][i] = Dict((:δ, [i]) => 1)
        ineq_consts["δ_lb"][i] = 0
        ineq_consts["δ_ub"][i] = -vls[i]["cap"]
    end

    kvl_coeffs = Dict()
    bals_coeffs = Dict()
    bals_consts = Dict()
    for t in 1:T
        for l in 1:nbranch
            v = branches[l]
            kvl_coeffs[l,t] = Dict((:f, [l,t]) => 1, (:θ, [v["f_bus"],t]) => -v["b"], (:θ, [v["t_bus"],t]) => v["b"])
        end

        for k in 1:nbus
            bals_coeffs[k,t] = Dict()
            for i in buses[k]["gens"]
                bals_coeffs[k,t][(:p, [i,t])] = 1
            end
            for i in buses[k]["edges_in"]
                bals_coeffs[k,t][(:f, [i,t])] = 1
            end
            for i in buses[k]["edges_out"]
                bals_coeffs[k,t][(:f, [i,t])] = -1
            end
            for v in vls_out[(k,t)]
                bals_coeffs[k,t][(:δ, [v])] = 1
            end
            for v in vls_in[(k,t)]
                bals_coeffs[k,t][(:δ, [v])] = -1
            end

            bals_consts[k,t] = length(buses[k]["loads"]) > 0 ? -sum(loads[i]["cap"][t] for i in buses[k]["loads"]) : 0
        end
    end

    n_ineq_constrs = n_s_lb + n_s_ub + n_s_ru + n_s_rd + n_f_lb + n_f_ub + n_δ_lb + n_δ_ub + n_d_lb
    ineq_names = Dict{Int64, String}()

    optimizer = optimizer_with_attributes(SCIP.Optimizer)
    m = Model(optimizer)

    @variable(m, p[1:ngen, 1:T])                      # supply
    @variable(m, f[1:nbranch, 1:T])                   # flow forward
    @variable(m, θ[1:nbus, 1:T])
    @variable(m, δ[1:nvls])
    @variable(m, s[1:n_ineq_constrs] >= 0)
    @variable(m, y[1:n_ineq_constrs], Bin)
    @variable(m, λ[1:n_ineq_constrs] >= 0)
    @variable(m, pai[1:nbus, 1:T])
    @variable(m, μ[1:nbranch, 1:T])
    @variable(m, d[1:nbus, 1:T]) ## NOTE: in this code d is NEGATIVE of net demand
    @variable(m, t) # This denotes the parameter for the box constraint now

    bals = Dict()
    ineq_constrs = Dict()
    grad_expr = Dict(i => AffExpr(0) for i in union(p, f, θ, δ))

    ct = 0
    for type in ["s_lb", "s_ub", "s_ru", "s_rd", "f_lb", "f_ub", "d_lb"]
        coeffs = ineq_coeffs[type]
        consts = ineq_consts[type]
        for (i,tt) in keys(coeffs)
            ct += 1
            ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i,tt]) + consts[i,tt] + s[ct] == 0)
            ineq_names[ct] = type * "[$(i),$(tt)]"
            for ((var, idx),num) in coeffs[i,tt]
                v = m[var][idx...]
                grad_expr[v] += λ[ct] * num
            end
        end
    end

    for type in ["δ_lb", "δ_ub"]
        coeffs = ineq_coeffs[type]
        consts = ineq_consts[type]
        for i in keys(coeffs)
            ct += 1
            ineq_constrs[ct] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in coeffs[i]) + consts[i] + s[ct] == 0)
            ineq_names[ct] = type * "[$(i)]"
            for ((var, idx),num) in coeffs[i]
                v = m[var][idx...]
                grad_expr[v] += λ[ct] * num
            end
        end
    end

    @constraint(m, sum(λ) == 1)
    # @constraint(m, sum(λ[1:n_ineq_constrs - n_δ_lb - n_δ_ub]) == 1)

    for i in 1:n_ineq_constrs
        @constraint(m, s[i] <= U * (1-y[i]))
        @constraint(m, λ[i] <= y[i])
    end

    for tt in 1:T
        for k in 1:nbranch
            @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in kvl_coeffs[k,tt]) == 0)
            for ((var, idx),num) in kvl_coeffs[k,tt]
                v = m[var][idx...]
                grad_expr[v] += μ[k, tt] * num
            end
        end

        # Node balance constraint
        for k in 1:nbus
            # bals[k,tt] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in bals_coeffs[k,tt]) + bals_consts[k,tt] == 0)
            bals[k,tt] = @constraint(m, sum(m[var][idx...] * num for ((var, idx),num) in bals_coeffs[k,tt]) + d[k,tt] == 0)
            for ((var, idx),num) in bals_coeffs[k,tt]
                v = m[var][idx...]
                grad_expr[v] += pai[k, tt] * num
            end
        end
    end

    for i in 1:nload
        j = loads[i]["bus"]
        nominal = [bals_consts[j,tt] for tt in 1:T]
        @constraint(m, (d[j,:] .- nominal)' * loads[i]["Vinv"] * (d[j,:] .- nominal) <= t)
    end

    for (_,expr) in grad_expr
        @constraint(m, expr == 0)
    end

    @objective(m, Min, t)

    # optimize!(m)

    return m, ineq_constrs, ineq_names
end