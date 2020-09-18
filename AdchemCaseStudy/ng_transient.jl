# Created by Sungho Shin (sungho.shin@wisc.edu) and Kaarthik Sundar (kaarthik@lanl.gov)

function ng_transient_6a(;Nt=1,Mt=1,Ms=1,sig=0.,seed=0,plasmo=false)
    Random.seed!(seed)
    data = parse_files(
        joinpath(@__DIR__,"data/case-6.m"),
        joinpath(@__DIR__,"data/time-series-case-6a.csv");
        total_time=86400*Nt,
        time_step=3600/Mt,
        spatial_discretization=1e4/Ms,
        additional_time=0)

    make_periodic!(data,Mt)
    sig!=0. && add_noise!(data,sig)
    
    gm = GasModels.build_ref(
        data,
        ref_extensions = [GasModels.ref_add_transient!, ref_fix_transient!])
    m, var, con, sol = plasmo ? post_transient_model_plasmo(gm) : post_transient_model(gm) 
    m.ext[:gm]=gm
    return m
end

function ng_transient_6b(;Nt=1,Mt=1,Ms=1,sig=0.,seed=0)
    data = parse_files(
        joinpath(@__DIR__,"data/case-6.m"),
        joinpath(@__DIR__,"data/time-series-case-6b.csv");
        total_time=86400*Nt,
        time_step=3600/Mt,
        spatial_discretization=1e4/Ms,
        additional_time=0)
    
    make_periodic!(data,Mt)
    sig!=0. && add_noise!(data,sig)
    
    gm = GasModels.build_ref(
        data,
        ref_extensions = [GasModels.ref_add_transient!, ref_fix_transient!])
    @time begin
    m, var, con, sol = plasmo ? post_transient_model_plasmo(gm) : post_transient_model(gm) 
    end
    m.ext[:gm]=gm
    return m
end


function make_periodic!(data,Mt)
    for i=1:length(data["nw"])
        for j=1:length(data["nw"]["1"]["transfer"])
            data["nw"]["$i"]["transfer"]["$j"]["withdrawal_max"] =
                data["nw"]["$(mod(i-1,24*Mt)+1)"]["transfer"]["$j"]["withdrawal_max"] 
        end
    end
end

function add_noise!(data,sig)
    for i=1:length(data["nw"])
        for j=1:length(data["nw"]["1"]["transfer"])
            data["nw"]["$i"]["transfer"]["$j"]["withdrawal_max"] *= 1+sig*randn()
        end
    end
end

function post_transient_model(ref)
    m = JuMP.Model()

    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    for (nw_id, nw) in ref[:nw]
        nw_var = var[:nw][nw_id] = Dict{Symbol,Any}()
        nw_con = con[:nw][nw_id] = Dict{Symbol,Any}()
        nw_sol = sol[:nw][nw_id] = Dict{Symbol,Any}()
    end

    time_points = collect(1:length(ref[:nw]))
    start_t = time_points[1]
    end_t = time_points[end]

    # variables for first n-1 time points and expressions 
    for n in time_points[1:end-1]
        nw = ref[:nw][n]

        var[:nw][n][:density] = JuMP.@variable(
            m,
            [i in keys(nw[:junction])],
            lower_bound = nw[:junction][i]["p_min"],
            upper_bound = nw[:junction][i]["p_max"],
            base_name = "density_$(n)"
        )

        var[:nw][n][:pressure] = var[:nw][n][:density]

        var[:nw][n][:pipe_flux_avg] = JuMP.@variable(
            m,
            [i in keys(nw[:pipe])],
            lower_bound = nw[:pipe][i]["flux_min"],
            upper_bound = nw[:pipe][i]["flux_max"],
            base_name = "phi_avg_$(n)"
        )


        var[:nw][n][:pipe_flux_neg] = JuMP.@variable(
            m,
            [i in keys(nw[:pipe])],
            lower_bound = nw[:pipe][i]["flux_min"],
            upper_bound = nw[:pipe][i]["flux_max"],
            base_name = "phi_neg_$(n)"
        )

        var[:nw][n][:pipe_flux_fr] = Dict()
        var[:nw][n][:pipe_flux_to] = Dict()
        for i in keys(nw[:pipe])
            var[:nw][n][:pipe_flux_fr][i] =
                var[:nw][n][:pipe_flux_avg][i] - var[:nw][n][:pipe_flux_neg][i]
            var[:nw][n][:pipe_flux_to][i] =
                var[:nw][n][:pipe_flux_avg][i] + var[:nw][n][:pipe_flux_neg][i]
        end

        var[:nw][n][:pipe_flow_avg] = Dict()
        var[:nw][n][:pipe_flow_neg] = Dict()
        var[:nw][n][:pipe_flow_fr] = Dict()
        var[:nw][n][:pipe_flow_to] = Dict()
        for (i, pipe) in nw[:pipe]
            var[:nw][n][:pipe_flow_avg][i] = var[:nw][n][:pipe_flux_avg][i] * pipe["area"]
            var[:nw][n][:pipe_flow_neg][i] = var[:nw][n][:pipe_flux_neg][i] * pipe["area"]
            var[:nw][n][:pipe_flow_fr][i] = var[:nw][n][:pipe_flux_fr][i] * pipe["area"]
            var[:nw][n][:pipe_flow_to][i] = var[:nw][n][:pipe_flux_to][i] * pipe["area"]
        end 

        var[:nw][n][:compressor_flow] = JuMP.@variable(
            m,
            [i in keys(nw[:compressor])],
            lower_bound = nw[:compressor][i]["flow_min"],
            upper_bound = nw[:compressor][i]["flow_max"],
            base_name = "f_$(n)"
        )

        var[:nw][n][:c_ratio] = JuMP.@variable(
            m,
            [i in keys(nw[:compressor])],
            lower_bound = nw[:compressor][i]["c_ratio_min"],
            upper_bound = nw[:compressor][i]["c_ratio_max"],
            base_name = "alpha_$(n)"
        )

        var[:nw][n][:fs] = JuMP.@variable(
            m,
            [i in keys(nw[:dispatchable_receipt])],
            lower_bound = nw[:receipt][i]["injection_min"],
            upper_bound = nw[:receipt][i]["injection_max"],
            base_name = "s_$(n)"
        )

        var[:nw][n][:fd] = JuMP.@variable(
            m,
            [i in keys(nw[:dispatchable_delivery])],
            lower_bound = nw[:delivery][i]["withdrawal_min"],
            upper_bound = nw[:delivery][i]["withdrawal_max"],
            base_name = "d_$(n)"
        )

        var[:nw][n][:fts] = JuMP.@variable(
            m, 
            [i in keys(nw[:dispatchable_transfer])],
            lower_bound = 0.0,
            upper_bound = max(0.0, -nw[:transfer][i]["withdrawal_min"]),
            base_name = "ts_$(n)"
        )

        var[:nw][n][:ftd] = JuMP.@variable(
            m, 
            [i in keys(nw[:dispatchable_transfer])],
            lower_bound = 0.0,
            upper_bound = max(nw[:transfer][i]["withdrawal_max"], 0.0),
            base_name = "ts_$(n)"
        )

        var[:nw][n][:ft] = Dict() 
        for i in keys(nw[:dispatchable_transfer])
            d = var[:nw][n][:ftd][i]
            s = var[:nw][n][:fts][i]
            var[:nw][n][:ft][i] = d - s
        end 
        

        var[:nw][n][:compressor_power] = Dict()
        for (i, compressor) in nw[:compressor]
            a = var[:nw][n][:c_ratio][i]
            f = var[:nw][n][:compressor_flow][i]
            m1 =
                (ref[:specific_heat_capacity_ratio] - 1) /
                ref[:specific_heat_capacity_ratio]
            W = 286.76 * ref[:temperature] / ref[:gas_specific_gravity] / m1
            var[:nw][n][:compressor_power][i] =
                JuMP.@NLexpression(m, W * abs(f) * (a^m1 - 1.0))
        end

        var[:nw][n][:net_nodal_injection] = Dict()
        var[:nw][n][:net_nodal_edge_out_flow] = Dict()
        for (i, junction) in nw[:junction]
            var[:nw][n][:net_nodal_injection][i] = 0
            for j in nw[:dispatchable_receipts_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] += var[:nw][n][:fs][j]
            end
            for j in nw[:dispatchable_deliveries_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -= var[:nw][n][:fd][j]
            end
            for j in nw[:dispatchable_transfers_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -= var[:nw][n][:ft][j]
            end
            for j in nw[:nondispatchable_receipts_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] += nw[:receipt][j]["injection_nominal"]
            end
            for j in nw[:nondispatchable_deliveries_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -=
                    nw[:delivery][j]["withdrawal_nominal"]
            end
            for j in nw[:nondispatchable_transfers_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -=
                    nw[:transfer][j]["withdrawal_nominal"]
            end
        end

        for (i, junction) in nw[:junction]
            var[:nw][n][:net_nodal_edge_out_flow][i] = 0
            for j in nw[:pipes_fr][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] +=
                    (var[:nw][n][:pipe_flux_fr][j] * nw[:pipe][j]["area"])
            end
            for j in nw[:compressors_fr][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] += var[:nw][n][:compressor_flow][j]
            end
            for j in nw[:pipes_to][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] -=
                    (var[:nw][n][:pipe_flux_to][j] * nw[:pipe][j]["area"])
            end
            for j in nw[:compressors_to][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] -= var[:nw][n][:compressor_flow][j]
            end
        end
    end

    # enforcing time-periodicity without adding additional variables (pointers)
    var[:nw][end_t][:density] = var[:nw][start_t][:density]
    var[:nw][end_t][:pipe_flux_avg] = var[:nw][start_t][:pipe_flux_avg]
    var[:nw][end_t][:pipe_flux_neg] = var[:nw][start_t][:pipe_flux_neg]
    var[:nw][end_t][:pipe_flux_fr] = var[:nw][start_t][:pipe_flux_fr]
    var[:nw][end_t][:pipe_flux_to] = var[:nw][start_t][:pipe_flux_to]
    var[:nw][end_t][:pipe_flow_avg] = var[:nw][start_t][:pipe_flow_avg]
    var[:nw][end_t][:pipe_flow_neg] = var[:nw][start_t][:pipe_flow_neg]
    var[:nw][end_t][:pipe_flow_fr] = var[:nw][start_t][:pipe_flow_fr]
    var[:nw][end_t][:pipe_flow_to] = var[:nw][start_t][:pipe_flow_to]
    var[:nw][end_t][:compressor_flow] = var[:nw][start_t][:compressor_flow]
    var[:nw][end_t][:c_ratio] = var[:nw][start_t][:c_ratio]
    var[:nw][end_t][:fs] = var[:nw][start_t][:fs]
    var[:nw][end_t][:fd] = var[:nw][start_t][:fd]
    var[:nw][end_t][:ft] = var[:nw][start_t][:ft]
    var[:nw][end_t][:compressor_power] = var[:nw][start_t][:compressor_power]
    var[:nw][end_t][:net_nodal_injection] = var[:nw][start_t][:net_nodal_injection]
    var[:nw][end_t][:net_nodal_edge_out_flow] = var[:nw][start_t][:net_nodal_edge_out_flow]

    # density derivatives
    for n in time_points[1:end-1]
        nw = ref[:nw][n]
        prev = n - 1
        (n == start_t) && (prev = time_points[end-1])
        var[:nw][n][:density_derivative] = Dict{Int,Any}()
        for (i, junction) in nw[:junction]
            if junction["junction_type"] == 1
                var[:nw][n][:density_derivative][i] = 0
            else
                var[:nw][n][:density_derivative][i] =
                    (var[:nw][n][:density][i] - var[:nw][prev][:density][i]) /
                    ref[:time_step]
            end
        end
    end

    # constraints 
    for n in time_points[1:end-1]
        nw = ref[:nw][n]

        con[:nw][n][:slack_density] = Dict()
        for (i, junction) in nw[:slack_junctions]
            rho = var[:nw][n][:density][i]
            rho_fixed = junction["p_nominal"]
            con[:nw][n][:slack_density][i] = JuMP.@constraint(m, rho == junction["p_nominal"])
        end

        con[:nw][n][:nodal_balance] = Dict()
        for (i, junction) in nw[:junction]
            net_injection = var[:nw][n][:net_nodal_injection][i]
            net_nodal_edge_out_flow = var[:nw][n][:net_nodal_edge_out_flow][i]
            con[:nw][n][:nodal_balance][i] =
                JuMP.@constraint(m, net_injection == net_nodal_edge_out_flow)
        end

        con[:nw][n][:compressor_flow_dir] = Dict()
        con[:nw][n][:compressor_boost] = Dict()
        con[:nw][n][:compressor_power] = Dict()
        for (i, compressor) in nw[:compressor]
            p_fr = var[:nw][n][:density][compressor["fr_junction"]]
            p_to = var[:nw][n][:density][compressor["to_junction"]]
            a = var[:nw][n][:c_ratio][i]
            f = var[:nw][n][:compressor_flow][i]
            power = var[:nw][n][:compressor_power][i]
            power_max = compressor["power_max"]
            con[:nw][n][:compressor_boost][i] = JuMP.@constraint(m, p_to == a * p_fr)
            con[:nw][n][:compressor_flow_dir][i] = JuMP.@constraint(m, f * (p_fr - p_to) <= 0)
            con[:nw][n][:compressor_power][i] = JuMP.@NLconstraint(m, power <= power_max)
        end

        con[:nw][n][:pipe_momentum_balance] = Dict()
        con[:nw][n][:pipe_mass_balance] = Dict()
        for (i, pipe) in nw[:pipe]
            p_fr = var[:nw][n][:density][pipe["fr_junction"]]
            p_to = var[:nw][n][:density][pipe["to_junction"]]
            resistance = pipe["resistance"]
            f = var[:nw][n][:pipe_flux_avg][i]
            con[:nw][n][:pipe_momentum_balance][i] =
                JuMP.@NLconstraint(m, p_fr^2 - p_to^2 - resistance * f * abs(f) == 0)

            p_fr_dot = var[:nw][n][:density_derivative][pipe["fr_junction"]]
            p_to_dot = var[:nw][n][:density_derivative][pipe["to_junction"]]
            L = pipe["length"]
            phi = var[:nw][n][:pipe_flux_neg][i]
            con[:nw][n][:pipe_mass_balance] =
                JuMP.@constraint(m, L * (p_fr_dot + p_to_dot) + 4 * phi == 0)
        end
    end

    econ_weight = ref[:economic_weighting]
    if econ_weight == 1.0
        load_shed_expression = 0
        for n in time_points[1:end-1]
            for (i, receipt) in ref[:nw][n][:dispatchable_receipt]
                load_shed_expression += (receipt["offer_price"] * var[:nw][n][:fs][i])
            end
            for (i, delivery) in ref[:nw][n][:dispatchable_delivery]
                load_shed_expression -= (delivery["bid_price"] * var[:nw][n][:fd][i])
            end
            for (i, transfer) in ref[:nw][n][:dispatchable_transfer]
                load_shed_expression += (
                    transfer["offer_price"] * var[:nw][n][:fts][i] -
                    transfer["bid_price"] * var[:nw][n][:ftd][i]
                )
            end
        end 
        JuMP.@objective(m, Min, load_shed_expression)
    elseif econ_weight == 0.0 
        compressor_power_expressions = []
        for n in time_points[1:end-1]
            for (i, compressor) in ref[:nw][n][:compressor]
                push!(compressor_power_expressions, var[:nw][n][:compressor_power][i])
            end
        end
        if length(compressor_power_expressions) != 0
            JuMP.@NLobjective(
                m,
                Min,
                sum(
                    compressor_power_expressions[i]
                    for i = 1:length(compressor_power_expressions)
                )
            )
        end
    else 
        load_shed_expressions = []
        compressor_power_expressions = []
        for n in time_points[1:end-1]
            for (i, receipt) in ref[:nw][n][:dispatchable_receipt]
                push!(
                    load_shed_expressions,
                    JuMP.@NLexpression(
                        m,
                        receipt["offer_price"] * var[:nw][n][:fs][i]
                    )
                )
            end
            for (i, delivery) in ref[:nw][n][:dispatchable_delivery]
                push!(
                    load_shed_expressions,
                    JuMP.@NLexpression(
                        m,
                        -delivery["bid_price"] * var[:nw][n][:fd][i]
                    )
                )
            end
            for (i, transfer) in ref[:nw][n][:dispatchable_transfer]
                push!(
                    load_shed_expressions,
                    JuMP.@NLexpression(
                        m,
                        transfer["offer_price"] * var[:nw][n][:fts][i] -
                        transfer["bid_price"] * var[:nw][n][:ftd][i]
                    )
                )
            end
            for (i, compressor) in ref[:nw][n][:compressor]
                push!(compressor_power_expressions, var[:nw][n][:compressor_power][i])
            end
        end
        if length(load_shed_expressions) != 0 && length(compressor_power_expressions) != 0
            JuMP.@NLobjective(
                m,
                Min,
                econ_weight *
                sum(load_shed_expressions[i] for i = 1:length(load_shed_expressions)) +
                (1 - econ_weight) *
                sum(compressor_power_expressions[i] for i = 1:length(compressor_power_expressions))
            )
        end
    end

    return m, var, con, sol
end


function post_transient_model_plasmo(ref)
    graph = Plasmo.OptiGraph()
    
    node= Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    edge= Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    for (nw_id, nw) in ref[:nw]
        var[:nw][nw_id] = Dict{Symbol,Any}()
        con[:nw][nw_id] = Dict{Symbol,Any}()
        sol[:nw][nw_id] = Dict{Symbol,Any}()
    end



    time_points = collect(1:length(ref[:nw]))
    start_t = time_points[1]
    end_t = time_points[end]

    # create nodes
    for n in time_points[1:end-1]
        node[:nw][n] = Plasmo.add_node!(graph)
    end
    for n in time_points[1:end-1]
        edge[:nw][n] = Plasmo.add_edge!(graph,[node[:nw][n],node[:nw][n!=1 ? n-1 : length(ref[:nw])-1]])
    end

    # variables for first n-1 time points and expressions 
    Threads.@threads for n in time_points[1:end-1]
        nw = ref[:nw][n]

        var[:nw][n][:density] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:junction])],
            lower_bound = nw[:junction][i]["p_min"],
            upper_bound = nw[:junction][i]["p_max"],
            base_name = "density_$(n)"
        )

        var[:nw][n][:pressure] = var[:nw][n][:density]

        var[:nw][n][:pipe_flux_avg] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:pipe])],
            lower_bound = nw[:pipe][i]["flux_min"],
            upper_bound = nw[:pipe][i]["flux_max"],
            base_name = "phi_avg_$(n)"
        )


        var[:nw][n][:pipe_flux_neg] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:pipe])],
            lower_bound = nw[:pipe][i]["flux_min"],
            upper_bound = nw[:pipe][i]["flux_max"],
            base_name = "phi_neg_$(n)"
        )

        var[:nw][n][:pipe_flux_fr] = Dict()
        var[:nw][n][:pipe_flux_to] = Dict()
        for i in keys(nw[:pipe])
            var[:nw][n][:pipe_flux_fr][i] =
                var[:nw][n][:pipe_flux_avg][i] - var[:nw][n][:pipe_flux_neg][i]
            var[:nw][n][:pipe_flux_to][i] =
                var[:nw][n][:pipe_flux_avg][i] + var[:nw][n][:pipe_flux_neg][i]
        end

        var[:nw][n][:pipe_flow_avg] = Dict()
        var[:nw][n][:pipe_flow_neg] = Dict()
        var[:nw][n][:pipe_flow_fr] = Dict()
        var[:nw][n][:pipe_flow_to] = Dict()
        for (i, pipe) in nw[:pipe]
            var[:nw][n][:pipe_flow_avg][i] = var[:nw][n][:pipe_flux_avg][i] * pipe["area"]
            var[:nw][n][:pipe_flow_neg][i] = var[:nw][n][:pipe_flux_neg][i] * pipe["area"]
            var[:nw][n][:pipe_flow_fr][i] = var[:nw][n][:pipe_flux_fr][i] * pipe["area"]
            var[:nw][n][:pipe_flow_to][i] = var[:nw][n][:pipe_flux_to][i] * pipe["area"]
        end 

        var[:nw][n][:compressor_flow] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:compressor])],
            lower_bound = nw[:compressor][i]["flow_min"],
            upper_bound = nw[:compressor][i]["flow_max"],
            base_name = "f_$(n)"
        )

        var[:nw][n][:c_ratio] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:compressor])],
            lower_bound = nw[:compressor][i]["c_ratio_min"],
            upper_bound = nw[:compressor][i]["c_ratio_max"],
            base_name = "alpha_$(n)"
        )

        var[:nw][n][:fs] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:dispatchable_receipt])],
            lower_bound = nw[:receipt][i]["injection_min"],
            upper_bound = nw[:receipt][i]["injection_max"],
            base_name = "s_$(n)"
        )

        var[:nw][n][:fd] = JuMP.@variable(
            node[:nw][n],
            [i in keys(nw[:dispatchable_delivery])],
            lower_bound = nw[:delivery][i]["withdrawal_min"],
            upper_bound = nw[:delivery][i]["withdrawal_max"],
            base_name = "d_$(n)"
        )

        var[:nw][n][:fts] = JuMP.@variable(
            node[:nw][n], 
            [i in keys(nw[:dispatchable_transfer])],
            lower_bound = 0.0,
            upper_bound = max(0.0, -nw[:transfer][i]["withdrawal_min"]),
            base_name = "ts_$(n)"
        )

        var[:nw][n][:ftd] = JuMP.@variable(
            node[:nw][n], 
            [i in keys(nw[:dispatchable_transfer])],
            lower_bound = 0.0,
            upper_bound = max(nw[:transfer][i]["withdrawal_max"], 0.0),
            base_name = "ts_$(n)"
        )

        var[:nw][n][:ft] = Dict() 
        for i in keys(nw[:dispatchable_transfer])
            d = var[:nw][n][:ftd][i]
            s = var[:nw][n][:fts][i]
            var[:nw][n][:ft][i] = d - s
        end 
        

        var[:nw][n][:compressor_power] = Dict()
        for (i, compressor) in nw[:compressor]
            a = var[:nw][n][:c_ratio][i]
            f = var[:nw][n][:compressor_flow][i]
            m1 =
                (ref[:specific_heat_capacity_ratio] - 1) /
                ref[:specific_heat_capacity_ratio]
            W = 286.76 * ref[:temperature] / ref[:gas_specific_gravity] / m1
            var[:nw][n][:compressor_power][i] =
                JuMP.@NLexpression(node[:nw][n].model, W * abs(f) * (a^m1 - 1.0))
        end

        var[:nw][n][:net_nodal_injection] = Dict()
        var[:nw][n][:net_nodal_edge_out_flow] = Dict()
        for (i, junction) in nw[:junction]
            var[:nw][n][:net_nodal_injection][i] = 0
            for j in nw[:dispatchable_receipts_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] += var[:nw][n][:fs][j]
            end
            for j in nw[:dispatchable_deliveries_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -= var[:nw][n][:fd][j]
            end
            for j in nw[:dispatchable_transfers_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -= var[:nw][n][:ft][j]
            end
            for j in nw[:nondispatchable_receipts_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] += nw[:receipt][j]["injection_nominal"]
            end
            for j in nw[:nondispatchable_deliveries_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -=
                    nw[:delivery][j]["withdrawal_nominal"]
            end
            for j in nw[:nondispatchable_transfers_in_junction][i]
                var[:nw][n][:net_nodal_injection][i] -=
                    nw[:transfer][j]["withdrawal_nominal"]
            end
        end

        for (i, junction) in nw[:junction]
            var[:nw][n][:net_nodal_edge_out_flow][i] = 0
            for j in nw[:pipes_fr][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] +=
                    (var[:nw][n][:pipe_flux_fr][j] * nw[:pipe][j]["area"])
            end
            for j in nw[:compressors_fr][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] += var[:nw][n][:compressor_flow][j]
            end
            for j in nw[:pipes_to][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] -=
                    (var[:nw][n][:pipe_flux_to][j] * nw[:pipe][j]["area"])
            end
            for j in nw[:compressors_to][i]
                var[:nw][n][:net_nodal_edge_out_flow][i] -= var[:nw][n][:compressor_flow][j]
            end
        end
    end

    # enforcing time-periodicity without adding additional variables (pointers)
    var[:nw][end_t][:density] = var[:nw][start_t][:density]
    var[:nw][end_t][:pipe_flux_avg] = var[:nw][start_t][:pipe_flux_avg]
    var[:nw][end_t][:pipe_flux_neg] = var[:nw][start_t][:pipe_flux_neg]
    var[:nw][end_t][:pipe_flux_fr] = var[:nw][start_t][:pipe_flux_fr]
    var[:nw][end_t][:pipe_flux_to] = var[:nw][start_t][:pipe_flux_to]
    var[:nw][end_t][:pipe_flow_avg] = var[:nw][start_t][:pipe_flow_avg]
    var[:nw][end_t][:pipe_flow_neg] = var[:nw][start_t][:pipe_flow_neg]
    var[:nw][end_t][:pipe_flow_fr] = var[:nw][start_t][:pipe_flow_fr]
    var[:nw][end_t][:pipe_flow_to] = var[:nw][start_t][:pipe_flow_to]
    var[:nw][end_t][:compressor_flow] = var[:nw][start_t][:compressor_flow]
    var[:nw][end_t][:c_ratio] = var[:nw][start_t][:c_ratio]
    var[:nw][end_t][:fs] = var[:nw][start_t][:fs]
    var[:nw][end_t][:fd] = var[:nw][start_t][:fd]
    var[:nw][end_t][:ft] = var[:nw][start_t][:ft]
    var[:nw][end_t][:compressor_power] = var[:nw][start_t][:compressor_power]
    var[:nw][end_t][:net_nodal_injection] = var[:nw][start_t][:net_nodal_injection]
    var[:nw][end_t][:net_nodal_edge_out_flow] = var[:nw][start_t][:net_nodal_edge_out_flow]

    # density derivatives
    Threads.@threads for n in time_points[1:end-1]
        nw = ref[:nw][n]
        prev = n - 1
        (n == start_t) && (prev = time_points[end-1])
        var[:nw][n][:density_derivative] = Dict{Int,Any}()
        for (i, junction) in nw[:junction]
            if junction["junction_type"] == 1
                var[:nw][n][:density_derivative][i] = 0
            else
                var[:nw][n][:density_derivative][i] =
                    (var[:nw][n][:density][i] - var[:nw][prev][:density][i]) /
                    ref[:time_step]
            end
        end
    end
    
    # constraints 
    Threads.@threads for n in time_points[1:end-1]
        nw = ref[:nw][n]

        con[:nw][n][:slack_density] = Dict()
        for (i, junction) in nw[:slack_junctions]
            rho = var[:nw][n][:density][i]
            rho_fixed = junction["p_nominal"]
            con[:nw][n][:slack_density][i] = JuMP.@constraint(node[:nw][n], rho == junction["p_nominal"])
        end
        
        con[:nw][n][:nodal_balance] = Dict()
        for (i, junction) in nw[:junction]
            net_injection = var[:nw][n][:net_nodal_injection][i]
            net_nodal_edge_out_flow = var[:nw][n][:net_nodal_edge_out_flow][i]
            con[:nw][n][:nodal_balance][i] =
                JuMP.@constraint(node[:nw][n], net_injection == net_nodal_edge_out_flow)
        end
        
        con[:nw][n][:compressor_flow_dir] = Dict()
        con[:nw][n][:compressor_boost] = Dict()
        con[:nw][n][:compressor_power] = Dict()
        for (i, compressor) in nw[:compressor]
            p_fr = var[:nw][n][:density][compressor["fr_junction"]]
            p_to = var[:nw][n][:density][compressor["to_junction"]]
            a = var[:nw][n][:c_ratio][i]
            f = var[:nw][n][:compressor_flow][i]
            power = var[:nw][n][:compressor_power][i]
            power_max = compressor["power_max"]
            con[:nw][n][:compressor_boost][i] = JuMP.@constraint(node[:nw][n], p_to == a * p_fr)
            con[:nw][n][:compressor_flow_dir][i] = JuMP.@constraint(node[:nw][n], f * (p_fr - p_to) <= 0)
            con[:nw][n][:compressor_power][i] = Plasmo.@NLnodeconstraint(node[:nw][n], power <= power_max)
        end

        con[:nw][n][:pipe_momentum_balance] = Dict()
        for (i, pipe) in nw[:pipe]
            p_fr = var[:nw][n][:density][pipe["fr_junction"]]
            p_to = var[:nw][n][:density][pipe["to_junction"]]
            resistance = pipe["resistance"]
            f = var[:nw][n][:pipe_flux_avg][i]
            con[:nw][n][:pipe_momentum_balance][i] =
                Plasmo.@NLnodeconstraint(node[:nw][n], p_fr^2 - p_to^2 - resistance * f * abs(f) == 0)
        end
    end

    Threads.@threads for n in time_points[1:end-1]
        nw = ref[:nw][n]
        con[:nw][n][:pipe_mass_balance] = Dict()
        for (i, pipe) in nw[:pipe]
            p_fr_dot = var[:nw][n][:density_derivative][pipe["fr_junction"]]
            p_to_dot = var[:nw][n][:density_derivative][pipe["to_junction"]]
            L = pipe["length"]
            phi = var[:nw][n][:pipe_flux_neg][i]
            
            con[:nw][n][:pipe_mass_balance][i] =
                Plasmo.@linkconstraint(edge[:nw][n], L * (p_fr_dot + p_to_dot) + 4 * phi == 0, attach = node[:nw][n])
        end
    end

    econ_weight = ref[:economic_weighting]
    load_shed_expressions = [[] for n in time_points[1:end-1]]
    compressor_power_expressions = [[] for n in time_points[1:end-1]]
    Threads.@threads for n in time_points[1:end-1]
        for (i, receipt) in ref[:nw][n][:dispatchable_receipt]
            push!(
                load_shed_expressions[n],
                JuMP.@NLexpression(
                    node[:nw][n].model,
                    receipt["offer_price"] * var[:nw][n][:fs][i]
                )
            )
        end
        for (i, delivery) in ref[:nw][n][:dispatchable_delivery]
            push!(
                load_shed_expressions[n],
                JuMP.@NLexpression(
                    node[:nw][n].model,
                    -delivery["bid_price"] * var[:nw][n][:fd][i]
                )
            )
        end
        for (i, transfer) in ref[:nw][n][:dispatchable_transfer]
            push!(
                load_shed_expressions[n],
                JuMP.@NLexpression(
                    node[:nw][n].model,
                    transfer["offer_price"] * var[:nw][n][:fts][i] -
                    transfer["bid_price"] * var[:nw][n][:ftd][i]
                )
            )
        end
        for (i, compressor) in ref[:nw][n][:compressor]
            push!(compressor_power_expressions[n], var[:nw][n][:compressor_power][i])
        end
    end

    # objective
    Threads.@threads for n in time_points[1:end-1]
        if length(load_shed_expressions) != 0 && length(compressor_power_expressions) != 0
            Plasmo.@NLnodeobjective(
                node[:nw][n],
                Min,
                econ_weight *
                sum(ex for ex in load_shed_expressions[n]) +
                (1 - econ_weight) *
                sum(ex for ex in compressor_power_expressions[n])
            )
        end
    end
    
    return graph, var, con, sol
end


function ref_fix_transient!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if InfrastructureModels.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]

        for (i, pipe) in nw_ref[:pipe]
            pd_min = pipe["pd_min"]
            pd_max = pipe["pd_max"]
            lambda = pipe["friction_factor"]
            L = pipe["length"] * ref[:base_length]
            D = pipe["diameter"]
            pipe["resistance"] = lambda * L / D
            w = 1 / pipe["resistance"]
            min_flux = pd_min < 0 ? -sqrt(w * abs(pd_min)) : sqrt(w * abs(pd_min))
            max_flux = pd_max < 0 ? -sqrt(w * abs(pd_max)) : sqrt(w * abs(pd_max))
            pipe["flux_min"] = min_flux
            pipe["flux_max"] = max_flux
            pipe["flow_min"] = pipe["flux_min"] * pipe["area"]
            pipe["flow_max"] = pipe["flux_max"] * pipe["area"]
        end

        for (i, pipe) in nw_ref[:compressor]
            pd_min = pipe["pd_min"]
            pd_max = pipe["pd_max"]
            lambda = pipe["friction_factor"]
            L = pipe["length"] * ref[:base_length]
            D = pipe["diameter"]
            pipe["resistance"] = lambda * L / D
            w = 1 / pipe["resistance"]
            min_flux = pd_min < 0 ? -sqrt(w * abs(pd_min)) : sqrt(w * abs(pd_min))
            max_flux = pd_max < 0 ? -sqrt(w * abs(pd_max)) : sqrt(w * abs(pd_max))
            pipe["flux_min"] = min_flux
            pipe["flux_max"] = max_flux
            if pipe["directionality"] == 1
                pipe["flux_min"] = 0.0
                min_flux = 0.0
            end
            pipe["flow_min"] = pipe["area"] * min_flux
            pipe["flow_max"] = pipe["area"] * max_flux
        end
    end
end

