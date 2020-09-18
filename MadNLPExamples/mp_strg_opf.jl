# Created by Sungho Shin (sungho.shin@wisc.edu) and Carleton Coffrin (cjc@lanl.gov)

import PowerModels: constraint_storage_losses

function mp_strg_opf_14_dc(;Nt=1,Mt=4,sig=0,seed=0,plasmo=false)
    Random.seed!(seed)
    case = PowerModels.parse_file(
        joinpath(@__DIR__,"data/pglib_opf_case14_ieee_mod.m"))
    profile = kron(ones(Nt),get_profile(0:1/Mt:24-1/Mt))
    mn_case = PowerModels.replicate(case, length(profile))
    make_periodic_add_noise!(mn_case,profile,sig)
    pm = PowerModels.instantiate_model(
        mn_case,DCPPowerModel,plasmo ? build_mn_opf_strg_nl_plasmo : build_mn_opf_strg_nl)
    pm.model.ext[:pm]=pm
    return pm.model
end

function mp_strg_opf_14_ac(;Nt=1,Mt=4,sig=0,seed=0,plasmo=false)
    Random.seed!(seed)
    case = PowerModels.parse_file(
        joinpath(@__DIR__,"data/pglib_opf_case14_ieee_mod.m"))
    profile = kron(ones(Nt),get_profile(0:1/Mt:24-1/Mt))
    mn_case = PowerModels.replicate(case, length(profile))
    
    make_periodic_add_noise!(mn_case,profile,sig)
    pm = PowerModels.instantiate_model(
        mn_case,ACPPowerModel,plasmo ? build_mn_opf_strg_nl_plasmo : build_mn_opf_strg_nl)
    pm.model.ext[:pm]=pm
    return pm.model
end

function mp_strg_opf_14_soc(;Nt=1,Mt=4,sig=0,seed=0,plasmo=false)
    Random.seed!(seed)
    case = PowerModels.parse_file(
        joinpath(@__DIR__,"data/pglib_opf_case14_ieee_mod.m"))
    profile = kron(ones(Nt),get_profile(0:1/Mt:24-1/Mt))
    mn_case = PowerModels.replicate(case, length(profile))
    make_periodic_add_noise!(mn_case,profile,sig)
    pm = PowerModels.instantiate_model(
        mn_case,SOCWRPowerModel, plasmo ? build_mn_opf_strg_nl_plasmo : build_mn_opf_strg_nl)
    pm.model.ext[:pm]=pm
    return pm.model
end

function make_periodic_add_noise!(mn_case,profile,sig)
    for (i,scalar) in enumerate(profile)
        network = mn_case["nw"]["$(i)"]
        network["scalar"] = scalar*(1+sig*randn())
        for (i,load) in network["load"]
            load["pd"] = network["scalar"]*load["pd"]
            load["qd"] = network["scalar"]*load["qd"]
        end
    end
end

function slice(mn_data, lookup)
    @assert PowerModels.InfrastructureModels.ismultinetwork(mn_data)
    vals = []
    nw_keys = sort(collect(keys(mn_data["nw"])), by=x->parse(Int, x))
    for n in nw_keys
        comp = mn_data["nw"][n]
        for lookup_val in lookup
            comp = comp[lookup_val]
        end
        push!(vals, comp)
    end
    return vals
end

# from RTS 96 paper
summer_wkdy_hour_scalar = [
    .64, .60, .58, .56, .56, .58, .64, .76, .87, .95, .99, 1.0, .99, 1.0, 1.0,
    .97, .96, .96, .93, .92, .92, .93, .87, .72, .64]
get_profile=scale(interpolate(summer_wkdy_hour_scalar,BSpline(Linear())),0:24)

summer_wkdy_15min_scalar = Float64[]

function build_mn_opf_strg_nl(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n)
        variable_storage_power(pm, nw=n)
        variable_branch_power(pm, nw=n)
        variable_dcline_power(pm, nw=n)

        constraint_model_voltage(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, :storage, nw=n)
            constraint_storage_losses(pm, i, nw=n)
            constraint_storage_thermal_limit(pm, i, nw=n)
        end

        for i in ids(pm, :branch, nw=n)
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)

            constraint_voltage_angle_difference(pm, i, nw=n)

            constraint_thermal_limit_from(pm, i, nw=n)
            constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in ids(pm, :dcline, nw=n)
            constraint_dcline(pm, i, nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    objective_min_fuel_and_flow_cost(pm)
end

function build_mn_opf_strg_nl_plasmo(pm::AbstractPowerModel)
    graph = OptiGraph()
    pms  = Dict(n=>shallow_copy(pm) for (n, network) in nws(pm))
    node = Dict(n=>@optinode(graph) for (n, network) in nws(pm))

    Threads.@threads for (n, network) in collect(nws(pm))
        pms[n].model = node[n].model
        variable_bus_voltage(pms[n], nw=n)
        variable_gen_power(pms[n], nw=n)
        variable_storage_power(pms[n], nw=n)
        variable_branch_power(pms[n], nw=n)
        variable_dcline_power(pms[n], nw=n)

        constraint_model_voltage(pms[n], nw=n)

        for i in ids(pms[n], :ref_buses, nw=n)
            constraint_theta_ref(pms[n], i, nw=n)
        end

        for i in ids(pms[n], :bus, nw=n)
            constraint_power_balance(pms[n], i, nw=n)
        end

        for i in ids(pms[n], :storage, nw=n)
            # constraint_storage_complementarity_nl(pms[n], i, nw=n)
            constraint_storage_losses(pms[n], i, nw=n)
            constraint_storage_thermal_limit(pms[n], i, nw=n)
        end

        for i in ids(pms[n], :branch, nw=n)
            constraint_ohms_yt_from(pms[n], i, nw=n)
            constraint_ohms_yt_to(pms[n], i, nw=n)

            constraint_voltage_angle_difference(pms[n], i, nw=n)

            constraint_thermal_limit_from(pms[n], i, nw=n)
            constraint_thermal_limit_to(pms[n], i, nw=n)
        end
        for i in ids(pm, :dcline, nw=n)
            constraint_dcline(pm, i, nw=n)
        end

        JuMP.@objective(node[n], Min,
                        sum(gen_cost(gen["cost"],var(pms[n],n,:pg,i),conductor_ids(pms[n], n))
                            for (i,gen) in network[:gen]))
    end

    n_1 = 1
    for (i,storage) in ref(pm, :storage, nw=n_1)
        JuMP.@constraint(
            node[n_1], var(pm, n_1, :se, i) - storage["energy"] ==
            ref(pm, n_1, :time_elapsed)*(storage["charge_efficiency"]*var(pm, n_1, :sc, i) -
                                         var(pm, n_1, :sd, i)/storage["discharge_efficiency"]))
    end

    for n_2 in 2:length(nws(pm))
        for (i,storage) in ref(pm, :storage, nw=n_2)
            Plasmo.@linkconstraint(
                graph, var(pm,n_2,:se,i) - var(pm,n_1,:se,i) ==
                ref(pm, n_1, :time_elapsed)*(storage["charge_efficiency"]*var(pm,n_2,:sc,i)-
                                             var(pm,n_2,:sd,i)/storage["discharge_efficiency"]),
                attach=node[n_1])
        end
        n_1 = n_2
    end

    pm.model = graph 
end

shallow_copy(pm) = typeof(pm)([getfield(pm,f) for f in fieldnames(typeof(pm))]...)
function gen_cost(cost,pgs,conductor_ids)
    pg = sum(pgs[c] for c in conductor_ids)
    if length(cost) == 1
        return cost[1]
    elseif length(cost) == 2
        return cost[1]*pg + cost[2]
    elseif length(cost) == 3
        return cost[1]*pg^2 + cost[2]*pg + cost[3]
    else
        return 0.0
    end
end
    

function constraint_storage_losses(pm::AbstractACPModel, n::Int, i, bus, r, x, p_loss, q_loss; conductors=[1])
    vm = var(pm, n, :vm, bus)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@NLconstraint(pm.model,sum(ps[c] for c in conductors) + (sd - sc)== 0)
    JuMP.@NLconstraint(pm.model,sum(qs[c] for c in conductors) == qsc)
end


