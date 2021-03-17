using Random,JuMP,PowerModels,Ipopt,DelimitedFiles,Plots,LightGraphs,LinearAlgebra
const PM = PowerModels
import PowerModels: build_opf, constraint_power_balance, constraint_voltage_angle_difference,
    constraint_ohms_yt_from, constraint_thermal_limit_from, constraint_ohms_yt_to, constraint_thermal_limit_to
import JuMP: dual

dual(x::Nothing) = 0

function build_opf(pm::AbstractPowerModel)
    variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)

    objective_min_fuel_and_flow_cost(pm)

    con(pm)[:va_ref]=Dict()
    con(pm)[:va_diff_min]=Dict()
    con(pm)[:va_diff_max]=Dict()
    con(pm)[:bal_p]=Dict()
    con(pm)[:bal_q]=Dict()
    con(pm)[:ohms_from_p]=Dict()
    con(pm)[:ohms_from_q]=Dict()
    con(pm)[:ohms_to_p]=Dict()
    con(pm)[:ohms_to_q]=Dict()
    con(pm)[:therm_from]=Dict()
    con(pm)[:therm_to]=Dict()
    
    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        con(pm,:va_ref)[i] = constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        con(pm,:bal_p)[i], con(pm,:bal_q)[i] = constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end
end

function constraint_thermal_limit_from(pm::AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    con(pm,:therm_from)[f_idx[1]] = JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2)
end

function constraint_thermal_limit_to(pm::AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    con(pm,:therm_to)[t_idx[1]] = JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2)
end

function constraint_ne_ohms_yt_from(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = var(pm, n, :p_ne, f_idx)
    q_fr  = var(pm, n, :q_ne, f_idx)
    vm_fr = var(pm, n,   :vm, f_bus)
    vm_to = var(pm, n,   :vm, t_bus)
    va_fr = var(pm, n,   :va, f_bus)
    va_to = var(pm, n,   :va, t_bus)
    z = var(pm, n, :branch_ne, i)
    
    con(pm,:ohms_from_p)[f_idx[1]] = JuMP.@NLconstraint(pm.model, p_fr == z*( (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
    con(pm,:ohms_from_q)[f_idx[1]] = JuMP.@NLconstraint(pm.model, q_fr == z*(-(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
end

function constraint_ne_ohms_yt_to(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_to = var(pm, n, :p_ne, t_idx)
    q_to = var(pm, n, :q_ne, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :branch_ne, i)

    con(pm,:ohms_to_p)[f_idx[1]] = JuMP.@NLconstraint(pm.model, p_to == z*( (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
    con(pm,:ohms_to_p)[f_idx[1]] = JuMP.@NLconstraint(pm.model, q_to == z*(-(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
end


function constraint_voltage_angle_difference(pm::AbstractPolarModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    con(pm,:va_diff_max)[f_idx[1]] = JuMP.@constraint(pm.model, va_fr - va_to <= angmax)
    con(pm,:va_diff_max)[f_idx[1]] = JuMP.@constraint(pm.model, va_fr - va_to >= angmin)
end


function constraint_power_balance(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = PM.var(pm, n, :vm, i)
    p    = PM.get(PM.var(pm, n),    :p, Dict()); PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = PM.get(PM.var(pm, n),    :q, Dict()); PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = PM.get(PM.var(pm, n),   :pg, Dict()); PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = PM.get(PM.var(pm, n),   :qg, Dict()); PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = PM.get(PM.var(pm, n),   :ps, Dict()); PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = PM.get(PM.var(pm, n),   :qs, Dict()); PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = PM.get(PM.var(pm, n),  :psw, Dict()); PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = PM.get(PM.var(pm, n),  :qsw, Dict()); PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = PM.get(PM.var(pm, n), :p_dc, Dict()); PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = PM.get(PM.var(pm, n), :q_dc, Dict()); PM._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    if !nl_form
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    else
        cstr_p = JuMP.@NLconstraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    end

    if !nl_form
        cstr_q = JuMP.@constraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    else
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    end

    return cstr_p,cstr_q
end

getsol(pm)=[vcat([value(var(pm,:va,i)),
                  value(var(pm,:vm,i)),
                  [value(var(pm,:pg,j)) for j in ref(pm,:bus_gens)[i]],
                  [value(var(pm,:qg,j)) for j in ref(pm,:bus_gens)[i]],
                  [value(var(pm,:p,a)) for a in ref(pm,:bus_arcs)[i]],
                  [value(var(pm,:q,a)) for a in ref(pm,:bus_arcs)[i]],
                  dual(UpperBoundRef(var(pm,:vm,i))),
                  dual(LowerBoundRef(var(pm,:vm,i))),
                  [dual(UpperBoundRef(var(pm,:pg,j))) for j in ref(pm,:bus_gens)[i]],
                  [dual(UpperBoundRef(var(pm,:qg,j))) for j in ref(pm,:bus_gens)[i]],
                  [dual(UpperBoundRef(var(pm,:p,a))) for a in ref(pm,:bus_arcs)[i]],
                  [dual(UpperBoundRef(var(pm,:q,a))) for a in ref(pm,:bus_arcs)[i]],
                  [dual(LowerBoundRef(var(pm,:pg,j))) for j in ref(pm,:bus_gens)[i]],
                  [dual(LowerBoundRef(var(pm,:qg,j))) for j in ref(pm,:bus_gens)[i]],
                  [dual(LowerBoundRef(var(pm,:p,a))) for a in ref(pm,:bus_arcs)[i]],
                  [dual(LowerBoundRef(var(pm,:q,a))) for a in ref(pm,:bus_arcs)[i]],
                  [dual(con(pm,:va_ref,j)) for j in keys(ref(pm,:ref_buses)) if i==j],
                  dual(con(pm,:bal_p,i)),
                  dual(con(pm,:bal_q,i)),
                  [dual(con(pm,:va_diff_max,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:va_diff_min,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:ohms_to_p,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:ohms_to_q,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:ohms_from_p,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:ohms_from_q,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:therm_from,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i],
                  [dual(con(pm,:therm_to,a[1])) for a in ref(pm,:bus_arcs)[i] if a[2]<i]
                  ]...)
            for i in sort(parse.(Int64,keys(pm.data["bus"])))]

Random.seed!(1)
fcolor = [0,0,1]
bcolor = [.95,.95,.95]
function coloring(x)
    x=min(1,x)
    arr = (1-x)*bcolor+x*fcolor
    return RGB(arr[1],arr[2],arr[3])
end

sig = 1e-3
n_sample = 30
n_grad = 4
casename = "case500_tamu"

dat_path = "data/pglib_opf_"*casename*".m"
pos_path = "data/pos_"*casename*".csv"

j = 1
params = [(1e6,10),(0,10),(1e6,0),(0,0)]

for cnt=1:length(params)
    (eta,b) = params[cnt]

    data = PowerModels.parse_file(dat_path)
    pos = readdlm(pos_path,',')
    buslist = sort(parse.(Int64,keys(data["bus"])))
    inv_buslist = Dict(buslist[i]=>i for i=1:length(buslist))

    g = Graph(length(data["bus"]))
    for br in values(data["branch"])
        add_edge!(g,br["t_bus"],br["f_bus"])
    end


    data["gen"]["0"]=Dict(
        "ncost"      => 0,
        "qc1max"     => 0.0,
        "pg"         => 0.0,
        "model"      => 2,
        "shutdown"   => 0.0,
        "startup"    => 0.0,
        "qc2max"     => 0.0,
        "ramp_agc"   => 0.0,
        "qg"         => 0.0,
        "gen_bus"    => j,
        "pmax"       => 0.0,
        "ramp_10"    => 0.0,
        "vg"         => 1.005,
        "mbase"      => 100.0,
        "source_id"  => Any["gen", 0],
        "pc2"        => 0.0,
        "index"      => 0,
        "cost"       => Float64[],
        "qmax"       => 0.0,
        "gen_status" => 1,
        "qmin"       => 0.0,
        "qc1min"     => 0.0,
        "qc2min"     => 0.0,
        "pc1"        => 0.0,
        "ramp_q"     => 0.0,
        "ramp_30"    => 0.0,
        "pmin"       => 0.0,
        "apf"        => 0.0)

    for gen in values(data["gen"])
        gen["pmax"]+=b
        gen["qmax"]+=b
        gen["pmin"]-=b
        gen["pmin"]-=b
    end

    pm = PowerModels.instantiate_model(data, ACPPowerModel, PowerModels.build_opf)
    set_optimizer(pm.model,Ipopt.Optimizer)
    set_optimizer_attribute(pm.model,"linear_solver","ma57")
    set_optimizer_attribute(pm.model,"tol",1e-10)
    set_optimizer_attribute(pm.model,"print_level",0)

    reg = sum((pm.var[:nw][0][:vm].-1).^2) +
        sum((pm.var[:nw][0][:va][src(e)] - pm.var[:nw][0][:va][dst(e)])^2
            for e in edges(g))
    @objective(pm.model,Min,objective_function(pm.model) + eta * reg)


    result=optimize_model!(pm)
    set_start_value.(all_variables(pm.model),value.(all_variables(pm.model)))
    sol_ref = getsol(pm)

    sols = []
    ds = []
    for i=1:n_sample
        println("$i")
        dp = sig*(rand(2).-.5)
        push!(ds,dp)

        set_upper_bound(pm.var[:nw][0][:pg][0],dp[1])
        set_upper_bound(pm.var[:nw][0][:qg][0],dp[2])
        set_lower_bound(pm.var[:nw][0][:pg][0],dp[1])
        set_lower_bound(pm.var[:nw][0][:qg][0],dp[2])
        
        result = optimize_model!(pm)
        
        set_start_value.(all_variables(pm.model),value.(all_variables(pm.model)))
        push!(sols,getsol(pm))
    end

    i=1
    C= maximum(hcat([norm.(sols[i].-sol_ref) for i=1:n_sample]...),dims=2)
    mc= maximum(C)
    if C[j]<=1e-10
        error("too small")
    end
    C = C/C[j]
    pyplot()
    p=scatter(leg=false,ticks=nothing,border=:none,size=(500,500));
    for br in values(data["branch"])
        plot!(p,
              range(pos[inv_buslist[br["t_bus"]],2];
                    stop = pos[inv_buslist[br["f_bus"]],2],
                    length=n_grad),
              range(pos[inv_buslist[br["t_bus"]],3];
                    stop = pos[inv_buslist[br["f_bus"]],3],
                    length=n_grad),
              lc= cgrad([coloring(C[inv_buslist[br["t_bus"]]]),
                         coloring(C[inv_buslist[br["f_bus"]]])]),
              line_z =1:n_grad)
    end
    for i=1:size(pos,1)
        scatter!(p,[pos[i,2]],[pos[i,3]],markersize=8,markercolor=coloring(C[i]),
                 markerstrokewidth=0);
    end
    scatter!(p,[pos[inv_buslist[j],2]],[pos[inv_buslist[j],3]],markersize=18,markercolor=false,markerstrokewidth=2,markerstrokecolor=:red);
    savefig(p, "fig/opf-hm-$cnt.pdf")
    pgfplotsx()
    q=scatter([length(a_star(g,i,j)) for i=1:nv(g)][C[:].>1e-10],C[:][C[:].>1e-10],
              size=(250,250), leg=false, framestyle=:box,
              yscale=:log10,
              ylims=(1e-5,10),
              markerstrokewidth=0,color=:black)
    savefig(q, "fig/opf-sc-$cnt.pdf")
end
