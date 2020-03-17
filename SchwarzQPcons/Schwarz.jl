using JuMP, LightGraphs, Random, LinearAlgebra, Printf, Metis, DelimitedFiles, PowerModels, Ipopt, SparseArrays

function get_data(path)
    data = PowerModels.parse_file(path)
    A = calc_susceptance_matrix(data)

    args = Dict()
    args[:N] = size(A.matrix,1)
    args[:y] = A.matrix; 
    args[:del] = spzeros(args[:N],args[:N])    
    args[:g] = Graph(args[:N])
    for e in values(data["branch"])
        i = A.bus_to_idx[e["f_bus"]]
        j = A.bus_to_idx[e["t_bus"]]
        add_edge!(args[:g],i,j)
        args[:del][i,j] = e["angmin"]
        args[:del][j,i] = e["angmin"]
    end

    args[:val] =-ones(args[:N])*Inf
    args[:vau] = ones(args[:N])*Inf
    args[:ref] =  A.bus_to_idx[reference_bus(data)["index"]]
    args[:val][args[:ref]] = 0; args[:vau][args[:ref]] = 0
    
    gens = bus_gen_lookup(data["gen"],data["bus"])
    args[:c1]  = [[1e6,-1e6] for i=1:args[:N]]
    args[:c2]  = [[.0,.0] for i=1:args[:N]]
    args[:sl]= [[.0,-1e2] for i=1:args[:N]]
    args[:su]= [[1e2,.0] for i=1:args[:N]]
    args[:ng] = 2*ones(Int64,args[:N])


    for i=1:args[:N]
        bus = A.idx_to_bus[i]
        for j=1:length(gens[bus])
            if length(gens[bus][j]["cost"])!=0
                push!(args[:c1][i],gens[bus][j]["cost"][1])
                push!(args[:c2][i],gens[bus][j]["cost"][2])
                push!(args[:sl][i],gens[bus][j]["pmin"])
                push!(args[:su][i],gens[bus][j]["pmax"])
                args[:ng][i] += 1
            end
        end
    end

    args[:Ng]=sum(args[:ng])

    args[:sd] = zeros(args[:N])
    for v in values(data["load"])
        args[:sd][A.bus_to_idx[v["load_bus"]]] = v["pd"]
    end

    return args
end


function closed_neighbors(g,V,d)
    V = collect(V)
    nbr = copy(V)
    newnbr = copy(V)
    oldnbr = []
    for k=1:d
        for i in newnbr
            union!(nbr, neighbors(g,i))
        end
        union!(oldnbr,newnbr)
        newnbr = setdiff(nbr,oldnbr)
    end
    return sort(nbr)
end

function tic()
    global tic_toc_time=-time()
    return 
end

function toc()
    tt = time()+tic_toc_time
    println("Elapased time = $tt sec")
    return tt
end

function dcopf!(m,V_om,V_om_compl,va,s,pg,sep,pg_exp,sep_exp,obj_exp,args;gamma=0)
    for i in V_om
        va[i] = @variable(m,upper_bound = args[:vau][i], lower_bound = args[:val][i])
        s[i] = @variable(m, [j=1:args[:ng][i]],upper_bound = args[:su][i][j], lower_bound = args[:sl][i][j])
    end
    for i in V_om_compl
        va[i] = @variable(m); fix(va[i],0)
        pg[i] = @variable(m); fix(pg[i],0)
        sep[i] = [@variable(m) for j in neighbors(args[:g],i)]; fix.(sep[i],0)
    end
    for i in [V_om;V_om_compl]
        pg_exp[i] = sum(s[i][j] for j=1:args[:ng][i]) - sum(args[:y][i,j]*(va[i]-va[j]) for j in neighbors(args[:g],i)) - args[:sd][i]
        sep_exp[i]= [va[i] - va[j] - args[:del][i,j] for j in neighbors(args[:g],i)]
        obj_exp[i]= sum(args[:c1][i][j]*s[i][j]+args[:c2][i][j]*s[i][j]^2 for j=1:args[:ng][i]) +
            1/4*gamma * sum((va[i]-va[j])^2 for j in neighbors(args[:g],i)) #+ 1/2*args[:rho] * va[i]^2
    end
    for i in V_om
        pg[i] = @constraint(m,pg_exp[i]==0)
    end
    for i in V_om
        sep[i]= @constraint(m,sep_exp[i].>=0)
    end
    @objective(m, Min, sum(obj_exp[i] for i in [V_om;V_om_compl])
               -sum(sep[i]'*sep_exp[i]+pg[i]*pg_exp[i] for i in V_om_compl))
    return m
end


function Schwarz_init(V0,V_om0,V_compl0,V_om_compl0,V_give0,gamma,args0)
    global args
    global V,V_om,V_compl,V_om_compl,V_give
    global m,va,s,pg,sep,pg_exp,sep_exp,obj_exp
    global va_val,s_val,pg_val,sep_val

    V=V0; V_om=V_om0; V_compl=V_compl0; V_om_compl=V_om_compl0; V_give=V_give0
    args=args0
    m = Model(with_optimizer(Ipopt.Optimizer,print_level=0,linear_solver="ma57",tol=1e-12))

    va = Array{Any,1}(zeros(args[:N]))
    s = Array{Any,1}([zeros(args[:ng][i]) for i=1:args[:N]])
    pg = Array{Any,1}(zeros(args[:N]))
    sep = Array{Any,1}(zeros(args[:N]))
    pg_exp = Array{Any,1}(zeros(args[:N]))
    sep_exp = Array{Any,1}(zeros(args[:N]))
    obj_exp = Array{Any,1}(zeros(args[:N]))

    va_val = Array{Any,1}(zeros(args[:N]))
    s_val = [zeros(args[:ng][i]) for i=1:args[:N]]
    pg_val = Array{Any,1}(zeros(args[:N]))
    sep_val = [zeros(length(neighbors(args[:g],i))) for i=1:args[:N]]
    
    dcopf!(m,V_om,V_om_compl,va,s,pg,sep,pg_exp,sep_exp,obj_exp,args,gamma=gamma)

    return true
end


function Schwarz_iter_update(va_compl,pg_compl,sep_compl)

    obj = sum(value.(obj_exp[V]))
    err_pr = maximum([0;va_val[V_compl].-va_compl]);
    err_du = max(maximum([0;pg_val[V_compl].-pg_compl]),0)
    slk = maximum(vcat([s_val[i][1] for i in V]...)) - minimum(sum(vcat([s_val[i][2] for i in V]...)))

    return obj,err_pr,err_du,slk
end

function Schwarz_iter_solve(va_om_compl,pg_om_compl,sep_om_compl,om)
    fix.(va[V_om_compl],va_om_compl)
    fix.(pg[V_om_compl],pg_om_compl)
    fix.(vcat(sep[V_om_compl]...),vcat(sep_om_compl...))

    optimize!(m)

    
    va_val[V] = value.(va[V])
    s_val[V]  = [value.(s[i]) for i in V]
    pg_val[V] = dual.(pg[V])
    sep_val[V]= [dual.(sep[i]) for i in V]

    if om!=0
        va_val[V_compl] = value.(va[V_compl])
        pg_val[V_compl] = dual.(pg[V_compl])
        sep_val[V_compl]= [dual.(sep[i]) for i in V_compl]
    else
        va_val[V_compl] = value.(va[V_compl])
        pg_val[V_compl] = value.(pg[V_compl])
        sep_val[V_compl]= [value.(sep[i]) for i in V_compl]
    end

    set_start_value.(all_variables(m), value.(all_variables(m)))
    
    return va_val[V_give],pg_val[V_give],sep_val[V_give]
end

Schwarz_get_value() = s_val[V]
cost(s,V,args) = sum(args[:c1][i][j]*s[i][j]+args[:c2][i][j]*s[i][j]^2 for i in V for j=1:args[:ng][i])
