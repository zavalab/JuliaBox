include("multigrid.jl")

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
n_proc = MPI.Comm_size(comm)

# warm-up                         
m=Model(solver=IpoptSolver(print_level=0))
@variable(m,x[1:2])
@objective(m,Min,x[1]^2+x[2]^2)
@NLconstraint(m,[1],x[1]^2/x[2]==1)
solve(m)

par = Dict()
if rank == 0
    results = []
end

for cs = 1:7
    result = Dict()
    if rank == 0
        if cs == 1
            data = load("jlds/pglib_opf_case500_tamu.jld")            
            par[:rho] = 1
            par[:scale]  = 1e6
            par[:divider]= 10
        elseif cs == 2
            data = load("jlds/pglib_opf_case2853_sdet.jld") #
            par[:rho] = 1
            par[:scale]  = 1e6
            par[:divider]= 10
        elseif cs == 3
            data = load("jlds/pglib_opf_case3120sp_k.jld")
            par[:rho] = 1
            par[:scale]  = 1e6
            par[:divider]= 20
        elseif cs == 4
            data = load("jlds/pglib_opf_case4661_sdet.jld") #
            par[:rho] = 1
            par[:scale]  = 2e6
            par[:divider]= 20
        elseif cs == 5
            data = load("jlds/pglib_opf_case6468_rte.jld") #
            par[:rho] = 1
            par[:scale]  = 2e6
            par[:divider]= 20
        elseif cs == 6
            data = load("jlds/pglib_opf_case6515_rte.jld")
            par[:rho] = 1
            par[:scale]  = 8e6
            par[:divider]= 20
        elseif cs == 7
            data = load("jlds/pglib_opf_case9241_pegase.jld")
            par[:rho] = 1
            par[:scale]  = 64e6
            par[:divider]= 20
        end          
        par[:n_part] = n_proc-1
        par[:price]  = 1e4
        par[:maxiter]= 1000
        
        net = get_net(data)

        W = get_partition(net,par[:n_part])
        W1= [closed_neighbors(net[:g],W[k],1) for k=1:par[:n_part]]
        Wz= sort(union([setdiff(W1[k],W[k]) for k=1:par[:n_part]]...))
        zmap=Dict(i=>[] for i in Wz)
        for i in Wz
            for k=1:par[:n_part]
                if i in W1[k]
                    push!(zmap[i],k)
                end
            end
        end
        Wy= [intersect(W1[k],Wz) for k=1:par[:n_part]]

        G = [vcat(net[:genmap][W[k]]...) for k=1:par[:n_part]]
        E = [[] for k=1:par[:n_part]];
        for e in edges(net[:g])
            for k=1:par[:n_part]
                if (src(e) in W1[k]) && (dst(e) in W1[k])
                    push!(E[k],e)
                end
            end
        end
        
        n=sum(length(W1[k]) for k=1:par[:n_part])
        p=sum(length(Wy[k]) for k=1:par[:n_part])

        submap,subimap,subnet = get_subpartition(net,W,Wz,par[:divider])
        subW = [uniq(subimap[W[k]]) for k=1:par[:n_part]]
        subW1 = [uniq(subimap[W1[k]]) for k=1:par[:n_part]]
        subWy = [uniq(subimap[Wy[k]]) for k=1:par[:n_part]]
        subWz = uniq(subimap[Wz])

        subE =[[] for k=1:par[:n_part]];
        for e in edges(subnet[:g])
            for k=1:par[:n_part]
                if (src(e) in subW1[k]) && (dst(e) in subW1[k])
                    push!(subE[k],e)
                end
            end
        end

        # Centralized
        m_cen = ac_opf(net,W,W1,G,E,par)
        @variable(m_cen,V_z[i in Wz], start = 1)
        @variable(m_cen,T_z[i in Wz], start = 0)
        @constraint(m_cen, V_cc[k in 1:par[:n_part],i in Wy[k]], V_z[i]-m_cen[:V][k,i]==0)
        @constraint(m_cen, T_cc[k in 1:par[:n_part],i in Wy[k]], T_z[i]-m_cen[:T][k,i]==0)
        m_cen.solver= IpoptSolver(linear_solver="ma57",print_level=0)
        
        time_cen=time_ns()
        solve(m_cen)
        time_cen=(time_ns()-time_cen)/1e9

        obj_cen = getobjectivevalue(m_cen) * par[:scale]
        @printf("Central time:      %11.6e\n",time_cen)
        @printf("Central objective: %11.6e\n",obj_cen)

        V_z_cen = Dict(i=>getvalue(m_cen[:V_z][i]) for i in Wz)
        T_z_cen = Dict(i=>getvalue(m_cen[:T_z][i]) for i in Wz)
        V_y_cen = [Dict(i=>getdual(m_cen[:V_cc][k,i]) for i in Wy[k]) for k = 1:par[:n_part]]
        T_y_cen = [Dict(i=>getdual(m_cen[:T_cc][k,i]) for i in Wy[k]) for k = 1:par[:n_part]]

        # Coarse
        m_crs = ac_opf(subnet,subW,subW1,G,subE,par)
        m_crs.solver= IpoptSolver(linear_solver="ma57",print_level=0)
        @variable(m_crs,V_z[i in subWz], start = 1)
        @variable(m_crs,T_z[i in subWz], start = 0)
        @constraint(m_crs, V_cc[k in 1:par[:n_part],i in subWy[k]], V_z[i]-m_crs[:V][k,i]==0)
        @constraint(m_crs, T_cc[k in 1:par[:n_part],i in subWy[k]], T_z[i]-m_crs[:T][k,i]==0)
        
        time_crs=time_ns()
        solve(m_crs)
        time_crs=(time_ns()-time_crs)/1e9
        obj_crs = getobjectivevalue(m_crs) * par[:scale]
        @printf("Coarse time:      %11.6e\n",time_crs)
        @printf("Coarse objective: %11.6e\n",obj_crs)


        starter = Dict(:V => zeros(net[:n]),
                       :T => zeros(net[:n]),
                       :P => zeros(net[:ngen]),
                       :Q => zeros(net[:ngen]),
                       :Ps => zeros(net[:n]),
                       :Qs => zeros(net[:n]),
                       :Tdu=> [Dict(i=>.0 for i in Wy[k]) for k=1:par[:n_part]],
                       :Vdu=> [Dict(i=>.0 for i in Wy[k]) for k=1:par[:n_part]])

        for k=1:par[:n_part]
            for i in W[k]
                starter[:V][i] = getvalue(m_crs[:V][k,subimap[i]])
                starter[:T][i] = getvalue(m_crs[:T][k,subimap[i]])
                starter[:Ps][i] = getvalue(m_crs[:Ps][k,subimap[i]])
                starter[:Qs][i] = getvalue(m_crs[:Qs][k,subimap[i]])
            end
            for i in G[k]
                starter[:P][i] = getvalue(m_crs[:P][k,i])
                starter[:Q][i] = getvalue(m_crs[:Q][k,i])
            end
            for i in Wy[k]
                starter[:Vdu][k][i] = getdual(m_crs[:V_cc][k,subimap[i]])
                starter[:Tdu][k][i] = getdual(m_crs[:T_cc][k,subimap[i]])
            end
        end

        # Warm
        m_warm = ac_opf(net,W,W1,G,E,par)
        @variable(m_warm,V_z[i in Wz], start = starter[:V][i])
        @variable(m_warm,T_z[i in Wz], start = starter[:T][i])
        @constraint(m_warm, V_cc[k in 1:par[:n_part],i in Wy[k]], V_z[i]-m_warm[:V][k,i]==0)
        @constraint(m_warm, T_cc[k in 1:par[:n_part],i in Wy[k]], T_z[i]-m_warm[:T][k,i]==0)

        m_warm = ac_opf(net,[1:net[:n]],[1:net[:n]],[1:net[:ngen]],[edges(net[:g])],par)
        m_warm.solver= IpoptSolver(linear_solver="ma57",print_level=0)
        
        time_warm=time_ns()
        solve(m_warm)
        time_warm=(time_ns()-time_warm)/1e9

        obj_warm = getobjectivevalue(m_warm) * par[:scale]
        @printf("Warm time:      %11.6e\n",time_warm)
        @printf("Warm objective: %11.6e\n",obj_warm)

        # Print graphs
        nlist, elist, subnlist, subelist = printgraph(par, net, subnet, W, subW, submap, subimap)
        writecsv("graph/nlist-$cs.csv",nlist)
        writecsv("graph/elist-$cs.csv",elist)
        writecsv("graph/subnlist-$cs.csv",subnlist)
        writecsv("graph/subelist-$cs.csv",subelist)
    else
        par = nothing
        net = nothing
        subnet = nothing
        W = nothing
        W1 = nothing
        Wy = nothing
        Wz = nothing
        E = nothing

        subW = nothing
        subW1 = nothing
        subWy = nothing
        subWz = nothing
        subE = nothing

        G = nothing
starter = nothing
end

net = MPI.bcast(net, 0, comm)
subnet = MPI.bcast(subnet, 0, comm)
par = MPI.bcast(par, 0, comm)

W = MPI.bcast(W, 0, comm)
W1 = MPI.bcast(W1, 0, comm)
Wy = MPI.bcast(Wy, 0, comm)
Wz = MPI.bcast(Wz, 0, comm)
E = MPI.bcast(E, 0, comm)

subW = MPI.bcast(subW, 0, comm)
subW1 = MPI.bcast(subW1, 0, comm)
subWy = MPI.bcast(subWy, 0, comm)
subWz = MPI.bcast(subWz, 0, comm)
subE = MPI.bcast(subE, 0, comm)

G = MPI.bcast(G, 0, comm)
starter = MPI.bcast(starter, 0, comm)

for mode = [:mul, :dec]
    if rank == 0
        tol_abs = 1e-3/2
        tol_rel = 1e-3/2

        tol_pr = 0
        tol_du = 0

        err_r = Inf
        err_s = Inf

        T_r =[Dict(i => .0 for i in Wy[k]) for k=1:par[:n_part]]
        V_r =[Dict(i => .0 for i in Wy[k]) for k=1:par[:n_part]]
        T_s = Dict(i => .0 for i in Wz)
        V_s = Dict(i => .0 for i in Wz)

        iter = 0
        t0=time_ns()
    else
        m =ac_opf(net,Dict(rank=>W[rank]),Dict(rank=>W1[rank]),Dict(rank=>G[rank]),Dict(rank=>E[rank]),par)
    end


    if mode == :dec
        rank == 0 && println("decentralized")
        T_x =[Dict(i => .0 for i in W1[k]) for k=1:par[:n_part]]
        V_x =[Dict(i =>1.0 for i in W1[k]) for k=1:par[:n_part]]

        T_y =[Dict(i => .0 for i in Wy[k]) for k=1:par[:n_part]]
        V_y =[Dict(i => .0 for i in Wy[k]) for k=1:par[:n_part]]

        T_z = Dict(i => .0 for i in Wz)
        V_z = Dict(i => .0 for i in Wz)
        
    elseif mode == :mul
        rank == 0 && println("multigrid")
        T_x =[Dict(i =>starter[:T][i] for i in W1[k]) for k=1:par[:n_part]]
        V_x =[Dict(i =>starter[:V][i] for i in W1[k]) for k=1:par[:n_part]]

        T_z = Dict(i =>starter[:T][i] for i in Wz)
        V_z = Dict(i =>starter[:V][i] for i in Wz)
        
        T_y =[Dict(i =>starter[:Tdu][k][i] for i in Wy[k]) for k=1:par[:n_part]]
        V_y =[Dict(i => starter[:Vdu][k][i] for i in Wy[k]) for k=1:par[:n_part]]
    end

    cont = true
    obj = zeros(par[:n_part])
    AL = zeros(par[:n_part])
    tol_pr = nothing
    tol_du = nothing
    if rank == 0
        lyap_save = []
        AL_save = []
        obj_save= []
        pr_save = []
        du_save = []
        time_save=[]
    end
    
    while cont        
        if rank != 0
            k = rank
            @objective(m,Min,
                       sum(net[:c][i,1] + net[:c][i,2]*m[:P][k,i] + net[:c][i,3]*m[:P][k,i]^2 for i in G[k])
                       +par[:price]/par[:scale]*sum(m[:Ps_abs][k,i] + m[:Qs_abs][k,i] for i in 1:net[:n] for i in W[k])
                       +sum(T_y[k][i]*(m[:T][k,i]-T_z[i]) +V_y[k][i]*(m[:V][k,i]-V_z[i]) for i in Wy[k])
                       +1/2*par[:rho]*sum((T_z[i]-m[:T][k,i])^2 + (V_z[i]-m[:V][k,i])^2 for i in Wy[k]))
            solve(m)
            for i in Wy[k]
                T_x[k][i] = getvalue(m[:T][rank,i])
                V_x[k][i] = getvalue(m[:V][rank,i])
            end
            AL[k] = getobjectivevalue(m) * par[:scale]
            obj[k] = AL[k]
                      -sum(T_y[k][i]*(T_x[k][i]-T_z[i]) +V_y[k][i]*(V_x[k][i]-V_z[i]) for i in Wy[k]) * par[:scale]
                      -1/2*par[:rho]*sum((T_z[i]-T_x[k][i])^2 + (V_z[i]-V_x[k][i])^2 for i in Wy[k]) * par[:scale]
        end

        for k=1:par[:n_part]
            T_x[k] = MPI.bcast(T_x[k], k, comm)
            V_x[k] = MPI.bcast(V_x[k], k, comm)
            obj[k] = MPI.bcast(obj[k], k, comm)
            AL[k] = MPI.bcast(AL[k], k, comm)
        end

        if rank == 0
            # Update s
            for i in Wz
                T_tmp = 0
                V_tmp = 0
                for k in zmap[i]
                    T_tmp += T_x[k][i]
                    V_tmp += V_x[k][i]
                end
                T_s[i] = T_tmp/length(zmap[i]) - T_z[i]
                V_s[i] = V_tmp/length(zmap[i]) - V_z[i]
            end

            # Update z
            for i in Wz
                T_z[i] += T_s[i]
                V_z[i] += V_s[i]
            end

            # Update r
            for k=1:par[:n_part]
                for i in Wy[k]
                    T_r[k][i] = T_x[k][i] - T_z[i]
                    V_r[k][i] = V_x[k][i] - V_z[i]
                end
            end

            # Update y
            for k=1:par[:n_part]
                for i in Wy[k]
                    T_y[k][i] += par[:rho] * T_r[k][i]
                    V_y[k][i] += par[:rho] * V_r[k][i]
                end
            end
            
            # Update error
            tmp = vcat([[collect(values(T_r[k]));collect(values(V_r[k]))] for k=1:par[:n_part]]...)
            err_r = norm(tmp)
            err_s = par[:rho]*norm([collect(values(T_s)) collect(values(V_s))])
            
            tmp = [vcat([collect(values(T_x[k])) for k=1:par[:n_part]]...);
                   vcat([collect(values(V_x[k])) for k=1:par[:n_part]]...)]
            tmp2= [collect(values(T_z));
                   collect(values(V_z))]
            tmp3= [vcat([collect(values(T_y[k])) for k=1:par[:n_part]]...);
                   vcat([collect(values(V_y[k])) for k=1:par[:n_part]]...)]
            
            tol_pr = sqrt(p) * tol_abs + tol_rel * max(norm(tmp),norm(tmp2))
            tol_du = sqrt(n) * tol_abs + tol_rel * norm(tmp3)
            @printf("%4i %11.6e %7.2e %7.2e %7.2e\n",iter,sum(obj),err_r/tol_pr,err_s/tol_du,par[:rho])
            iter += 1
            cont = ((err_r > tol_pr) || (err_s > tol_du)) && (iter < par[:maxiter])

            # Lyapunov function
            lyap = sum((1/par[:rho])*(V_y[k][i] - V_y_cen[k][i])^2
                       +(1/par[:rho])*(T_y[k][i] - T_y_cen[k][i])^2
                       +par[:rho]*(V_z[i] - V_z_cen[i])^2
                       +par[:rho]*(T_z[i] - T_z_cen[i])^2 for k in 1:par[:n_part] for i in Wy[k])

            # Varying penalty parameter
            # if err_r > 10 * err_s
            #     par[:rho] *= 2
            # elseif err_r < err_s / 10
            #     par[:rho] /= 2
            # end

            push!(lyap_save,lyap)
            push!(AL_save,sum(AL))
            push!(obj_save,sum(obj))
            push!(pr_save,err_r)
            push!(du_save,err_s)
            push!(time_save,(time_ns()-t0)/1e9)
        end
        # Bcast cont
        for k=1:par[:n_part]
            T_y[k] = MPI.bcast(T_y[k], 0, comm)
            V_y[k] = MPI.bcast(V_y[k], 0, comm)
        end
        T_z = MPI.bcast(T_z, 0, comm)
        V_z = MPI.bcast(V_z, 0, comm)
        cont = MPI.bcast(cont, 0, comm)
        # par[:rho] = MPI.bcast(par[:rho], 0, comm)
    end                         # while
if rank == 0
    @printf("Elpased time is: %5.2f seconds\n", (time_ns()-t0)/1e9)

    if mode == :dec
        result[:dec_lyap] = lyap_save
        result[:dec_AL] = AL_save
        result[:dec_obj]= obj_save
        result[:dec_pr] = pr_save
        result[:dec_du] = du_save
        result[:dec_time] = time_save
        result[:dec_tol_pr] = tol_pr
        result[:dec_tol_du] = tol_du
    elseif mode == :mul
        result[:mul_lyap] = lyap_save
        result[:mul_AL] = AL_save
        result[:mul_obj]= obj_save
        result[:mul_pr] = pr_save
        result[:mul_du] = du_save
        result[:mul_time] = time_save
        result[:mul_tol_pr] = tol_pr
        result[:mul_tol_du] = tol_du
    end
end
end                             # for mul dec
if rank == 0
    result[:cen_n_var]=MathProgBase.numvar(m_cen)
    result[:cen_n_constr]=MathProgBase.numconstr(m_cen)
    result[:cen_time]=time_cen
    result[:cen_obj] =getobjectivevalue(m_cen) * par[:scale]
    result[:crs_obj] =getobjectivevalue(m_crs) * par[:scale]
    
    result[:crs_n_var]=MathProgBase.numvar(m_crs)
    result[:crs_n_constr]=MathProgBase.numconstr(m_crs)
    result[:crs_time]=time_crs

    result[:warm_time]=time_warm

    writecsv("output/dec-$cs.csv",[result[:dec_pr] result[:dec_du] result[:dec_obj] result[:dec_AL] result[:dec_lyap] result[:dec_time]])
    writecsv("output/mul-$cs.csv",[result[:mul_pr] result[:mul_du] result[:mul_obj] result[:mul_AL] result[:mul_lyap] result[:mul_time]])
    writecsv("output/stat-$cs.csv",[result[:cen_time], result[:crs_time], result[:warm_time],
                                    result[:cen_n_var],result[:cen_n_constr],
                                    result[:crs_n_var],result[:crs_n_constr],
                                    result[:dec_tol_pr],result[:dec_tol_du],
                                    result[:mul_tol_pr],result[:mul_tol_du],
                                    result[:cen_obj],result[:crs_obj],par[:rho]*par[:scale]])
end
end                             # for cs
MPI.Finalize()
