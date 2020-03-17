# Written by Sungho Shin (sungho.shin@wisc.edu), Mihai Anitescu (anitescu@mcs.anl.gov), and Victor Zavala (victor.zavala@wisc.edu)

using Distributed
K = 16
if nprocs() < K+1
    addprocs(-nprocs()+K+1)
end
@sync @everywhere include("Schwarz.jl")
Random.seed!(1)


# path = "../../pglib-opf/pglib_opf_case300_ieee.m"
# path = "../../pglib-opf/pglib_opf_case1354_pegase.m"
# path = "../../pglib-opf/pglib_opf_case2000_tamu.m"
path = "../../pglib-opf/pglib_opf_case9241_pegase.m"
# path = "../../pglib-opf/pglib_opf_case10000_tamu.m"
# path = "../../pglib-opf/pglib_opf_case13659_pegase.m"


args = get_data(path)

gammas = [1e5]
oms = [1 5 10]
max_iter = 1000
tol_pr = 1e-2; tol_du = 1e2


for cnt_gamma=1:length(gammas)
    gamma = gammas[cnt_gamma]
    # non-regularized centralized solution
    m_zr = Model(with_optimizer(Ipopt.Optimizer,linear_solver="ma57",tol=1e-12))
    # m_zr = Model(with_optimizer(Gurobi.Optimizer))
    va_zr = Array{Any,1}(zeros(args[:N]))
    s_zr = Array{Any,1}([zeros(args[:ng][i]) for i=1:args[:N]])
    pg_zr = Array{Any,1}(zeros(args[:N]))
    sep_zr = Array{Any,1}(zeros(args[:N]))
    pg_exp_zr = Array{Any,1}(zeros(args[:N]))
    sep_exp_zr = Array{Any,1}(zeros(args[:N]))
    obj_exp_zr = Array{Any,1}(zeros(args[:N]))

    dcopf!(m_zr,1:args[:N],[],va_zr,s_zr,pg_zr,sep_zr,
           pg_exp_zr,sep_exp_zr,obj_exp_zr,args,gamma=0)
    tic()
    optimize!(m_zr)
    tt_zr = toc()
    obj_zr = objective_value(m_zr)

    va_val_zr = value.(va_zr)
    s_val_zr = [value.(s_zr[i]) for i=1:args[:N]]
    pg_val_zr = dual.(pg_zr)
    sep_val_zr = [dual.(sep_zr[i]) for i=1:args[:N]]


    # centralized solution
    m_cen = Model(with_optimizer(Ipopt.Optimizer,linear_solver="ma57"))
    # m_cen = Model(with_optimizer(Gurobi.Optimizer))
    va_cen = Array{Any,1}(zeros(args[:N]))
    s_cen = Array{Any,1}([zeros(args[:ng][i]) for i=1:args[:N]])
    pg_cen = Array{Any,1}(zeros(args[:N]))
    sep_cen = Array{Any,1}(zeros(args[:N]))
    pg_exp_cen = Array{Any,1}(zeros(args[:N]))
    sep_exp_cen = Array{Any,1}(zeros(args[:N]))
    obj_exp_cen = Array{Any,1}(zeros(args[:N]))


    dcopf!(m_cen,1:args[:N],[],va_cen,s_cen,pg_cen,sep_cen,
           pg_exp_cen,sep_exp_cen,obj_exp_cen,args,gamma=gamma)
    tic()
    optimize!(m_cen)
    tt_cen = toc()
    obj_cen = objective_value(m_cen)

    va_val_cen = value.(va_cen)
    s_val_cen = [value.(s_cen[i]) for i=1:args[:N]]
    pg_val_cen = dual.(pg_cen)
    sep_val_cen = [dual.(sep_cen[i]) for i=1:args[:N]]

    for cnt_om=1:length(oms)
        # decentralized solution
        om = oms[cnt_om]

        part = Metis.partition(args[:g],K,alg=:KWAY)
        Vs   = [findall(part.==k) for k=1:K]
        V_compls = [setdiff(closed_neighbors(args[:g],Vs[k],1),Vs[k]) for k=1:K]
        V_oms = [closed_neighbors(args[:g],Vs[k],om) for k=1:K]
        V_om_compls= [setdiff(closed_neighbors(args[:g],V_oms[k],1),V_oms[k]) for k=1:K]
        V_compls_union=union([V_compls;V_om_compls]...)
        V_gives = [intersect(V_compls_union,Vs[k]) for k=1:K]

        va = Array{Any,1}(zeros(args[:N]))
        s = Array{Any,1}([zeros(args[:ng][i]) for i=1:args[:N]])
        pg = Array{Any,1}(zeros(args[:N]))
        sep = Array{Any,1}(zeros(args[:N]))

        va_val = Array{Any,1}(zeros(args[:N]))
        s_val = [zeros(args[:ng][i]) for i=1:args[:N]]
        pg_val = Array{Any,1}(zeros(args[:N]))
        sep_val = [zeros(length(neighbors(args[:g],i))) for i=1:args[:N]]

        r=[remotecall(Schwarz_init,workers()[k],Vs[k],V_oms[k],V_compls[k],V_om_compls[k],V_gives[k],gamma,args) for k=1:K]
        fetch.(r)
        tbl=[]

        tt =-time()
        for iter =0:max_iter
            # tic()

            r_solve = [remotecall(Schwarz_iter_solve,workers()[k],
                                  va_val[V_om_compls[k]],pg_val[V_om_compls[k]],sep_val[V_om_compls[k]],om) for k=1:K]
            for k=1:K
                va_val[V_gives[k]],pg_val[V_gives[k]],sep_val[V_gives[k]] = fetch(r_solve[k])
            end
            
            r_update = [remotecall(Schwarz_iter_update,workers()[k],
                                   va_val[V_compls[k]],pg_val[V_compls[k]],sep_val[V_compls[k]]) for k=1:K]
            obj = 0; err_pr = 0; err_du=0; slk = 0;
            for k=1:K
                obj_k,err_pr_k,err_du_k,slk_k = fetch(r_update[k]);
                obj += obj_k;
                slk = max(slk,slk_k)
                err_pr = max(err_pr,err_pr_k);
                err_du = max(err_du,err_du_k);
            end
            if mod(iter,10)==0
                @printf("iter   objective     err_pr   err_du   slk      time\n")
            end
            @printf("%4i  %11.7e %6.2e %6.2e %6.2e %6.2e\n",iter,obj,err_pr,err_du,slk,tt+time())
            push!(tbl ,[iter,obj,err_pr,err_du,slk,tt+time()])
            if err_pr < tol_pr && err_du < tol_du
                break
            end
        end
        tt+=time()
        for k=1:K
            s_val[Vs[k]]=remotecall_fetch(Schwarz_get_value,workers()[k])
        end
        c = cost(s_val,1:args[:N],args)
        c_zr = cost(s_val_zr,1:args[:N],args)

        tbl = hcat(tbl...)'
        writedlm("output/tbl-$cnt_gamma-$cnt_om.csv", tbl, ',')
        writedlm("output/etc-$cnt_gamma-$cnt_om.csv", [c,c_zr,tt,tt_cen,obj_cen], ',')
    end
end

run(`python3 plotter.py`)
