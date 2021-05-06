using SimpleNL, SimpleNLUtils, SimpleSchwarz, SimpleADMM, Ipopt, Random

# include("case_study_include.jl")

for (name,compute_vals,perturb,N,K,NS,j,tol,omega,omegas,rho,rhos,etabs,etablabels,etabs_eds,etabs_eds_labels) = [
    (:hehnandrea,:compute_vals_hehnandrea,:perturb_hehnandrea,
     12000,20,40,25,1e-5,
     1.,[0.1,0.5,1.0],
     1.,[1e-2,1.,1e2],
     [(1,1),(.1,1e0),(1,.1),(.1,.1)],
     ["(1,1)","(0.1,1)","(1,0.1)","(0.1,0.1)"],
     [(1,1),(.1,1e0),(1,.1),(.1,.1)],
     ["(1,1)","(0.1,1)","(1,0.1)","(0.1,0.1)"])
    (:storage,:compute_vals_storage,:perturb_storage,
     12,4,8,9,1e-5,
     1.2,[.3,.6,1.2],
     1.,[1e-2,1.,1e2],
     [(1,1),(.01,1e0),(1,.01),(.01,.01)],
     ["(1,1)","(0.01,1)","(1,0.01)","(0.01,0.01)"],
     [(1,1),(.01,1e0),(1,.01),(.01,.01)],
     ["(1,1)","(0.01,1)","(1,0.01)","(0.01,0.01)"])
    (:thinplate,:compute_vals_thinplate,:perturb_thinplate,
     400,5,9,41,1e-5,
     2.,[0.1,0.5,1.0],
     1.,[.5,1.,2.],
     [(1,1),(.1,1e0),(1,.1),(.1,.1)],
     ["(1,1)","(0.1,1)","(1,0.1)","(0.1,0.1)"],
     [(1,1),(.1,1),(1,.1),(.1,.1)],
     ["(1,1)","(0.1,1)","(1,0.1)","(0.1,0.1)"])
    (:acopf,:compute_vals_acopf,:perturb_acopf,
     0,20,1,105,1e-5,
     1.,[0.1,0.3,1.0],
     1.,[1e-2,1.,1e2],
     [(1e6,1),(1e2,1),(1e6,1e-1),(1e2,1e-1)],
     ["(10^6,1)","(10^2,1)","(10^6,0.1)","(10^2,0.1)"],
     [(1e2,1),(1e0,1),(1e2,5e-2),(1e0,5e-2)],
     ["(10^2,1)","(1,1)","(10^2,0.05)","(1,0.05)"])
]

    get_model = @eval $name
    compute_vals = @eval $compute_vals
    perturb = @eval $perturb

    # case-by-case fix
    edsopt = name == :hehnandrea ? [:dt=>.1] : []
    rho_admm() = name == :acopf ? 1e4 : 1 
    
    # # exponential decay of sensitivity
    # println("** EDS")
    # dist = nothing
    # outputs = []
    # for (eta,b) in etabs_eds
    #     m_ref = get_model(Ipopt.Optimizer,NS,2;eta=eta,b=b,edsopt...);
    #     m_per = get_model(Ipopt.Optimizer,NS,2;eta=eta,b=b,edsopt...);
    #     optimize!(m_ref)
    #     vals = zeros(nv(m_ref[:g]))
    #     Random.seed!(13)
    #     for i=1:30
    #         perturb(m_ref,m_per,j,1e-3)
    #         optimize!(m_per)
    #         vals.=max.(vals,compute_vals(m_ref,m_per))
    #     end
    #     vals./= maximum(vals)
    #     pyplot()
    #     plt = plot_graph(m_ref[:g],m_ref[:pos],vals,j)
    #     savefig(plt,"../fig/$name-graph-$eta-$b.pdf")
    #     gr()

    #     dist = gdistances(m_ref[:g],j)
    #     push!(outputs,vals)
    # end
    # plt = plot_eds(dist,outputs,"(\\eta,b)",etabs_eds_labels)
    # savefig(plt,"../fig/$name-eds.pdf")

    # effect of conditioning
    println("** effect of conditioning")
    outputs = []
    for (eta,b) in etabs
        m_schwarz  = get_model(SimpleSchwarz.Optimizer,N,K;eta=eta,b=b,rho=rho,omega=omega,maxiter=40,tol=tol,save_output=true) # ,subproblem_option=Dict(:print_level=>5)
        SimpleSchwarz.instantiate!(m_schwarz)
        SimpleSchwarz.optimize!(m_schwarz)
        push!(outputs,m_schwarz[:output])
    end
    plt1, plt2 =plot_err_profile(outputs,"(\\eta,b)",etablabels)
    savefig(plt1,"../fig/$name-etab.pdf")
    savefig(plt2,"../fig/$name-etab-t.pdf")

    # effect of overlap
    println("** effect of overlap")
    outputs = []
    for omega in omegas
        m_schwarz  = get_model(SimpleSchwarz.Optimizer,N,K;rho=rho,omega=omega,maxiter=40,tol=tol,save_output=true)
        SimpleSchwarz.instantiate!(m_schwarz)
        SimpleSchwarz.optimize!(m_schwarz)
        push!(outputs,m_schwarz[:output])
    end

    plt1, plt2 =plot_err_profile(outputs,"\\widetilde{\\omega}",omegas)
    savefig(plt1,"../fig/$name-overlap.pdf")
    savefig(plt2,"../fig/$name-overlap-t.pdf")
    
    # effect of penalty
    println("** effect of penalty")
    outputs = []
    for rho in rhos
        m_schwarz  = get_model(SimpleSchwarz.Optimizer,N,K;rho=rho,omega=omega,maxiter=40,tol=tol,save_output=true)
        SimpleSchwarz.instantiate!(m_schwarz)
        SimpleSchwarz.optimize!(m_schwarz)
        push!(outputs,m_schwarz[:output])
    end

    plt1, plt2 =plot_err_profile(outputs,"\\mu",rhos)
    savefig(plt1,"../fig/$name-penalty.pdf")
    savefig(plt2,"../fig/$name-penalty-t.pdf")
    
    # benchmark
    println("** benchmark against ADMM and Ipopt")
    m_ref = get_model(Ipopt.Optimizer,N,K)
    instantiate!(m_ref)
    GC.gc(); GC.enable(false)
    time_ref = @elapsed begin
        optimize!(m_ref)
    end
    GC.enable(true); GC.gc()
    err_ref = SimpleNLUtils.KKTErrorEvaluator(m_ref)(m_ref.x,m_ref.l,m_ref.gl)

    m_schwarz  = get_model(SimpleSchwarz.Optimizer,N,K;rho=rho,omega=omega,maxiter=40,tol=tol,save_output=true)
    SimpleSchwarz.instantiate!(m_schwarz)
    GC.gc(); GC.enable(false)
    SimpleSchwarz.optimize!(m_schwarz)
    GC.enable(true); GC.gc()

    m_admm = get_model(SimpleADMM.Optimizer,N,K;print_level=0,rho=rho_admm(),maxiter=40,tol=tol,save_output=true)
    instantiate!(m_admm)
    GC.gc(); GC.enable(false)
    @time optimize!(m_admm)
    GC.enable(true); GC.gc()

    plt1,plt2 = plot_benchmark(err_ref,time_ref,m_schwarz[:output],m_admm[:output])
    savefig(plt1,"../fig/$name-benchmark.pdf")
    savefig(plt2,"../fig/$name-benchmark-t.pdf")
end
