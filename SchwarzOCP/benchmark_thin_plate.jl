# benchmark - thin plate

(K,T,L,dt,dx) = (1:20,3600*24,1,10.,.1) 
rhos = [1, 10, 100]

N = round(Int,T/dt)
n = round(Int,L/dx)

m_ref = thinplate(Ipopt.Optimizer,n,N,K;dt=dt,dx=dx)
instantiate!(m_ref)
GC.gc(); GC.enable(false)
time_ref = @elapsed begin
    @time optimize!(m_ref)
end
GC.enable(true); GC.gc()
set_KKT_error_evaluator!(m_ref)
err_ref = get_KKT_error(m_ref)

err_schwarzs = []
time_schwarzs = []
for rho in rhos
    m = thinplate(Ipopt.Optimizer,n,N,K;dt=dt,dx=dx,print_level=0)
    schwarz = SchwarzModel(m;rho=1,omega=1.)
    set_KKT_error_evaluator!(m)
    GC.gc(); GC.enable(false)
    start = time()
    err_schwarz = Float64[get_KKT_error(m)]
    time_schwarz = Float64[start]
    @time optimize!(schwarz;tol=0,maxiter=30,optional=schwarz->(push!(err_schwarz,get_KKT_error(schwarz.model)),push!(time_schwarz,time())))
    time_schwarz.-= start
    GC.enable(true); GC.gc()
    push!(err_schwarzs,err_schwarz)
    push!(time_schwarzs,time_schwarz)
end

err_admms = []
time_admms = []
for rho in rhos
    m = thinplate(Ipopt.Optimizer,n,N,K;dt=dt,dx=dx,print_level=0)
    admm = ADMMModel(m;rho=rho)
    set_KKT_error_evaluator!(m)
    GC.gc(); GC.enable(false)
    start = time()
    err_admm = Float64[get_KKT_error(m)]
    time_admm = Float64[start]
    @time optimize!(admm;tol=0,maxiter=30,optional=admm->(push!(err_admm,get_KKT_error(admm.model)),push!(time_admm,time())))
    time_admm.-= start
    GC.enable(true); GC.gc()
    
    push!(err_admms,err_admm)
    push!(time_admms,time_admm)
end

plt = plot_benchmark(err_ref,time_ref,err_schwarzs,time_schwarzs,err_admms,time_admms,rhos)
savefig(plt,"fig/pde.pdf")

