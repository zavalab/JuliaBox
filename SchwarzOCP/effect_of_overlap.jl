# quadrotor effect of overlap
dt= 0.005
K = 1:20
T = 60
N = round(Int,T/dt)
omegas = [.3,.5,1.]

err_schwarz=[]
ms = Dict()
for omega in omegas
    m = hehnandrea(Ipopt.Optimizer,N,K;print_level=0,dt=dt)
    ms[omega] = m
    schwarz = SchwarzModel(m;rho=1,omega=omega)
    set_KKT_error_evaluator!(m)
    start = time()
    push!(err_schwarz,Float64[get_KKT_error(m)])
    @time optimize!(schwarz;tol=1e-6,optional=schwarz->push!(err_schwarz[end],get_KKT_error(schwarz.model)))
end

plt = plot_err_profile(err_schwarz,omegas)
savefig(plt,"fig/err_profile.pdf")
