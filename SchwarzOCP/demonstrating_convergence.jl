# demonstrating conergence
dt= 0.005
K = 1:3
T = 60
N = round(Int,T/dt)

m_ref = hehnandrea(Ipopt.Optimizer,N,K;dt=dt)
optimize!(m_ref)

m = hehnandrea(Ipopt.Optimizer,N,K;print_level=0,dt=dt)
schwarz = SchwarzModel(m;rho=1,omega=.1);

iters = [[] for k in K]
optimize!(schwarz;tol=1e-6,
          optional=schwarz->[push!(iters[k],copy(schwarz.submodels[k].model.x)) for k in K]);

idx1 = []
idx2 = []
idx3 = []
for k in K
    push!(idx1,indexin(SimpleNLModels.index.(m[:x][:,1]),m[:Ws][k]))
    idx1[k]=idx1[k][idx1[k].!=nothing]
    push!(idx2,indexin(SimpleNLModels.index.(m[:x][:,3]),m[:Ws][k]))
    idx2[k]=idx2[k][idx2[k].!=nothing]
    push!(idx3,indexin(SimpleNLModels.index.(m[:x][:,5]),m[:Ws][k]))
    idx3[k]=idx3[k][idx3[k].!=nothing]
end
coords = [[(view(itr,idx1[k]),view(itr,idx2[k]),view(itr,idx3[k])) for itr in iters[k]] for k in K]

x_ref = value.(m_ref[:x][:,1])
y_ref = value.(m_ref[:x][:,3])
z_ref = value.(m_ref[:x][:,5])

plts = plot_demonstration(coords,x_ref,y_ref,z_ref)

for i in eachindex(plts)
    savefig(plts[i], "fig/iter-$i.pdf")
end
