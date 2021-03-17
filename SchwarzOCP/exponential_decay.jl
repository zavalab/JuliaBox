# exponential decay of sensitivity
dt= 0.005
K = 1:20
T = 60
N = round(Int,T/dt)
nsam = 30
per_std = .5

m_ref = hehnandrea(Ipopt.Optimizer,N,K;dt=dt)
optimize!(m_ref)

x_ref = value.(m_ref[:x][:,1])
y_ref = value.(m_ref[:x][:,3])
z_ref = value.(m_ref[:x][:,5])
lx_ref = dual.(m_ref[:c][:,1])
ly_ref = dual.(m_ref[:c][:,3])
lz_ref = dual.(m_ref[:c][:,5])
x_pers = []
y_pers = []
z_pers = []
lx_pers = []
ly_pers = []
lz_pers = []

for i=1:nsam
    x0_new = [randn(),0,randn(),0,randn(),0,0,0,0] * .3
    setvalue.(m_ref[:x0],x0_new)
    optimize!(m_ref)
    push!(x_pers, value.(m_ref[:x][:,1]))
    push!(y_pers, value.(m_ref[:x][:,3]))
    push!(z_pers, value.(m_ref[:x][:,5]))
    push!(lx_pers, dual.(m_ref[:c][:,1]))
    push!(ly_pers, dual.(m_ref[:c][:,3]))
    push!(lz_pers, dual.(m_ref[:c][:,5]))
end

plt1,plt2 = plot_eds(x_ref,y_ref,z_ref,x_pers,y_pers,z_pers,lx_ref,ly_ref,lz_ref,lx_pers,ly_pers,lz_pers)
savefig(plt1,"fig/eds-x.pdf")
savefig(plt2,"fig/eds-l.pdf")
