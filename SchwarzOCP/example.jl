srand(1)

# Set parallel
ncores=16
N = 240
K = 4   # number of subproblems
dt= 1/10
tol=1e-6
mu=10
sig=1e-1

x0 = [1,0,1,0,1,0,0,0,0]
xN = [0,0,0,0,0,0,0,0,0]

xref = zeros(N,9)

for k=1:Int(N/30)
    xref[30*(k-1)+1:30*k,[1,3,5]] .= ones(30)*randn(3)'
end

uref = [9.8,0,0,0]    
uref = kron(ones(N),uref')

is_progress=true
progress=[]

if nprocs() < ncores+1
    addprocs(-nprocs()+ncores+1)
end
@everywhere include("quadrotor.jl")

# ######################################################################
# Central

# Get Param
srand(1)
par = getpar(xref,uref,x0,xN,K=K,N=N,dt=dt,tol=tol,mu=mu)
par[:solver]=IpoptSolver(print_level=0,linear_solver="ma27",tol=par[:tol])

# Solve central
t_cen,x_cen,l_cen,m_cen,status_cen = do_centralized(par)

######################################################################
# Decentral
err_saves = []; ts=[]
Oms = [3 6 9]
for Om in Oms
    srand(1)
    par = getpar(xref,uref,x0,xN,K=K,N=N,Om=Om,dt=dt,tol=tol,mu=mu)
    par[:solver]=IpoptSolver(print_level=0,linear_solver="ma27",tol=1e-10)


    # Initialize decentral
    x = ones(par[:N]) * par[:x0]'
    l = ones(par[:N]) * par[:lN]'
    u = ones(par[:N]) * par[:uN]'

    err_pr = Inf
    err_du = Inf

    r = [remotecall(do_init,workers()[k],par,k) for k=1:par[:K]]
    fetch.(r)

    # Iteration
    err_save = []
    t =-time()
    cnt=0
    # for fixedcounter = 1:20
    while err_pr > par[:tol] || err_du > par[:tol]
        prf=x[par[:n1][2:end],:]; duf=l[par[:n1][2:end],:]
        r = [remotecall(do_iter,workers()[k],
                        x[par[:n1Om][k],:],x[par[:n2Om][k],:],
                        l[par[:n2Om][k],:],u[par[:n2Om][k],:]) for k=1:par[:K]]
        for k=1:par[:K]
            x[par[:ind][k],:],l[par[:ind][k],:],u[par[:ind][k],:],prfk,dufk = fetch(r[k])
            if k!=par[:K]
                prf[k,:]-=prfk; duf[k,:]-=dufk;
            end
        end
        err_pr=norm(prf[:],Inf); err_du=norm(duf[:],Inf)
        @printf("%4i %7.2e %7.2e\n",cnt,err_pr,err_du)
        push!(err_save,[err_pr err_du])
        cnt+=1
    end
    t+= time()

    println("central: $t_cen decentral: $t")
    push!(err_saves,vcat(err_save...))
    push!(ts,t)
end

writecsv("ts.csv",ts)
plot_err(par,err_saves,Oms=Oms)



# Get Param
srand(1)
par = getpar(xref,uref,x0,xN,K=1,N=Int(N/K),dt=dt,tol=tol,mu=mu)
par[:solver]=IpoptSolver(print_level=0,linear_solver="ma27")

t_cen,x_cen,l_cen,m_cen,status_cen = do_centralized(par)

# Sensitivity test
x_pers=[];l_pers=[];u_pers=[];pars=[]; objs=[]
for i=1:30
    println(i)
    par[:x0]=x0 + sig*randn(par[:nx])
    par[:xN]=xN + sig*randn(par[:nx])
    par[:uN]= sig*randn(par[:nu])
    par[:lN]= 10*sig*randn(par[:nx])

    do_setpar(m_cen,par[:x0],par[:xN],par[:lN],par[:uN])
    status = solve(m_cen); status!=:Optimal && continue
    
    x = getvalue(m_cen[:x][:,:])
    u = getvalue(m_cen[:u][:,:])
    l = getdual(m_cen[:l][:,:])
    
    push!(pars,par)
    push!(x_pers,x[1:end-1,:])
    push!(l_pers,l)
    push!(u_pers,u)
    push!(objs,getobjectivevalue(m_cen))
end
outs=Dict(:par=>par,:x_cen=>x_cen[1:end-1,:],:l_cen=>l_cen,:x_pers=>x_pers,:u_pers=>u_pers,:l_pers=>l_pers)
plot_cen_T(outs[:par],outs[:x_cen],outs[:l_cen],
         outs[:x_pers],outs[:l_pers],ks=[1,3,5])



# # scp -r shin79@baptiste.che.wisc.edu:/home/shin79/argonne-etc/quadrotor/*.pdf .
