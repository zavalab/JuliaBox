using Plots, LaTeXStrings, PyCall
PyDict(pyimport("matplotlib")["rcParams"])["text.usetex"] = [true]
PyDict(pyimport("matplotlib")["rcParams"])["font.family"] = ["serif"]
pyplot()

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

freq = 30
for k=1:Int(N/freq)
    xref[freq*(k-1)+1:freq*k,[1,3,5]] .= ones(freq)*randn(3)'
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
    x[end,:] = par[:xN]
    
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

mkr = [:diamond,:rect,:circle]
clr = [RGB(1,0,0),RGB(0,1,0),RGB(0,0,1)]

p=plot(size=(400,200),box=true,grid=true,fontfamily="serif",
       xlim=(0,size(err_saves[1],1)),
       yscale=:log,
       ylabel=L"\textrm{primal error}",
       xlabel=L"\textrm{iteration steps}");
for k=1:length(err_saves)
    plot!(p,err_saves[k][:,1],marker=mkr[k],color=clr[k],linewidth=.5,markersize=4,markerstrokecolor=clr[k],markercolor=:transparent,label=LaTeXString("\$\\omega=$(Oms[k])\$"))
end
plot!(p,par[:tol]*ones(size(err_saves[1],1)+1),color=:black,linestyle=:dash,linewidth=.5,label="")
savefig(p,"inf-x.pdf")

p=plot(size=(400,200),box=true,grid=true,fontfamily="serif",
       xlim=(0,size(err_saves[1],1)),
       yscale=:log,
       ylabel=L"\textrm{dual error}",
       xlabel=L"\textrm{iteration steps}");
for k=1:length(err_saves)
    plot!(p,err_saves[k][:,2],marker=mkr[k],color=clr[k],linewidth=.5,markersize=4,markerstrokecolor=clr[k],markercolor=:transparent,label=LaTeXString("\$\\omega=$(Oms[k])\$"))
end
plot!(p,par[:tol]*ones(size(err_saves[1],1)+1),color=:black,linestyle=:dash,linewidth=.5,label="")
savefig(p,"inf-l.pdf")


######################################################################
# one-by-one

maxiter = 6

err_saves = []; ts=[]
Om = 6
srand(1)
par = getpar(xref,uref,x0,xN,K=K,N=N,Om=Om,dt=dt,tol=tol,mu=mu)
par[:solver]=IpoptSolver(print_level=0,linear_solver="ma27",tol=1e-10)


# Initialize decentral
x = ones(par[:N]) * par[:x0]'
l = ones(par[:N]) * par[:lN]'
u = ones(par[:N]) * par[:uN]'
x[end,:] = par[:xN]

r = [remotecall(do_init,workers()[k],par,k) for k=1:par[:K]]
fetch.(r)

# Iteration
trajs=[]
cnt=0
# for fixedcounter = 1:20
for cnt=1:maxiter
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
    # @printf("%4i %7.2e %7.2e\n",cnt,err_pr,err_du)
    
    traj=[remotecall_fetch(get_sol,workers()[k]) for k=1:par[:K]]
    push!(trajs,traj)
end

for cnt = 1:maxiter
    clr = [RGB(1,0,0),RGB(0,1,0),RGB(0,0,1),RGB(.5,.5,.5)]
    mkr = [:diamond,:rect,:circle,:dtriangle]
    p=plot(size=(500,500),color=:black,fontfamily="serif",
           leg=false, box=true,xlabel=L"$X$",ylabel=L"$Y$",zlabel=L"$Z$");
    for k=1:par[:K]
        x = trajs[cnt][k][1]
        plot!(p,
              x[:,1],x[:,3],x[:,5],
              color=clr[k],linewidth=.5,
              markershape=mkr[k],markersize=4,
              markerstrokecolor=clr[k],markercolor=:transparent)
    end
    plot!(p,x_cen[:,1],x_cen[:,3],x_cen[:,5],color=:black);
    plot!(p,[x_cen[1,1]],[x_cen[1,3]],[x_cen[1,5]],
          marker=:circle,markersize=8,markerstrokewidth=1,markercolor=:transparent);
    plot!(p,[x_cen[end,1]],[x_cen[end,3]],[x_cen[end,5]],
          marker=:diamond,markersize=8,markerstrokewidth=1,markercolor=:transparent);
    savefig(p,"iter-x-$cnt.pdf")

    p=plot(size=(500,500),color=:black,fontfamily="serif",
           leg=false, box=true,xlabel=L"$\lambda_X$",ylabel=L"$\lambda_Y$",zlabel=L"$\lambda_Z$");
    for k=1:par[:K]
        l = trajs[cnt][k][2]
        plot!(p,
              l[:,1],l[:,3],l[:,5],
              color=clr[k],linewidth=.5,
              markershape=mkr[k],markersize=4,
              markerstrokecolor=clr[k],markercolor=:transparent)
    end
    plot!(p,l_cen[:,1],l_cen[:,3],l_cen[:,5],color=:black);
    plot!(p,[l_cen[1,1]],[l_cen[1,3]],[l_cen[1,5]],
          marker=:circle,markersize=8,markerstrokewidth=1,markercolor=:transparent);
    plot!(p,[l_cen[end,1]],[l_cen[end,3]],[l_cen[end,5]],
          marker=:diamond,markersize=8,markerstrokewidth=1,markercolor=:transparent);
    savefig(p,"iter-l-$cnt.pdf")
end

######################################################################
# sensitivity study
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
    status = solve(m_cen); 
    
    x = getvalue(m_cen[:x][:,:])
    u = getvalue(m_cen[:u][:,:])
    l = getdual(m_cen[:l][:,:])
    
    push!(pars,par)
    push!(x_pers,x[1:end-1,:])
    push!(l_pers,l)
    push!(u_pers,u)
    push!(objs,getobjectivevalue(m_cen))
end


p=plot(size=(500,500),color=:black,fontfamily="serif",grid=true,
       leg=false, box=true,xlabel=L"$X$",ylabel=L"$Y$",zlabel=L"$Z$");
for x_per in x_pers
    plot!(p,
          x_per[:,1],x_per[:,3],x_per[:,5],
          color=:lightblue,linewidth=.5)
end
plot!(p,x_cen[:,1],x_cen[:,3],x_cen[:,5],color=:black);
plot!(p,[x_cen[1,1]],[x_cen[1,3]],[x_cen[1,5]],
      marker=:circle,markersize=8,markerstrokewidth=1,markercolor=:transparent);
plot!(p,[x_cen[end,1]],[x_cen[end,3]],[x_cen[end,5]],
      marker=:diamond,markersize=8,markerstrokewidth=1,markercolor=:transparent);
savefig(p,"x.pdf")

p=plot(size=(500,500),color=:black,fontfamily="serif",
       leg=false, box=true,xlabel=L"$\lambda_X$",ylabel=L"$\lambda_Y$",zlabel=L"$\lambda_Z$");
for l_per in l_pers
    plot!(p,
          l_per[:,1],l_per[:,3],l_per[:,5],
          color=:lightblue,linewidth=.5)
end
plot!(p,l_cen[:,1],l_cen[:,3],l_cen[:,5],color=:black);
plot!(p,[l_cen[1,1]],[l_cen[1,3]],[l_cen[1,5]],
      marker=:circle,markersize=8,markerstrokewidth=1,markercolor=:transparent);
plot!(p,[l_cen[end,1]],[l_cen[end,3]],[l_cen[end,5]],
      marker=:diamond,markersize=8,markerstrokewidth=1,markercolor=:transparent);
savefig(p,"l.pdf")
