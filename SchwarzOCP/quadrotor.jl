using JuMP, Ipopt, PyPlot

function P(par,n1,n2)
    m = Model(solver=par[:solver])
    @variable(m, par[:xlb][j] <= x[i in n1:n2+1, j in 1:par[:nx]] <= par[:xub][j],
              start = 0)

    @variable(m, par[:ulb][j] <= u[i in n1:n2-1, j in 1:par[:nu]] <= par[:uub][j],
              start = 0)
    
    m[:x0] = @NLparameter(m, x0[1:par[:nx]] == 0)
    m[:xN] = @NLparameter(m, xN[1:par[:nx]] == 0)
    m[:lN] = @NLparameter(m, lN[1:par[:nx]] == 0)
    m[:uN] = @NLparameter(m, uN[1:par[:nu]] == 0)

    @constraintref l[n1:n2,1:par[:nx]] ; m[:l]=l

    for i in n1+1:n2+1
        if i==n2+1
            dyn=par[:f](m,x[i-1,:],uN)
        else
            dyn=par[:f](m,x[i-1,:],u[i-1,:])
        end
        for j in 1:par[:nx]
            m[:l][i-1,j] = @NLconstraint(m, x[i,j] == x[i-1,j] + dyn[j] * par[:dt])
        end
    end
    m[:l][n1,:] =  @NLconstraint(m, [j in 1:par[:nx]], x[n1,j] == x0[j])
    
    @NLobjective(m, Min,
                 sum(par[:Q][j]*(x[i,j]^2/2-par[:xref][i,j]*x[i,j])
                     for i in n1:n2 for j=1:par[:nx])
                 + sum(par[:R][j]*(u[i,j]^2/2-par[:uref][i,j]*u[i,j])
                       for i in n1:n2-1 for j=1:par[:nu])
                 + sum(par[:R][j]*(uN[j]^2/2-par[:uref][n2,j]*uN[j])
                       for j=1:par[:nu])
                 + sum(par[:mu]*(x[n2,j]^2/2-x[n2,j]*xN[j]) for j=1:par[:nx])
                 + sum(lN[j] * x[n2+1,j] for j=1:par[:nx]))
    
    return m
end


function getpar(xref,uref,x0,xN;seed=1,N=100,K=4,Om=1,dt=1/60,
                solver=IpoptSolver(print_level=0,linear_solver="ma57"),
                tol=1e-4,mu=1e-3)
    par = Dict()

    par[:g] = 9.8
        
    par[:nx]= 9
    par[:nu]= 4
    par[:Om]= Om
    par[:N] = N
    par[:K] = K
    par[:M] = Int64(N/K)
    par[:n1] = collect(0:par[:M]:par[:N])[1:end-1]+1
    par[:n2] = collect(0:par[:M]:par[:N])[2:end]+1
    par[:n2][end] -= 1
    par[:n1Om]= max.(par[:n1]-par[:Om],1)
    par[:n2Om]= min.(par[:n2]+par[:Om],par[:N])
    par[:mu]=mu


    par[:xlb] = [-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf]
    par[:xub] = [ Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
    par[:ulb] = -Inf*ones(4)
    par[:uub] = Inf*ones(4)

    function f(m,x,u)
        return [
            @NLexpression(m,x[2]) 
            @NLexpression(m,u[1]*cos(x[7])*sin(x[8])*cos(x[9])+u[1]*sin(x[7])*sin(x[9]))
            @NLexpression(m,x[4])
            @NLexpression(m,u[1]*cos(x[7])*sin(x[8])*sin(x[9])-u[1]*sin(x[7])*cos(x[9]))
            @NLexpression(m,x[6])
            @NLexpression(m,u[1]*cos(x[7])*cos(x[8])-par[:g])
            @NLexpression(m,u[2]*cos(x[7])/cos(x[8])+u[3]*sin(x[7])/cos(x[8]))
            @NLexpression(m,-u[2]*sin(x[7])+u[3]*cos(x[7]))
            @NLexpression(m,u[2]*cos(x[7])*tan(x[8])+u[3]*sin(x[7])*tan(x[8])+u[4])
        ]
    end
    
    par[:f] = f
    par[:dt]= dt
    
    par[:solver]=solver
    par[:tol]=tol
        
    par[:Q] = [10,1,2,1,10,1,1,1,1]
    par[:R] = [1,1,1,1]
    
    
    par[:xref] = xref
    par[:uref] = uref

    par[:x0] = x0
    par[:xN] = xN
    par[:lN] = zeros(9)
    par[:uN] = zeros(par[:nu])
    
    ind_gl = union(par[:n1Om][2:end],par[:n1][2:end],par[:n2Om][1:end-1])
    par[:ind] = [par[:n1][k]:par[:n2][k]-1 for k=1:par[:K]]
    
    return par
end


function do_centralized(par)
    m_cen = P(par,1,par[:N])

    do_setpar(m_cen,par[:x0],par[:xN],par[:lN],par[:uN])
    
    t_cen =-time()
    status=solve(m_cen)
    t_cen+= time()
    
    x_cen = getvalue(m_cen[:x][:,:])
    l_cen = getdual(m_cen[:l][:,:])
    
    return t_cen,x_cen,l_cen,m_cen,status
end


function do_init(par0,k0)
    global k = k0
    global par = par0
    global m = P(par,par[:n1Om][k],par[:n2Om][k])
    return 0
end


function do_iter(x0,xN,lN,uN)
    
    do_setpar(m,x0,xN,lN,uN)
    status = solve(m); status!=:Optimal && error("Suboptimal solution detected at $k")
    
    prf=nothing; duf=nothing
    xk=hcat([getvalue(m[:x][i,:]) for i in par[:ind][k]]...)'
    uk=hcat([getvalue(m[:u][i,:]) for i in par[:ind][k]]...)'
    lk=hcat([getdual(m[:l][i,:]) for i in par[:ind][k]]...)'
    if k!=par[:K]
        prf = getvalue(m[:x][par[:n2][k],:])
        duf = getdual(m[:l][par[:n2][k],:])
    end
    return xk, lk, uk, prf, duf
end

function do_setpar(m,x0,xN,lN,uN)
    setvalue.(m[:x0],x0)
    setvalue.(m[:xN],xN)
    setvalue.(m[:lN],lN)
    setvalue.(m[:uN],uN)
end


function plot_cen_T(par,x_cen,l_cen,x_pers,l_pers;ks=[1,2,3])
    ioff()
    rc("text",usetex="True")
    matplotlib["rcParams"][:update](Dict("font.size" => 8))
    
    lst = ["D","x","o"]
    clr = [(1, 0, 0),(0, 0, 1),(0, .9, 0)]

    tgrid1 =(0:(par[:N]-1))*par[:dt]
    tgrid2 =(0:(par[:N]-1))*par[:dt]

    ul = [L"T_c\,(\mathrm{K})",L"F\,(\mathrm{m^3/min})"]

    fg=figure(figsize=(3.25,3.9))

    for k=1:3
        subplot(3,1,k)
        for x_per in x_pers
            plot(tgrid1,x_per[:,ks[k]],color=(.8,.8,1),linewidth=.5)
        end
        plot(tgrid1,x_cen[:,ks[k]],color="black",linewidth=1)
        ylabel(LaTeXString("\${x}[$(ks[k])]\$"));
        grid("on",color="lightgray");
        xlim([tgrid1[1],tgrid1[end]+par[:dt]])
    end    
    xlabel(L"\textrm{time (sec)}");
    tight_layout()
    savefig("x.pdf",bbox_inches="tight")
    PyPlot.close_figs()

    fg=figure(figsize=(3.25,3.9))
    for k=1:3
        subplot(3,1,k)
        for l_per in l_pers
            plot(tgrid1,l_per[:,ks[k]],color=(.8,.8,1),linewidth=.5)
        end
        plot(tgrid1,l_cen[:,ks[k]],color="black",linewidth=1)
        ylabel(LaTeXString("\${\\lambda}[$(ks[k])]\$"));
        grid("on",color="lightgray");
        xlim([tgrid1[1],tgrid1[end]+par[:dt]])
    end    
    xlabel(L"\textrm{time (sec)}");
    tight_layout()
    savefig("l.pdf",bbox_inches="tight")
    PyPlot.close_figs()
end

function plot_err(par,err_saves;Oms=[1,2,3])
    ioff();
    rc("text",usetex="True")
    matplotlib["rcParams"][:update](Dict("font.size" => 8))
    lst = ["D-","x-","o-"]
    clr = [(1, 0, 0),(0, 0, 1),(0, .9, 0)]

    subplots(2,1,sharex=true,figsize=(3.25,2.6)); 
    subplot(2,1,1)
    for k=1:length(err_saves)
        plot(err_saves[k][:,1],color=clr[k],linewidth=.5,lst[k],markersize=4,markerfacecolor="none",label=LaTeXString("\$\\omega=$(Oms[k])\$"))
    end
    len=size(err_saves[1],1)
    plot(par[:tol]*ones(len+1),color="black","--",linewidth=.5)
    yscale("log")
    ylabel(L"\textrm{primal error}");
    grid("on",color="lightgray");
    xlim([0,len]);

    subplot(2,1,2)
    for k=1:length(err_saves)
        plot(err_saves[k][:,2],color=clr[k],linewidth=.5,lst[k],markersize=4,markerfacecolor="none",label=LaTeXString("\$\\omega=$(Oms[k])\$"))
    end
    ylabel(L"\textrm{dual error}");
    plot(par[:tol]*ones(len+1),color="black","--",linewidth=.5)
    grid("on",color="lightgray");
    
    xlabel(L"\textrm{iteration steps}");
    yscale("log")
    xlim([0,len]);
    legend();
    tight_layout()
    savefig("inf.pdf",bbox_inches="tight")
    close()
end
