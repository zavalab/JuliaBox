module SchwarzOCPExamples

using SimpleNLModels
import SimpleNLModels: variable, Model, Constraint
using SparseArrays
using LinearAlgebra
using Plots, LaTeXStrings; gr()

function variable(m::Model,k;opt...)
    push!(m[:Vs][k],num_variables(m)+1)
    variable(m;opt...)
end

function Model(optimizer,K;opt...)
    m = Model(optimizer;opt...)
    m[:Vs]= Dict(k=>Int[] for k in K)
    return m
end

function hehnandrea(
    optimizer,N,K; x0_val = [0,0,0,0,0,0,0,0,0], dt = .01,
    Q = [1,0,1,0,1,0,1,1,1], Qf= [1,0,1,0,1,0,1,1,1]/dt, R = ones(4)/10,
    d_func = (i,j)->(j==1 ? 1*sin(2*pi/N*i) : 0)+(j==3 ? 2*sin(4*pi/N*i) : 0)+(j==5 ? 2*i/N : 0),
    opt...)
    
    n = 9
    p = 4    

    m = SimpleNLModels.Model(optimizer,K;opt...)
    
    x=[variable(m,floor(Int,(i-1)/(N+1)*length(K))+1) for i=1:N+1,j=1:n]
    u=[variable(m,floor(Int,(i-1)/(N+1)*length(K))+1) for i=1:N,j=1:p]

    x0 = [parameter(m,x0_val[i]) for i=1:n]
    d = [parameter(m,d_func(i,j)) for i=1:N+1,j=1:n]
         
    c = Array{Constraint,2}(undef,N+1,n)
    
    for j=1:n
        c[1,j] = constraint(m,x[1,j]-x0[j])
    end


    for i=1:N
        c[i+1,1] = constraint(m, -x[i+1,1] + x[i,1] + (x[i,2])*dt)
        c[i+1,2] = constraint(m, -x[i+1,2] + x[i,2] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*cos(x[i,9])+u[i,1]*sin(x[i,7])*sin(x[i,9]))*dt)
        c[i+1,3] = constraint(m, -x[i+1,3] + x[i,3] + (x[i,4])*dt)
        c[i+1,4] = constraint(m, -x[i+1,4] + x[i,4] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*sin(x[i,9])-u[i,1]*sin(x[i,7])*cos(x[i,9]))*dt)
        c[i+1,5] = constraint(m, -x[i+1,5] + x[i,5] + (x[i,6])*dt)
        c[i+1,6] = constraint(m, -x[i+1,6] + x[i,6] + (u[i,1]*cos(x[i,7])*cos(x[i,8])-9.8)*dt)
        c[i+1,7] = constraint(m, -x[i+1,7] + x[i,7] + (u[i,2]*cos(x[i,7])/cos(x[i,8])+u[i,3]*sin(x[i,7])/cos(x[i,8]))*dt)
        c[i+1,8] = constraint(m, -x[i+1,8] + x[i,8] + (-u[i,2]*sin(x[i,7])+u[i,3]*cos(x[i,7]))*dt)
        c[i+1,9] = constraint(m, -x[i+1,9] + x[i,9] + (u[i,2]*cos(x[i,7])*tan(x[i,8])+u[i,3]*sin(x[i,7])*tan(x[i,8])+u[i,4])*dt)
    end
    
    for i=1:N
        for j=1:n
            objective(m,.5*Q[j]*(x[i,j]-d[i,j])^2)
        end
        for j=1:p
            objective(m,.5*R[j]*(u[i,j]^2))
        end
    end

    for j=1:n
        objective(m,.5*Qf[j]*(x[N+1,j]-d[N+1,j])^2)
    end

    m[:x] = x
    m[:u] = u
    m[:x0] = x0
    m[:d] = d
    m[:c] = c
    
    return m
end

function thinplate(
    optimizer,n,N,K;
    x0 = ones(n,n)*298, dt = .01, dx = 1.,
    Q = ones(n,n), Qf= ones(n,n)/dt, R = ones(n,n)/10,
    d = (i,j,k)->100*sin(2*pi*(4*i/N-12*j/n-12*k/n)) + 400,
    opt...)

    kappa = 400 # thermal conductivity of copper, W/(m-K)
    rho = 8960 # density of copper, kg/m^3
    specificHeat = 386 # specific heat of copper, J/(kg-K)
    thick = .01 # plate thickness in meters
    stefanBoltz = 5.670373e-8 # Stefan-Boltzmann constant, W/(m^2-K^4)
    hCoeff = 1 # Convection coefficient, W/(m^2-K)
    Ta = 300 # The ambient temperature is assumed to be 300 degrees-Kelvin.
    emiss = .5 # emissivity of the plate surface
    
    m = SimpleNLModels.Model(optimizer,K;opt...)
    
    x=Dict((i,j,k) => (j ==0 || j== n+1 || j== n+1 || j== n+1) ? .0 : variable(m,floor(Int,(i-1)/(N+1)*length(K))+1)
           for i=1:N+1,j=0:n+1,k=0:n+1)
    u=Dict((i,j,k) => (j ==0 || j== n+1) ? .0 : variable(m,floor(Int,(i-1)/(N+1)*length(K))+1)
           for i=1:N,j=0:n+1,k=0:n+1)
    
    for j=1:n
        for k=1:n
            constraint(m,x[1,j,k]-x0[j,k])
        end
    end
    
    for i=1:N
        for j=1:n
            for k=1:n
                constraint(m, x[i+1,j,k] ==  x[i,j,k] + (1/rho/specificHeat/thick)*(kappa*thick*(-4*x[i,j,k] - x[i,j,k-1] - x[i,j,k+1] - x[i,j-1,k] - x[i,j+1,k])/dx^2  -.25* u[i,j,k] + 2*hCoeff*(x[i,j,k]-Ta) + 2*emiss*stefanBoltz*(x[i,j,k]^4-Ta^4))*dt)
            end
        end
    end

    for k=1:n
        for j=1:n
            for i=1:N
                objective(m,.5*Q[j,k]*(x[i,j,k]-d(i,j,k))^2)
                objective(m,.5*R[j,k]*(u[i,j,k]^2))
            end
            objective(m,.5*Qf[j,k]*(x[N+1,j,k]-d(N+1,j,k))^2)
        end
    end
    
    m[:x] = x
    m[:u] = u
    
    return m
end

mutable struct KKTErrorEvaluator
    f
    c
    jac
    jac_view
    
    grad!
    con!
    jac!
end
    
function set_KKT_error_evaluator!(m::SimpleNLModels.Model)
    _grad! = SimpleNLModels.eval_grad!(m.objs,m.p)
    _con! = SimpleNLModels.eval_con!(m.cons,m.p)
    _jac!,_jac_sparsity!,nnz_jac = SimpleNLModels.eval_jac!(m.cons,m.p)
    
    I = Vector{Int}(undef,nnz_jac)
    J = Vector{Int}(undef,nnz_jac)
    _jac_sparsity!(I,J)

    f = similar(m.x)
    c = similar(m.l)
    
    tmp = sparse(J,I,collect(1:nnz_jac),num_variables(m),num_constraints(m))
    jac = SparseMatrixCSC{Float64,Int64}(tmp.m,tmp.n,tmp.colptr,tmp.rowval,Vector{Float64}(undef,nnz_jac))

    idx = Vector{Int}(undef,length(tmp.nzval))
    for i=eachindex(tmp.nzval)
        idx[tmp.nzval[i]] = i
    end
    jac_view = view(jac.nzval,idx)
    
    m[:KKTErrorEvaluator] = KKTErrorEvaluator(f,c,jac,jac_view,_grad!,_con!,_jac!)
end

function get_KKT_error(m::SimpleNLModels.Model)
    evaluator = m[:KKTErrorEvaluator]

    evaluator.grad!(evaluator.f,m.x)
    evaluator.con!(evaluator.c,m.x)
    evaluator.jac!(evaluator.jac_view,m.x)

    mul!(evaluator.f,evaluator.jac,m.l,1.,1.)
    evaluator.c .-= m.gl

    return max(norm(evaluator.f,Inf),norm(evaluator.c,Inf))
end


function plot_benchmark(err_ref,time_ref,err_schwarz,time_schwarz,err_admm,time_admm)
    plt = plot(;xlabel="Wall Time (sec)",ylabel="KKT Error",framestyle=:box,yscale=:log10,xlim=(0,min(time_schwarz[end],time_admm[end])),legend=:right,fontfamily="Computer Modern");
    plot!(plt,time_schwarz,err_schwarz,label="Schwarz",markershape=:diamond);
    plot!(plt,time_admm,err_admm,label="ADMM",markershape=:circle);
    vline!(plt,[time_ref],linestyle=:dash,color=:black,label="Ipopt Time")
    hline!(plt,[err_ref],linestyle=:dot,color=:black,label="Ipopt KKT error")
end


function plot_eds(x_ref,y_ref,z_ref,x_pers,y_pers,z_pers,
                  lx_ref,ly_ref,lz_ref,lx_pers,ly_pers,lz_pers)
    plt1 = plot(;xlabel=L"$x$",ylabel=L"$y$",zlabel=L"$z$",fontfamily="Computer Modern");
    for i in eachindex(x_pers)
        plot!(plt1,x_pers[i],y_pers[i],z_pers[i],color=:lightblue,
              label= i==1 ? "Perturbed Solutions" : :none)
        plot!(plt1,(x_pers[i][1],y_pers[i][1],z_pers[i][1]),
              marker=:diamond,label=:none,color=:lightblue)
        plot!(plt1,(x_pers[i][end],y_pers[i][end],z_pers[i][end]),
              marker=:circle,label=:none,color=:lightblue)
    end
    plot!(plt1,x_ref,y_ref,z_ref,label="Reference solution",color=:black)
    plot!(plt1,(x_ref[1],y_ref[1],z_ref[1]),marker=:diamond,label=:none,color=:black)
    plot!(plt1,(x_ref[end],y_ref[end],z_ref[end]),marker=:circle,label=:none,color=:black)

    plt2 = plot(;xlabel=L"$\lambda_x$",ylabel=L"$\lambda_y$",zlabel=L"$\lambda_z$",
                fontfamily="Computer Modern");
    for i in eachindex(lx_pers)
        plot!(plt2,lx_pers[i],ly_pers[i],lz_pers[i],color=:lightblue,
              label= i==1 ? "Perturbed Solutions" : :none)
        plot!(plt2,(lx_pers[i][1],ly_pers[i][1],lz_pers[i][1]),
              marker=:diamond,label=:none,color=:lightblue)
        plot!(plt2,(lx_pers[i][end],ly_pers[i][end],lz_pers[i][end]),
              marker=:circle,label=:none,color=:lightblue)
    end
    plot!(plt2,lx_ref,ly_ref,lz_ref,label="Reference solution",color=:black)
    plot!(plt2,(lx_ref[1],ly_ref[1],lz_ref[1]),marker=:diamond,label=:none,color=:black)
    plot!(plt2,(lx_ref[end],ly_ref[end],lz_ref[end]),marker=:circle,label=:none,color=:black)
    return plt1,plt2
end

function plot_demonstration(coords,x_ref,y_ref,z_ref)

    markers = [:diamond,:circle,:square]
    plts = []
    for i in eachindex(coords[1])
        plt = plot(;xlabel=L"$x$",ylabel=L"$y$",zlabel=L"$z$",fontfamily="Computer Modern");
        for k in eachindex(coords)
            (x,y,z) = coords[k][i]
            scatter!(plt,x[1:25:end],y[1:25:end],z[1:25:end],label="Subproblem #$k",marker=markers[k])
        end
        plot!(plt,x_ref,y_ref,z_ref,color=:black,label="Full solution")
        push!(plts,plt)
    end
    return plts
end

function plot_err_profile(err_schwarz,omegas)
    k = 0
    markers = [:diamond,:circle,:square]
    plt = plot(;xlabel="Iteration Steps",ylabel="KKT Error",framestyle=:box,
               yscale=:log10,fontfamily="Computer Modern",
               xlim=(0,maximum(length(err) for err in err_schwarz)),
               ylim=(maximum(minimum(err) for err in err_schwarz),
                     maximum(maximum(err) for err in err_schwarz)*2));
    for i in eachindex(err_schwarz)
        plot!(0:length(err_schwarz[i])-1,err_schwarz[i],label=latexstring("\$\\widetilde{\\tau}= $(omegas[i]) \$"),marker=markers[i])
    end
    return plt
end

export hehnandrea, thinplate, set_KKT_error_evaluator!, get_KKT_error, plot_benchmark, plot_eds, savefig, plot_demonstration, plot_err_profile

end # module
