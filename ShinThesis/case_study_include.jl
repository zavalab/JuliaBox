using LightGraphs,Plots,PowerModels,Metis,LaTeXStrings,LinearAlgebra,DelimitedFiles,Random; PowerModels.silence()

function add_branch!(g,G,B,i,j,gg,bb)
    add_edge!(g,i,j)
    G[i,j]=G[j,i]=gg
    B[i,j]=B[j,i]=bb
end

function get_network(str,K)
    dat = build_ref(PowerModels.parse_file(str))[:it][:pm][:nw][0]

    g = Graph(length(dat[:bus]))
    G = Dict{Tuple{Int,Int},Float64}()
    B = Dict{Tuple{Int,Int},Float64}()

    cnt = 0
    buslist = Int[]
    buslistinv = Dict{Int,Int}()
    mapper = Dict(1=>1,2=>2,3=>3,9=>4)
    
    for (i,bus) in dat[:bus]
        cnt += 1
        push!(buslist,i)
        buslistinv[i] = cnt
    end
    for (i,br) in dat[:branch]
        gg,bb = calc_branch_y(br)
        add_branch!(g,G,B,buslistinv[br["f_bus"]],buslistinv[br["t_bus"]],gg,bb)
    end
    
    return g,G,B,Metis.partition(g,K)
end

function acopf(optimizer,N,K; eta = 1e6, b= 1, opt...)

    Random.seed!(0)
    g,G,B,part = get_network(N == 0 ? "pglib_opf_case30000_goc.m" : "pglib_opf_case500_tamu.m",K)
    
    m = SimpleNL.Model(optimizer,1:K;opt...)
    va = Dict(i=>variable(m,part[i];start=0) for i in vertices(g))
    vm = Dict(i=>variable(m,part[i];start=1) for i in vertices(g))
    pg = Dict(i=>variable(m,part[i];start=0) for i in vertices(g))
    qg = Dict(i=>variable(m,part[i];start=0) for i in vertices(g))

    pd = Dict(i=>parameter(m,randn()) for i in vertices(g))
    qd = Dict(i=>parameter(m,.25*randn()) for i in vertices(g))
    vas = Dict(i=>parameter(m,0.1*randn()) for i in vertices(g))
    vms = Dict(i=>parameter(m,1+0.01*randn()) for i in vertices(g))
    ps = Dict(i=>parameter(m,0.1*randn()) for i in vertices(g))
    qs = Dict(i=>parameter(m,1+0.01*randn()) for i in vertices(g))
    
    bp = [constraint(m,sum( G[i,j]*vm[i]^2 + (-G[i,j])*(vm[i]*vm[j]*cos(va[i]-va[j])) + (-B[i,j])*(vm[i]*vm[j]*sin(va[i]-va[j])) for j in neighbors(g,i))== b*pg[i] + pd[i]) for i in vertices(g)]
    bq = [constraint(m,sum(-B[i,j]*vm[i]^2 - (-B[i,j])*(vm[i]*vm[j]*cos(va[i]-va[j])) + (-G[i,j])*(vm[i]*vm[j]*sin(va[i]-va[j])) for j in neighbors(g,i))== b*qg[i] + qd[i]) for i in vertices(g)]
    
    objective(m,sum((pg[i]-ps[i])^2 + (qg[i]-qs[i])^2 + eta*(vm[i]-vms[i])^2+eta*(va[i]-vas[i])^2  for i in vertices(g)))

    m[:vas] = vas
    m[:vms] = vms
    m[:va] = va
    m[:vm] = vm
    m[:pg] = pg
    m[:qg] = qg
    
    m[:ps] = ps
    m[:qs] = qs
    m[:ps] = ps
    m[:qs] = qs
    m[:pd] = pd
    m[:qd] = qd

    m[:bp] = bp
    m[:bq] = bq

    m[:g] = g
    m[:pos] = N == 0 ? nothing : readdlm("pos.csv",',')[:,2:3]
    
    return m
end

function storage(
    optimizer,N,K;
    NS = 3, x0 = 0, eta = 1., b = 1.,
    opt...)

    g = Graph()
    add_vertex!(g)
    dd  = [0.]
    pr = [1.]
    ds = [-1. 0. 1.]
    part = [1]
    part_cnt = 1
    changing = [K]
    pos= [[.0,.0]]
    
    for t in 1:N-1
        ps = nv(g).-(NS^(t-1):-1:1).+1
        for p in ps
            for j=1:NS
                add_vertex!(g)
                add_edge!(g,p,nv(g))
                push!(dd,ds[j])
                push!(pr,1/NS^t)
                push!(part,t+1 in changing ? part_cnt+=1 : part[p])
                push!(pos,[t,pos[p][2]+NS^(-t+.0)*((NS+1)/2-j)])
            end
        end
    end
    
    pos = hcat(pos...)'
    
    m = SimpleNL.Model(optimizer,1:maximum(part);opt...)
    x = [variable(m,part[i]) for i in 1:nv(g)]
    u = [variable(m,part[i]) for i in 1:nv(g)]

    d = [parameter(m,dd[j]) for j=1:nv(g)]

    objective(m, sum(pr[i]*(1/2*eta*x[i]^2+1/2*u[i]^2) for i=1:nv(g)))
    l = Array{Any,1}(nothing,nv(g))
    l[1] = constraint(m, x[1] == x0)
    for i=2:nv(g)
        l[i] = constraint(m, x[i] == x[inneighbors(g,i)[1]] + b*u[inneighbors(g,i)[1]]  + d[i])
    end

    m[:x] = x
    m[:u] = u
    m[:d] = d
    m[:g] = g
    m[:l] = l
    m[:pos] = pos
    
    return m
end

function thinplate(
    optimizer,N,K;
    x0 = ones(N,N)*298, dx = .1, eta =1, b= 1,
    d_func = (i,j)->sin(2*pi*(12*i/N-12*j/N)),
    opt...)
    Q = ones(N,N)*eta/dx^2
    R = ones(N,N)/10/dx^2

    kappa = 400. # thermal conductivity of copper, W/(m-K)
    rho = 8960. # density of copper, kg/m^3
    specificHeat = 386. # specific heat of copper, J/(kg-K)
    thick = .01 # plate thickness in meters
    stefanBoltz = 5.670373e-8 # Stefan-Boltzmann constant, W/(m^2-K^4)
    hCoeff = 1. # Convection coefficient, W/(m^2-K)
    Ta = 300. # The ambient temperature is assumed to be 300 degrees-Kelvin.
    emiss = .5 # emissivity of the plate surface

    k(i,j) = floor(Int,(i-1)/(N+1)*K)*K+floor(Int,(j-1)/(N+1)*K)+1
    
    m = SimpleNL.Model(optimizer,1:K^2;opt...)
    x=Dict((i,j) => (i ==0 || i== N+1 || j== 0 || j== N+1) ? Ta : variable(m,k(i,j)) for i=0:N+1,j=0:N+1)
    u=[variable(m,k(i,j)) for i=1:N,j=1:N]
    d=[parameter(m,d_func(i,j)) for i=1:N,j=1:N]
    l=[constraint(m, kappa*thick*(4*x[i,j] - x[i,j-1] - x[i,j+1] - x[i-1,j] - x[i+1,j])/dx^2  -.25* u[i,j]*b + 2*hCoeff*(x[i,j]-Ta) + 2*emiss*stefanBoltz*(x[i,j]^4-Ta^4) + d[i,j]) for i=1:N, j=1:N]

    objective(m,sum(sum(.5*Q[i,j]*(x[i,j]-Ta)^2 + .5*R[i,j]*(u[i,j]^2) for j=1:N) for i=1:N))

    g = Graph()
    pos = []
    for i in 1:N
        for j in 1:N
            add_vertex!(g)
            add_edge!(g,N*(i-1)+j,N*(i-2)+j)
            j!= 1 && add_edge!(g,N*(i-1)+j,N*(i-1)+j-1)
            push!(pos,[i,j])
        end
    end
    pos = hcat(pos...)'
    
    m[:x] = [x[i,j] for i=1:N, j=1:N]
    m[:u] = u
    m[:d] = d
    m[:l] = l
    m[:g] = g
    m[:pos] = pos
    
    return m
end

function hehnandrea(
    optimizer,N,K; x0_val = [0,0,0,0,0,0,0,0,0], dt = .005, eta = 1,b = 1,
    d_func = (i,j)->(j==1 ? 1*sin(2*pi/N*i) : 0)+(j==3 ? 2*sin(4*pi/N*i) : 0)+(j==5 ? 2*i/N : 0),
    opt...)
    
    Q = [1,0,1,0,1,0,1,1,1]*eta
    Qf= [1,0,1,0,1,0,1,1,1]/dt*eta
    R = ones(4)/10
    
    n = 9
    p = 4    
    k(i,j) = floor(Int,(i-1)/(N+1)*K)+1
    
    m = SimpleNL.Model(optimizer,1:K;opt...)
    x=[variable(m,k(i,j)) for i=1:N+1,j=1:n]
    u=[variable(m,k(i,j)) for i=1:N,j=1:p]

    x0 = [parameter(m,x0_val[i]) for i=1:n]
    d = [parameter(m,d_func(i,j)) for i=1:N+1,j=1:n]
         
    c = Array{SimpleNL.Constraint,2}(undef,N+1,n)
    
    for j=1:n
        c[1,j] = constraint(m,x[1,j] == x0[j])
    end

    for i=1:N
        c[i+1,1] = constraint(m, -x[i+1,1] + x[i,1] + (x[i,2])*dt)
        c[i+1,2] = constraint(m, -x[i+1,2] + x[i,2] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*cos(x[i,9])+u[i,1]*sin(x[i,7])*sin(x[i,9]))*dt)
        c[i+1,3] = constraint(m, -x[i+1,3] + x[i,3] + (x[i,4])*dt)
        c[i+1,4] = constraint(m, -x[i+1,4] + x[i,4] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*sin(x[i,9])-u[i,1]*sin(x[i,7])*cos(x[i,9]))*dt)
        c[i+1,5] = constraint(m, -x[i+1,5] + x[i,5] + (x[i,6])*dt)
        c[i+1,6] = constraint(m, -x[i+1,6] + x[i,6] + (u[i,1]*cos(x[i,7])*cos(x[i,8])-9.8)*dt)
        c[i+1,7] = constraint(m, -x[i+1,7] + x[i,7] + (b*u[i,2]*cos(x[i,7])/cos(x[i,8])+u[i,3]*sin(x[i,7])/cos(x[i,8]))*dt)
        c[i+1,8] = constraint(m, -x[i+1,8] + x[i,8] + (-b*u[i,2]*sin(x[i,7])+u[i,3]*cos(x[i,7]))*dt)
        c[i+1,9] = constraint(m, -x[i+1,9] + x[i,9] + (b*u[i,2]*cos(x[i,7])*tan(x[i,8])+u[i,3]*sin(x[i,7])*tan(x[i,8])+u[i,4])*dt)
    end
    
    objective(m,sum(.5*Q[j]*(x[i,j]-d[i,j])^2 for i=1:N for j=1:n) + sum(.5*R[j]*(u[i,j]^2) for i=1:N for j=1:p)+ sum(.5*Qf[j]*(x[N+1,j]-d[N+1,j])^2 for j=1:n))

    g = Graph(N)
    for i=1:N
        add_edge!(g,i,i+1)
    end
    
    m[:x] = x
    m[:u] = u
    m[:x0] = x0
    m[:d] = d
    m[:c] = c
    m[:g] = g
    m[:pos] = hcat([[i,0] for i=1:N]...)'
    
    return m
end

const markers = [:diamond,:circle,:square,:utriangle]

function plot_err_profile(err_input,str,lbls)
    k = 0
    errs = [get_errs_times(err)[1] for err in err_input]
    times = [get_errs_times(err)[2] for err in err_input]
    
    plt1 = plot(;xlabel="Iteration Steps",ylabel="KKT Error",framestyle=:box,
               yscale=:log10,fontfamily="Computer Modern",
               xlim=(0,maximum(length(err) for err in errs)),
               ylim=(minimum(minimum(err) for err in errs),
                     maximum(maximum(err)*1e3 for err in errs)*2));
    for i in eachindex(errs)
        plot!(0:length(errs[i])-1,errs[i],label=latexstring("\$ $(str)= $(lbls[i]) \$"),marker=markers[i])
    end

    plt2 = plot(;xlabel="Wall Time (sec)",ylabel="KKT Error",framestyle=:box,
                yscale=:log10,fontfamily="Computer Modern",
                xlim=(0,maximum(maximum(time) for time in times)),
                ylim=(minimum(minimum(err) for err in errs),
                      maximum(maximum(err)*1e3 for err in errs)*2));
    for i in eachindex(errs)
        plot!(times[i],errs[i],label=latexstring("\$ $(str)= $(lbls[i]) \$"),marker=markers[i])
    end

    return plt1,plt2
end

function plot_benchmark(err_ref,time_ref,output_schwarz,output_admm)
    err_schwarz,time_schwarz = get_errs_times(output_schwarz)
    err_admm,time_admm = get_errs_times(output_admm)
    
    plt1 = plot(;xlabel="Iteration Steps",ylabel="KKT Error",framestyle=:box,yscale=:log10,xlim=(0,max(length(err_schwarz),length(err_admm))),legend=:topright,fontfamily="Computer Modern");
    plot!(plt1,0:length(err_schwarz)-1,err_schwarz,label="Schwarz",markershape=:diamond);
    plot!(plt1,0:length(err_admm)-1,err_admm,label="ADMM",markershape=:circle);
    # hline!(plt1,[err_ref],linestyle=:dot,color=:black,label="Ipopt KKT error")

    plt2 = plot(;xlabel="Wall Time (sec)",ylabel="KKT Error",framestyle=:box,yscale=:log10,xlim=(0,1.1*max(time_ref,time_schwarz[end],time_admm[end])),legend=:topright,fontfamily="Computer Modern");
    plot!(plt2,time_schwarz,err_schwarz,label="Schwarz",markershape=:diamond);
    plot!(plt2,time_admm,err_admm,label="ADMM",markershape=:circle);
    vline!(plt2,[time_ref],linestyle=:dash,color=:black,label="Ipopt Time")
    # hline!(plt2,[err_ref],linestyle=:dot,color=:black,label="Ipopt KKT error")

    return plt1,plt2
end

function get_errs_times(errtimes)
    errs = Vector{Float64}(undef,length(errtimes))
    times= Vector{Float64}(undef,length(errtimes))

    for i in eachindex(errtimes)
        errs[i],times[i] = errtimes[i]
    end
    
    return errs,times
end

function coloring(x;fcolor = [0,0,1],bcolor = [.95,.95,.95])
    x=min(1,x)
    arr = (1-x)*bcolor+x*fcolor
    return RGB(arr[1],arr[2],arr[3])
end

function plot_graph(g,pos,vals,j)
    plt = plot(;grid=false,legend=false,axis=([], false),size=(500,500))
    plot!(plt,[0,0],[-1e-6,1e-6],linealpha=0)
    for e in edges(g)
        plot!(plt,[pos[e.src,1],pos[e.dst,1]],[pos[e.src,2],pos[e.dst,2]],linecolor = coloring(0.5*(vals[e.src]+vals[e.dst])))
    end
    scatter!(plt,pos[:,1],pos[:,2],markersize=4,markercolor=coloring.(vals),markerstrokewidth=0)
    scatter!(plt,[pos[j,1]],[pos[j,2]],markersize=8,markerstrokewidth=2,markerstrokecolor=:red,markercolor=:transparent)
end

function compute_vals_storage(m_ref,m_per)
    [norm([value(m_ref[:x][i])-value(m_per[:x][i]),dual(m_ref[:l][i])-dual(m_per[:l][i])]) for i in 1:nv(m_ref[:g])]
end
function compute_vals_thinplate(m_ref,m_per)
    [norm([
        value(m_ref[:x][i])-value(m_per[:x][i]),
        value(m_ref[:u][i])-value(m_per[:u][i]),
        dual(m_ref[:l][i])-dual(m_per[:l][i])
    ]) for i in 1:nv(m_ref[:g])]
end
function compute_vals_acopf(m_ref,m_per)
    [norm([
        value(m_ref[:va][i])-value(m_per[:va][i]),
        value(m_ref[:vm][i])-value(m_per[:vm][i]),
        value(m_ref[:pg][i])-value(m_per[:pg][i]),
        value(m_ref[:qg][i])-value(m_per[:qg][i]),
        dual(m_ref[:bp][i])-dual(m_per[:bp][i]),
        dual(m_ref[:bq][i])-dual(m_per[:bq][i]),
    ]) for i in 1:nv(m_ref[:g])]
end

function perturb_storage(m_ref,m_per,j,sig)
    setvalue(m_per[:d][j],value(m_ref[:d][j]) + sig*randn())
end
function perturb_thinplate(m_ref,m_per,j,sig)
    setvalue(m_per[:d][j],value(m_ref[:d][j]) + sig*randn())
end
function perturb_acopf(m_ref,m_per,j,sig)
    setvalue(m_per[:vas][j],value(m_ref[:vas][j]) + sig*randn())
    setvalue(m_per[:vms][j],value(m_ref[:vms][j]) + sig*randn())
    setvalue(m_per[:pd][j],value(m_ref[:pd][j]) + sig*randn())
    setvalue(m_per[:qd][j],value(m_ref[:qd][j]) + sig*randn())
    setvalue(m_per[:ps][j],value(m_ref[:ps][j]) + sig*randn())
    setvalue(m_per[:qs][j],value(m_ref[:qs][j]) + sig*randn())
end

function compute_vals_hehnandrea(m_ref,m_per)
    val = [norm([
        value.(m_ref[:x][i,:])-value.(m_per[:x][i,:]);
        value(m_ref[:u][i,1])-value(m_per[:u][i,1]);
        dual(m_ref[:c][i,1])-dual(m_per[:c][i,1])
    ]) for i in 1:nv(m_ref[:g])]
    return val
end
function perturb_hehnandrea(m_ref,m_per,j,sig)
    setvalue.(m_per[:d][j,:],value.(m_ref[:d][j,:]) .+ sig*randn(size(m_ref[:d],2)))
end


function plot_eds(dist,outputs,str,lbls)
    plt = plot(;xlabel=L"d_G(i,j)",ylabel=L"\overline{C}_{ij}",framestyle=:box,yscale=:log10,fontfamily="Computer Modern");
    for i in eachindex(outputs)
        scatter!(plt,dist,outputs[i],label=latexstring("\$ $(str)= $(lbls[i]) \$"),marker=markers[i])
    end
    return plt
end

