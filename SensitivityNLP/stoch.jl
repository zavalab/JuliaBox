using Random, JuMP, Ipopt, Plots, LightGraphs, LinearAlgebra

fcolor = [0,0,1]
bcolor = [.95,.95,.95]

function coloring(x)
    x=min(1,x)
    arr = (1-x)*bcolor+x*fcolor
    return RGB(arr[1],arr[2],arr[3])
end

Random.seed!(1)

norms(A;p=2) = [norm(A[i,:],p) for i=1:size(A,1)]
sig = .001
n_sample = 30
n_grad = 8

T = 6;  # Number of stages
NS= 3; # Number of scenarios of loads at each stage
j = 3;
params = [(1,1), (1e-2,1), (1,1e-2), (1e-2,1e-2)]

g = Graph()
add_vertex!(g)
d  = [0.]
pr = [1.]
ds = [-1. 0. 1.]
pos= [[.0,.0]]
for t in 1:T-1
    ps = nv(g).-(NS^(t-1):-1:1).+1
    for p in ps
        for j=1:NS
            add_vertex!(g)
            add_edge!(g,p,nv(g))
            push!(d,ds[j])
            push!(pr,1/NS^t)
            push!(pos,[t,pos[p][2]+NS^(-t+.0)*((NS+1)/2-j)])
        end
    end
end
pos = hcat(pos...)'



for cnt = 1:length(params)
    eta,b = params[cnt]
    m = Model(Ipopt.Optimizer)
    @variable(m, x[1:nv(g)])
    @variable(m, u[1:nv(g)])
    @variable(m, v[1:nv(g)])

    @NLparameter(m,ee[1:nv(g)]==0)
    @NLparameter(m,dd[j=1:nv(g)]==d[j])
    @NLparameter(m,pi[1:nv(g)]==0)

    @NLobjective(m, Min, sum(pr[i]*(1/2*eta*x[i]^2+1/2*u[i]^2+pi[i]*v[i]) for i=1:nv(g)))
    l = Array{Any,1}(nothing,nv(g))
    l[1]     = @NLconstraint(m, x[1] == ee[1])
    for i=2:nv(g)
        l[i] = @NLconstraint(m, x[i] == x[inneighbors(g,i)[1]] + b*u[inneighbors(g,i)[1]] + ee[i])
    end
    mu= Array{Any,1}(nothing,nv(g))
    for i=1:nv(g)
        mu[i] = @NLconstraint(m, v[i] == u[i] + dd[i])
    end

    optimize!(m)
    set_start_value.(all_variables(m),value.(all_variables(m)))
    sol_ref = [value.(x) value.(u) value.(v) dual.(l) dual.(mu)]

    sols = []
    ps   = []
    for i=1:n_sample
        p = sig*(rand(3).-.5)
        set_value(pi[j],p[1])
        set_value(dd[j],d[j]+p[2])
        set_value(ee[j],p[3])
        
        optimize!(m)
        Int(termination_status(m)) in [1,4]||error("suboptimal") 
        set_start_value.(all_variables(m),value.(all_variables(m)))
        push!(sols,[value.(x) value.(u) value.(v) dual.(l) dual.(mu)])
        push!(ps,  p)
    end
    C = maximum(hcat([norms(sols[i].-sol_ref)/norm(ps[i]) for i=1:n_sample]...),dims=2)
    C = C./C[j]

    pyplot()
    p=scatter(leg=false,ticks=nothing,border=:none,size=(500,500));
    for e in edges(g)
        plot!(p,
              range(pos[src(e),1];
                    stop = pos[dst(e),1],
                    length=n_grad),
              range(pos[src(e),2];
                    stop = pos[dst(e),2],
                    length=n_grad),
              lc= cgrad([coloring(C[src(e)]),
                         coloring(C[dst(e)])]),
              line_z =1:n_grad);
    end
    for i=1:size(pos,1)
        scatter!(p,[pos[i,1]],[pos[i,2]],markersize=8,markercolor=coloring(C[i]),markerstrokecolor=coloring(C[i]),markerstrokewidth=0);
    end
    scatter!(p,[pos[j,1]],[pos[j,2]],markersize=18,markercolor=:transparent,markerstrokewidth=2,markerstrokecolor=:red);
    savefig(p, "fig/stoch-hm-$cnt.pdf")

    pgfplotsx()
    q=scatter([length(a_star(g,i,j)) for i=1:nv(g)][C[:].>1e-10],C[:][C[:].>1e-10],
              size=(300,300), leg=false, framestyle=:box,
              yscale=:log10,
              ylims=(1e-3,10),
              markerstrokewidth=0,color=:black)
    savefig(q, "fig/stoch-sc-$cnt.pdf")
end
