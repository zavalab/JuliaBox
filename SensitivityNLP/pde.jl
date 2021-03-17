using Random, JuMP, Ipopt, Plots, LightGraphs, LinearAlgebra

Random.seed!(1)

n_grad = 4
fcolor = [0,0,1]
bcolor = [.95,.95,.95]
function coloring(x)
    x=min(1,x)
    arr = (1-x)*bcolor+x*fcolor
    return RGB(arr[1],arr[2],arr[3])
end
norms(A;p=2) = sqrt.(abs.(sum(A.^2)))
params = [(1,1), (1e-2,1), (1,1e-2), (1e-2,1e-2)]
(i0,j0) = (5,5)

sig = .001
n_sample = 30
NX = 9
NY = 9

k=400;
tz=.01;
hc=1;
eps=.5;
sb=5.680373e-8
Ta=300

grid = hcat([[(i,j) for i=1:NX] for j=1:NY]...)
nbr = hcat([[intersect([(i-1,j),(i+1,j),(i,j-1),(i,j+1)],grid)  for i=1:NX] for j=1:NY]...)

for cnt = 1:length(params)
    global sols
    (eta,b)=params[cnt]

    m=Model(Ipopt.Optimizer)
    @variable(m,x[1:NX,1:NY])
    @variable(m,u[1:NX,1:NY])
    @NLparameter(m,d[1:NX,1:NY]==0)
    @NLparameter(m,xref[1:NX,1:NY]==300)
    
    l=@NLconstraint(m,[i=1:NX,j=1:NY],
                    sum(x[i,j]-x[i1,j1] for (i1,j1) in nbr[i,j])==
                    -2hc/k/tz*(x[i,j]-300)
                    -2eps*sb/k/tz*(x[i,j]^4-300^4)
                    +(b*u[i,j] + d[i,j])/k/tz)
    @NLobjective(m,Min,sum(eta*(x[i,j]-xref[i,j])^2+u[i,j]^2 for i=1:NX for j=1:NY))
    optimize!(m)
    set_start_value.(all_variables(m),value.(all_variables(m)))
    sol_ref = [value.(x),value.(u),dual.(l)]

    sols = []
    ps   = []


    for i=1:n_sample
        p = sig*(rand(3).-.5)
        set_value(d[i0,j0],p[1])
        set_value(xref[i0,j0],300+p[2])
        optimize!(m)
        push!(sols,[value.(x),value.(u),dual.(l)])
        push!(ps,  p)
    end
    C = maximum(cat([norms(sols[i].-sol_ref)/norm(ps[i]) for i=1:n_sample]...,dims=3),dims=3)[:,:,1]
    C = C./C[i0,j0]

    pyplot()
    p=scatter(leg=false,ticks=nothing,framestyle=:none,size=(500,500));
    for i in 1:NX
        for j in 1:NY
            if i!=NX
                plot!(p,
                      range(i;
                            stop = i+1,
                            length=n_grad),
                      range(j;
                            stop = j,
                            length=n_grad),
                      lc= cgrad([coloring(C[i,j]),
                                 coloring(C[i+1,j])]),
                      line_z =1:n_grad);
            end
            if j!=NY
                plot!(p,
                      range(i;
                            stop = i,
                            length=n_grad),
                      range(j;
                            stop = j+1,
                            length=n_grad),
                      lc= cgrad([coloring(C[i,j]),
                                 coloring(C[i,j+1])]),
                      line_z =1:n_grad);
            end
        end
    end
    for i=1:NX
        for j=1:NY
            scatter!(p,[i],[j],markersize=8,markercolor=coloring(C[i,j]),markerstrokecolor=coloring(C[i,j]),markerstrokewidth=0);
        end
    end
    scatter!(p,[i0],[j0],markersize=18,markercolor=:transparent,markerstrokewidth=2,markerstrokecolor=:red);
    savefig(p, "fig/pde-hm-$cnt.pdf")

    pgfplotsx()
    q=scatter([abs(i-i0)+abs(j-j0) for j=1:NY for i=1:NX][C[:].>1e-10],C[:][C[:].>1e-10],
              size=(300,300), leg=false, framestyle=:box,
              yscale=:log10,
              ylims=(1e-3,10),
              markerstrokewidth=0,color=:black)
    savefig(q, "fig/pde-sc-$cnt.pdf")    
end
