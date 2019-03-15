using JLD, LightGraphs, JuMP, Ipopt, MPI

function get_net(data)
    mva = 100

    # bus
    n   = size(data["bus"],1)
    bus = 1:n
    bustype = Int64.(data["bus"][:,2])
    rbus = find(bustype.==3)
    Pd  = data["bus"][:,3]/mva
    Qd  = data["bus"][:,4]/mva
    Gs  = data["bus"][:,5]/mva
    Bs  = data["bus"][:,5]/mva
    zone  = data["bus"][:,10]
    Vmax  = data["bus"][:,12]
    Vmin  = data["bus"][:,13]
    map = Dict(Int64.(data["bus"][i,1]) => i for i=1:n)
    mapper(i) = map[i]

    # gen:
    genbus= mapper.(Int64.(data["gen"][:,1]))
    Qmax  = data["gen"][:,4]/mva
    Qmin  = data["gen"][:,5]/mva
    Pmax  = data["gen"][:,9]/mva
    Pmin  = data["gen"][:,10]/mva
    ngen = length(genbus)

    # gencost:
    c = data["gencost"][:,end:-1:5].*[1 mva mva^2]

    # branch:
    fbus = mapper.(Int64.(data["branch"][:,1]))
    tbus = mapper.(Int64.(data["branch"][:,2]))
    l = length(fbus)
    r = data["branch"][:,3]
    x = data["branch"][:,4]
    b = data["branch"][:,5]
    tau  = data["branch"][:,9]
    shift= deg2rad.(data["branch"][:,10])
    angmin = deg2rad.(data["branch"][:,12])
    angmax = deg2rad.(data["branch"][:,13])

    genmap= [[] for i=1:n]

    for i=1:ngen
        push!(genmap[genbus[i]],i)
    end



    for i=1:l
        if tau[i]==.0
            tau[i]=1.0
        end
    end


    g = Graph(n)
    for i=1:l
        add_edge!(g,fbus[i],tbus[i])
    end
    ys = 1./(r + x * im)

    Yff = (ys + im*b/2)./(tau).^2
    Yft = -ys.*(1./tau./exp.(-im*shift))
    Ytf = -ys.*(1./tau./exp.(im*shift))
    Ytt = (ys + im*b/2)

    Cf = spzeros(l,n)
    Ct = spzeros(l,n)
    for i=1:l
        j=fbus[i];k=tbus[i]
        Cf[i,j]=1
        Ct[i,k]=1
    end

    Ysh = sparse(1:n,1:n,Gs/mva + Bs*im/mva)
    Y = Cf'*spdiagm(Yff)*Cf + Cf'*spdiagm(Yft)*Ct + Ct'*spdiagm(Ytf)*Cf + Ct'*spdiagm(Ytt)*Ct + Ysh

    net = Dict()
    net[:rbus] = rbus
    net[:G] = real(Y)
    net[:B] = imag(Y)
    net[:n] = n
    net[:g] = g
    net[:ngen] = ngen
    net[:Vmax]= Vmax
    net[:Vmin]= Vmin
    net[:Pmax]= Pmax
    net[:Pmin]= Pmin
    net[:Qmax]= Qmax
    net[:Qmin]= Qmin
    net[:Pd]= Pd
    net[:Qd]= Qd
    net[:angmin] = angmin
    net[:angmax] = angmax
    net[:c] = c / par[:scale]
    net[:genmap] = genmap
    
    return net
end

function ac_opf(net,W,W1,G,E,par;
                starter=Dict(:T=>[0 for i=1:net[:n]],
                             :V=>[1 for i=1:net[:n]],
                             :P=>[(net[:Pmin][i]+net[:Pmax][i])/2 for i=1:net[:ngen]],
                             :Q=>[(net[:Qmin][i]+net[:Qmax][i])/2 for i=1:net[:ngen]],
                             :Ps=>[0 for i=1:net[:n]],
                             :Qs=>[0 for i=1:net[:n]]))
    K = keys(W)
    m = Model(solver=IpoptSolver(linear_solver="ma57",print_level=0))
    @variable(m, net[:Pmin][i] <= P[k in K,i in G[k]] <=net[:Pmax][i],start = starter[:P][i])
    @variable(m, net[:Qmin][i] <= Q[k in K,i in G[k]] <=net[:Qmax][i],start = starter[:Q][i])
    # @variable(m, Q[k in K,i in G[k]],start = starter[:Q][i])
    @variable(m, T[k in K,i in W1[k]], start = starter[:T][i])
    @variable(m, net[:Vmin][i]<=V[k in K,i in W1[k]]<=net[:Vmax][i],
              start = starter[:V][i])

    @variable(m, Ps[k in K,i in W[k]], start = starter[:Ps][i])
    @variable(m, Qs[k in K,i in W[k]], start = starter[:Qs][i])

    @variable(m, Ps_abs[k in K,i in W[k]], start = abs(starter[:Ps][i]))
    @variable(m, Qs_abs[k in K,i in W[k]], start = abs(starter[:Qs][i]))
    
    @constraint(m, [k in K,i in W[k]], Ps_abs[k,i]>=Ps[k,i])
    @constraint(m, [k in K,i in W[k]], Ps_abs[k,i]>=-Ps[k,i])
    @constraint(m, [k in K,i in W[k]], Qs_abs[k,i]>=Qs[k,i])
    @constraint(m, [k in K,i in W[k]], Qs_abs[k,i]>=-Qs[k,i])
    
    
    @constraint(m, [k in K,i in intersect(W[k],net[:rbus])], T[k,i] == 0)
    @constraint(m, [k in K,e in E[k]], -pi/6<=T[k,src(e)]-T[k,dst(e)]<=pi/6)
    @NLconstraint(m, [k in K,i in W[k]], sum(P[k,j] for j in net[:genmap][i]) + Ps[k,i] - net[:Pd][i] == net[:G][i,i]*V[k,i]^2 + sum( V[k,i]*V[k,j]*( net[:G][i,j] * cos(T[k,i]-T[k,j]) + net[:B][i,j] * sin(T[k,i]-T[k,j])) for j in neighbors(net[:g],i)))
    @NLconstraint(m, [k in K,i in W[k]], sum(Q[k,j] for j in net[:genmap][i]) + Qs[k,i] - net[:Qd][i] ==-net[:B][i,i]*V[k,i]^2 + sum( V[k,i]*V[k,j]*( net[:G][i,j] * sin(T[k,i]-T[k,j]) - net[:B][i,j] * cos(T[k,i]-T[k,j])) for j in neighbors(net[:g],i)))

    @objective(m, Min, sum(net[:c][i,1] + net[:c][i,2]*P[k,i] + net[:c][i,3]*P[k,i]^2 for k in K for i in G[k]) 
               +par[:price]/par[:scale]*sum(Ps_abs[k,i] + Qs_abs[k,i] for k in K for i in W[k]))
    return m
end

function get_partition(net,n_part)
    run(`touch 0.graph`);
    f=open("0.graph","w");
    @printf(f,"%i %i 001\n",nv(net[:g]),ne(net[:g]));
    for i=1:nv(net[:g])
        for j in neighbors(net[:g],i)
            @printf(f,"%i %i ",j,ceil(Int64,sqrt(1000*(net[:B][i,j]^2+net[:G][i,j]^2))))
        end
        @printf(f,"\n")
    end
    close(f)
    run(`gpmetis 0.graph $n_part`)
    part = Int64.(readdlm("0.graph.part.$n_part"))

    run(`rm 0.graph`)
    run(`rm 0.graph.part.$n_part`)
    return [find(part.==i-1) for i=1:n_part]
end

function closed_neighbors(g,V,n)
    nbr = Set([])
    newnbr = Set(V)
    for k=1:n
        for i in newnbr
            union!(nbr, neighbors(g,i))
        end
        newnbr = setdiff(nbr,newnbr)
    end
    return sort(union(V,nbr))
end

function get_subpartition(net,W,Wz,divider)
    submap = []
    # Vmap = nothing
    for k=1:length(W)
        n_subp = ceil(Int64,length(W[k])/divider)
        subg,Vmap = induced_subgraph(net[:g],W[k])        
        run(`touch 1.graph`);
        f=open("1.graph","w");
        @printf(f,"%i %i 001\n",nv(subg),ne(subg));
        for i=1:nv(subg)
            for j in neighbors(subg,i)
                @printf(f,"%i %i ",j,ceil(Int64,sqrt(1000*(net[:B][Vmap[i],Vmap[j]]^2+net[:G][Vmap[i],Vmap[j]]^2))))
            end
            @printf(f,"\n")
        end
        close(f)
        run(`gpmetis 1.graph $n_subp`)
        part = Int64.(readdlm("1.graph.part.$n_subp"))
        run(`rm 1.graph`)
        run(`rm 1.graph.part.$n_subp`)
        for i=1:n_subp
            sub=Vmap[find(part.==i-1)]
            if length(sub)!=0
                push!(submap,sub)
            end
        end
    end

    # Special handling of coupling variables
    tmp = []
    for kk=1:length(submap)
        tmp2=setdiff(submap[kk],Wz)
        if length(tmp2)!=0
            push!(tmp,tmp2)
        end
    end
    tmp = submap
    submap = [submap;[[Wz[i]] for i=1:length(Wz)]]
    subimap= zeros(Int64,net[:n])
    cnt=0
    for kk=1:length(submap)
        for i=1:length(submap[kk])
            subimap[submap[kk][i]]=kk
        end
    end

    subnet=Dict()
    subnet[:g]=Graph(length(submap))
    subnet[:n]=length(submap)
    subnet[:B]=spzeros(subnet[:n],subnet[:n])
    subnet[:G]=spzeros(subnet[:n],subnet[:n])
    subnet[:c]=net[:c]
    subnet[:Pd]=[sum(net[:Pd][submap[i]]) for i=1:subnet[:n]]
    subnet[:Qd]=[sum(net[:Qd][submap[i]]) for i=1:subnet[:n]]
    subnet[:Qmax]=net[:Qmax]
    subnet[:Qmin]=net[:Qmin]
    subnet[:Pmax]=net[:Pmax]
    subnet[:Pmin]=net[:Pmin]
    subnet[:ngen]=net[:ngen]
    subnet[:genmap]=[vcat(net[:genmap][submap[i]]...) for i=1:subnet[:n]]
    subnet[:Vmax]=[maximum(net[:Vmax][submap[i]]) for i=1:subnet[:n]]
    subnet[:Vmin]=[maximum(net[:Vmin][submap[i]]) for i=1:subnet[:n]]
    subnet[:rbus]=subimap[net[:rbus]]

    for i in 1:net[:n]
        for j in [i;neighbors(net[:g],i)]
            ic=subimap[i]
            jc=subimap[j]
            if ic!=jc
                add_edge!(subnet[:g],ic,jc)
            end
            subnet[:B][ic,jc]+=net[:B][i,j]
            subnet[:G][ic,jc]+=net[:G][i,j]
        end
    end

    return submap,subimap,subnet
end

function uniq(v)
    return sort(collect(keys(Set(v).dict)))
end

function printgraph(par, net, subnet, W, subW, submap, subimap)
    Wmap = zeros(Int64,net[:n])
    subWmap = zeros(Int64,subnet[:n])
    for k in 1:par[:n_part]
        for i in W[k]
            Wmap[i] = k
        end
    end
    for k in 1:par[:n_part]
        for i in subW[k]
            subWmap[i] = k
        end
    end
    
    nlist = [["Id" "Partition" "Subpartition"]]
    elist = [["Source" "Target"]]
    for k in 1:subnet[:n]
        for i in submap[k]
            push!(nlist,[string(i) string(Wmap[i]) string(subimap[i])])
        end
    end
    for e in edges(net[:g])
        push!(elist,[string(src(e)) string(dst(e))])
    end

    subnlist = [["Id" "Partition" "Subpartition"]]
    subelist = [["Source" "Target"]]
    for i in 1:subnet[:n]
        push!(subnlist,[string(i) string(subWmap[i]) string(i)])
    end
    for e in edges(subnet[:g])
        push!(subelist,[string(src(e)) string(dst(e))])
    end
    return nlist, elist, subnlist, subelist
end
