## Market coordination case study
## Market analysis of the rock river area

println("----------Reading Model-----------")
tic()
include("market_model.jl"); # Reading model file
toc()
println("----------Begin solving-----------")

solve(m) # Solve model (may get unbounded result when lambda is small)

solve(m) # Solve again (solve again may avoid that)

# Analyze results
# dual variable - clearing price
π = Dict((NODES[1],PRODS[1]) => 0.5);
for i in NODES
    for pr in PRODS
        π[(i,pr)] = -getdual(balance[i,pr]);
    end
end

# count players
ns = sum(getvalue(sup[ss]) >= 0.01 for ss in SUPS)
nd = sum(getvalue(dem[dd]) >= 0.01 for dd in DEMS)
ntr = sum(getvalue(f[i,j,pr]) >= 0.01 for i in NODES for j in NODES for pr in PRODS)
ntp = 0;
for tp in TECH_PRVD
    if getvalue(m[:x][tp_site[tp],tech_refprod[tp_tech[tp]],tp_tech[tp]]) <= -0.01
        ntp = ntp + 1;
    end
end

# calculate market players' profits
ϕs = zeros(length(SUPS),1);
ϕs = Dict(zip(SUPS,ϕs));
for ss in SUPS
    ϕs[ss] = getvalue(m[:sup][ss])*sum(π[i,pr] - sup_bid[ss] for i in NODES for pr in PRODS
                    if sup_node[ss] == i && sup_prod[ss] == pr);
end

ϕd = zeros(length(DEMS),1);
ϕd = Dict(zip(DEMS,ϕd));
for dd in DEMS
    ϕd[dd] = getvalue(m[:dem][dd])*sum(dem_bid[dd] - π[i,pr] for i in NODES for pr in PRODS
                    if dem_node[dd] == i && dem_prod[dd] == pr);
end

ϕl = Dict((NODES[1], NODES[2],PRODS[1]) => 0.5)
for i in NODES
    for j in NODES
        for pr in PRODS
            ϕl[i,j,pr] = getvalue(m[:f][i,j,pr])*sum(π[j,pr] - π[i,pr]
                            - prod_trans[pr]*distance[i,j]);
        end
    end
end

ϕt = zeros(length(TECH_PRVD),1);
ϕt = Dict(zip(TECH_PRVD,ϕt));
for tp in TECH_PRVD
    ϕt[tp] = sum(π[tp_site[tp],pr]*getvalue(m[:x][tp_site[tp],pr,tp_tech[tp]]) for pr in PRODS) +
                    getvalue(m[:x][tp_site[tp],tech_refprod[tp_tech[tp]],tp_tech[tp]])*tech_bid[tp_tech[tp]];
end

# print main results
println("\n----------Print results----------")
println("Total social welfare (million \$\year): \t\t\t\t\t", getvalue(m[:swf]));
println("Total revenue (million \$/year): \t\t\t\t\t", getvalue(m[:demrevn]));
println("Total supply cost (million \$/year): \t\t\t\t\t", getvalue(m[:supcost]));
println("Total transportation cost (million \$/year): \t\t\t\t", getvalue(m[:transcost]));
println("Total technology cost (million \$/year): \t\t\t\t", getvalue(m[:opcost]));
println("Total supplier profit (million \$/year): \t\t\t\t", sum(ϕs[ss] for ss in SUPS));
println("Total customer profit (million \$/year): \t\t\t\t", sum(ϕd[dd] for dd in DEMS));
println("Total transportation provider profit (million \$/year): \t\t\t", sum(ϕl[i,j,pr] for i in NODES for j in NODES for pr in PRODS));
println("Total technology provider profit (million \$/year): \t\t\t", sum(ϕt[tp] for tp in TECH_PRVD));
println("Percentage of suppliers involved (\%): \t\t\t\t\t", ns/length(SUPS)*100);
println("Percentage of customers involved (\%): \t\t\t\t\t", nd/length(DEMS)*100);
println("Percentage of transportation provider involved (\%): \t\t\t", ntr/(length(NODES)^2*(length(PRODS)-1))*100);
println("Percentage of technology provider involved (\%): \t\t\t", ntp/length(TECH_PRVD)*100);
println("Ratio of applied P and uptaken P by crops: \t\t\t\t", getvalue(sum(Pn) - Pn["n1371"] - Pn["n1372"])/(1388328.73*0.000454));
