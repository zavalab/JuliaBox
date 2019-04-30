## Market coordination case study
## Market analysis of the rock river area

println("----------Printing result summary----------");
open("summary_results.csv", "a") do rp
    println(rp,"Solving time",",", "$(Dates.year(now()))-$(Dates.month(now()))-$(Dates.day(now())) $(Dates.hour(now())):$(Dates.minute(now()))");
    println(rp,"Total social welfare (million \$\year)",",",getvalue(m[:swf]));
    println(rp,"Total revenue (million \$/year)",",", getvalue(m[:demrevn]));
    println(rp,"Total supply cost (million \$/year)",",", getvalue(m[:supcost]));
    println(rp,"Total transportation cost (million \$/year)",",", getvalue(m[:transcost]));
    println(rp,"Total technology cost (million \$/year)",",", getvalue(m[:opcost]));
    println(rp,"Total supplier profit (million \$/year)",",", sum(ϕs[ss] for ss in SUPS));
    println(rp,"Total customer profit (million \$/year)",",", sum(ϕd[dd] for dd in DEMS));
    println(rp,"Total transportation provider profit (million \$/year)",",", sum(ϕl[i,j,pr] for i in NODES for j in NODES for pr in PRODS));
    println(rp,"Total technology provider profit (million \$/year)",",", sum(ϕt[tp] for tp in TECH_PRVD));
    println(rp,"Percentage of suppliers involved (\%)",",", ns/length(SUPS)*100);
    println(rp,"Percentage of customers involved (\%)",",", nd/length(DEMS)*100);
    println(rp,"Percentage of transportation provider involved (\%)",",", ntr/(length(NODES)^2*(length(PRODS)-1))*100);
    println(rp,"Percentage of technology provider involved (\%)",",", ntp/length(TECH_PRVD)*100);
    println(rp,"Ratio of applied P and uptaken P by crops",",", getvalue(sum(Pn) - Pn["n1371"] - Pn["n1372"])/(1388328.73*0.000454));
    println(rp,"\n");
end

println("----------Printing network flows----------");
for pp in PRODS
    open("flow_results_"*"$(pp)"*".csv", "w") do ff
        print(ff,",")
        for j in NODES # Prints the header with node index
            print(ff,j,",")
        end
        println(ff,"") # Used to enter next line
        for j in NODES # Prints the 1st row entry i.e. the sender node
            print(ff,j,",")
            for k in NODES # Prints the flow value from node j to node k with product p
                print(ff,getvalue(f[j,k,pp]),",")
            end
            println(ff)
        end
    end
end

println("----------Printing clearing prices----------");
open("clearing_prices.csv", "w") do cp
    print(cp,",")
    for pr in PRODS # Prints the header with product index
        print(cp,pr,",")
    end
    println(cp,"") # Used to enter next line
    for i in NODES
        print(cp,i,",") # Prints the node locations
        for pr in PRODS # Prints the clearing price of product pr at node i
            print(cp,π[i,pr],",")
        end
        println(cp)
    end
end

println("----------Printing market players profits----------");
open("market_player.csv","w") do mp # We can divide this to several files when there are too many players
    print(mp,"Market player", ",", "Profit (million \$ /year)", ",", "If involved?");
    println(mp)
    for ss in SUPS
        print(mp, ss, ",", ϕs[ss], ",");
        if getvalue(m[:sup][ss]) >= 0.01
            print(mp, "Yes");
        else
            print(mp, "No");
        end
        println(mp);
    end
    for dd in DEMS
        print(mp, dd, ",", ϕd[dd], ",");
        if getvalue(m[:dem][dd]) >= 0.01
            print(mp, "Yes");
        else
            print(mp, "No");
        end
        println(mp);
    end
    for tp in TECH_PRVD
        print(mp, tp, ",", ϕt[tp], ",");
        if getvalue(m[:x][tp_site[tp],tech_refprod[tp_tech[tp]],tp_tech[tp]]) <= -0.01
            print(mp, "Yes");
        else
            print(mp, "No");
        end
        println(mp);
    end
    print(mp, "transportation provider", ",", sum(ϕl[i,j,pr] for i in NODES for j in NODES for pr in PRODS), ",");
    print(mp, ntr/(length(NODES)^2*(length(PRODS)-1))*100);
    println(mp);
end

println("----------Printing supplies and demands----------");
open("supply_results.csv","w") do sp
    print(sp,"supply",",","node",",","prod","bid","value",",","cap");
    println(sp)
    for ss in SUPS
        print(sp, ss, ",", sup_node[ss], ",", sup_prod[ss], ",", sup_bid[ss], ",", getvalue(m[:sup][ss]),",", sup_cap[ss]);
        println(sp)
    end
end
open("demand_results.csv","w") do dp
    print(dp,"demand",",","node",",","prod","bid","value",",","cap");
    println(dp)
    for dd in DEMS
        print(dp, dd, ",", dem_node[dd], ",", dem_prod[dd], ",", dem_bid[dd], ",", getvalue(m[:dem][dd]),",",dem_cap[dd]);
        println(dp)
    end
end

open("P_results.csv","w") do PP
    print(PP, "#node","P-limit","P amount");
    println(PP);
    for n in NODES
        print(PP,n,",",P_limit[n],",",getvalue(Pn[n]));
        println(PP);
    end
end
