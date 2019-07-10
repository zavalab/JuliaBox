println("-------------------------------------------\n");

# Read the model
println("Reading model, ");
@time(include("SC_OPT.jl"));

# SOLVE
println("-------------------------------------------");
println("SOLVER BEGINS");
println("-------------------------------------------");


solve(m)

println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------");

println("SUMMARY");

println("Objective Value (\$):", getobjectivevalue(m))
println("Investment Cost (\$): ", getvalue(Cinv))
println("Transportation Cost (\$): ", getvalue(Ctrans))
println("Operating Cost (\$)): ", getvalue(Cop))
println("Total Revenue (\$): ", getvalue(Revenue))
println("Total Technologies Installed: ", 0)
println("Total Pellet Generated (t): ", getvalue(sum(d[n,pr,t] for n in NODES for pr in ["p10";"p11";"p12"] for t in TIME)))
println("Total Struvite Generated (t): ", getvalue(sum(d[n,"p13",t] for n in NODES for t in TIME)))
println("P from fertilizers (t): ", getvalue(sumferP));
println("N from fertilizers (t): ", getvalue(sumferN));
println("Total excess P (t): ", getvalue(sumTP));
println("Total excess N (t): ", getvalue(sumTN));

println("-------------------------------------------");
println("PRINTING");
println("-------------------------------------------");

open("results_summary.csv","w") do rp
    println(rp, "Objective Value (\$):",",", getobjectivevalue(m))
    println(rp, "Investment Cost (\$): ",",", getvalue(Cinv))
    println(rp, "Transportation Cost (\$): ",",", getvalue(Ctrans))
    println(rp, "Operating Cost (\$)): ",",", getvalue(Cop))
    println(rp, "Total Revenue (\$): ",",", getvalue(Revenue))
    println(rp, "Total Technologies Installed: ",",", 0)
    println(rp, "Total Pellet Generated (t): ",",", getvalue(sum(d[n,pr,t] for n in NODES for pr in ["p10";"p11";"p12"] for t in TIME)))
    println(rp, "Total Struvite Generated (t): ",",", getvalue(sum(d[n,"p13",t] for n in NODES for t in TIME)))
    println(rp, "P from fertilizers (t): ",",", getvalue(sumferP));
    println(rp, "N from fertilizers (t): ",",", getvalue(sumferN));
    println(rp, "Total excess P (t): ",",", getvalue(sumTP));
    println(rp, "Total excess N (t): ",",", getvalue(sumTN));
end

open("technology_sites.csv", "w") do lp
	println(lp, "#Node Number", ",", "Latitude",",","Longitude",",", "Technology", ",", "Capacity")
	for n in NODES_C
		for t in TECHS
			#if round(Int, getvalue(y[n,t])) == 1
			#println(lp, n, ",", node_lat[n], ",", node_long[n],",",t, ",", tech_cap[t])
			#end
		end
	end
end

open("demand_results.csv","w") do dp
    println(dp,"node",",","product",",","amount(tonne/week)",",","time");
    for n in NODES
        for pr in PRODS
            for t in TIME
                if getvalue(d[n,pr,t]) >= 0.001
                    println(dp, n,",", pr,",",getvalue(d[n,pr,t]),",",t);
                end
            end
        end
    end
end

open("inventory_results.csv","w") do ip
	println(ip, "node",",", "product",",","amount(tonne/week)",",","time");
	for n in NODES
		for pr in PRODS
			for t in TIME
				if getvalue(I[n,pr,t]) >= 0.001
					println(ip, n, ",",pr,",",getvalue(I[n,pr,t]),",",t);
				end
			end 
		end
	end
end




open("P_results.csv","w") do PP
    open("N_results.csv","w") do NN
        print(PP, "#node/time");
        print(NN, "#node/time");
        for t in TIME
            print(PP,",",t);
            print(NN,",",t);
        end
        println(PP);
        println(NN);
        for n in NODES
            print(PP, n);
            print(NN, n);
            for t in TIME
                print(PP,",",getvalue(TP[n,t]));
                print(NN,",",getvalue(TN[n,t]));
            end
            println(PP);
            println(NN);
        end
    end
end

for pp in PRODS
    open("flow_results_"*"$(pp)"*".csv","w") do ff
	print(ff,",")
	for j in NODES
	    print(ff,j,",")
	end
	println(ff,"")
	for j in NODES 
	    print(ff,j,",")
	    for k in NODES
		print(ff,getvalue(f[j,k,pp,8]),",")
	    end
	    println(ff)
	end
    end 
end
