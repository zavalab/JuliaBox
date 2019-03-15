println("----------Printing P balances----------");
open("./Summary_Output/p_balances.csv", "w") do cp
    print(cp,"node", ",", "Latitude", ",", "Longitude", ",", "P difference", ",", "Waste Demand")
    println(cp,"") # Used to enter next line
    for i in PNODES
        print(cp, i, ",")
        print(cp,node_lat[i],",") # Prints the node locations
        print(cp,node_long[i],",")
        print(cp, getvalue(t[i]),",")
        print(cp, getvalue(d[i,"p2"]))
        println(cp)
    end
end

open("./Summary_Output/demand_values.csv", "w") do dp
print(dp, "demand",",","node",",", "product",",","value")
              for i in DEMS
                  print(dp, i, ",")
                  print(dp, dem_node[i], ",")
                  print(dp, dem_prod[i], ",")
                  print(dp, getvalue(d[dem_node[i], dem_prod[i]]))
                  println(dp)
              end
end

println("----------Printing network flows----------");
for pp in PRODS
    open("./Summary_Output/flow_results_"*"$(pp)"*".csv", "w") do ff
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
