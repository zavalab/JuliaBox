LAND_dup = node_matrix[182:1348,1];
open("ag_node_matrix.csv","w") do node
    open("ag_demand_matrix.csv","w") do dem
        open("ag_Plimit_matrix.csv","w") do PP
            open("ag_Nlimit_matrix.csv","w") do NN
                open("ag_node_area_crop.csv","w") do crop
                println(node,"#ag_node",",","latitude",",","longitude",",","contains:");
                println(dem,"#dem",",","node",",","product",",","cap(tonne/week)",",","price(\$/tonne)");
                println(crop,"#ag_node",",","area (acre)",",","crop");
                print(PP, "#node/time");
                print(NN, "#node/time");
                for t in 1:15
                    print(PP,",",t);
                    print(NN,",",t);
                end
                println(PP);
                println(NN);
                for n in FARM1
                    print(PP, n);
                    print(NN, n);
                    for t in 1:15
                        print(PP,",",0);
                        print(NN,",",0);
                    end
                    println(PP);
                    println(NN);
                    println(node, n,",", node_lat[n],",", node_long[n])
                    println(crop, n,",", 0, ",", 0);
                end

                global nn = 0;
                global ii = 0;

                for i in 1:length(LAND)
                    if in(LAND[i],LAND_dup)
                        ag = [LAND[i]];
                        index = Getindex(LAND_dup,LAND[i]);
                        delete_index = [index];
                        nn = nn + 1;
                        for j in 1:length(LAND_dup)
                            if D[LAND[i],LAND_dup[j]] <= 5 && croptype[LAND[i]] == croptype[LAND_dup[j]] && LAND[i] != LAND_dup[j] && length(ag) < 20
                                ag = push!(ag,LAND_dup[j]);
                                delete_index = push!(delete_index,j);
                            end
                        end
                        #splice!(LAND_dup,index);
                        for jj in 1:length(delete_index)
                            splice!(LAND_dup,delete_index[jj]-(jj-1));
                        end
                        #println(delete_index)

                        lat = sum(node_lat[n] for n in ag)/length(ag);
                        long = sum(node_long[n] for n in ag)/length(ag);
                        area = sum(landarea[n] for n in ag);
                        # print node
                        #println(node,"ag_node$nn",",",lat,",",long,",",ag);
                        println(node,"ag_node$nn",",",lat,",",long);
                        println(crop,"ag_node$nn",",",area,",",croptype[LAND[i]]);
                        # print demand
                        for pr in ["p1";"p2";"p3";"p4";"p5";"p6";"p7";"p8";"p9";"p14";"p15";"p16"];
                            price = 0*(pr == "p1" || pr == "p2" || pr == "p3") + 0.000002*(pr == "p4" || pr == "p5" || pr == "p6") + 0.000005*(pr == "p7" || pr == "p8" || pr == "p9") + 0.001*(pr == "p14" || pr == "p15" || pr =="p6");
                            ii = ii + 1;
                            println
                            println(dem,"d$ii",",","ag_node$nn",",",pr,",",1e30,",",",",price);
                            #println(dem,",",",",pr,",",1e30,",",",",price);
                        end
                        # print nutrient
                        print(PP, "ag_node$nn");
                        print(NN, "ag_node$nn");
                        for t in 1:15
                            P = sum(Plimit[n,t] for n in ag);
                            N = sum(Nlimit[n,t] for n in ag);
                            print(PP,",",P);
                            print(NN,",",N)
                        end
                        println(PP);
                        println(NN);
                    end
                end

            for pr in ["p7";"p8";"p9";"p10";"p11";"p12";"p13"];
                price = 0.00002*(pr == "p7" || pr == "p8" || pr == "p9") + 0.00035*(pr == "p10" || pr == "p11" || pr == "p12") + 800*(pr == "p13");
                cap = 1e30 * (pr == "p10" || pr == "p11" || pr == "p12" || pr == "p13") + 1000*(pr == "p7" || pr == "p8" || pr == "p9");
                println(dem,"d$ii",",","company",",",pr,",",cap,",",",",price);
                ii = ii + 1;
            end

            for n in Other
                print(PP, n);
                print(NN, n);
                for t in 1:15
                    print(PP,",",0);
                    print(NN,",",0);
                end
                println(PP);
                println(NN);
                println(node, n, ",", node_lat[n],",", node_long[n]);
                println(crop, n, ",", 0,",",0);
            end

        end
     end
  end
end
end
