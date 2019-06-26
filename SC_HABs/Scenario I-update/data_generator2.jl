node_matrix       = readdlm("node_matrix.csv",',');
nutrient_info     = readdlm("nutrients_analysis.csv",',');

FARM1 = node_matrix[1:181,1];
LAND = node_matrix[182:1348,1];
Other = node_matrix[1349:end,1];

croptype = Dict(zip(LAND, nutrient_info[:,2]));
Pneed    = Dict(zip(LAND, nutrient_info[:,5]));
Nneed    = Dict(zip(LAND, nutrient_info[:,8]));

open("Plimit_matrix.csv","w") do PP
    open("Nlimit_matrix.csv","w") do NN
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
        end
        for i in LAND
            print(PP, i);
            print(NN, i);
            for t in 1:15
                if croptype[i] == 1
                    if 3 <= t <= 15
                        print(PP,",",2*Pneed[i]*0.453592);   # convert lb to kg
                        print(NN,",",2*Nneed[i]*0.453592);
                    else
                        print(PP,",",0);
                        print(NN,",",0);
                    end
                elseif croptype[i] == 36
                    if 1 <= t <= 9
                        print(PP,",",2*Pneed[i]*0.453592);
                        print(NN,",",2*Nneed[i]*0.453592);
                    else
                        print(PP,",",0);
                        print(NN,",",0);
                    end
                elseif croptype[i] == 176
                    if 1 <= t <= 8
                        print(PP,",",2*Pneed[i]*0.453592);
                        print(NN,",",2*Nneed[i]*0.453592);
                    else
                        print(PP,",",0);
                        print(NN,",",0);
                    end
                elseif croptype[i] == 5
                    if 4 <= t <= 11
                        print(PP,",",2*Pneed[i]*0.453592);
                        print(NN,",",2*Nneed[i]*0.453592);
                    else
                        print(PP,",",0);
                        print(NN,",",0);
                    end
                elseif croptype[i] == 24
                    if 1 <= t <= 9
                        print(PP,",",2*Pneed[i]*0.453592);
                        print(NN,",",2*Nneed[i]*0.453592);
                    else
                        print(PP,",",0);
                        print(NN,",",0);
                    end
                end
            end
            println(PP);
            println(NN);
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
        end
    end
end
