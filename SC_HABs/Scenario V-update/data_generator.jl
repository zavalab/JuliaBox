node_matrix       = readdlm("node_matrix.csv",',');

LAND = node_matrix[182:1348,1];

open("demand_matrix.csv","w") do dd
    println(dd,"#dem",",","node",",","product",",","cap(tonne/week)",",","time",",","price(\$/tonne)");
    i = 1;
    for t in 1:30
        for pr in ["p1";"p2";"p3";"p4";"p5";"p6";"p7";"p8";"p9";"p14";"p15";"p16"];
            for n in LAND
                price = 0*(pr == "p1" || pr == "p2" || pr == "p3") + 0.000002*(pr == "p4" || pr == "p5" || pr == "p6") + 0.000005*(pr == "p7" || pr == "p8" || pr == "p9") + 0.001*(pr == "p14" || pr == "p15" || pr =="p6");
                println(dd,"d$i",",",n,",",pr,",",1e30,",",t,",",price);
                i = i + 1;
            end
        end
        for pr in ["p7";"p8";"p9";"p10";"p11";"p12";"p13"];
            price = 0.00002*(pr == "p7" || pr == "p8" || pr == "p9") + 0.00035*(pr == "p10" || pr == "p11" || pr == "p12") + 800*(pr == "p13");
            cap = 1e30 * (pr == "p10" || pr == "p11" || pr == "p12" || pr == "p13") + 1000*(pr == "p7" || pr == "p8" || pr == "p9");
            println(dd,"d$i",",","company",",",pr,",",cap,",",t,",",price);
            i = i + 1;
        end
    end
end
