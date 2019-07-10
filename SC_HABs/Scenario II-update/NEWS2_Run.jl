## Reading data

# Set of nutrients and nutrient forms
E = ["P";"N"];
F = ["DI";"DO"];


# data needed
TIME = 1:15;
A = 529.81;                                 # watershed area (km2)
hw_frem = Dict(zip(E, [0.95;0.85]));        # fraction of E that is removed by WWTP
I = 0;                                      # population in the study area - no WWTP in the region
WShw_E = Dict(zip(E, 2*[0.0105;0.0875]));     # nutrients in human waste (kg/cap/2week)
maxhw_frem_N = 0.9                          # maximum hw_frem["N"] that is observed in the country

rainfall = readdlm("rainfall.csv",',');     # rainfall matrix [mmH2O]

node_matrix = readdlm("ag_node_matrix.csv",',');    # read node matrix
crop_area   = readdlm("ag_node_area_crop.csv",',');
NODES = node_matrix[:,1];                           # nodes
area  = Dict(zip(NODES, crop_area[:,2]));           # in acre;
crop  = Dict(zip(NODES, crop_area[:,3]));


Presult     = readdlm("P_results.csv",',');
Nresult     = readdlm("N_results.csv",',');

open("NEWS2_results.csv","w") do pr
    open("daily_runoff.csv","w") do rr
        global Q_stor0 = 0;                         # initial runoff  stored in last week [mmH2O]
        println(pr,"#TIME",",","N(kg)",",","P(kg)",",","Rnat(mH2O)",",","DIN (P)",",","DIN (NP)",",","DON (P)",",","DON (NP)",",","DIP (P)",",","DIP (NP)",",","DOP (P)",",","DOP (NP)");

for t in TIME
#t = 1;
    ## NEWS2 Nutrient Fate and Transport Model (dissolved part)
    ##Coded by Yicheng Hu 2018-07


    ## Point source part
    RSpnt_E = Dict(zip(E, zeros(length(E),1))); # point source emission of nutrient E (kg/2week)
    for e in E
        RSpnt_E[e] = (1-hw_frem[e])*I*WShw_E[e];
    end

    FE_pntF = Dict(("DI","N") => 0.485 + 0.255*(hw_frem["N"]/ maxhw_frem_N), ("DO","N") => 0.14, ("DI","P") => 1, ("DO","P") => 0.01);
    RSpnt_F = Dict((F[1],E[1]) => 0.1);         # point source emission of nutrient form F (kg/2week)
    for f in F
        for e in E
            RSpnt_F[(f,e)] = FE_pntF[f,e]*RSpnt_E[e];
        end
    end

    ## Non-point source part
    # rainfall and runoff
    # input data, define variables
    R_day   = rainfall[(t-1)*14+1:(t-1)*14+14,2];  # rainfall data in each day [mmH2O]
    S       = zeros(14,1);                       # rainfall retention parameter [mmH2O]
    CN      = 70*ones(14,1);                     # curve number for the day
    Qp_surf = zeros(14,1);                       # generated runoff in each day [mmH2O]
    Q_surf  = zeros(14,1);                       # runoff that reached river in each day [mmH2O]
    Q_stor  = zeros(14,1);                       # runoff store in each day [mmH2O]
    surlag  = 3.64                              # surface runoff lag coefficient (using global average) [day]
    tc      = zeros(14,1);                       # time of concentration [hr]
    l       = 109361;                           # we use 1/3 of the Yahara river  [ft]
    Y       = (318-236)/55.66/10                # average slope of Yahara river [%]
    for i in 1:14
        S[i]        = 25.4*(1000/CN[i] - 10);
        tc[i]       = l^0.8*(S[i]/25.4+1)^0.7/(1140*Y^0.5)/24;          # convert to [day]
        Qp_surf[i]  = (R_day[i] - 0.2*S[i])^2/(R_day[i] + 0.8*S[i]);
        if i == 1
            Q_surf[i] = (Qp_surf[i] + Q_stor0)*(1-exp(-surlag/tc[i]));
            Q_stor[i] = Q_stor0 + Qp_surf[i] - Q_surf[i];
        else
            Q_surf[i] = (Qp_surf[i] + Q_stor[i-1])*(1-exp(-surlag/tc[i]));
            Q_stor[i] = Q_stor[i-1] + Qp_surf[i] - Q_surf[i];
        end
        println(rr,Q_surf[i]);
    end
    Q_stor0 = Q_stor[end];
    Rnat    = sum(Q_surf)*52.1/1000/2           # annualize the runoff [mH2O]

    # the inexplicit non-point source
    f_F = Dict(("DI","N") => Rnat^1, ("DO","N") => Rnat^0.95, ("DI","P") => (1+(Rnat/0.85)^(-2))^-1, ("DO","P") => Rnat^0.95);
    EC_F = Dict(("DI","N") => 0, ("DO","N") => 280, ("DI","P") => 26, ("DO","P") => 15);
    RSdif_ecF = Dict((F[1],E[1]) => 0.1);      # emission to rivers from inexplicit non-point sources (kg/year)
    for e in E
        for f in F
            RSdif_ecF[(f,e)] = A*f_F[f,e]*EC_F[f,e]/52.1*2;  # (kg/2week)
        end
    end


    # the explicit non-point source
    # human activity part
    WSdif_feE = Dict((NODES[1], E[1]) => 0);          # non point source due to fertilizers in each node
    WSdif_maE = Dict((NODES[1], E[1]) => 0.1);          # non point source due to animal manure
    WSdif_fixantE = Dict((NODES[1], E[1]) => 0.1);      # non point source due to fixation by crops
    WSdif_depantE = Dict((NODES[1], E[1]) => 0.1);      # non point source due to deposition by crops
    WSdif_exE = Dict((NODES[1], E[1]) => 0);          # withdraw of E in crops and animal grazing [kg E/week]
    for i in 1:length(NODES)
        for e in 1:length(E)
            WSdif_feE[(NODES[i],E[e])] = 0;
            WSdif_maE[(NODES[i],E[e])] = (E[e] == "P")*Presult[i,t+1]*1000 + (E[e] == "N")*Nresult[i,t+1]*1000;
            WSdif_fixantE[(NODES[i],E[e])] = (E[e] == "P")*0 + (E[e] == "N")*area[NODES[i]]*((crop[NODES[i]] == 5)*25/20.71*2/2.471054 + (crop[NODES[i]] == 1)*5/24.14*2/2.471054 +(crop[NODES[i]] == 36)*5/17.43*2/2.471054 +(crop[NODES[i]] == 176)*5/16.14*2/2.471054 +(crop[NODES[i]] == 24)*5/17.57*2/2.471054*0.7) ;  # 25kg/ha/year, 20.71 week, double week, 1 ha = 2.47 ac
            WSdif_depantE[(NODES[i],E[e])] = (E[e] == "P")*0 + (E[e] == "N")*area[NODES[i]]*0.04795*14/2.471054;
            WSdif_exE[(NODES[i],E[e])] = 0;
        end
    end

    WSdif_antE = Dict(zip(E, zeros(2,1)));      # diffuse source due to human activities [kg/2week]
    for e in E
        WSdif_antE[e] = 52.1/2*(sum(max(WSdif_feE[n,e] + WSdif_maE[n,e] + WSdif_fixantE[n,e] + WSdif_depantE[n,e] - WSdif_exE[n,e], 0) for n in NODES) - max(WSdif_feE[NODES[270],e] + WSdif_maE[NODES[270],e] + WSdif_fixantE[NODES[270],e] + WSdif_depantE[NODES[270],e] - WSdif_exE[NODES[270],e], 0)); #annualize
    end


    e_F = Dict(("DI","N") => 0.94, ("DO","N") => 0.01, ("DI","P") => 0.29, ("DO","P") => 0.01);
    FE_wsF = Dict((F[1],E[1]) => 0.1);
    for f in F
        for e in E
            FE_wsF[(f,e)] =  f_F[f,e]* e_F[f,e];
        end
    end

    # natural activity part
    WSdif_natE = Dict(zip(E,zeros(2,1)));       # diffuse source due to natural activities [kg/2week]
    for e in E
        WSdif_natE[e] = 1/6*WSdif_antE[e];      # assume the natural nutrient release is 1/6 of the anthropogenic part (cite report)
    end
    e_natF = Dict(("DI","N") => 0.1, ("DO","N") => 0.01, ("DI","P") => 0.29, ("DO","P") => 0.01);
    FE_wsnatF = Dict((F[1],E[1]) => 0.1);
    for f in F
        for e in E
            FE_wsnatF[(f,e)] =  f_F[f,e]* e_natF[f,e];
        end
    end

    RSdif_explF = Dict((F[1],E[1]) => 0.1);    # emission to rivers from explicit non-point sources (kg/week)
    for e in E
        for f in F
            RSdif_explF[(f,e)] = (FE_wsnatF[(f,e)]*WSdif_natE[e] + FE_wsF[(f,e)]*WSdif_antE[e])/52.1*2;
        end
    end

    Yield = Dict((F[1],E[1]) => 0.1);           # total emission to rivers (kg/2week)
    for e in E
        for f in F
            #Yield[f,e] = RSdif_explF[(f,e)] + RSpnt_F[(f,e)];
            Yield[f,e] = RSdif_ecF[(f,e)] + RSdif_explF[(f,e)];
        end
    end

    println(pr, t,",", Yield["DI","N"]+Yield["DO","N"],",",Yield["DI","P"]+Yield["DO","P"],",",Rnat,",",RSpnt_F["DI","N"],",",RSdif_explF["DI","N"],",",RSpnt_F["DO","N"],",",RSdif_explF["DO","N"],",",RSpnt_F["DI","P"],",",RSdif_explF["DI","P"],",",RSpnt_F["DO","P"],",",RSdif_explF["DO","P"]);
end
end
end
